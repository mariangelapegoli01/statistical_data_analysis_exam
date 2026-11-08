library(tidyverse)
library(mice)
library(glmnet)
library(survival)

set.seed(42)

B                   <- 200    # Number of bootstrap resamples for internal validation
CNA_AS_PREDICTORS <- TRUE     # FALSE = CNAs (except del_17p) excluded as predictors in MICE;
                              # TRUE  = only predefined CNA variables are excluded

dataset <- read_rds("Clinical_Genomic_Dataset_Pre_Imputation.rds")

vars_genomic <- grep("^(del_|amp_)", colnames(dataset), value = TRUE)  
vars_clinical <- setdiff(colnames(dataset),
                         c("PFS_I_MONTHS", "PFS_I_EVENT", "MPC",
                           vars_genomic))

cat("Clinical variables (", length(vars_clinical), "):",
    paste(vars_clinical, collapse = ", "), "\n")
cat("Genomic variables (", length(vars_genomic), "):",
    paste(vars_genomic, collapse = ", "), "\n")

all_vars <- c(vars_clinical, vars_genomic)

alphas <- c(Integrated_Lasso = 1.0, Integrated_EN = 0.5, Integrated_Ridge = 0.0)
nm     <- names(alphas)

# Define CNA exclusion rules for the MICE imputation model
cna_cols     <- setdiff(vars_genomic, "del_17p")
cols_no_pred <- if (CNA_AS_PREDICTORS) c("amp_7p", "amp_5p", "amp_18p", "amp_3p") else cna_cols

# Construct the design matrix used for penalized Cox regression
make_x <- function(df) {
  vars_ok <- intersect(all_vars, colnames(df))
  model.matrix(reformulate(vars_ok), data = df)[, -1, drop = FALSE]
}

# Align predictors across bootstrap samples
align_x <- function(x_raw, master) {
  out <- matrix(0, nrow = nrow(x_raw), ncol = length(master),
                dimnames = list(NULL, master))
  comm <- intersect(colnames(x_raw), master)
  if (length(comm) > 0) {
    out[, comm] <- x_raw[, comm, drop = FALSE]
  }
  out
}

# Compute Harrell's concordance index (C-index) from the model linear predictor.
cindex <- function(lp, time, event) {
  if (sd(lp) == 0) return(0.5)
  d <- data.frame(t = time, e = event, inv = -lp)
  concordance(Surv(t, e) ~ inv, data = d)$concordance
}

# Estimate the calibration slope using the linear predictor
slope_fun <- function(lp, time, event) {
  if (sd(lp) == 0) return(NA_real_)
  d <- data.frame(t = time, e = event, lp = lp)
  unname(coef(coxph(Surv(t, e) ~ lp, data = d)))
}

# Fit cross-validated penalized Cox model
fit_penalty <- function(x, t, e, foldid, alphas) {
  y <- Surv(t, e)
  lapply(alphas, function(a) {
    cv <- cv.glmnet(x, y, family = "cox", alpha = a,
                    foldid = foldid, cox.ties = "breslow")
    setNames(as.numeric(coef(cv, s = "lambda.min")), colnames(x))
  })
}

# Match coefficients to the reference predictor set
to_master <- function(b, master) {
  out  <- setNames(numeric(length(master)), master)
  comm <- intersect(names(b), master)
  out[comm] <- b[comm]
  out
}

# Construct the predictor matrix used for MICE imputation.
build_pred <- function(dat, cols_no_pred = character(0)) {
  pm <- make.predictorMatrix(dat)
  pm[, "PFS_I_MONTHS"] <- 0                  # Raw survival time excluded as a direct predictor
  pm[, c("chaz", "PFS_I_EVENT")] <- 1        # Auxiliary Nelson-Aalen hazard and event status included
  
  if (length(cols_no_pred) > 0) {
    cols_excl <- intersect(cols_no_pred, colnames(pm))
    pm[, cols_excl] <- 0                     # Deactivate specified CNAs from imputing other variables
  }
  
  diag(pm) <- 0
  pm[c("PFS_I_MONTHS", "PFS_I_EVENT", "chaz"), ] <- 0
  pm
}

# Impute missing values with MICE
impute_and_build <- function(dat, cols_no_pred, seed) {
  dat$chaz <- nelsonaalen(dat, PFS_I_MONTHS, PFS_I_EVENT)
  imp <- mice(dat, m = 1, predictorMatrix = build_pred(dat, cols_no_pred),
              printFlag = FALSE, seed = seed)
  d <- complete(imp, 1)
  list(x = make_x(d), t = d$PFS_I_MONTHS, e = d$PFS_I_EVENT)
}

cat("\n--- DATA PREPARATION AND IMPUTATION ---\n")

# Replace non-positive survival times with a small positive value
dataset <- dataset %>% mutate(PFS_I_MONTHS = ifelse(PFS_I_MONTHS <= 0, 0.1, PFS_I_MONTHS))

dati  <- dataset %>% select(-any_of("MPC"))
n_pat <- nrow(dati)

# Generate the imputed datasets from the original data
ORIG   <- impute_and_build(dati, cols_no_pred, seed = 123)

# Define reference master feature space from the first imputed original dataset
master <- colnames(ORIG$x)
p_cols <- length(master)

cat("Patients:", n_pat, "| Events:", sum(dati$PFS_I_EVENT), "| Columns:", p_cols, "\n")

set.seed(1)
foldid_orig <- sample(rep(1:10, length.out = n_pat))

# Fit penalized models on the imputed original data
beta_orig <- fit_penalty(ORIG$x, ORIG$t, ORIG$e, foldid_orig, alphas)

#Coefficients and variable-selection results
sel_orig  <- lapply(beta_orig, function(b) as.numeric(b != 0))

# Compute linear predictor, apparent discrimination (C-index) and calibration slope 
lp_orig        <- lapply(beta_orig, function(b) as.numeric(ORIG$x %*% b))
apparent_c     <- sapply(lp_orig, cindex,    time = ORIG$t, event = ORIG$e)
apparent_slope <- sapply(lp_orig, slope_fun, time = ORIG$t, event = ORIG$e)

cat("\n---  BOOTSTRAP VALIDATION ---\n")
mk <- function() matrix(NA_real_, B, length(nm), dimnames = list(NULL, nm))
app_mat       <- mk()
test_mat      <- mk()
slope_mat     <- mk()   # Test* calibration slope (used to estimate the uniform shrinkage factor)
slope_app_mat <- mk()   # Apparent* calibration slope

# Store bootstrap coefficient
boot_betas <- lapply(setNames(nm, nm), function(s)
  matrix(NA_real_, length(master), B, dimnames = list(master, NULL)))

for (b in seq_len(B)) {
  if (b %% 25 == 0 || b == 1 || b == B) {
    cat("Bootstrap iteration", b, "of", B, "\n")
  }
  
  # Resample incomplete data with replacement
  set.seed(10000 + b)
  idx       <- sample(seq_len(n_pat), replace = TRUE)
  foldid_b  <- sample(rep(1:10, length.out = n_pat))[idx]
  data_boot <- dati[idx, ]
  
  # Re-run imputation inside the bootstrap sample
  BOOT <- impute_and_build(data_boot, cols_no_pred, seed = 1000 + b)
  
  # Fit models on the bootstrap data
  x_b    <- align_x(BOOT$x, master)
  beta_b <- fit_penalty(x_b, BOOT$t, BOOT$e, foldid_b, alphas)
  
  # Evaluate Apparent* (bootstrap data) and Test* (original data) performance
  for (s in nm) {
    beta <- beta_b[[s]]
    boot_betas[[s]][, b] <- beta
    
    lp_app  <- as.numeric(x_b %*% beta)
    lp_test <- as.numeric(ORIG$x %*% beta)
    
    app_mat[b, s]       <- cindex(lp_app,  BOOT$t, BOOT$e)
    test_mat[b, s]      <- cindex(lp_test, ORIG$t, ORIG$e)
    slope_app_mat[b, s] <- slope_fun(lp_app,  BOOT$t, BOOT$e)
    slope_mat[b, s]     <- slope_fun(lp_test, ORIG$t, ORIG$e)
  }
}

# Calculate optimism in C-index
optimism_mat <- app_mat - test_mat

# Calculate optimism-corrected calibration slope
slope_ott  <- colMeans(slope_app_mat - slope_mat, na.rm = TRUE)
slope_corr <- apparent_slope[nm] - slope_ott

cat("\n--- RESULTS SUMMARY ---\n")

final_summary <- data.frame(
  Set             = nm,
  Colonne_p       = p_cols,
  C_Apparente     = apparent_c[nm],
  Slope_Apparente = apparent_slope[nm],
  Ottimismo_Medio = colMeans(optimism_mat, na.rm = TRUE),
  C_Validato      = apparent_c[nm] - colMeans(optimism_mat, na.rm = TRUE),
  Shrinkage_s     = colMeans(slope_mat, na.rm = TRUE),
  Slope_App_Boot  = colMeans(slope_app_mat, na.rm = TRUE),
  Slope_Ottimismo = slope_ott[nm],
  Slope_Validata  = slope_corr[nm],
  row.names = NULL
) %>% mutate(across(where(is.numeric) & !c(Colonne_p), ~ round(., 3)))

cat("\n--- PERFORMANCE SUMMARY (PENALTY COMPARISON) ---\n")
print(final_summary, row.names = FALSE)

# Generate a coefficient summary table
make_coef_table <- function(set) {
  cols_set <- master
  sh       <- mean(slope_mat[, set], na.rm = TRUE)
  sh_ok    <- ifelse(is.na(sh), 1, sh)
  
  b       <- beta_orig[[set]][cols_set]
  b_mat   <- boot_betas[[set]][cols_set, , drop = FALSE]
  
  # Variability conditional on variable selection
  # (computed only among non-zero coefficients)
  sd_cond <- apply(b_mat, 1, function(x) {
    nz <- x[!is.na(x) & x != 0]
    if (length(nz) > 1) sd(nz) else NA_real_
  })
  
  data.frame(Variabile    = cols_set,
             Beta_Orig    = b,
             HR     = exp(b),
             Media_Boot   = apply(b_mat, 1, mean, na.rm = TRUE),
             SD_Boot      = apply(b_mat, 1, sd,   na.rm = TRUE),
             SD_Boot_Cond = sd_cond,
             Sel_Orig = as.numeric(b != 0),
             Boot_Sel_Pct = rowMeans(b_mat != 0, na.rm = TRUE) * 100,
             row.names    = NULL) %>%
    filter(Beta_Orig != 0 | Boot_Sel_Pct > 0) %>%
    arrange(desc(Boot_Sel_Pct), desc(abs(Beta_Orig))) %>%
    mutate(across(c(Beta_Orig, HR, Media_Boot, SD_Boot, SD_Boot_Cond), ~ round(., 3)),
           Boot_Sel_Pct = round(Boot_Sel_Pct, 1))
}

# Create coefficient summary tables for all penalization strategies
coef_tables <- setNames(lapply(nm, make_coef_table), nm)

for (s in nm) {
  cat("\n--- COEFFICIENTS & BOOTSTRAP SD:", s, "---\n")
  print(coef_tables[[s]], row.names = FALSE)
}

cat("\n--- PROPORTIONAL HAZARDS ASSUMPTION TEST ---\n")

d <- ORIG
ph_lp <- ph_vars <- fit_lp_list <- zph_lp_list <- list()

for (s in nm) {
  b        <- beta_orig[[s]]
  sel_vars <- names(b[b != 0])
  
  # A. Linear predictor
  lp     <- as.numeric(d$x %*% b)
  fit_lp <- coxph(Surv(d$t, d$e) ~ lp)
  zph_lp <- cox.zph(fit_lp)
  fit_lp_list[[s]] <- fit_lp
  zph_lp_list[[s]] <- zph_lp
  ph_lp[[s]] <- data.frame(Set = s,
                           Chisq   = zph_lp$table["lp", "chisq"],
                           P_value = zph_lp$table["lp", "p"])
  
  # B. Selected covariates
  df_cox <- data.frame(t = d$t, e = d$e, d$x[, sel_vars, drop = FALSE])
  tab    <- cox.zph(coxph(Surv(t, e) ~ ., data = df_cox))$table
  p_cov  <- tab[sel_vars, "p"]
  ph_vars[[s]] <- data.frame(
    Set        = s,
    Variabile  = c(sel_vars, "GLOBAL", "LINEAR_PREDICTOR"),
    P_value    = c(p_cov, tab["GLOBAL", "p"], zph_lp$table["lp", "p"]),
    P_value_BH = c(p.adjust(p_cov, method = "BH"), NA, NA),
    row.names  = NULL
  )
}

cat("\n--- PH TEST ON THE LINEAR PREDICTOR ---\n")
print(bind_rows(ph_lp) %>% mutate(across(where(is.numeric), ~ round(., 4))), row.names = FALSE)

cat("\n--- PH TEST BY COVARIATE (BH-adjusted within model), GLOBAL, AND LP ---\n")
print(bind_rows(ph_vars) %>% mutate(across(c(P_value, P_value_BH), ~ round(., 4))), row.names = FALSE)

# Scaled Schoenfeld residual plots for the linear predictor
par(mfrow = c(1, length(nm)), mar = c(4.5, 4.5, 3, 1))
for (s in nm) {
  plot(zph_lp_list[[s]],
       main = paste0(s, "\n(p = ", round(ph_lp[[s]]$P_value, 3), ")"),
       ylab = "Beta(t) for Linear Predictor", xlab = "Time (Months)")
  abline(h = 0, col = "gray60", lty = 3)
  abline(h = coef(fit_lp_list[[s]]), col = "red", lty = 2, lwd = 1.5)
}
par(mfrow = c(1, 1))

cat("\n---  DIAGNOSTIC CHECKS ---\n")
cat("NA count in slope_mat (test*):\n");         print(colSums(is.na(slope_mat)))
cat("NA count in slope_app_mat (apparent*):\n"); print(colSums(is.na(slope_app_mat)))
