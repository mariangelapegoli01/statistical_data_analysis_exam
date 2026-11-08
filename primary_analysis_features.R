library(tidyverse)
library(mice)
library(glmnet)
library(survival)

set.seed(42)

B                   <- 200     # Number of Bootstrap Samples
CNA_AS_PREDICTORS   <- TRUE    # FALSE = all CNAs (except del_17p) excluded as predictors in MICE; 
                               #TRUE = only selected CNAs are excluded as predictors 

dataset <- read_rds("Clinical_Genomic_Dataset_Pre_Imputation.rds")

vars_genomic <- grep("^(del_|amp_)", colnames(dataset), value = TRUE)  

vars_clinical <- setdiff(colnames(dataset),
                         c("PFS_I_MONTHS", "PFS_I_EVENT", "MPC",
                           vars_genomic))

cat("Clinical variables (", length(vars_clinical), "):",
    paste(vars_clinical, collapse = ", "), "\n")
cat("Genomic variables (", length(vars_genomic), "):",
    paste(vars_genomic, collapse = ", "), "\n")

feature_sets <- list(
  Clinical   = vars_clinical,
  Genomic    = vars_genomic,
  Integrated = c(vars_clinical, vars_genomic)
)
all_vars <- unique(unlist(feature_sets))

# Define CNA exclusion rules for the MICE imputation model (retaining del_17p)
cna_cols     <- setdiff(vars_genomic, "del_17p")
cols_no_pred <- if (CNA_AS_PREDICTORS) c("amp_7p", "amp_5p", "amp_18p", "amp_3p") else cna_cols

# Construct the design matrix used for penalized Cox regression
make_x <- function(df) {
  vars_ok <- intersect(all_vars, colnames(df))
  mm <- model.matrix(reformulate(vars_ok), data = df)
  x  <- mm[, -1, drop = FALSE]
  attr(x, "var") <- vars_ok[attr(mm, "assign")[-1]]
  x
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
fit_all_specs <- function(x, y, foldid, set_cols, alpha = 1) {
  lapply(set_cols, function(cols) {
    cv <- cv.glmnet(x[, cols, drop = FALSE], y, family = "cox",
                    alpha = alpha, foldid = foldid, cox.ties = "breslow")
    b <- setNames(numeric(ncol(x)), colnames(x))
    b[cols] <- as.numeric(coef(cv, s = "lambda.min"))
    b
  })
}

# Construct the predictor matrix used for MICE imputation.
build_pred <- function(dat, cols_no_pred = character(0)) {
  pm <- make.predictorMatrix(dat)
  pm[, "PFS_I_MONTHS"] <- 0             # Raw survival time excluded as a direct predictor                 
  pm[, c("chaz", "PFS_I_EVENT")] <- 1   # Auxiliary Nelson-Aalen hazard and event status included       
  
  if (length(cols_no_pred) > 0) {
    pm[, cols_no_pred] <- 0             # Deactivate specified CNAs from imputing other variables             
  }
  
  diag(pm) <- 0
  pm[c("PFS_I_MONTHS", "PFS_I_EVENT", "chaz"), ] <- 0
  pm
}

# Impute missing values with MICE
impute_and_build <- function(dat, cols_no_pred, seed, label = "MICE", verbose_log = TRUE) {
  dat$chaz <- nelsonaalen(dat, PFS_I_MONTHS, PFS_I_EVENT)
  imp <- mice(dat, m = 1, predictorMatrix = build_pred(dat, cols_no_pred),
              printFlag = FALSE, seed = seed)
  
  if (verbose_log) {
    if (!is.null(imp$loggedEvents) && nrow(imp$loggedEvents) > 0) {
      cat("\n[", label, "] MICE ha registrato", nrow(imp$loggedEvents), "logged events:\n")
      print(imp$loggedEvents)
    } else cat("\n[", label, "] Nessun logged event registrato da MICE.\n")
  }
  
  d <- complete(imp, 1)
  list(x = make_x(d), t = d$PFS_I_MONTHS, e = d$PFS_I_EVENT)
}

# Replace non-positive survival times with a small positive value
dataset <- dataset %>% mutate(PFS_I_MONTHS = ifelse(PFS_I_MONTHS <= 0, 0.1, PFS_I_MONTHS))

dati  <- dataset %>% select(-MPC)
n_pat <- nrow(dati)

# Generate the original imputed dataset
ORIG <- impute_and_build(dati, cols_no_pred, seed = 123, label = "ORIGINAL DATA")
logged_orig <- attr(ORIG, "loggedEvents")

# Define reference master feature space from the first imputed original dataset
master        <- colnames(ORIG$x)
var_of_master <- attr(ORIG$x, "var")
set_cols      <- lapply(feature_sets, function(vs) master[var_of_master %in% vs])

print(sapply(set_cols, length))
cat("Patients:", n_pat, "| Events:", sum(dati$PFS_I_EVENT), "\n")

nm <- names(set_cols)

set.seed(1)
foldid_orig <- sample(rep(1:10, length.out = n_pat))

#Coefficients and variable-selection results
beta_orig <- fit_all_specs(ORIG$x, Surv(ORIG$t, ORIG$e), foldid_orig, set_cols)
lp_orig        <- lapply(beta_orig, function(b) as.numeric(ORIG$x %*% b))
apparent_c     <- sapply(lp_orig, cindex,    time = ORIG$t, event = ORIG$e)
apparent_slope <- sapply(lp_orig, slope_fun, time = ORIG$t, event = ORIG$e)

# Compute apparent discrimination and calibration
for (s in nm) {
  lp <- as.numeric(ORIG$x %*% beta_orig[[s]])
  apparent_c[s]     <- cindex(lp, ORIG$t, ORIG$e)
  apparent_slope[s] <- slope_fun(lp, ORIG$t, ORIG$e)
}

cat("\n--- BOOTSTRAP VALIDATION ---\n")
mk <- function() matrix(NA_real_, B, length(nm), dimnames = list(NULL, nm))
app_mat <- mk(); test_mat <- mk(); slope_mat <- mk(); slope_app_mat <- mk()

# Store bootstrap coefficient
boot_betas <- lapply(setNames(nm, nm), function(s)
  matrix(0, length(master), B, dimnames = list(master, NULL)))

for (b in seq_len(B)) {
  cat("Bootstrap iteration", b, "of", B, "\n")
  
  # Resample incomplete data with replacement
  set.seed(10000 + b)
  idx <- sample(seq_len(n_pat), replace = TRUE)
  foldid_b <- sample(rep(1:10, length.out = n_pat))[idx]
  data_boot <- dati[idx, ]
  
  # Perform a single imputation within the bootstrap sample
  BOOT <- impute_and_build(data_boot, cols_no_pred,
                           seed = 1000 + b,
                           label = paste("BOOT", b), verbose_log = FALSE)
  
  # Fit penalized Cox models in the bootstrap sample
  x_b    <- align_x(BOOT$x, master)
  beta_b <- fit_all_specs(x_b, Surv(BOOT$t, BOOT$e), foldid_b, set_cols)
  
  
  # Evaluate bootstrap apparent performance and test performance
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

# Bootstrap estimate of optimism (C-index)
optimism_mat <- app_mat - test_mat

# Optimism-corrected calibration slope
slope_ott  <- colMeans(slope_app_mat - slope_mat, na.rm = TRUE)
slope_corr <- apparent_slope[nm] - slope_ott

# Summarize apparent and optimism-corrected model performance
final_summary <- data.frame(
  Set             = nm,
  C_Apparente     = apparent_c[nm],
  Slope_Apparente = apparent_slope[nm],
  Ottimismo_Medio = colMeans(optimism_mat, na.rm = TRUE),
  C_Validato      = apparent_c[nm] - colMeans(optimism_mat, na.rm = TRUE),
  Shrinkage_s     = colMeans(slope_mat, na.rm = TRUE),
  Slope_App_Boot  = colMeans(slope_app_mat, na.rm = TRUE),
  Slope_Ottimismo = slope_ott[nm],
  Slope_Validata  = slope_corr[nm],
  row.names = NULL
) %>% mutate(across(where(is.numeric), ~ round(., 3)))
print(final_summary, row.names = FALSE)

# Generate a coefficient summary table
make_coef_table <- function(set) {
  cols_set <- set_cols[[set]]
  b        <- beta_orig[[set]][cols_set]
  b_mat    <- boot_betas[[set]][cols_set, , drop = FALSE]
  
  data.frame(Variabile    = cols_set,
             Beta_Orig    = b,
             HR     = exp(b),
             Media_Boot   = rowMeans(b_mat),
             SD_Boot      = apply(b_mat, 1, sd),
             SD_Boot_Cond = apply(b_mat, 1, function(x) { nz <- x[x != 0]; if (length(nz) > 1) sd(nz) else NA_real_ }),
             Sel_Orig     = as.integer(b != 0),
             Boot_Sel_Pct = rowMeans(b_mat != 0) * 100,
             row.names    = NULL) %>%
    filter(Beta_Orig != 0 | Boot_Sel_Pct > 0) %>%
    arrange(desc(Boot_Sel_Pct), desc(abs(Beta_Orig))) %>%
    mutate(across(c(Beta_Orig, HR, Media_Boot, SD_Boot, SD_Boot_Cond), ~ round(., 3)),
           Boot_Sel_Pct = round(Boot_Sel_Pct, 1))
}


coef_tables <- setNames(lapply(nm, make_coef_table), nm)

for (s in nm) {
  cat("\n--- COEFFICIENTS & BOOTSTRAP SD:", s, "---\n")
  print(coef_tables[[s]], row.names = FALSE)
}

cat("\n--- PROPORTIONAL HAZARDS ASSUMPTION TEST ---\n")

# PH test on the linear predictor and on the selected covariates (unpenalized refit)
ph_check <- function(s, b, d) {
  
  sel <- names(b)[b != 0]
  
  # A. Linear predictor
  lp     <- as.numeric(d$x %*% b)
  fit_lp <- coxph(Surv(d$t, d$e) ~ lp)
  zph_lp <- cox.zph(fit_lp)
  
  # B. Selected covariates and global test
  tab  <- cox.zph(coxph(Surv(t, e) ~ ., data = data.frame(t = d$t, e = d$e, d$x[, sel, drop = FALSE])))$table
  cov  <- rownames(tab) != "GLOBAL"
  p    <- tab[cov, "p"]
  glob <- if ("GLOBAL" %in% rownames(tab)) tab["GLOBAL", "p"] else NA_real_
  
  list(
    lp     = data.frame(Set = s, Chisq = zph_lp$table["lp", "chisq"], P_value = zph_lp$table["lp", "p"]),
    vars   = data.frame(Set        = s,
                        Variabile  = c(sel, "GLOBAL"),
                        P_value    = c(p, glob),
                        P_value_BH = c(p.adjust(p, method = "BH"), NA)),
    zph_lp = zph_lp,
    fit_lp = fit_lp
  )
}

ph_res <- lapply(names(set_cols), function(s) ph_check(s, beta_orig[[s]], ORIG)) %>%
  setNames(names(set_cols))

round_df <- function(d) mutate(d, across(where(is.numeric), ~ round(., 4)))

cat("\n--- PH ASSUMPTION TEST ON THE LINEAR PREDICTOR ---\n")
print(round_df(bind_rows(lapply(ph_res, `[[`, "lp"))), row.names = FALSE)

cat("\n--- PH ASSUMPTION TEST BY COVARIATE AND GLOBAL ---\n")
print(round_df(bind_rows(lapply(ph_res, `[[`, "vars"))), row.names = FALSE)

# Scaled Schoenfeld residuals for the linear predictor
par(mfrow = c(1, length(ph_res)), mar = c(4.5, 4.5, 3, 1))
for (r in ph_res) {
  plot(r$zph_lp, main = sprintf("%s Model\n(p = %.3f)", r$lp$Set, r$lp$P_value),
       ylab = "Beta(t) for Linear Predictor", xlab = "Time (Months)")
  abline(h = 0, col = "gray60", lty = 3)
  abline(h = coef(r$fit_lp), col = "red", lty = 2, lwd = 1.5)
}
par(mfrow = c(1, 1))
# dev.off()

# Final consistency and diagnostic checks
cat("NA count in slope_mat (test*):\n");         print(colSums(is.na(slope_mat)))
cat("NA count in slope_app_mat (apparent*):\n"); print(colSums(is.na(slope_app_mat)))

#Correlation Matrix
COR_MODELS    <- "Integrated" 
COR_SEL_TRESH <- 0.10           # Minimum bootstrap selection frequency
COR_METHOD    <- "pearson"

# Clustered heatmap 
plot_cor <- function(cm, title) {
  n <- ncol(cm)
  corrplot::corrplot(cm, method = "color", type = "upper", order = "hclust",
                     hclust.method = "average", addCoef.col = "black",
                     number.cex = max(0.35, 0.9 - n * 0.015),
                     tl.cex     = max(0.5, 1.0 - n * 0.012),
                     tl.col = "black", tl.srt = 45, diag = FALSE,
                     title = title, mar = c(0, 0, 2, 0))
}

for (mod in intersect(COR_MODELS, nm)) {
  b_mat  <- t(boot_betas[[mod]][set_cols[[mod]], , drop = FALSE])
  active <- colnames(b_mat)[colMeans(b_mat != 0) >= COR_SEL_TRESH]
  active <- active[apply(b_mat[, active, drop = FALSE], 2, sd) > 0]
  
  if (length(active) < 2) {
    message("Model ", mod, ": fewer than 2 variables selected in >= ",
            COR_SEL_TRESH * 100, "% of the bootstraps. Correlation plots skipped.")
    next
  }
  
  x_mod <- ORIG$x[, active, drop = FALSE]
  x_mod <- x_mod[, apply(x_mod, 2, sd) > 0, drop = FALSE]
  
  size <- max(8, length(active) * 0.55)
  pdf(paste0("cor_", mod, ".pdf"), width = size, height = size)
  plot_cor(cor(b_mat[, active], method = COR_METHOD),
           paste0("Bootstrap coefficient correlations (", mod, ")"))
  if (ncol(x_mod) >= 2) {
    plot_cor(cor(x_mod, method = COR_METHOD),
             paste0("Predictor correlations (", mod, ")"))
  }
  dev.off()
  cat("Saved: cor_", mod, ".pdf\n", sep = "")
}
