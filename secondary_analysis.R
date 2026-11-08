library(tidyverse)
library(mice)
library(glmnet)
library(survival)

set.seed(42)        
B                   <- 200     # Number of bootstrap resamples
CNA_AS_PREDICTORS <- FALSE     # FALSE = CNAs (except del_17p) excluded as predictors in MICE; 
                               #TRUE  = only predefined CNA variables are excluded
dataset <- read_rds("Imaging_Dataset_Pre_Imputation.rds")

vars_genomic <- grep("^(del_|amp_)", colnames(dataset), value = TRUE)  

vars_imaging_candidate <- c("PM_SUVmax", "BM_Suv_max", "EM_SUV_max", "FL_SUV_max",
                            "gDNA_tumor_fraction", "cfDNA_tumor_fraction")
vars_imaging <- intersect(vars_imaging_candidate, colnames(dataset))

vars_clinical <- setdiff(colnames(dataset),
                         c("PFS_I_MONTHS", "PFS_I_EVENT", "MPC",
                           vars_genomic, vars_imaging))

dataset <- dataset %>%
  select(any_of(c("PFS_I_MONTHS", "PFS_I_EVENT", "MPC",
                  vars_clinical, vars_genomic, vars_imaging)))

cat("Clinical variables (", length(vars_clinical), "):",
    paste(vars_clinical, collapse = ", "), "\n")
cat("Genomic variables (", length(vars_genomic), "):",
    paste(vars_genomic, collapse = ", "), "\n")
cat("Imaging/Liquid Biopsy variables (", length(vars_imaging), "):",
    ifelse(length(vars_imaging) > 0, paste(vars_imaging, collapse = ", "), "None"), "\n")


feature_sets <- list(
  Clinical        = vars_clinical,
  Clin_Gen        = c(vars_clinical, vars_genomic),
  Imaging_Only    = vars_imaging,
  Clin_Imag       = c(vars_clinical, vars_imaging),
  Full_Integrated = c(vars_clinical, vars_genomic, vars_imaging)
)


all_vars <- unique(unlist(feature_sets))

# Define CNA exclusion rules for the MICE imputation model
# del_17p is always retained as a predictor.
cna_cols     <- setdiff(vars_genomic, "del_17p")
cols_no_pred <- if (CNA_AS_PREDICTORS) c("amp_7p", "amp_5p", "amp_18p", "amp_3p") else cna_cols

# Construct the design matrix used for penalized Cox regression
make_x <- function(df, target_vars = all_vars) {
  vars_ok <- intersect(target_vars, colnames(df))
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
  if (sd(lp) == 0) return(0.5)
  d <- data.frame(t = time, e = event, inv = -lp)
  concordance(Surv(t, e) ~ inv, data = d)$concordance
}

# Estimate the calibration slope
slope_fun <- function(lp, time, event) {
  if (sd(lp) == 0) return(NA_real_)
  d <- data.frame(t = time, e = event, lp = lp)
  unname(coef(coxph(Surv(t, e) ~ lp, data = d)))
}

# Fit cross-validated penalized Cox models
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
    cols_excl <- intersect(cols_no_pred, colnames(pm))
    pm[, cols_excl] <- 0                  
  }
  
  # Prevent high-missingness imaging variables (~67% NA) from imputing baseline clinical variables
  if (length(vars_imaging) > 0) {
    pm[intersect(vars_clinical, rownames(pm)), intersect(vars_imaging, colnames(pm))] <- 0
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

# Replace non-positive survival times with a small positive value
dataset <- dataset %>% mutate(PFS_I_MONTHS = ifelse(PFS_I_MONTHS <= 0, 0.1, PFS_I_MONTHS))

dati  <- dataset %>% select(-any_of("MPC"))
n_pat <- nrow(dati)

# Generate original imputed dataset
ORIG <- impute_and_build(dati, cols_no_pred, seed = 123)

# Define reference master feature space from the first imputed original dataset
master        <- colnames(ORIG$x)
var_of_master <- attr(ORIG$x, "var")
set_cols      <- lapply(feature_sets, function(vs) master[var_of_master %in% vs])

cat("\nExpanded design matrix columns per model:\n")
print(sapply(set_cols, length))
cat("Total Patients:", n_pat, "| Total Events:", sum(dati$PFS_I_EVENT), "\n")

nm <- names(set_cols)

set.seed(1)
foldid_orig <- sample(rep(1:10, length.out = n_pat))

# Fit the penalized Cox models on the original imputed dataset
beta_orig <- fit_all_specs(ORIG$x, Surv(ORIG$t, ORIG$e), foldid_orig, set_cols)

# Estimate linear predictor, apparent discrimination and calibration.
lp_orig        <- lapply(beta_orig, function(b) as.numeric(ORIG$x %*% b))
apparent_c     <- sapply(lp_orig, cindex,    time = ORIG$t, event = ORIG$e)
apparent_slope <- sapply(lp_orig, slope_fun, time = ORIG$t, event = ORIG$e)

cat("\nApparent C-index:\n");           print(round(apparent_c, 3))
cat("Apparent calibration slope:\n"); print(round(apparent_slope, 3))

cat("\n--- BOOTSTRAP VALIDATION (IMPUTED PIPELINE) ---\n")
mk <- function() matrix(NA_real_, B, length(nm), dimnames = list(NULL, nm))
app_mat <- mk(); test_mat <- mk()
slope_mat     <- mk()   # Test* calibration slope (bootstrap model evaluated on the original dataset)
slope_app_mat <- mk()   # Apparent* calibration slope (bootstrap model evaluated on its own bootstrap sample)

# Store bootstrap coefficient
boot_betas <- lapply(setNames(nm, nm), function(s)
  matrix(NA_real_, length(master), B, dimnames = list(master, NULL)))

for (b in seq_len(B)) {
  if (b %% 25 == 0 || b == 1 || b == B) {
    cat("Bootstrap iteration", b, "of", B, "\n")
  }
  
  # Resample incomplete data with replacement
  set.seed(10000 + b)
  idx <- sample(seq_len(n_pat), replace = TRUE)
  foldid_b <- sample(rep(1:10, length.out = n_pat))[idx]
  
  data_boot <- dati[idx, ]
  
  # Perform a single imputation within the bootstrap sample
  BOOT   <- impute_and_build(data_boot, cols_no_pred, seed = 1000 + b)
  
  
  # Fit models on bootstrap sample
  x_b    <- align_x(BOOT$x, master)
  beta_b <- fit_all_specs(x_b, Surv(BOOT$t, BOOT$e), foldid_b, set_cols)
  
  # Evaluate Apparent* (on BOOT) and Test* (on ORIG) performance (C-index and slope)
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
  Colonne_p       = sapply(set_cols, length)[nm],
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

cat("\n--- PERFORMANCE SUMMARY (IMPUTED MODELS) ---\n")
print(final_summary, row.names = FALSE)

make_coef_table <- function(set) {
  cols_set <- set_cols[[set]]
  sh       <- mean(slope_mat[, set], na.rm = TRUE)
  b        <- beta_orig[[set]][cols_set]
  b_mat    <- boot_betas[[set]][cols_set, , drop = FALSE]
  
  data.frame(Variabile    = cols_set,
             Beta_Orig    = b,
             HR     = exp(b),
             Media_Boot   = rowMeans(b_mat, na.rm = TRUE),
             SD_Boot      = apply(b_mat, 1, sd, na.rm = TRUE),
             SD_Boot_Cond = apply(b_mat, 1, function(x) {
               nz <- x[!is.na(x) & x != 0]
               if (length(nz) > 1) sd(nz) else NA_real_
             }),
             Sel_Orig     = as.integer(b != 0),
             Boot_Sel_Pct = rowMeans(b_mat != 0, na.rm = TRUE) * 100,
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

cat("\n==============================================================================\n")
cat(" MICE vs COMPLETE-CASE BOOTSTRAP VALIDATION:\n")
cat("==============================================================================\n")

# Complete-case bootstrap validation
run_cc_validation <- function(df_raw, target_vars, master_cols, B_iter = 200, n_folds = 5, label = "CC") {
  
  # Select complete cases
  df_cc <- df_raw %>% drop_na(all_of(c("PFS_I_MONTHS", "PFS_I_EVENT", target_vars)))
  n_cc  <- nrow(df_cc)
  ev_cc <- sum(df_cc$PFS_I_EVENT)
  
  cat(sprintf("[%s] Complete-Case Cohort: n = %d patients | Events = %d (%.1f%%)\n",
              label, n_cc, ev_cc, 100 * ev_cc / n_cc))
  
  x_cc_orig <- align_x(make_x(df_cc, target_vars), master)[, master_cols, drop = FALSE]
  y_cc_orig <- Surv(df_cc$PFS_I_MONTHS, df_cc$PFS_I_EVENT)
  
  # Fit the complete-case model
  set.seed(2)
  foldid_cc_orig <- sample(rep(seq_len(n_folds), length.out = n_cc))
  cv_cc_orig     <- cv.glmnet(x_cc_orig, y_cc_orig, family = "cox",
                              alpha = 1, foldid = foldid_cc_orig, cox.ties = "breslow", lambda.min.ratio = 0.01)
  
  
  # Store bootstrap coefficients and selection frequencies
  beta_cc_orig        <- as.numeric(coef(cv_cc_orig, s = "lambda.min"))
  names(beta_cc_orig) <- master_cols
  
  # Estimate apparent discrimination and calibration.
  lp_cc_orig   <- as.numeric(x_cc_orig %*% beta_cc_orig)
  c_app_cc     <- cindex(lp_cc_orig, df_cc$PFS_I_MONTHS, df_cc$PFS_I_EVENT)
  slope_app_cc <- slope_fun(lp_cc_orig, df_cc$PFS_I_MONTHS, df_cc$PFS_I_EVENT)
  
  app_cc            <- rep(NA_real_, B_iter)
  test_cc           <- rep(NA_real_, B_iter)
  slope_cc          <- rep(NA_real_, B_iter)   # slope  test*
  slope_app_boot_cc <- rep(NA_real_, B_iter)   # slope apparent*
  
  # Store bootstrap selection frequencies and coefficient estimates.
  boot_betas_cc <- matrix(NA_real_, length(master_cols), B_iter, dimnames = list(master_cols, NULL))
  
  for (b in seq_len(B_iter)) {
    # Generate a bootstrap sample from the complete-case cohort.
    set.seed(20000 + b)
    idx_cc       <- sample(seq_len(n_cc), replace = TRUE)
    foldid_b_cc  <- sample(rep(seq_len(n_folds), length.out = n_cc))[idx_cc]
    dati_boot_cc <- df_cc[idx_cc, ]
    
    x_boot_cc <- align_x(make_x(dati_boot_cc, target_vars), master)[, master_cols, drop = FALSE]
    y_boot_cc <- Surv(dati_boot_cc$PFS_I_MONTHS, dati_boot_cc$PFS_I_EVENT)
    
    # Refit the penalized Cox model within the bootstrap sample.
    cv_b <- tryCatch(
      cv.glmnet(x_boot_cc, y_boot_cc, family = "cox", alpha = 1,
                foldid = foldid_b_cc, cox.ties = "breslow", lambda.min.ratio = 0.01),
      error = function(e) NULL
    )
    if (is.null(cv_b)) next
    
    # Extract and store coefficient
    beta_b <- setNames(as.numeric(coef(cv_b, s = "lambda.min")), master_cols)
    boot_betas_cc[, b] <- beta_b
    
    # Calculate linear predictors in the bootstrap and original cohorts
    lp_app  <- as.numeric(x_boot_cc %*% beta_b)
    lp_test <- as.numeric(x_cc_orig %*% beta_b)
    
    # Evaluate Apparent* (bootstrap sample) and Test* (original cohort)
    # discrimination and calibration.
    app_cc[b]            <- cindex(lp_app, dati_boot_cc$PFS_I_MONTHS, dati_boot_cc$PFS_I_EVENT)
    test_cc[b]           <- cindex(lp_test, df_cc$PFS_I_MONTHS, df_cc$PFS_I_EVENT)
    slope_app_boot_cc[b] <- slope_fun(lp_app, dati_boot_cc$PFS_I_MONTHS, dati_boot_cc$PFS_I_EVENT)
    slope_cc[b]          <- slope_fun(lp_test, df_cc$PFS_I_MONTHS, df_cc$PFS_I_EVENT)
  }
  
  # Calculate optimism-corrected performance and shrinkage
  opt_cc   <- mean(app_cc - test_cc, na.rm = TRUE)
  c_val_cc <- c_app_cc - opt_cc
  sh_cc    <- mean(slope_cc, na.rm = TRUE)
  sh_safe  <- ifelse(is.na(sh_cc), 1, sh_cc)
  
  # Bootstrap-based optimism correction for calibration
  slope_app_boot_mean <- mean(slope_app_boot_cc, na.rm = TRUE)
  slope_ott_cc        <- mean(slope_app_boot_cc - slope_cc, na.rm = TRUE)
  slope_corr_cc       <- slope_app_cc - slope_ott_cc
  
  # Summarize bootstrap mean and sd
  mean_boot_cc <- apply(boot_betas_cc, 1, mean, na.rm = TRUE)
  sd_boot_cc   <- apply(boot_betas_cc, 1, sd, na.rm = TRUE)
  
  # SD among bootstrap estimates with non-zero coefficients
  sd_cond_cc   <- apply(boot_betas_cc, 1, function(x) {
    nz <- x[x != 0 & !is.na(x)]
    if (length(nz) > 1) sd(nz) else NA_real_
  })
  
  # Generate the coefficient summary table
  coef_tab <- data.frame(
    Variabile    = master_cols,
    Beta_Orig    = beta_cc_orig,
    HR    = exp(beta_cc_orig),
    Media_Boot   = mean_boot_cc,
    SD_Boot      = sd_boot_cc,
    SD_Boot_Cond = sd_cond_cc,
    Sel_Orig     = as.numeric(beta_cc_orig != 0),
    Boot_Sel_Pct = rowMeans(boot_betas_cc != 0, na.rm = TRUE) * 100,
    row.names    = NULL
  ) %>%
    filter(Beta_Orig != 0 | Boot_Sel_Pct > 0) %>%
    arrange(desc(Boot_Sel_Pct), desc(abs(Beta_Orig))) %>%
    mutate(across(c(Beta_Orig, HR, Media_Boot, SD_Boot, SD_Boot_Cond), ~ round(., 3)),
           Boot_Sel_Pct = round(Boot_Sel_Pct, 1))
  
  # Summarize complete-case model performance
  summary_row <- data.frame(
    Set                = label,
    Strategia          = "Complete-case",
    N                  = n_cc,
    Eventi             = ev_cc,
    Colonne_p          = length(master_cols),
    Eventi_per_colonna = ev_cc / length(master_cols),
    C_Apparente        = c_app_cc,
    Slope_Apparente    = slope_app_cc,
    Ottimismo_Medio    = opt_cc,
    C_Validato         = c_val_cc,
    Shrinkage_s        = sh_cc,
    Slope_App_Boot     = slope_app_boot_mean,
    Slope_Ottimismo    = slope_ott_cc,
    Slope_Validata     = slope_corr_cc,
    row.names          = NULL
  )
  
  list(summary = summary_row, coef_table = coef_tab)
}

# Complete-case validation of the Clinical-Imaging model.
res_cc_clin_imag <- run_cc_validation(
  df_raw      = dati,
  target_vars = c(vars_clinical, vars_imaging),
  master_cols = set_cols[["Clin_Imag"]],
  B_iter      = B,
  n_folds     = 5,
  label       = "Clin_Imag"
)

# Complete-case validation of the Imaging-Only model.
res_cc_imag_only <- run_cc_validation(
  df_raw      = dati,
  target_vars = vars_imaging,
  master_cols = set_cols[["Imaging_Only"]],
  B_iter      = B,
  n_folds     = 5,
  label       = "Imaging_Only"
)

# Retrieve MICE results for comparison
imp_row <- function(set) {
  p_set <- length(set_cols[[set]])
  ev    <- sum(dati$PFS_I_EVENT)
  data.frame(
    Set                = set,
    Strategia          = "MICE",
    N                  = n_pat,
    Eventi             = ev,
    Colonne_p          = p_set,
    Eventi_per_colonna = ev / p_set,
    C_Apparente        = unname(apparent_c[set]),
    Slope_Apparente    = unname(apparent_slope[set]),
    Ottimismo_Medio    = mean(optimism_mat[, set], na.rm = TRUE),
    C_Validato         = unname(apparent_c[set]) - mean(optimism_mat[, set], na.rm = TRUE),
    Shrinkage_s        = mean(slope_mat[, set], na.rm = TRUE),
    Slope_App_Boot     = mean(slope_app_mat[, set], na.rm = TRUE),
    Slope_Ottimismo    = unname(slope_ott[set]),
    Slope_Validata     = unname(slope_corr[set]),
    row.names          = NULL
  )
}

comp_summary <- bind_rows(
  imp_row("Clin_Imag"),    res_cc_clin_imag$summary,
  imp_row("Imaging_Only"), res_cc_imag_only$summary
) %>% mutate(across(where(is.numeric) & !c(N, Eventi, Colonne_p), ~ round(., 3)))

cat("\n---  MICE vs COMPLETE-CASE: ---\n")
print(comp_summary, row.names = FALSE)

cat("\n--- COEFFICIENTS & BOOTSTRAP SD: Clin_Imag COMPLETE-CASE ---\n")
print(res_cc_clin_imag$coef_table, row.names = FALSE)

cat("\n--- COEFFICIENTS & BOOTSTRAP SD: Imaging_Only COMPLETE-CASE ---\n")
print(res_cc_imag_only$coef_table, row.names = FALSE)

cat("\n--- PROPORTIONAL HAZARDS ASSUMPTION TEST ---\n")
d <- ORIG
ph_lp <- ph_vars <- fit_lp_list <- zph_lp_list <- list()

for (s in nm) {
  b        <- beta_orig[[s]]
  sel_vars <- names(b[b != 0])
  if (length(sel_vars) == 0) next
  
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
# pdf("immagini/schoenfeld_lp.pdf", width = 12, height = 4)
par(mfrow = c(1, length(zph_lp_list)), mar = c(4.5, 4.5, 3, 1))
for (s in names(zph_lp_list)) {
  plot(zph_lp_list[[s]],
       main = paste0(s, " Model\n(p = ", round(ph_lp[[s]]$P_value, 3), ")"),
       ylab = "Beta(t) for Linear Predictor", xlab = "Time (Months)")
  abline(h = 0, col = "gray60", lty = 3)
  abline(h = coef(fit_lp_list[[s]]), col = "red", lty = 2, lwd = 1.5)
}

cat("\n--- 9 DIAGNOSTIC CHECKS ---\n")
cat("NA count in slope_mat (test*):\n");         print(colSums(is.na(slope_mat)))
cat("NA count in slope_app_mat (apparent*):\n"); print(colSums(is.na(slope_app_mat)))