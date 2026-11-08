library(tidyverse)
library(mice)
library(glmnet)
library(survival)

set.seed(42)

B                   <- 200     # Number of bootstrap replicates
CNA_AS_PREDICTORS <- TRUE      # FALSE = CNAs (except del_17p) excluded as predictors in MICE;
                               #TRUE  = only predefined CNA variables are excluded

dataset <- read_rds("Clinical_Genomic_Dataset_Pre_Imputation.rds")

vars_genomic <- grep("^(del_|amp_)", colnames(dataset), value = TRUE)

vars_clinical <- setdiff(colnames(dataset),
                         c("PFS_I_MONTHS", "PFS_I_EVENT", "MPC", vars_genomic))

cat("Clinical variables (", length(vars_clinical), "):",
    paste(vars_clinical, collapse = ", "), "\n")
cat("Genomic variables (", length(vars_genomic), "):",
    paste(vars_genomic, collapse = ", "), "\n")

all_vars <- c(vars_clinical, vars_genomic)

# Define CNA exclusion rules for the MICE imputation model (retaining del_17p)
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
  if (length(comm) > 0) out[, comm] <- x_raw[, comm, drop = FALSE]
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
  if (any(!is.finite(lp)) || sd(lp) == 0) return(NA_real_)
  d <- data.frame(t = time, e = event, lp = lp)
  unname(coef(coxph(Surv(t, e) ~ lp, data = d)))
}

# Fit a cross-validated Cox-LASSO model
fit_lasso_cox <- function(x, t, e, foldid, alpha = 1) {
  cv <- cv.glmnet(x, Surv(t, e), family = "cox",
                  alpha = alpha, foldid = foldid, cox.ties = "breslow")
  
  setNames(as.numeric(coef(cv, s = "lambda.min")), colnames(x))
}

# Evaluate apparent* and test* performance
eval_model <- function(be, x_app, t_app, e_app, x_test, t_test, e_test) {
  lp_a <- as.numeric(x_app %*% be)
  lp_t <- as.numeric(x_test %*% be)
  c(app     = cindex(lp_a, t_app, e_app),
    test    = cindex(lp_t, t_test, e_test),
    sl_app  = slope_fun(lp_a, t_app, e_app),
    sl_test = slope_fun(lp_t, t_test, e_test))
}

# ---- MICE ---------------------------------------------------------------------

# Construct the predictor matrix used for MICE imputation.
build_pred <- function(dat, cols_no_pred = character(0)) {
  pm <- make.predictorMatrix(dat)
  pm[, "PFS_I_MONTHS"] <- 0
  pm[, c("chaz", "PFS_I_EVENT")] <- 1
  if (length(cols_no_pred) > 0) pm[, intersect(cols_no_pred, colnames(pm))] <- 0
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

# ---- MEAN SUBSTITUTION --------------------------------------------------------

# Estimate mean/mode substitution parameters
make_impute_params <- function(df, vars) {
  lapply(df[vars], function(x) {
    if (is.numeric(x)) mean(x, na.rm = TRUE) else names(which.max(table(x)))
  })
}

# Apply mean/mode imputation parameters to missing values.
apply_impute <- function(df, params) {
  for (v in names(params)) {
    p <- params[[v]]
    if (is.logical(df[[v]])) p <- as.logical(p)
    df[[v]][is.na(df[[v]])] <- p
  }
  df
}

# Apply mean/mode imputation and build the design matrix.
mean_impute_and_build <- function(dat, params) {
  d <- apply_impute(dat, params)
  list(x = make_x(d), t = d$PFS_I_MONTHS, e = d$PFS_I_EVENT)
}

# Summarize performance metrics.
perf_row <- function(label, n, ev, p, c_app, sl_app, r) {
  opt   <- mean(r[, "app"] - r[, "test"], na.rm = TRUE)
  sl_ot <- mean(r[, "sl_app"] - r[, "sl_test"], na.rm = TRUE)
  data.frame(
    Modello            = label,
    N                  = n,
    Eventi             = ev,
    Eventi_per_colonna = ev / p,
    C_Apparente        = c_app,
    Slope_Apparente    = sl_app,
    Ottimismo_Medio    = opt,
    C_Validato         = c_app - opt,
    Shrinkage_s        = mean(r[, "sl_test"], na.rm = TRUE),   # media test* slopes
    Slope_App_Boot     = mean(r[, "sl_app"], na.rm = TRUE),
    Slope_Ottimismo    = sl_ot,
    Slope_Validata     = sl_app - sl_ot,                        # optimism corrected slope
    Boot_Validi_C      = sum(!is.na(r[, "app"] - r[, "test"])),
    Boot_Validi_Slope  = sum(!is.na(r[, "sl_app"] - r[, "sl_test"])),
    row.names = NULL
  )
}

# Generate a coefficient summary table
coef_table <- function(cols, b, bm) {
  b     <- b[cols]
  bm    <- bm[cols, , drop = FALSE]
  
  sd_cond <- apply(bm, 1, function(x) {
    nz <- x[!is.na(x) & x != 0]
    if (length(nz) > 1) sd(nz) else NA_real_
  })
  
  data.frame(Variabile    = cols,
             Beta_Orig    = b,
             HR     = exp(b),
             Media_Boot   = apply(bm, 1, mean, na.rm = TRUE),
             SD_Boot      = apply(bm, 1, sd, na.rm = TRUE),
             SD_Boot_Cond = sd_cond,
             Sel_Orig     = as.numeric(b != 0),
             Boot_Sel_Pct = rowMeans(bm != 0, na.rm = TRUE) * 100,
             row.names    = NULL) %>%
    filter(Beta_Orig != 0 | Boot_Sel_Pct > 0) %>%
    arrange(desc(Boot_Sel_Pct), desc(abs(Beta_Orig))) %>%
    mutate(across(c(Beta_Orig, HR, Media_Boot, SD_Boot, SD_Boot_Cond), ~ round(., 3)),
           Boot_Sel_Pct = round(Boot_Sel_Pct, 1))
}

# Replace non-positive survival times with a small positive value
dataset <- dataset %>% mutate(PFS_I_MONTHS = ifelse(PFS_I_MONTHS <= 0, 0.1, PFS_I_MONTHS))

dati  <- dataset %>% select(-MPC)
n_pat <- nrow(dati)
ev_all <- sum(dati$PFS_I_EVENT)

# --- MICE (single imputation)
ORIG   <- impute_and_build(dati, cols_no_pred, seed = 123)
master        <- colnames(ORIG$x)
p_cols <- length(master)

# --- Mean substitution
vars_imp    <- intersect(all_vars, colnames(dati))
par_orig    <- make_impute_params(dati, vars_imp)
ORIG_MEAN   <- mean_impute_and_build(dati, par_orig)
ORIG_MEAN$x <- align_x(ORIG_MEAN$x, master)

# --- Complete-case cohort: patients with all variables and outcome observed
S      <- complete.cases(dati[, c("PFS_I_MONTHS", "PFS_I_EVENT", vars_imp)])
df_cc  <- dati[S, ]
n_cc   <- nrow(df_cc)
ev_cc  <- sum(df_cc$PFS_I_EVENT)
x_cc   <- align_x(make_x(df_cc), master)
t_cc   <- df_cc$PFS_I_MONTHS
e_cc   <- df_cc$PFS_I_EVENT

cat("\nDesign matrix columns:", p_cols, "\n")
cat(sprintf("All patients     : n = %d | events = %d (%.1f%%) | events per column = %.1f\n",
            n_pat, ev_all, 100 * ev_all / n_pat, ev_all / p_cols))
cat(sprintf("Complete-case    : n = %d (%.1f%%) | events = %d (%.1f%%) | events per column = %.1f\n",
            n_cc, 100 * n_cc / n_pat, ev_cc, 100 * ev_cc / n_cc, ev_cc / p_cols))

set.seed(1)
foldid_orig    <- sample(rep(1:10, length.out = n_pat))
set.seed(2)
foldid_cc_orig <- sample(rep(1:10, length.out = n_cc))

beta_mice <- fit_lasso_cox(ORIG$x,      ORIG$t,      ORIG$e,      foldid_orig)
beta_mean <- fit_lasso_cox(ORIG_MEAN$x, ORIG_MEAN$t, ORIG_MEAN$e, foldid_orig)
beta_cc   <- fit_lasso_cox(x_cc,        t_cc,        e_cc,        foldid_cc_orig)

# Compute apparent performance
app_orig <- function(be, x, t, e) {
  lp <- as.numeric(x %*% be)
  c(c = cindex(lp, t, e), slope = slope_fun(lp, t, e))
}

ap_mice <- app_orig(beta_mice, ORIG$x,      ORIG$t,      ORIG$e)
ap_mean <- app_orig(beta_mean, ORIG_MEAN$x, ORIG_MEAN$t, ORIG_MEAN$e)
ap_cc   <- app_orig(beta_cc,   x_cc,        t_cc,        e_cc)

cat("\n--- BOOTSTRAP VALIDATION ---\n")

meths <- c("MICE", "Mean", "CC")

# Bootstrap performance metrics
res <- lapply(setNames(meths, meths), function(m)
  matrix(NA_real_, B, 4, dimnames = list(NULL, c("app", "test", "sl_app", "sl_test"))))

# Bootstrap coefficient estimates
boot_betas <- lapply(setNames(meths, meths), function(m)
  matrix(NA_real_, length(master), B, dimnames = list(master, NULL)))

for (b in seq_len(B)) {
  if (b %% 25 == 0 || b == 1 || b == B) cat("Bootstrap iteration", b, "of", B, "\n")
  
  set.seed(10000 + b)
  idx      <- sample(seq_len(n_pat), replace = TRUE)
  foldid_b <- sample(rep(1:10, length.out = n_pat))[idx]
  set.seed(20000 + b)
  idx_cc   <- sample(seq_len(n_cc), replace = TRUE)
  foldid_cc_b <- sample(rep(1:10, length.out = n_cc))[idx_cc]
  
  data_boot <- dati[idx, ]
  
  # --- MICE ---------------------------------------------------------------
  BOOT <- impute_and_build(data_boot, cols_no_pred, seed = 1000 + b)
  x_boot <- align_x(BOOT$x, master)
  be <- fit_lasso_cox(x_boot, BOOT$t, BOOT$e, foldid_b)
  boot_betas$MICE[, b] <- be
  res$MICE[b, ]        <- eval_model(be, x_boot, BOOT$t, BOOT$e, ORIG$x, ORIG$t, ORIG$e)
  
  
  # --- Mean substitution --------------------------------------------------
  BOOT_M   <- mean_impute_and_build(data_boot, make_impute_params(data_boot, vars_imp))
  x_boot_m <- align_x(BOOT_M$x, master)
  be <- fit_lasso_cox(x_boot_m, BOOT_M$t, BOOT_M$e, foldid_b)
  boot_betas$Mean[, b] <- be
  res$Mean[b, ] <- eval_model(be, x_boot_m, BOOT_M$t, BOOT_M$e, ORIG_MEAN$x, ORIG_MEAN$t, ORIG_MEAN$e)
  
  
  # --- Complete-case ------------------------------------------------------
  dat_cc_b <- df_cc[idx_cc, ]
  x_cc_b   <- tryCatch(align_x(make_x(dat_cc_b), master), error = function(err) NULL)
  if (!is.null(x_cc_b)) {
    t_b <- dat_cc_b$PFS_I_MONTHS
    e_b <- dat_cc_b$PFS_I_EVENT
    be  <- fit_lasso_cox(x_cc_b, t_b, e_b, foldid_cc_b)
    boot_betas$CC[, b] <- be
    res$CC[b, ]        <- eval_model(be, x_cc_b, t_b, e_b, x_cc, t_cc, e_cc)
  }
}

# Optimism-corrected performance of the three strategies
final_summary <- bind_rows(
  perf_row("MICE (all patients)",              n_pat, ev_all, p_cols, ap_mice["c"], ap_mice["slope"], res$MICE),
  perf_row("Mean substitution (all patients)", n_pat, ev_all, p_cols, ap_mean["c"], ap_mean["slope"], res$Mean),
  perf_row("Complete-case",                    n_cc,  ev_cc,  p_cols, ap_cc["c"],   ap_cc["slope"],   res$CC)
) %>% mutate(across(where(is.numeric) & !c(N, Eventi, Boot_Validi_C, Boot_Validi_Slope),
                    ~ round(., 3)))

cat("\n--- INTEGRATED MODEL: MICE vs MEAN SUBSTITUTION vs COMPLETE-CASE ---\n")
cat("Shrinkage_s = mean test* slope; Slope_Validata = apparent slope - mean(apparent* - test* slope).\n")
print(final_summary, row.names = FALSE)

# Coefficient tables with bootstrap SD
sh_mice <- mean(res$MICE[, "sl_test"], na.rm = TRUE)
sh_mean <- mean(res$Mean[, "sl_test"], na.rm = TRUE)
sh_cc   <- mean(res$CC[,   "sl_test"], na.rm = TRUE)

tab_mice <- coef_table(master, beta_mice, boot_betas$MICE)
tab_mean <- coef_table(master, beta_mean, boot_betas$Mean)
tab_cc   <- coef_table(master, beta_cc,   boot_betas$CC)

cat("\n--- COEFFICIENTS (MICE) ---\n");              print(tab_mice, row.names = FALSE)
cat("\n--- COEFFICIENTS (MEAN SUBSTITUTION) ---\n"); print(tab_mean, row.names = FALSE)
cat("\n--- COEFFICIENTS (COMPLETE-CASE) ---\n");     print(tab_cc,   row.names = FALSE)

cat("\n--- PROPORTIONAL HAZARDS ASSUMPTION TEST (Integrated) ---\n")

ph_check <- function(label, beta, x, t, e) {
  sel <- names(beta)[beta != 0]
  
  # PH test on the linear predictor
  lp     <- as.numeric(x %*% beta)
  fit_lp <- coxph(Surv(t, e) ~ lp)
  zph_lp <- cox.zph(fit_lp)
  
  # PH test on selected covariates (unpenalized refit)
  df_cox <- data.frame(t = t, e = e, x[, sel, drop = FALSE])
  tab    <- cox.zph(coxph(Surv(t, e) ~ ., data = df_cox))$table
  cov    <- rownames(tab) != "GLOBAL"
  
  list(
    lp     = data.frame(Modello = label,
                        Chisq   = zph_lp$table["lp", "chisq"],
                        P_value = zph_lp$table["lp", "p"]),
    vars   = data.frame(Modello    = label,
                        Variabile  = c(sel, "GLOBAL"),
                        P_value    = c(tab[cov, "p"], tab["GLOBAL", "p"]),
                        P_value_BH = c(p.adjust(tab[cov, "p"], method = "BH"), NA)),
    zph_lp = zph_lp,
    fit_lp = fit_lp
  )
}

ph_res <- list(
  MICE = ph_check("MICE",              beta_mice, ORIG$x,      ORIG$t,      ORIG$e),
  Mean = ph_check("Mean substitution", beta_mean, ORIG_MEAN$x, ORIG_MEAN$t, ORIG_MEAN$e),
  CC   = ph_check("Complete-case",     beta_cc,   x_cc,        t_cc,        e_cc)
)

round_df <- function(d) mutate(d, across(where(is.numeric), ~ round(., 4)))

cat("\n--- PH test on the linear predictor ---\n")
print(round_df(bind_rows(lapply(ph_res, `[[`, "lp"))), row.names = FALSE)

for (r in ph_res) {
  cat("\n--- PH test by covariate (", r$vars$Modello[1], ") ---\n", sep = "")
  print(round_df(r$vars), row.names = FALSE)
}

# Schoenfeld residual plots for the linear predictor
par(mfrow = c(1, length(ph_res)), mar = c(4.5, 4.5, 3, 1))
for (r in ph_res) {
  plot(r$zph_lp, main = sprintf("%s\n(p = %.3f)", r$lp$Modello, r$lp$P_value),
       ylab = "Beta(t) for Linear Predictor", xlab = "Time (Months)")
  abline(h = 0, col = "gray60", lty = 3)
  abline(h = coef(r$fit_lp), col = "red", lty = 2, lwd = 1.5)
}
par(mfrow = c(1, 1))

cat("\n--- DIAGNOSTIC CHECKS ---\n")

diag_tab <- data.frame(
  Method        = meths,
  NA_slope_app  = sapply(res, function(r) sum(is.na(r[, "sl_app"]))),
  NA_slope_test = sapply(res, function(r) sum(is.na(r[, "sl_test"]))),
  row.names     = NULL
)
print(diag_tab, row.names = FALSE)

cat("Selected variables in the original fit (MICE / Mean / CC):",
    sum(beta_mice != 0), "/", sum(beta_mean != 0), "/", sum(beta_cc != 0), "\n")