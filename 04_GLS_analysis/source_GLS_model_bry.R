## =========================================================
## 2. Helpers
## =========================================================

safe_scale <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  sdx <- sd(x, na.rm = TRUE)
  if (is.na(sdx) || sdx == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / sdx)
}

make_rhs <- function(terms) {
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
}

build_formula_string <- function(resp, t_col,
                                 env_terms = character(0),
                                 pulse_terms = character(0),
                                 predictor_term = NULL) {
  rhs_terms <- c(
    t_col,
    paste0("I(", t_col, "^2)"),
    env_terms,
    pulse_terms,
    predictor_term
  )
  paste0(resp, " ~ ", make_rhs(rhs_terms))
}

circular_shift <- function(x, k) {
  n <- length(x)
  if (n <= 1) return(x)
  k <- k %% n
  if (k == 0) return(x)
  c(x[(k + 1):n], x[1:k])
}

get_valid_shifts <- function(n, min_shift = 10) {
  if (n <= 2) return(integer(0))

  if (n > (2 * min_shift)) {
    return(seq.int(min_shift, n - min_shift))
  } else {
    return(seq_len(n - 1))
  }
}

fit_model_safe <- function(formula_str, data, t_col) {
  fm <- as.formula(formula_str)

  fit1 <- tryCatch(
    gls(
      fm,
      data = data,
      correlation = corAR1(form = as.formula(paste0("~", t_col))),
      method = "ML",
      na.action = na.omit
    ),
    error = function(e) NULL
  )
  if (!is.null(fit1)) return(list(model = fit1, fit_method = "gls_ar1"))

  fit2 <- tryCatch(
    gls(
      fm,
      data = data,
      method = "ML",
      na.action = na.omit
    ),
    error = function(e) NULL
  )
  if (!is.null(fit2)) return(list(model = fit2, fit_method = "gls_no_ar1"))

  fit3 <- tryCatch(
    lm(fm, data = data),
    error = function(e) NULL
  )
  if (!is.null(fit3)) return(list(model = fit3, fit_method = "lm"))

  list(model = NULL, fit_method = NA_character_)
}

fit_model_fixed_method <- function(formula_str, data, t_col, method) {
  fm <- as.formula(formula_str)

  if (identical(method, "gls_ar1")) {
    fit <- tryCatch(
      gls(
        fm,
        data = data,
        correlation = corAR1(form = as.formula(paste0("~", t_col))),
        method = "ML",
        na.action = na.omit
      ),
      error = function(e) NULL
    )
    return(list(model = fit, fit_method = if (!is.null(fit)) "gls_ar1" else NA_character_))
  }

  if (identical(method, "gls_no_ar1")) {
    fit <- tryCatch(
      gls(
        fm,
        data = data,
        method = "ML",
        na.action = na.omit
      ),
      error = function(e) NULL
    )
    return(list(model = fit, fit_method = if (!is.null(fit)) "gls_no_ar1" else NA_character_))
  }

  if (identical(method, "lm")) {
    fit <- tryCatch(
      lm(fm, data = data),
      error = function(e) NULL
    )
    return(list(model = fit, fit_method = if (!is.null(fit)) "lm" else NA_character_))
  }

  list(model = NULL, fit_method = NA_character_)
}

extract_term_stats <- function(fit_obj, term_name) {
  fit <- fit_obj$model

  if (is.null(fit)) {
    return(list(
      estimate   = NA_real_,
      t_value    = NA_real_,
      p_value    = NA_real_,
      AIC        = NA_real_,
      BIC        = NA_real_,
      logLik     = NA_real_,
      fit_method = fit_obj$fit_method
    ))
  }

  if (inherits(fit, "gls")) {
    tab <- summary(fit)$tTable
    est <- if (term_name %in% rownames(tab)) tab[term_name, "Value"]   else NA_real_
    tv  <- if (term_name %in% rownames(tab)) tab[term_name, "t-value"] else NA_real_
    pv  <- if (term_name %in% rownames(tab)) tab[term_name, "p-value"] else NA_real_
  } else if (inherits(fit, "lm")) {
    tab <- summary(fit)$coefficients
    est <- if (term_name %in% rownames(tab)) tab[term_name, "Estimate"] else NA_real_
    tv  <- if (term_name %in% rownames(tab)) tab[term_name, "t value"]  else NA_real_
    pv  <- if (term_name %in% rownames(tab)) tab[term_name, "Pr(>|t|)"] else NA_real_
  } else {
    est <- tv <- pv <- NA_real_
  }

  list(
    estimate   = est,
    t_value    = tv,
    p_value    = pv,
    AIC        = tryCatch(AIC(fit), error = function(e) NA_real_),
    BIC        = tryCatch(BIC(fit), error = function(e) NA_real_),
    logLik     = tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_),
    fit_method = fit_obj$fit_method
  )
}

calc_partial_r2_from_t <- function(t_value, df_resid) {
  if (is.na(t_value) || is.na(df_resid) || df_resid <= 0) return(NA_real_)
  (t_value^2) / (t_value^2 + df_resid)
}

safe_cor_stats <- function(x, y, method = "pearson") {
  out <- tryCatch(
    suppressWarnings(cor.test(x, y, method = method)),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(list(estimate = NA_real_, p_value = NA_real_))
  }
  list(
    estimate = unname(out$estimate),
    p_value  = out$p.value
  )
}

check_regular_time <- function(x, expected_step = 1, tol = 1e-8) {
  dx <- abs(diff(x))
  all(abs(dx - expected_step) < tol)
}

run_null_test_gls_sync <- function(df_common, resp_var_z, pred_col_z,
                                   env_terms, pulse_terms, t_col,
                                   n_iter = 1000, min_shift = 10) {
  form_str <- build_formula_string(
    resp = resp_var_z,
    t_col = t_col,
    env_terms = env_terms,
    pulse_terms = pulse_terms,
    predictor_term = ".pred_tmp"
  )

  df0 <- df_common
  df0$.pred_tmp <- df0[[pred_col_z]]

  if (!check_regular_time(df0[[t_col]], expected_step = 1)) {
    warning("Time series used in null test is not strictly regular after drop_na(); circular shift is being applied over observed rows, not exact 1-Ma positions.")
  }

  real_fit    <- fit_model_safe(form_str, df0, t_col)
  real_stats  <- extract_term_stats(real_fit, ".pred_tmp")
  real_t      <- real_stats$t_value
  real_method <- real_fit$fit_method

  if (is.na(real_t) || is.na(real_method)) {
    return(list(
      p_val = NA_real_,
      real_t = NA_real_,
      valid_n = 0,
      fit_method = NA_character_
    ))
  }

  n <- nrow(df0)
  valid_k <- get_valid_shifts(n, min_shift = min_shift)
  if (length(valid_k) == 0) {
    return(list(
      p_val = NA_real_,
      real_t = real_t,
      valid_n = 0,
      fit_method = real_method
    ))
  }

  null_t <- rep(NA_real_, n_iter)

  for (i in seq_len(n_iter)) {
    k <- sample(valid_k, 1)
    dfi <- df0
    dfi$.pred_tmp <- circular_shift(dfi$.pred_tmp, k)

    fit_i <- fit_model_fixed_method(form_str, dfi, t_col, real_method)
    st_i  <- extract_term_stats(fit_i, ".pred_tmp")
    null_t[i] <- st_i$t_value
  }

  valid_idx <- which(!is.na(null_t))
  valid_n <- length(valid_idx)

  if (valid_n == 0) {
    return(list(
      p_val = NA_real_,
      real_t = real_t,
      valid_n = 0,
      fit_method = real_method
    ))
  }

  b <- sum(abs(null_t[valid_idx]) >= abs(real_t))
  p_null <- (b + 1) / (valid_n + 1)

  list(
    p_val = p_null,
    real_t = real_t,
    valid_n = valid_n,
    fit_method = real_method
  )
}

run_null_test_lm_sync <- function(df_common, resp_var_z, pred_col_z,
                                  env_terms, pulse_terms, t_col,
                                  n_iter = 1000, min_shift = 10) {
  form_str <- build_formula_string(
    resp = resp_var_z,
    t_col = t_col,
    env_terms = env_terms,
    pulse_terms = pulse_terms,
    predictor_term = ".pred_tmp"
  )

  df0 <- df_common
  df0$.pred_tmp <- df0[[pred_col_z]]

  if (!check_regular_time(df0[[t_col]], expected_step = 1)) {
    warning("Time series used in null test is not strictly regular after drop_na(); circular shift is being applied over observed rows, not exact 1-Ma positions.")
  }

  real_fit   <- fit_model_fixed_method(form_str, df0, t_col, method = "lm")
  real_stats <- extract_term_stats(real_fit, ".pred_tmp")
  real_t     <- real_stats$t_value

  if (is.na(real_t)) {
    return(list(
      p_val = NA_real_,
      real_t = NA_real_,
      valid_n = 0,
      fit_method = "lm"
    ))
  }

  n <- nrow(df0)
  valid_k <- get_valid_shifts(n, min_shift = min_shift)
  if (length(valid_k) == 0) {
    return(list(
      p_val = NA_real_,
      real_t = real_t,
      valid_n = 0,
      fit_method = "lm"
    ))
  }

  null_t <- rep(NA_real_, n_iter)

  for (i in seq_len(n_iter)) {
    k <- sample(valid_k, 1)
    dfi <- df0
    dfi$.pred_tmp <- circular_shift(dfi$.pred_tmp, k)

    fit_i <- fit_model_fixed_method(form_str, dfi, t_col, method = "lm")
    st_i  <- extract_term_stats(fit_i, ".pred_tmp")
    null_t[i] <- st_i$t_value
  }

  valid_idx <- which(!is.na(null_t))
  valid_n <- length(valid_idx)

  if (valid_n == 0) {
    return(list(
      p_val = NA_real_,
      real_t = real_t,
      valid_n = 0,
      fit_method = "lm"
    ))
  }

  b <- sum(abs(null_t[valid_idx]) >= abs(real_t))
  p_null <- (b + 1) / (valid_n + 1)

  list(
    p_val = p_null,
    real_t = real_t,
    valid_n = valid_n,
    fit_method = "lm"
  )
}

## =========================================================
## 3. Read data, raw only
## =========================================================

dat_raw <- read.csv(infile, check.names = FALSE)

response_vars_final <- response_vars_user
predictor_var_final <- predictor_var_user
analysis_type <- "raw"

numeric_cols_pre <- unique(c(
  time_col,
  unname(response_vars_final),
  predictor_var_final,
  env_vars_final,
  pulse_vars
))
numeric_cols_pre <- intersect(numeric_cols_pre, names(dat_raw))

dat <- dat_raw
for (v in numeric_cols_pre) {
  dat[[v]] <- as.numeric(dat[[v]])
}

dat <- dat %>% arrange(desc(.data[[time_col]]))

if (any(is.na(dat[[time_col]]))) {
  stop("time_col contains NA after numeric conversion.")
}

if (any(duplicated(dat[[time_col]]))) {
  stop(time_col, " has duplicated values. Please fix before running.")
}

time_steps <- abs(diff(dat[[time_col]]))
cat("Time step summary (Ma):\n")
print(summary(time_steps))
if (sd(time_steps, na.rm = TRUE) > 1e-3) {
  warning("Time steps are irregular. Interpret AR1 and circular-shift results with extra caution.")
}

## =========================================================
## 4. Check required columns
## =========================================================

needed_cols <- unique(c(
  time_col,
  unname(response_vars_final),
  predictor_var_final,
  env_vars_final,
  pulse_vars
))

missing_cols <- setdiff(needed_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
}

## =========================================================
## 5. Standardization, raw only
## =========================================================

vars_to_scale <- unique(c(
  unname(response_vars_final),
  predictor_var_final,
  env_vars_final
))

for (v in vars_to_scale) {
  dat[[paste0(v, "_z")]] <- safe_scale(dat[[v]])
}

for (pv in pulse_vars) {
  dat[[pv]] <- as.numeric(dat[[pv]])
}

predictor_var_z <- paste0(predictor_var_final, "_z")
env_terms_z     <- if (length(env_vars_final) > 0) paste0(env_vars_final, "_z") else character(0)

## =========================================================
## 6. Main synchronous analysis
## =========================================================

all_model_rows    <- list()
sync_support_rows <- list()
all_fits          <- list()

set.seed(SEED_NULL)

for (r_name in names(response_vars_final)) {
  cat("\nRunning raw synchronous models for:", r_name, "...\n")

  resp_col   <- response_vars_final[[r_name]]
  resp_var_z <- paste0(resp_col, "_z")

  req_cols <- unique(c(
    time_col,
    resp_var_z,
    predictor_var_z,
    env_terms_z,
    pulse_vars
  ))

  df_common <- dat %>% tidyr::drop_na(dplyr::all_of(req_cols))
  cat("  Common N rows:", nrow(df_common), "\n")

  if (nrow(df_common) < 8) {
    warning("Too few complete rows for ", r_name, "; skipping.")
    next
  }

  form_m0 <- build_formula_string(
    resp = resp_var_z,
    t_col = time_col,
    env_terms = character(0),
    pulse_terms = character(0),
    predictor_term = NULL
  )

  form_m1 <- build_formula_string(
    resp = resp_var_z,
    t_col = time_col,
    env_terms = env_terms_z,
    pulse_terms = pulse_vars,
    predictor_term = NULL
  )

  form_m2 <- build_formula_string(
    resp = resp_var_z,
    t_col = time_col,
    env_terms = env_terms_z,
    pulse_terms = pulse_vars,
    predictor_term = predictor_var_z
  )

  fit_m0 <- fit_model_safe(form_m0, df_common, time_col)
  fit_m1 <- fit_model_safe(form_m1, df_common, time_col)
  fit_m2 <- fit_model_safe(form_m2, df_common, time_col)

  if (STORE_FITS) {
    all_fits[[r_name]] <- list(
      M0_time = fit_m0,
      M1_env  = fit_m1,
      M2_sync = fit_m2
    )
  }

  st_m0 <- extract_term_stats(fit_m0, predictor_var_z)
  st_m1 <- extract_term_stats(fit_m1, predictor_var_z)
  st_m2 <- extract_term_stats(fit_m2, predictor_var_z)

  model_df <- bind_rows(
    data.frame(
      data_type = analysis_type,
      response = r_name,
      response_col = resp_col,
      predictor_label = predictor_label_user,
      model_name = "M0_time",
      n = nrow(df_common),
      fit_method = st_m0$fit_method,
      pred_term = NA_character_,
      pred_estimate = NA_real_,
      pred_t = NA_real_,
      pred_p = NA_real_,
      AIC = st_m0$AIC,
      BIC = st_m0$BIC,
      logLik = st_m0$logLik
    ),
    data.frame(
      data_type = analysis_type,
      response = r_name,
      response_col = resp_col,
      predictor_label = predictor_label_user,
      model_name = "M1_env",
      n = nrow(df_common),
      fit_method = st_m1$fit_method,
      pred_term = NA_character_,
      pred_estimate = NA_real_,
      pred_t = NA_real_,
      pred_p = NA_real_,
      AIC = st_m1$AIC,
      BIC = st_m1$BIC,
      logLik = st_m1$logLik
    ),
    data.frame(
      data_type = analysis_type,
      response = r_name,
      response_col = resp_col,
      predictor_label = predictor_label_user,
      model_name = "M2_sync",
      n = nrow(df_common),
      fit_method = st_m2$fit_method,
      pred_term = predictor_var_z,
      pred_estimate = st_m2$estimate,
      pred_t = st_m2$t_value,
      pred_p = st_m2$p_value,
      AIC = st_m2$AIC,
      BIC = st_m2$BIC,
      logLik = st_m2$logLik
    )
  )

  if (all(is.na(model_df$AIC))) {
    model_df$deltaAIC <- NA_real_
  } else {
    min_aic <- min(model_df$AIC, na.rm = TRUE)
    model_df$deltaAIC <- model_df$AIC - min_aic
  }

  all_model_rows[[r_name]] <- model_df

  pearson_test  <- safe_cor_stats(df_common[[resp_var_z]], df_common[[predictor_var_z]], method = "pearson")
  spearman_test <- safe_cor_stats(df_common[[resp_var_z]], df_common[[predictor_var_z]], method = "spearman")

  lm_full <- tryCatch(lm(as.formula(form_m2), data = df_common), error = function(e) NULL)
  if (!is.null(lm_full)) {
    lm_tab <- summary(lm_full)$coefficients
    lm_beta <- if (predictor_var_z %in% rownames(lm_tab)) lm_tab[predictor_var_z, "Estimate"] else NA_real_
    lm_t    <- if (predictor_var_z %in% rownames(lm_tab)) lm_tab[predictor_var_z, "t value"]  else NA_real_
    lm_p    <- if (predictor_var_z %in% rownames(lm_tab)) lm_tab[predictor_var_z, "Pr(>|t|)"] else NA_real_
    lm_partial_r2 <- calc_partial_r2_from_t(lm_t, df.residual(lm_full))
  } else {
    lm_beta <- lm_t <- lm_p <- lm_partial_r2 <- NA_real_
  }

  null_gls <- run_null_test_gls_sync(
    df_common   = df_common,
    resp_var_z  = resp_var_z,
    pred_col_z  = predictor_var_z,
    env_terms   = env_terms_z,
    pulse_terms = pulse_vars,
    t_col       = time_col,
    n_iter      = N_NULL,
    min_shift   = MIN_SHIFT
  )

  null_lm <- run_null_test_lm_sync(
    df_common   = df_common,
    resp_var_z  = resp_var_z,
    pred_col_z  = predictor_var_z,
    env_terms   = env_terms_z,
    pulse_terms = pulse_vars,
    t_col       = time_col,
    n_iter      = N_NULL,
    min_shift   = MIN_SHIFT
  )

  aic_vals <- c(
    M0_time = st_m0$AIC,
    M1_env  = st_m1$AIC,
    M2_sync = st_m2$AIC
  )

  if (all(is.na(aic_vals))) {
    best_model <- NA_character_
    deltaAIC_vs_best <- NA_real_
  } else {
    aic_vals_tmp <- aic_vals
    aic_vals_tmp[is.na(aic_vals_tmp)] <- Inf
    best_model <- names(which.min(aic_vals_tmp))[1]
    deltaAIC_vs_best <- st_m2$AIC - min(aic_vals, na.rm = TRUE)
  }

  support_score <- sum(c(
    !is.na(st_m2$estimate),
    !is.na(st_m2$p_value) && st_m2$p_value < 0.05,
    !is.na(st_m2$AIC) && !is.na(st_m1$AIC) && (st_m2$AIC < st_m1$AIC),
    !is.na(null_gls$p_val) && null_gls$p_val < 0.05,
    !is.na(null_lm$p_val)  && null_lm$p_val  < 0.05
  ))

  sync_support_rows[[r_name]] <- data.frame(
    data_type = analysis_type,
    response = r_name,
    response_col = resp_col,
    predictor_col = predictor_var_final,
    predictor_label = predictor_label_user,
    n = nrow(df_common),

    pearson_r = pearson_test$estimate,
    pearson_p = pearson_test$p_value,
    spearman_rho = spearman_test$estimate,
    spearman_p = spearman_test$p_value,

    model_name = "M2_sync",
    fit_method_sync = st_m2$fit_method,
    pred_term = predictor_var_z,
    beta_sync = st_m2$estimate,
    t_sync = st_m2$t_value,
    p_sync = st_m2$p_value,

    AIC_M0 = st_m0$AIC,
    AIC_M1 = st_m1$AIC,
    AIC_M2 = st_m2$AIC,
    deltaAIC_vs_M0 = st_m2$AIC - st_m0$AIC,
    deltaAIC_vs_M1 = st_m2$AIC - st_m1$AIC,
    deltaAIC_vs_best = deltaAIC_vs_best,
    best_model = best_model,

    null_P_value = null_gls$p_val,
    gls_real_t = null_gls$real_t,
    gls_valid_n = null_gls$valid_n,
    null_fit_method = null_gls$fit_method,

    lm_beta = lm_beta,
    lm_t = lm_t,
    lm_p = lm_p,
    partial_R2_lm = lm_partial_r2,
    lm_null_P = null_lm$p_val,
    lm_real_t = null_lm$real_t,
    lm_valid_n = null_lm$valid_n,
    lm_null_fit_method = null_lm$fit_method,

    support_score = support_score
  )
}

## =========================================================
## 7. Combine and export
## =========================================================

all_models   <- bind_rows(all_model_rows)
sync_support <- bind_rows(sync_support_rows)

if (nrow(all_models) == 0 || nrow(sync_support) == 0) {
  stop("No valid results were produced.")
}

write.csv(all_models,   paste0(output_prefix, "_all_models.csv"), row.names = FALSE)
write.csv(sync_support, paste0(output_prefix, "_sync_support_table.csv"), row.names = FALSE)

if (STORE_FITS) {
  saveRDS(all_fits, paste0(output_prefix, "_all_fits.rds"))
}

## =========================================================
## 8. Optional quick plots
## =========================================================

plot_df1 <- all_models %>%
  mutate(model_name = factor(model_name, levels = c("M0_time", "M1_env", "M2_sync")))

p1 <- ggplot(plot_df1, aes(x = model_name, y = deltaAIC, fill = response)) +
  geom_col(position = position_dodge(width = 0.8)) +
  theme_bw(base_size = 12) +
  labs(x = "Model", y = expression(Delta*AIC), title = "Raw 0-lag synchronous model comparison")

ggsave(paste0(output_prefix, "_plot_modelBar_sync.pdf"), p1, width = 8, height = 5)

plot_df2 <- sync_support %>%
  arrange(beta_sync) %>%
  mutate(response = factor(response, levels = response))

p2 <- ggplot(plot_df2, aes(x = response, y = beta_sync, fill = response)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Standardized beta (M2_sync)", title = "Raw synchronous association effect sizes") +
  theme(legend.position = "none")

ggsave(paste0(output_prefix, "_plot_sync_beta.pdf"), p2, width = 7, height = 5)

cat("\nDone!\n")
cat("Main outputs:\n")
cat("  - ", paste0(output_prefix, "_all_models.csv"), "\n", sep = "")
cat("  - ", paste0(output_prefix, "_sync_support_table.csv"), "\n", sep = "")
if (STORE_FITS) {
  cat("  - ", paste0(output_prefix, "_all_fits.rds"), "\n", sep = "")
}
cat("  - ", paste0(output_prefix, "_plot_modelBar_sync.pdf"), "\n", sep = "")
cat("  - ", paste0(output_prefix, "_plot_sync_beta.pdf"), "\n", sep = "")
