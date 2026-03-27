# =============================================================================
# 04_instrument_module.R
# Rigobon (2003) IV construction + diagnostics
#
# INTEGRATION: Source this at the end of 03b_create_balanced_dataset.R
# It reads balanced_data and PARAM_VOL_MEASURE from your environment,
# adds vol_hat and regime columns, and overwrites balanced_data + the .rds
#
# USAGE (add to bottom of 03b):
#   source(here::here("R", "04_instrument_module.R"))
# =============================================================================

library(car)
library(lmtest)

# -----------------------------------------------------------------------------
# CONFIG: BoE regime shift dates
# Adjust to match your identified dates
# -----------------------------------------------------------------------------

BOE_REGIME_DATES <- list(
  corridor_end     = as.Date("2009-03-05"),  # QE1 onset / floor system begins
  qe_stable_start  = as.Date("2012-07-01"),  # Stable floor, pre-Brexit vote
  pre_sonia_end    = as.Date("2018-04-23"),  # SONIA reform
  post_truss_start = as.Date("2022-09-01")   # Post-Truss gilt stress
)

# -----------------------------------------------------------------------------
# STEP 1: Add regime column
# -----------------------------------------------------------------------------

make_regimes <- function(df, dates = BOE_REGIME_DATES) {
  df %>%
    mutate(
      regime = case_when(
        date < dates$corridor_end     ~ "1_corridor",
        date < dates$qe_stable_start  ~ "2_floor_qe1",
        date < dates$pre_sonia_end    ~ "3_floor_stable",
        date < dates$post_truss_start ~ "4_post_sonia",
        TRUE                          ~ "5_post_truss"
      ),
      regime = factor(regime)
    )
}

# -----------------------------------------------------------------------------
# STEP 2: First stage — Rigobon instrument
# vol_col: the centred volatility column name (string), e.g. "centered_sonia_20d_sd"
# Adds: vol_hat (fitted values), vol_resid (residuals)
# -----------------------------------------------------------------------------

make_instrument <- function(df, vol_col) {
  
  stopifnot(
    "regime" %in% names(df),
    vol_col  %in% names(df)
  )
  
  fml <- as.formula(paste(vol_col, "~ regime"))
  first_stage   <- lm(fml, data = df)
  df$vol_hat    <- fitted(first_stage)
  df$vol_resid  <- residuals(first_stage)
  
  attr(df, "first_stage") <- first_stage
  attr(df, "vol_col")     <- vol_col
  
  message(sprintf("✓ vol_hat constructed from: %s ~ regime", vol_col))
  return(df)
}

# -----------------------------------------------------------------------------
# STEP 3: Diagnostics
# -----------------------------------------------------------------------------

check_instrument <- function(df) {
  
  vol_col     <- attr(df, "vol_col")
  first_stage <- attr(df, "first_stage")
  
  # Refit if attributes lost through dplyr operations
  if (is.null(first_stage) || is.null(vol_col)) {
    stop("Attributes lost — rerun make_instrument(df, vol_col) first.")
  }
  
  cat("\n========================================================\n")
  cat("  INSTRUMENT DIAGNOSTICS: Rigobon Regime-Based IV\n")
  cat(sprintf("  Volatility measure: %s\n", vol_col))
  cat("========================================================\n\n")
  
  # 1. Within-regime descriptives
  cat("── 1. Within-regime volatility ──────────────────────────\n")
  df %>%
    group_by(regime) %>%
    summarise(
      n    = n(),
      mean = round(mean(.data[[vol_col]], na.rm = TRUE), 5),
      sd   = round(sd(.data[[vol_col]],   na.rm = TRUE), 5),
      .groups = "drop"
    ) %>%
    print()
  cat("\n")
  
  # 2. Levene test — key validity check
  cat("── 2. Levene test (H0: equal variances across regimes) ──\n")
  lev_fml <- as.formula(paste(vol_col, "~ regime"))
  print(leveneTest(lev_fml, data = df))
  cat("  ✓ Want: p < 0.05\n\n")
  
  # 3. First-stage F-stat
  cat("── 3. First-stage F-statistic ───────────────────────────\n")
  fs  <- summary(first_stage)
  fst <- fs$fstatistic
  fsp <- pf(fst[1], fst[2], fst[3], lower.tail = FALSE)
  cat(sprintf("  F(%g, %g) = %.2f,  p = %.4f\n", fst[2], fst[3], fst[1], fsp))
  cat("  ✓ Want: F > 10\n\n")
  
  # 4. R-squared
  cat("── 4. First-stage R² ────────────────────────────────────\n")
  cat(sprintf("  R² = %.3f,  Adj. R² = %.3f\n\n", fs$r.squared, fs$adj.r.squared))
  
  # 5. Breusch-Pagan — confirms heteroskedasticity motivating Rigobon
  cat("── 5. Breusch-Pagan (vol ~ time trend) ──────────────────\n")
  bp_fml <- as.formula(paste(vol_col, "~ as.numeric(date)"))
  print(bptest(lm(bp_fml, data = df)))
  cat("  ✓ Want: p < 0.05 (confirms raw vol is heteroskedastic)\n\n")
  
  # 6. Diagnostic plot
  cat("── 6. Plot ──────────────────────────────────────────────\n")
  regime_bands <- df %>%
    group_by(regime) %>%
    summarise(xmin = min(date), xmax = max(date), .groups = "drop")
  
  p <- ggplot(df, aes(x = date)) +
    geom_rect(
      data = regime_bands,
      aes(xmin = xmin, xmax = xmax, fill = regime),
      ymin = -Inf, ymax = Inf, alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_line(aes(y = .data[[vol_col]]), colour = "grey40", linewidth = 0.5) +
    geom_line(aes(y = vol_hat), colour = "firebrick",
              linewidth = 1, linetype = "dashed") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title    = "SONIA Volatility: Actual vs Instrumented (Rigobon)",
      subtitle = sprintf("Actual = %s  |  Dashed = vol_hat from regime first stage", vol_col),
      x = NULL, y = "Centred volatility", fill = "Regime"
    ) +
    theme_minimal(base_size = 12)
  
  print(p)
  
  cat("\n========================================================\n")
  cat("  Done. Instrument is valid if F > 10 and Levene p < 0.05\n")
  cat("========================================================\n\n")
  
  invisible(df)
}

# -----------------------------------------------------------------------------
# STEP 4: Add instrument interaction term (vol_hat × shock)
# This is what enters the IV formula in the LP as the excluded instrument
# -----------------------------------------------------------------------------

make_iv_interaction <- function(df, shock_col = "target") {
  
  stopifnot("vol_hat" %in% names(df), shock_col %in% names(df))
  
  df <- df %>%
    mutate(
      vol_hat_x_shock     = vol_hat * .data[[shock_col]],
      lag1_vol_hat_x_shock = lag(vol_hat_x_shock, 1),
      lag2_vol_hat_x_shock = lag(vol_hat_x_shock, 2)
    )
  
  message("✓ vol_hat_x_shock added (instrument interaction).")
  return(df)
}

# =============================================================================
# RUN — executes when sourced from 03b
# Expects balanced_data, PARAM_VOL_MEASURE, project_paths in environment
# =============================================================================

message("\n--- Running instrument module ---")

balanced_data <- balanced_data %>%
  make_regimes() %>%
  make_instrument(vol_col = vol_col) %>%   # vol_col already defined in 03b
  make_iv_interaction(shock_col = "target")

check_instrument(balanced_data)

# Overwrite processed dataset so LP script picks up vol_hat automatically
readr::write_rds(
  balanced_data,
  fs::path(project_paths[["data_processed"]], "monthly_all_series_WIDE.rds")
)

message("✓ balanced_data updated with regime, vol_hat, vol_hat_x_shock and saved.")
