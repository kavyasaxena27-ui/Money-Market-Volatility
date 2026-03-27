# =============================================================================
# 04_instrument_module.R
# Rigobon (2003) / Lewbel (2012) IV construction
#
# Follows Section 2.3.3 of Holm-Hadulla & Pool (ECB WP 3048):
#
#   Eq (2): sigma_{t-1}     ~ X_{t-p}  -> residuals e1
#   Eq (3): sigma_{t-1}*S_t ~ X_{t-p}  -> residuals e2
#   Instruments: iv1_k = (Z_k - Zbar_k) * e1   [for sigma]
#                iv2_k = (Z_k - Zbar_k) * e2   [for sigma*S]
#   k = 1:4 regime dummies -> 8 instruments total
#
# SOURCE: at top of 04_local_projections.R
# CALL:   after df is built (post-centering, pre-LP loop)
#
#   df <- df %>%
#     make_instrument(vol_col = vol_col, shock_col = "target",
#                     fs_controls = fs_controls_iv)
#   check_instrument(df)
# =============================================================================

library(car)
library(lmtest)

# -----------------------------------------------------------------------------
# Regime dates — update to your identified BoE shift dates
# -----------------------------------------------------------------------------

BOE_REGIME_DATES <- list(
  r1_end = as.Date("2009-03-05"),   # Corridor -> floor (QE onset)
  r2_end = as.Date("2012-02-01"),   # Post-GFC high vol -> ample reserves
  r3_end = as.Date("2016-08-04"),   # Ample -> pre-SONIA reform
  r4_end = as.Date("2022-09-01")    # Post-SONIA -> post-Truss
)

# -----------------------------------------------------------------------------
# make_instrument()
#
# Regime dummies and interaction are created inline here — no separate function.
# Regime dummies also enter X_{t-p} controls per Lewbel (2012).
#
# INPUTS:
#   df          - data frame, must already have vol_col, shock_col, fs_controls
#   vol_col     - string, centred vol column e.g. "centered_sonia_20d_sd"
#   shock_col   - string, e.g. "target"
#   fs_controls - character vector matching X_{t-p} in the LP
#   dates       - regime date list (defaults to BOE_REGIME_DATES above)
#
# OUTPUTS: df with added columns:
#   regime, dummy_r1:r4         - regime factor and shift dummies
#   vol_x_shock                 - sigma * S (LHS of eq 3)
#   fs1_resid, fs2_resid        - residuals from eq (2) and (3)
#   iv1_dummy_r1:r4             - (Z_k - Zbar) * e1
#   iv2_dummy_r1:r4             - (Z_k - Zbar) * e2
# -----------------------------------------------------------------------------

make_instrument <- function(df, vol_col, shock_col = "target",
                            fs_controls, dates = BOE_REGIME_DATES) {

  stopifnot(
    vol_col   %in% names(df),
    shock_col %in% names(df),
    all(fs_controls %in% names(df))
  )

  # -- Regime dummies (inlined) -----------------------------------------------
  df <- df %>%
    mutate(
      regime   = case_when(
        date <  dates$r1_end ~ "1_corridor",
        date <  dates$r2_end ~ "2_post_gfc",
        date <  dates$r3_end ~ "3_ample",
        date <  dates$r4_end ~ "4_post_sonia",
        TRUE                 ~ "5_post_truss"
      ),
      regime   = factor(regime),
      dummy_r1 = as.integer(date >= dates$r1_end),
      dummy_r2 = as.integer(date >= dates$r2_end),
      dummy_r3 = as.integer(date >= dates$r3_end),
      dummy_r4 = as.integer(date >= dates$r4_end),
      # LHS of equation (3)
      vol_x_shock = .data[[vol_col]] * .data[[shock_col]]
    )

  # -- First-stage controls: X_{t-p} + regime dummies (Lewbel 2012) -----------
  regime_dummies  <- c("dummy_r1", "dummy_r2", "dummy_r3", "dummy_r4")
  all_fs_controls <- unique(c(fs_controls, regime_dummies))
  rhs             <- paste(all_fs_controls, collapse = " + ")

  # -- Equation (2): sigma_{t-1} ~ X_{t-p} ------------------------------------
  fml2     <- as.formula(paste(vol_col, "~", rhs))
  fs1      <- lm(fml2, data = df, na.action = na.exclude)
  df$fs1_resid <- residuals(fs1)

  # -- Equation (3): sigma_{t-1}*S_t ~ X_{t-p} --------------------------------
  fml3     <- as.formula(paste("vol_x_shock ~", rhs))
  fs2      <- lm(fml3, data = df, na.action = na.exclude)
  df$fs2_resid <- residuals(fs2)

  # -- Instruments: (Z_k - Zbar_k) * e_hat ------------------------------------
  df <- df %>%
    mutate(
      across(
        all_of(regime_dummies),
        list(
          iv1 = \(z) (z - mean(z, na.rm = TRUE)) * fs1_resid,
          iv2 = \(z) (z - mean(z, na.rm = TRUE)) * fs2_resid
        ),
        .names = "{.fn}_{.col}"
      )
    )

  attr(df, "fs1")       <- fs1
  attr(df, "fs2")       <- fs2
  attr(df, "vol_col")   <- vol_col
  attr(df, "shock_col") <- shock_col

  message(sprintf("checkmark: Instruments built for [%s] and [%s x %s]",
                  vol_col, vol_col, shock_col))
  message("  iv1_dummy_r1:r4 -> use for ", vol_col)
  message("  iv2_dummy_r1:r4 -> use for inter_term (vol x shock)")

  return(df)
}

# -----------------------------------------------------------------------------
# check_instrument()
# Replicates Tables 2 & 3 of the paper diagnostics
# -----------------------------------------------------------------------------

check_instrument <- function(df) {

  vol_col <- attr(df, "vol_col")
  fs1     <- attr(df, "fs1")
  fs2     <- attr(df, "fs2")

  if (is.null(fs1) || is.null(vol_col))
    stop("Attributes lost — rerun make_instrument() first.")

  cat("\n================================================================\n")
  cat("  INSTRUMENT DIAGNOSTICS\n")
  cat(sprintf("  Vol measure: %s\n", vol_col))
  cat("================================================================\n\n")

  # 1. Within-regime descriptives
  cat("-- 1. Within-regime volatility (cf. Table 5) -------------------\n")
  df %>%
    group_by(regime) %>%
    summarise(
      n    = n(),
      mean = round(mean(.data[[vol_col]], na.rm = TRUE), 5),
      sd   = round(sd(.data[[vol_col]],   na.rm = TRUE), 5),
      var  = round(var(.data[[vol_col]],   na.rm = TRUE), 5),
      .groups = "drop"
    ) %>% print()
  cat("\n")

  # 2. Breusch-Pagan on both first stages (replicates Table 3)
  cat("-- 2. Breusch-Pagan on first stages (cf. Table 3) -------------\n")
  cat("   Want: p < 0.05 for BOTH -> instruments not weak\n\n")
  bp1 <- bptest(fs1)
  bp2 <- bptest(fs2)
  cat(sprintf("   Eq (2) [sigma ~ X]:    chi2(%g) = %.2f,  p = %.4f\n",
              bp1$parameter, bp1$statistic, bp1$p.value))
  cat(sprintf("   Eq (3) [sigma*S ~ X]:  chi2(%g) = %.2f,  p = %.4f\n\n",
              bp2$parameter, bp2$statistic, bp2$p.value))

  # 3. First-stage F-statistics
  cat("-- 3. First-stage F-statistics ---------------------------------\n")
  cat("   Want: F > 10\n\n")
  for (tag in c("Eq (2)", "Eq (3)")) {
    fs  <- if (tag == "Eq (2)") fs1 else fs2
    fst <- summary(fs)$fstatistic
    fsp <- pf(fst[1], fst[2], fst[3], lower.tail = FALSE)
    cat(sprintf("   %s:  F(%g, %g) = %.2f,  p = %.4f\n",
                tag, fst[2], fst[3], fst[1], fsp))
  }
  cat("\n")

  # 4. R-squared
  cat("-- 4. First-stage R-squared ------------------------------------\n")
  cat(sprintf("   Eq (2):  R2 = %.3f\n", summary(fs1)$r.squared))
  cat(sprintf("   Eq (3):  R2 = %.3f\n\n", summary(fs2)$r.squared))

  # 5. Levene test
  cat("-- 5. Levene test (variance differs across regimes) -----------\n")
  cat("   Want: p < 0.05\n")
  lev_fml <- as.formula(paste(vol_col, "~ regime"))
  print(leveneTest(lev_fml, data = df))
  cat("\n")

  # 6. Diagnostic plot
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
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "SONIA/RONIA Volatility across BoE regimes",
      x = NULL, y = vol_col, fill = "Regime"
    ) +
    theme_minimal(base_size = 12)

  print(p)

  cat("================================================================\n")
  cat("  NOTE: Pagan-Hall test (cf. Table 2) run AFTER LP estimation.\n")
  cat("================================================================\n\n")

  invisible(df)
}
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
