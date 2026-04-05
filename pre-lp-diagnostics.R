# =============================================================================
# 02b_pre_lp_diagnostics.R
# Pre-LP diagnostic tests
# 1. Normality of key variables
# 2. Structural break tests on vol series
# 3. Granger causality: vol -> outcomes
# 4. Cross-correlation plots: vol vs outcomes
#
# SOURCE: Run after 03b_create_balanced_dataset.R and instrument module
# INPUT:  df (with regime dummies and centered_vol already constructed)
# =============================================================================

source(here::here("R", "00_setup.R"))
library(strucchange)   # Bai-Perron / Chow
library(vars)          # Granger causality
library(tseries)       # Jarque-Bera
library(patchwork)

# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------

# Variables to test
outcome_vars <- c(
  "log_ecy2_dp_m", "core", "household_mortgage_rates",
  "gb2yt", "gb5yt", "gb10yt"
)

vol_var <- "centered_vol"   # update if using sonia vs ronia

# Load data — assumes df already exists in environment from 04_local_projections
# If running standalone, load here:
# df <- readr::read_rds(path(project_paths[["data_processed"]], "monthly_all_series_WIDE.rds"))

# -----------------------------------------------------------------------------
# 1. NORMALITY — Jarque-Bera test + QQ plots
# -----------------------------------------------------------------------------

cat("\n================================================================\n")
cat("  1. NORMALITY TESTS (Jarque-Bera)\n")
cat("  H0: variable is normally distributed\n")
cat("  Want: context-dependent — flag severe departures\n")
cat("================================================================\n\n")

jb_results <- c(vol_var, outcome_vars) %>%
  map_dfr(function(v) {
    x <- df[[v]]
    x <- x[!is.na(x)]
    jb <- jarque.bera.test(x)
    tibble(
      variable  = v,
      statistic = round(jb$statistic, 3),
      p_value   = round(jb$p.value, 4),
      skewness  = round(moments::skewness(x), 3),
      kurtosis  = round(moments::kurtosis(x), 3)
    )
  })

print(jb_results)

# QQ plots
qq_plots <- c(vol_var, outcome_vars) %>%
  map(function(v) {
    ggplot(df %>% filter(!is.na(.data[[v]])),
           aes(sample = .data[[v]])) +
      stat_qq(size = 0.8, colour = "steelblue") +
      stat_qq_line(colour = "firebrick") +
      labs(title = v, x = "Theoretical", y = "Sample") +
      theme_minimal(base_size = 10)
  })

wrap_plots(qq_plots, ncol = 3) +
  plot_annotation(title = "QQ Plots — key variables")

# -----------------------------------------------------------------------------
# 2. STRUCTURAL BREAK TESTS on vol series
# Bai-Perron multiple breakpoint test
# Want: breaks close to your regime dates
# -----------------------------------------------------------------------------

cat("\n================================================================\n")
cat("  2. STRUCTURAL BREAK TESTS (Bai-Perron)\n")
cat("  H0: no structural break\n")
cat("  Want: breaks near your regime dates\n")
cat("================================================================\n\n")

vol_ts <- df %>%
  filter(!is.na(.data[[vol_var]])) %>%
  arrange(date) %>%
  pull(.data[[vol_var]])

# Bai-Perron test — up to 5 breaks
bp_test <- breakpoints(vol_ts ~ 1, breaks = 5)
print(summary(bp_test))

# Map breakpoint indices back to dates
bp_dates <- df %>%
  filter(!is.na(.data[[vol_var]])) %>%
  arrange(date) %>%
  slice(bp_test$breakpoints) %>%
  pull(date)

cat("\nEstimated breakpoint dates:\n")
print(bp_dates)
cat("\nYour regime dates for comparison:\n")
print(unlist(BOE_REGIME_DATES))

# Plot vol with breakpoints and regime boundaries
regime_bands <- df %>%
  group_by(regime) %>%
  summarise(xmin = min(date), xmax = max(date), .groups = "drop")

ggplot(df %>% filter(!is.na(.data[[vol_var]])),
       aes(x = date, y = .data[[vol_var]])) +
  geom_rect(
    data = regime_bands,
    aes(xmin = xmin, xmax = xmax, fill = regime),
    ymin = -Inf, ymax = Inf, alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_line(colour = "grey30", linewidth = 0.5) +
  geom_vline(xintercept = as.numeric(bp_dates),
             colour = "firebrick", linetype = "dashed", linewidth = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Vol series with Bai-Perron breakpoints (red dashed) and regime shading",
    x = NULL, y = vol_var, fill = "Regime"
  ) +
  theme_minimal(base_size = 12)

# -----------------------------------------------------------------------------
# 3. GRANGER CAUSALITY: vol -> outcomes
# Tests whether lagged vol helps predict each outcome
# Want: rejection (p < 0.05) — prior plausibility for the LP
# -----------------------------------------------------------------------------

cat("\n================================================================\n")
cat("  3. GRANGER CAUSALITY: vol -> outcomes\n")
cat("  H0: vol does not Granger-cause outcome\n")
cat("  Want: p < 0.05 for key outcomes\n")
cat("================================================================\n\n")

granger_results <- outcome_vars %>%
  map_dfr(function(v) {
    # Need complete cases for both variables
    test_df <- df %>%
      select(all_of(c(vol_var, v))) %>%
      na.omit()

    # VAR with 2 lags then Granger test
    tryCatch({
      var_fit <- VAR(test_df, p = 2, type = "const")
      gt <- causality(var_fit, cause = vol_var)$Granger
      tibble(
        outcome   = v,
        statistic = round(gt$statistic, 3),
        p_value   = round(gt$p.value, 4),
        df1       = gt$parameter[1],
        df2       = gt$parameter[2]
      )
    }, error = function(e) {
      tibble(outcome = v, statistic = NA, p_value = NA, df1 = NA, df2 = NA)
    })
  })

print(granger_results)

# -----------------------------------------------------------------------------
# 4. CROSS-CORRELATION PLOTS: vol vs outcomes
# Shows lead-lag relationship between vol and each outcome
# Useful descriptive evidence for the paper
# -----------------------------------------------------------------------------

cat("\n================================================================\n")
cat("  4. CROSS-CORRELATION PLOTS\n")
cat("================================================================\n\n")

ccf_plots <- outcome_vars %>%
  map(function(v) {
    ccf_data <- df %>%
      select(all_of(c(vol_var, v))) %>%
      na.omit()

    ccf_result <- ccf(
      ccf_data[[vol_var]],
      ccf_data[[v]],
      lag.max = 24,
      plot    = FALSE
    )

    tibble(
      lag = ccf_result$lag[,,1],
      acf = ccf_result$acf[,,1]
    ) %>%
      ggplot(aes(x = lag, y = acf)) +
      geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
      geom_hline(
        yintercept = c(-1, 1) * qnorm(0.975) / sqrt(nrow(ccf_data)),
        colour = "steelblue", linetype = "dashed", linewidth = 0.5
      ) +
      geom_segment(aes(xend = lag, yend = 0), colour = "grey40") +
      geom_point(size = 1.5, colour = "steelblue") +
      labs(
        title = glue("CCF: {vol_var} vs {v}"),
        x = "Lag (months)", y = "Correlation"
      ) +
      theme_minimal(base_size = 10)
  })

wrap_plots(ccf_plots, ncol = 3) +
  plot_annotation(
    title    = glue("Cross-correlations: {vol_var} vs outcome variables"),
    subtitle = "Dashed lines = 95% confidence bounds"
  )

cat("\n================================================================\n")
cat("  Diagnostics complete.\n")
cat("================================================================\n\n")
