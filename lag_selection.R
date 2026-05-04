# =============================================================================
# 04c_lag_selection.R
# Lag selection diagnostics for Local Projections
#
# Methods:
#   1. VAR lag selection (AIC, BIC, HQ, FPE) — joint system
#   2. LP-specific lag selection per outcome — AIC/BIC at h=0
#   3. ACF of residuals — visual check per outcome
#
# USAGE:
#   source(here::here("R", "04c_lag_selection.R"))
#   lag_results <- run_lag_selection(df, outcomes, common_controls)
#   print(lag_results$summary)
#
# SOURCE: After df is built and instruments constructed in 04_local_projections.R
# =============================================================================

source(here::here("R", "00_setup.R"))
library(vars)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. VAR lag selection — joint system
# -----------------------------------------------------------------------------

var_lag_selection <- function(df, outcomes, shock = "target",
                               vol = "centered_vol", max_lags = 12) {

  cat("\n================================================================\n")
  cat("  1. VAR LAG SELECTION\n")
  cat("  AIC tends to over-select, BIC tends to under-select\n")
  cat("  Use BIC as baseline, AIC as robustness\n")
  cat("================================================================\n\n")

  var_vars <- as.data.frame(df) %>%
    dplyr::select(all_of(c(outcomes, shock, vol))) %>%
    na.omit()

  sel <- VARselect(var_vars, lag.max = max_lags, type = "const")

  cat("Optimal lags by criterion:\n")
  print(sel$selection)
  cat("\nCriteria values across lags:\n")
  print(round(sel$criteria, 2))

  return(sel)
}

# -----------------------------------------------------------------------------
# 2. LP-specific lag selection per outcome
# Runs LP at h=0 with 1:max_lags outcome lags, compares AIC/BIC
# -----------------------------------------------------------------------------

lp_lag_select_outcome <- function(df, outcome, common_controls,
                                   max_lags = 12, cumulative = TRUE) {

  work_data <- as.data.frame(df) %>%
    arrange(date) %>%
    mutate(
      inter_term = centered_vol * target,
      lhs_col    = if (cumulative) {
        lead(.data[[outcome]], 0) - lag(.data[[outcome]], 1)
      } else {
        lead(.data[[outcome]], 0)
      }
    )

  map_dfr(1:max_lags, function(l) {

    # Outcome lags — differenced if cumulative, levels if not
    if (cumulative) {
      outcome_lags <- glue("diff_{outcome}_lag{1:l}")
    } else {
      outcome_lags <- glue("lag{1:l}_{outcome}")
    }

    controls <- c(
      "target", "centered_vol", "inter_term",
      outcome_lags,
      common_controls
    )

    # Keep only columns that exist in data
    controls <- controls[controls %in% names(work_data)]

    rhs <- paste(controls, collapse = " + ")
    fml <- as.formula(glue("lhs_col ~ {rhs}"))

    tryCatch({
      m <- lm(fml, data = work_data, na.action = na.omit)
      tibble(
        outcome   = outcome,
        n_lags    = l,
        aic       = AIC(m),
        bic       = BIC(m),
        adj_r2    = round(summary(m)$adj.r.squared, 4),
        n_obs     = nobs(m),
        n_params  = length(coef(m))
      )
    }, error = function(e) {
      tibble(outcome = outcome, n_lags = l,
             aic = NA, bic = NA, adj_r2 = NA,
             n_obs = NA, n_params = NA)
    })
  })
}

# -----------------------------------------------------------------------------
# 3. ACF of LP residuals per outcome and lag choice
# -----------------------------------------------------------------------------

acf_lp_residuals <- function(df, outcome, n_lags, common_controls,
                              cumulative = TRUE, max_acf_lags = 24) {

  work_data <- as.data.frame(df) %>%
    arrange(date) %>%
    mutate(
      inter_term = centered_vol * target,
      lhs_col    = if (cumulative) {
        lead(.data[[outcome]], 0) - lag(.data[[outcome]], 1)
      } else {
        lead(.data[[outcome]], 0)
      }
    )

  if (cumulative) {
    outcome_lags <- glue("diff_{outcome}_lag{1:n_lags}")
  } else {
    outcome_lags <- glue("lag{1:n_lags}_{outcome}")
  }

  controls <- c("target", "centered_vol", "inter_term",
                outcome_lags, common_controls)
  controls <- controls[controls %in% names(work_data)]

  rhs <- paste(controls, collapse = " + ")
  fml <- as.formula(glue("lhs_col ~ {rhs}"))

  m   <- lm(fml, data = work_data, na.action = na.omit)
  res <- residuals(m)

  acf_data <- acf(res, lag.max = max_acf_lags, plot = FALSE)
  ci       <- qnorm(0.975) / sqrt(length(res))

  tibble(
    lag = acf_data$lag[,,1],
    acf = acf_data$acf[,,1]
  ) %>%
    ggplot(aes(x = lag, y = acf)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    geom_hline(yintercept = c(-ci, ci),
               colour = "steelblue", linetype = "dashed") +
    geom_segment(aes(xend = lag, yend = 0), colour = "grey40") +
    geom_point(size = 1.5, colour = "steelblue") +
    labs(
      title    = glue("{outcome} — {n_lags} lags"),
      subtitle = "Dashed = 95% CI  |  Spikes outside = remaining autocorrelation",
      x = "Lag", y = "ACF"
    ) +
    theme_minimal(base_size = 10)
}

# -----------------------------------------------------------------------------
# MAIN: run_lag_selection()
# Loops over all outcomes, returns summary table and plots
#
# INPUTS:
#   df              — data frame
#   outcomes        — character vector of outcome names
#   common_controls — character vector of common controls
#   cumulative_outcomes — which outcomes use cumulative spec
#   max_lags        — maximum lags to consider (default 12)
# -----------------------------------------------------------------------------

run_lag_selection <- function(df, outcomes, common_controls,
                               cumulative_outcomes = NULL,
                               max_lags = 12) {

  if (is.null(cumulative_outcomes)) {
    # Default: log level outcomes are cumulative, rates are not
    cumulative_outcomes <- outcomes[grepl("log_", outcomes)]
  }

  cat("\n================================================================\n")
  cat("  LP LAG SELECTION DIAGNOSTICS\n")
  cat(glue("  Outcomes: {paste(outcomes, collapse = ', ')}\n"))
  cat(glue("  Max lags: {max_lags}\n"))
  cat("================================================================\n\n")

  # -- Step 1: VAR lag selection --------------------------------------------
  var_sel <- var_lag_selection(df, outcomes,
                                max_lags = max_lags)

  # -- Step 2: LP lag selection per outcome ---------------------------------
  cat("\n-- LP-specific lag selection (AIC/BIC at h=0) ------------------\n\n")

  lp_results <- outcomes %>%
    map_dfr(function(v) {
      cumulative <- v %in% cumulative_outcomes
      lp_lag_select_outcome(df, v, common_controls,
                             max_lags   = max_lags,
                             cumulative = cumulative)
    })

  # Summary: optimal lags per outcome
  summary_tbl <- lp_results %>%
    group_by(outcome) %>%
    summarise(
      aic_opt    = n_lags[which.min(aic)],
      bic_opt    = n_lags[which.min(bic)],
      adj_r2_at_bic = adj_r2[which.min(bic)],
      n_obs      = max(n_obs, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      recommendation = pmax(aic_opt, bic_opt),  # conservative: take max
      cumulative     = outcome %in% cumulative_outcomes
    )

  cat("Optimal lags per outcome:\n\n")
  print(summary_tbl)

  # -- Step 3: AIC/BIC plots per outcome ------------------------------------
  cat("\n-- Generating AIC/BIC plots... ----------------------------------\n")

  aic_bic_plots <- outcomes %>%
    map(function(v) {
      lp_results %>%
        filter(outcome == v) %>%
        pivot_longer(cols = c(aic, bic),
                     names_to = "criterion",
                     values_to = "value") %>%
        ggplot(aes(x = n_lags, y = value, colour = criterion)) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 1.5) +
        scale_colour_manual(values = c(aic = "#0072B2", bic = "#E69F00")) +
        labs(title = v, x = "Number of lags", y = "Information criterion",
             colour = NULL) +
        theme_minimal(base_size = 10) +
        theme(legend.position = "bottom")
    })

  wrap_plots(aic_bic_plots, ncol = 3) +
    plot_annotation(title = "AIC/BIC by lag length per outcome",
                    subtitle = "Blue = AIC, Orange = BIC | Lower is better") %>%
    print()

  # -- Step 4: ACF plots at recommended lag length --------------------------
  cat("\n-- Generating ACF plots at recommended lags... ------------------\n")

  acf_plots <- outcomes %>%
    map(function(v) {
      cumulative  <- v %in% cumulative_outcomes
      rec_lags    <- summary_tbl$bic_opt[summary_tbl$outcome == v]
      acf_lp_residuals(df, v, rec_lags, common_controls,
                        cumulative = cumulative)
    })

  wrap_plots(acf_plots, ncol = 3) +
    plot_annotation(
      title    = "Residual ACF at BIC-optimal lag length",
      subtitle = "Spikes outside dashed lines = remaining autocorrelation → add more lags"
    ) %>%
    print()

  cat("\n================================================================\n")
  cat("  RECOMMENDATION SUMMARY\n")
  cat("  Use 'recommendation' column as baseline lag length\n")
  cat("  Run robustness with ±2 lags\n")
  cat("================================================================\n\n")

  print(summary_tbl %>%
          dplyr::select(outcome, bic_opt, aic_opt,
                        recommendation, cumulative))

  invisible(list(
    var_selection = var_sel,
    lp_results    = lp_results,
    summary       = summary_tbl
  ))
}
