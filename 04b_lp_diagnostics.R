# =============================================================================
# 04b_lp_diagnostics.R
# Two things:
#
# 1. inspect_lp_spec() — prints exact formula, data snapshot, summary stats
#    and sample info for a given outcome and horizon. Run BEFORE spec_grid
#    to verify the specification is correct.
#
# 2. get_controls() — corrected version that uses DIFFERENCED lags of the
#    outcome variable on the RHS, consistent with the cumulative LHS.
#    Other controls remain in levels.
#
# LITERATURE:
#   Jordà (2005) AER — baseline LP specification
#   Ramey (2016) Handbook of Macroeconomics — LP specification choices
#   Plagborg-Møller & Wolf (2021) Econometrica — stationarity of LP controls;
#     shows lagged dependent variable controls should match LHS transformation
#   Stock & Watson (2018) — I(1) variables in LP
#
# SOURCE: at top of 04_local_projections.R, after 04_instrument_module.R
# =============================================================================

# -----------------------------------------------------------------------------
# get_controls() — CORRECTED
#
# For cumulative LP (LHS = y_{t+h} - y_{t-1}):
#   - Lags of outcome enter as DIFFERENCES: diff_y_{t-1}, diff_y_{t-2}
#     (Plagborg-Møller & Wolf 2021 — controls must be stationary)
#   - All other controls remain in LEVELS
#     (standard in LP literature — absorbing confounders not structural)
#
# For non-cumulative LP (LHS = y_{t+h}, stationary variable):
#   - Lags of outcome enter in LEVELS as normal
#
# INPUTS:
#   outcome_name — string, e.g. "log_gdp_sa"
#   cumulative   — logical, TRUE if LHS is cumulative change
#
# The differenced lag columns must exist in df, named:
#   diff_{outcome}_lag1, diff_{outcome}_lag2
# These are created in the mutate block below (add_outcome_lags())
# -----------------------------------------------------------------------------

get_controls <- function(outcome_name, cumulative = TRUE) {
  if (cumulative) {
    # Differenced lags of outcome (stationary, matches cumulative LHS)
    outcome_lags <- c(
      glue("diff_{outcome_name}_lag1"),
      glue("diff_{outcome_name}_lag2")
    )
  } else {
    # Level lags (outcome already stationary)
    outcome_lags <- glue("lag{1:2}_{outcome_name}")
  }
  c(outcome_lags, common_controls)
}

# -----------------------------------------------------------------------------
# add_outcome_lags()
#
# Creates differenced lag columns for each outcome in df.
# Call this once after df is built, before the LP loop.
#
# Adds columns: diff_{outcome}_lag1, diff_{outcome}_lag2
# for each outcome in the outcomes vector
# -----------------------------------------------------------------------------

add_outcome_lags <- function(df, outcomes) {
  df <- as.data.frame(df)
  for (v in outcomes) {
    df <- df %>%
      mutate(
        !!glue("diff_{v}_lag1") := lag(.data[[v]], 1) - lag(.data[[v]], 2),
        !!glue("diff_{v}_lag2") := lag(.data[[v]], 2) - lag(.data[[v]], 3)
      )
  }
  message("checkmark: Differenced outcome lags added for: ",
          paste(outcomes, collapse = ", "))
  return(df)
}

# -----------------------------------------------------------------------------
# inspect_lp_spec()
#
# Prints full specification and data snapshot for one outcome at one horizon.
# Run before spec_grid to verify everything looks right.
#
# INPUTS:
#   df          — data frame (will be coerced from panel if needed)
#   outcome     — string, outcome variable name
#   shock       — string, shock variable name (default "target")
#   h           — integer, horizon to inspect (default 0)
#   cumulative  — logical (default TRUE)
#   n_rows      — number of rows to print in snapshot (default 6)
# -----------------------------------------------------------------------------

inspect_lp_spec <- function(df, outcome, shock = "target",
                             h = 0, cumulative = TRUE, n_rows = 6) {

  df <- as.data.frame(df) %>% arrange(date)

  controls   <- get_controls(outcome, cumulative)
  inter_name <- "inter_term"

  # Build LHS
  work_data <- df %>%
    mutate(
      inter_term = centered_vol * .data[[shock]],
      lhs_col    = if (cumulative) {
        lead(.data[[outcome]], h) - lag(.data[[outcome]], 1)
      } else {
        lead(.data[[outcome]], h)
      }
    )

  # OLS formula (IV structure same but with | endog ~ instr)
  rhs_ols <- paste(
    c(shock, "centered_vol", inter_name, controls),
    collapse = " + "
  )
  fml_ols <- glue("lhs_col ~ {rhs_ols}")

  rhs_iv_exog <- paste(c(shock, controls), collapse = " + ")
  fml_iv  <- glue(
    "lhs_col ~ {rhs_iv_exog} | ",
    "centered_vol + inter_term ~ ",
    "iv1_dummy_r1 + iv1_dummy_r2 + iv1_dummy_r3 + iv1_dummy_r4 + ",
    "iv2_dummy_r1 + iv2_dummy_r2 + iv2_dummy_r3 + iv2_dummy_r4"
  )

  cat("\n================================================================\n")
  cat(glue("  LP SPECIFICATION INSPECTOR\n"))
  cat(glue("  Outcome: {outcome}  |  Horizon: h = {h}  |  ",
           "Cumulative: {cumulative}\n"))
  cat("================================================================\n\n")

  # 1. Formula
  cat("-- 1. OLS formula ------------------------------------------\n")
  cat(fml_ols, "\n\n")
  cat("-- 1b. IV formula ------------------------------------------\n")
  cat(fml_iv, "\n\n")

  # 2. LHS description
  cat("-- 2. LHS variable -----------------------------------------\n")
  if (cumulative) {
    cat(glue("  {outcome}_{{t+{h}}} - {outcome}_{{t-1}}\n"))
    cat("  (cumulative change — stationary per Plagborg-Møller & Wolf 2021)\n\n")
  } else {
    cat(glue("  {outcome}_{{t+{h}}}  (level — variable already stationary)\n\n"))
  }

  # 3. RHS controls breakdown
  cat("-- 3. RHS controls -----------------------------------------\n")
  cat("  Shock:              ", shock, "\n")
  cat("  Vol (endogenous):   centered_vol\n")
  cat("  Interaction:        inter_term = centered_vol *", shock, "\n")
  cat("  Outcome lags:      ", paste(get_controls(outcome, cumulative)[1:2],
                                     collapse = ", "), "\n")
  cat("  (differenced lags — stationary, Plagborg-Møller & Wolf 2021)\n")
  cat("  Other controls:    ", paste(common_controls, collapse = "\n               "),
      "\n")
  cat("  (levels — standard LP practice, Ramey 2016)\n\n")

  # 4. Data snapshot
  cat("-- 4. Data snapshot (first", n_rows, "rows) ----------------------\n")
  snap_cols <- c("date", "lhs_col", shock, "centered_vol", "inter_term",
                 get_controls(outcome, cumulative)[1:2])
  snap_cols <- snap_cols[snap_cols %in% names(work_data)]

  work_data %>%
    dplyr::select(all_of(snap_cols)) %>%
    head(n_rows) %>%
    print()
  cat("\n")

  # 5. Summary stats
  cat("-- 5. Summary stats of key variables -----------------------\n")
  work_data %>%
    dplyr::select(lhs_col, all_of(shock), centered_vol, inter_term) %>%
    summary() %>%
    print()
  cat("\n")

  # 6. Sample info
  complete_rows <- work_data %>%
    dplyr::select(lhs_col, all_of(shock), centered_vol,
                  all_of(get_controls(outcome, cumulative))) %>%
    complete.cases() %>%
    sum()

  date_range <- work_data %>%
    filter(!is.na(lhs_col)) %>%
    summarise(start = min(date), end = max(date))

  cat("-- 6. Sample info ------------------------------------------\n")
  cat(sprintf("  Complete cases:  %d / %d rows\n",
              complete_rows, nrow(work_data)))
  cat(sprintf("  Date range:      %s to %s\n\n",
              date_range$start, date_range$end))

  # 7. Check all control columns exist
  cat("-- 7. Control column check ---------------------------------\n")
  missing_cols <- controls[!controls %in% names(df)]
  if (length(missing_cols) == 0) {
    cat("  checkmark: All control columns found in df\n\n")
  } else {
    cat("  WARNING: Missing columns:\n")
    cat(" ", paste(missing_cols, collapse = "\n  "), "\n\n")
  }

  cat("================================================================\n\n")

  invisible(work_data)
}

# =============================================================================
# USAGE — add these lines in 04_local_projections.R after df is built:
#
#   source(here::here("R", "04b_lp_diagnostics.R"))
#
#   # Add differenced outcome lags to df
#   df <- add_outcome_lags(df, outcomes)
#
#   # Inspect specification before running
#   inspect_lp_spec(df, outcome = outcomes[1], h = 0,  cumulative = TRUE)
#   inspect_lp_spec(df, outcome = outcomes[1], h = 12, cumulative = TRUE)
# =============================================================================
