# =============================================================================
# 04_local_projections.R
#
# State-dependent local projections with Rigobon/Lewbel IV
# Follows Holm-Hadulla & Pool (ECB WP 3048), adapted for UK / BoE
#
# IV NOTE: Currently 2SLS via fixest::feols
#          Switch to GMM later for final paper version
#
# LITERATURE:
#   Jordà (2005) AER — baseline LP specification
#   Ramey (2016) Handbook of Macroeconomics — LP specification choices
#   Plagborg-Møller & Wolf (2021) Econometrica — stationarity of LP controls
#   Rigobon (2003) RES — heteroskedasticity-based IV
#   Lewbel (2012) JBES — instrument construction
#
# =============================================================================
#
# TO-DO LIST:
#
#   🔴 1. FIX PATH VARIABLE
#         lag1_path / lag2_path has incomplete coverage -> 78 missing obs
#         Options:
#           a) Replace with BoE balance sheet size (fuller coverage)
#           b) Use path shocks only where available and trim sample
#           c) Drop from fs_controls_iv if no good alternative
#         Check coverage: df %>% summarise(n_missing = sum(is.na(path)))
#
#   🔴 2. RUN IV SPEC
#         Once path is fixed, flip spec_grid iv to c(FALSE, TRUE)
#         Check: Sargan test should improve with clean data
#         Check: BP tests still reject on both first stages
#
#   🟡 3. COMPARE OLS VS IV IRFs SIDE BY SIDE
#         Key question: does dampening effect on GDP sharpen under IV?
#         Expected: yes — OLS attenuates delta downward
#         Plot both on same axes for each outcome
#         Code hint:
#           ols_grid <- use_grid %>% filter(!iv)
#           iv_grid  <- use_grid %>% filter(iv)
#           # then combine fig_data and plot with linetype = iv
#
#   🟡 4. SWITCH TO GMM
#         Replace feols 2SLS with ivreg/estimatr GMM
#         Fixes: vcov non-PSD warnings, closer to paper's exact estimator
#         Add Newey-West bandwidth explicitly (paper uses NW SEs)
#
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here::here("R", "04_instrument_module.R"))
source(here::here("R", "04b_lp_diagnostics.R"))   # corrected get_controls,
                                                    # add_outcome_lags,
                                                    # inspect_lp_spec
library(fixest)
library(patchwork)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

PARAM_SHOCKS          <- c("target")
PARAM_HORIZONS        <- 24
PARAM_CONSTANT_SAMPLE <- TRUE
PARAM_VOL_MEASURE     <- "rolling_sonia"

outcomes <- c(
  "log_gdp_sa",
  "log_core_cpi",
  "household_mortgage_rates_sa",
  "gb2yt", "gb5yt", "gb10yt"
)

# -----------------------------------------------------------------------------
# Load data and build df
# -----------------------------------------------------------------------------

df_raw <- readr::read_rds(
  path(project_paths[["data_processed"]], "monthly_all_series_WIDE.rds")
) %>%
  arrange(date)

vol_col <- switch(
  PARAM_VOL_MEASURE,
  rolling_sonia = "centered_sonia_20d_sd",
  garch_sonia   = "centered_sonia_garch_vol",
  rolling_ronia = "centered_ronia_20d_sd",
  garch_ronia   = "centered_ronia_garch_vol"
)

vol_means <- df_raw %>%
  filter(date >= PARAM_START_DATE) %>%
  summarise(across(
    any_of(c("sonia_20d_sd", "ronia_20d_sd",
             "sonia_garch_vol", "ronia_garch_vol")),
    \(x) mean(x, na.rm = TRUE)
  ))

df <- df_raw %>%
  mutate(
    across(
      any_of(names(vol_means)),
      \(x) x - vol_means[[cur_column()]],
      .names = "centered_{.col}"
    ),
    centered_vol      = .data[[vol_col]],
    lag1_centered_vol = lag(centered_vol, 1),
    lag2_centered_vol = lag(centered_vol, 2),
    vol_mps           = centered_vol * target,
    lag1_vol_mps      = lag(vol_mps, 1),
    lag2_vol_mps      = lag(vol_mps, 2),
    lag1_target       = lag(target, 1),
    lag2_target       = lag(target, 2),
    id                = 1L
  ) %>%
  filter(date >= PARAM_START_DATE)

# -----------------------------------------------------------------------------
# Controls
# NOTE: get_controls() comes from 04b_lp_diagnostics.R
#       It uses differenced lags of outcome on RHS (Plagborg-Møller & Wolf 2021)
#       common_controls must be defined here before get_controls() is called
# -----------------------------------------------------------------------------

# Instrument names — 8 total (4 per endogenous variable)
iv_names <- c(
  paste0("iv1_dummy_r", 1:4),
  paste0("iv2_dummy_r", 1:4)
)

# First-stage controls X_{t-p} — regime dummies added automatically inside
# make_instrument()
# 🔴 TODO: replace lag1_path/lag2_path once coverage issue is fixed
fs_controls_iv <- c(
  "lag1_target", "lag2_target",
  "lag1_centered_vol", "lag2_centered_vol",
  "lag1_vol_mps", "lag2_vol_mps",
  "lag1_effective_exchange_rate",
  "lag1_log_import_price_index_rolling",
  "lag1_path", "lag2_path"
)

# LP controls — common across all outcomes
# Outcome-specific differenced lags added by get_controls() from 04b
common_controls <- c(
  "lag1_target", "lag2_target",
  "lag1_centered_vol", "lag2_centered_vol",
  "lag1_effective_exchange_rate",
  "lag1_log_import_price_index_rolling",
  "lag1_vol_mps", "lag2_vol_mps",
  "dummy_r1", "dummy_r2", "dummy_r3", "dummy_r4"
)

# -----------------------------------------------------------------------------
# Build instruments
# -----------------------------------------------------------------------------

df <- df %>%
  make_instrument(
    vol_col     = vol_col,
    shock_col   = "target",
    fs_controls = fs_controls_iv
  )

check_instrument(df)

# Add differenced outcome lags (needed by get_controls cumulative = TRUE)
df <- add_outcome_lags(df, outcomes)

# -----------------------------------------------------------------------------
# Inspect specification before running — change outcome/h as needed
# -----------------------------------------------------------------------------

inspect_lp_spec(df, outcome = outcomes[1], h = 0,  cumulative = TRUE)
inspect_lp_spec(df, outcome = outcomes[1], h = 12, cumulative = TRUE)

# -----------------------------------------------------------------------------
# LP estimation
# iv = TRUE  -> 2SLS (paper baseline)
# iv = FALSE -> OLS (robustness / comparison)
# -----------------------------------------------------------------------------

estimate_lp_state_dependent <- function(
    data,
    outcome,
    shock,
    interaction_var,
    controls,
    horizons   = 12,
    cumulative = TRUE,
    iv         = FALSE
) {

  work_data <- as.data.frame(data) %>%
    arrange(date) %>%
    mutate(inter_term = .data[[shock]] * .data[[interaction_var]])

  map(0:horizons, function(h) {

    horizon_data <- work_data %>%
      mutate(
        lhs_col = if (cumulative) {
          lead(.data[[outcome]], h) - lag(.data[[outcome]], 1)
        } else {
          lead(.data[[outcome]], h)
        }
      )

    vc <- if (h == 0) "hetero" else fixest::NW(h + 1)

    if (iv) {
      # 2SLS — two endogenous regressors, 8 instruments
      # Endogenous: interaction_var (sigma), inter_term (sigma*S)
      # Instruments: iv1_dummy_r1:r4, iv2_dummy_r1:r4
      exog  <- glue_collapse(c(shock, controls), sep = " + ")
      endog <- paste(interaction_var, "inter_term", sep = " + ")
      instr <- glue_collapse(iv_names, sep = " + ")

      fml <- as.formula(glue(
        "lhs_col ~ {exog} | {endog} ~ {instr}"
      ))

    } else {
      # OLS
      rhs <- glue_collapse(
        c(shock, interaction_var, "inter_term", controls),
        sep = " + "
      )
      fml <- as.formula(glue("lhs_col ~ {rhs}"))
    }

    m <- fixest::feols(
      fml,
      data     = horizon_data,
      panel.id = ~id + date,
      vcov     = vc,
      warn     = FALSE,
      notes    = FALSE
    )

    b <- coef(m)
    v <- vcov(m)

    # Coefficient names differ between OLS and IV in fixest
    inter_name <- if (iv) "fit_inter_term"   else "inter_term"
    vol_name   <- if (iv) glue("fit_{interaction_var}") else interaction_var

    beta     <- b[[shock]]
    delta    <- b[[inter_name]]
    se_beta  <- sqrt(v[shock, shock])
    se_delta <- sqrt(v[inter_name, inter_name])
    cov_bd   <- v[shock, inter_name]

    tibble(
      horizon        = h,
      beta           = beta,
      delta          = delta,
      se_beta        = se_beta,
      se_delta       = se_delta,
      cov_beta_delta = cov_bd
    )

  }) %>%
    list_rbind()
}

# -----------------------------------------------------------------------------
# IRF calculation
# -----------------------------------------------------------------------------

calculate_state_irfs <- function(
    results,
    data,
    interaction_var,
    sd_mult = 1,
    ci      = 0.90,
    scale   = 1,
    labels  = c(mean = "Mean state", hi = "High Vol (+1SD)")
) {
  z      <- qnorm((1 + ci) / 2)
  mean_x <- mean(data[[interaction_var]], na.rm = TRUE)
  sd_x   <- sd(data[[interaction_var]],  na.rm = TRUE)

  eval_pts <- tibble(
    curve = unname(labels),
    x_val = c(mean_x, mean_x + sd_mult * sd_x)
  )

  results %>%
    crossing(eval_pts) %>%
    mutate(
      est = beta + delta * x_val,
      se  = sqrt(se_beta^2 + x_val^2 * se_delta^2 + 2 * x_val * cov_beta_delta),
      lo  = est - z * se,
      hi  = est + z * se,
      across(c(est, lo, hi), \(v) v * scale)
    ) %>%
    select(-x_val, -se)
}

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

plot_irf <- function(fig_data, title, subtitle, ylab = "Response") {
  cols <- c("Mean state" = "#0072B2", "High Vol (+1SD)" = "#E69F00")

  ggplot(fig_data, aes(x = horizon, y = est, colour = curve, fill = curve)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    theme_minimal() +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Months since shock",
      y        = ylab
    ) +
    theme(legend.position = "bottom", legend.title = element_blank())
}

# -----------------------------------------------------------------------------
# Wrapper
# -----------------------------------------------------------------------------

run_analysis <- function(
    data,
    outcome,
    shock,
    interaction_var,
    horizons        = 24,
    ci              = 0.90,
    sd_mult         = 1,
    scale           = 1,
    cumulative      = TRUE,
    constant_sample = TRUE,
    iv              = FALSE
) {

  reg_data <- if (constant_sample) {
    as.data.frame(data) %>%
      filter(!is.na(lead(.data[[outcome]], horizons)))
  } else {
    as.data.frame(data)
  }

  res <- estimate_lp_state_dependent(
    data            = reg_data,
    outcome         = outcome,
    shock           = shock,
    interaction_var = interaction_var,
    controls        = get_controls(outcome, cumulative),  # corrected version
    horizons        = horizons,
    cumulative      = cumulative,
    iv              = iv
  )

  fig_data <- calculate_state_irfs(
    results         = res,
    data            = reg_data,
    interaction_var = interaction_var,
    sd_mult         = sd_mult,
    ci              = ci,
    scale           = scale
  )

  plot <- plot_irf(
    fig_data = fig_data,
    title    = glue("{outcome} response to {shock}"),
    subtitle = glue("State: {interaction_var} (+{sd_mult} SD)  |  IV: {iv}")
  )

  list(results = res, fig_data = fig_data, plot = plot)
}

# -----------------------------------------------------------------------------
# Spec grid
# All I(1) outcomes use cumulative = TRUE (Plagborg-Møller & Wolf 2021)
# iv = FALSE for now — flip to c(FALSE, TRUE) once path issue is fixed (TODO 1)
# -----------------------------------------------------------------------------

spec_grid <- crossing(
  outcome     = outcomes,
  shock       = PARAM_SHOCKS,
  interaction = "centered_vol",
  iv          = FALSE        # 🔴 TODO 2: change to c(FALSE, TRUE) after fixing path
) %>%
  mutate(cumulative = TRUE)  # all outcomes are I(1) — use cumulative LP

use_grid <- spec_grid %>%
  mutate(
    results = pmap(
      list(outcome, shock, interaction, cumulative, iv),
      function(out, shk, int, cum, iv) {
        run_analysis(
          data            = df,
          outcome         = out,
          shock           = shk,
          interaction_var = int,
          horizons        = PARAM_HORIZONS,
          cumulative      = cum,
          constant_sample = PARAM_CONSTANT_SAMPLE,
          iv              = iv
        )
      },
      .progress = "Running LPs..."
    )
  )

use_grid$results %>%
  map("plot") %>%
  wrap_plots(ncol = 3)

# -----------------------------------------------------------------------------
# TODO 3: OLS vs IV comparison plot (run after flipping iv above)
# -----------------------------------------------------------------------------
# ols_grid <- use_grid %>% filter(!iv)
# iv_grid  <- use_grid %>% filter(iv)
# -- combine fig_data with iv label and plot with linetype = iv
