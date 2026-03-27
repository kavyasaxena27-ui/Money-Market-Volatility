# =============================================================================
# 04_local_projections.R  [IV-updated sections only]
#
# Changes from original:
#   1. estimate_lp_state_dependent() gains iv = TRUE/FALSE argument
#   2. When iv = TRUE, feols uses IV syntax: endog ~ | inter_term ~ vol_hat_x_shock
#   3. run_analysis() passes iv flag through
#   4. spec_grid gains iv column
#
# Everything else (controls, calculate_state_irfs, plot_irf, run_analysis
# wrapper) is UNCHANGED — drop these replacements in directly.
# =============================================================================

source(here::here("R", "00_setup.R"))
library(fixest)
library(patchwork)

# Small helper — unchanged
get_controls <- function(outcome_name) {
  c(glue("lag{1:2}_{outcome_name}"), common_controls)
}

# -----------------------------------------------------------------------------
# LP estimation — IV-extended
# New argument: iv (logical, default TRUE)
#   TRUE  -> instruments inter_term with vol_hat_x_shock (Rigobon IV)
#   FALSE -> OLS as before (use for robustness comparison)
# -----------------------------------------------------------------------------

estimate_lp_state_dependent <- function(
    data,
    outcome,
    shock,
    interaction_var,
    controls,
    horizons   = 12,
    cumulative = FALSE,
    iv         = TRUE          # <-- NEW
) {
  
  work_data <- data %>%
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
      # ------------------------------------------------------------------
      # IV specification (Rigobon)
      # Endogenous: inter_term  (centered_vol × shock)
      # Instrument: vol_hat_x_shock  (vol_hat × shock)
      # Exogenous:  shock, all controls (including lagged inter_term)
      # ------------------------------------------------------------------
      rhs_exog <- glue_collapse(
        c(shock, controls),
        sep = " + "
      )
      
      fml <- as.formula(glue(
        "lhs_col ~ {rhs_exog} | inter_term ~ vol_hat_x_shock"
      ))
      
    } else {
      # ------------------------------------------------------------------
      # OLS specification (original)
      # ------------------------------------------------------------------
      rhs <- glue_collapse(
        c(shock, interaction_var, "inter_term", controls),
        sep = " + "
      )
      fml <- as.formula(glue("lhs_col ~ {rhs}"))
    }
    
    m <- fixest::feols(
      fml,
      data  = horizon_data,
      vcov  = vc,
      warn  = FALSE,
      notes = FALSE
    )
    
    b   <- coef(m)
    v   <- vcov(m)
    
    beta        <- b[[shock]]
    delta       <- b[["inter_term"]]
    se_beta     <- sqrt(v[shock, shock])
    se_delta    <- sqrt(v["inter_term", "inter_term"])
    cov_bd      <- v[shock, "inter_term"]
    
    tibble(
      horizon      = h,
      beta         = beta,
      delta        = delta,
      se_beta      = se_beta,
      se_delta     = se_delta,
      cov_beta_delta = cov_bd
    )
    
  }) %>%
    list_rbind()
}

# calculate_state_irfs — UNCHANGED (paste your original here)
calculate_state_irfs <- function(
    results,
    data,
    interaction_var,
    sd_mult = 1,
    ci      = 0.90,
    scale   = 1,
    labels  = c(mean = "Mean state", hi = "High Vol (+1SD)")
) {
  z <- qnorm((1 + ci) / 2)
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

# plot_irf — UNCHANGED (paste your original here)
plot_irf <- function(fig_data, title, subtitle, ylab = "Response") {
  cols <- c("Mean state" = "#0072B2", "High Vol (+1SD)" = "#E69F00")
  
  ggplot(fig_data, aes(x = horizon, y = est, colour = curve, fill = curve)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    theme_minimal() +
    labs(title = title, subtitle = subtitle,
         x = "Months since shock", y = ylab) +
    theme(legend.position = "bottom", legend.title = element_blank())
}

# -----------------------------------------------------------------------------
# run_analysis wrapper — gains iv argument, passes through
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
    cumulative      = FALSE,
    constant_sample = TRUE,
    iv              = TRUE     # <-- NEW
) {
  
  reg_data <- if (constant_sample) {
    data %>% filter(!is.na(lead(.data[[outcome]], horizons)))
  } else {
    data
  }
  
  res <- estimate_lp_state_dependent(
    data            = reg_data,
    outcome         = outcome,
    shock           = shock,
    interaction_var = interaction_var,
    controls        = get_controls(outcome),
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
# Manual parameters — unchanged except iv added to spec_grid
# -----------------------------------------------------------------------------

PARAM_SHOCKS         <- c("target")
PARAM_HORIZONS       <- 24
PARAM_CONSTANT_SAMPLE <- TRUE
PARAM_VOL_MEASURE    <- "garch_ronia"

common_controls <- c(
  "lag1_target", "lag2_target",
  "lag1_centered_vol", "lag2_centered_vol",
  "lag1_effective_exchange_rate",
  "lag1_log_import_price_index_rolling",
  "lag1_vol_mps", "lag2_vol_mps"
)

outcomes <- c("log_ecy2_dp_m", "core", "household_mortgage_rates",
              "gb2yt", "gb5yt", "gb10yt")

# Load data
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
  summarise(across(any_of(c(
    "sonia_20d_sd", "ronia_20d_sd",
    "sonia_garch_vol", "ronia_garch_vol"
  )), \(x) mean(x, na.rm = TRUE)))

df <- df_raw %>%
  mutate(
    across(any_of(names(vol_means)),
           \(x) x - vol_means[[cur_column()]],
           .names = "centered_{.col}"),
    centered_vol    = .data[[vol_col]],
    lag1_centered_vol = lag(centered_vol, 1),
    lag2_centered_vol = lag(centered_vol, 2),
    vol_mps         = centered_vol * target,
    lag1_vol_mps    = lag(vol_mps, 1),
    lag2_vol_mps    = lag(vol_mps, 2),
    # IV interaction — uses vol_hat from instrument module
    inter_term      = centered_vol * target
  ) %>%
  filter(date >= PARAM_START_DATE) %>%
  mutate(id = 1) %>%
  panel(panel.id = ~ id + date)

# -----------------------------------------------------------------------------
# Spec grid — now includes iv = TRUE (baseline) and iv = FALSE (robustness)
# -----------------------------------------------------------------------------

spec_grid <- crossing(
  outcome     = outcomes,
  shock       = PARAM_SHOCKS,
  interaction = "centered_vol",
  iv          = c(TRUE, FALSE)    # <-- run both; compare as robustness check
) %>%
  mutate(cumulative = outcome == "log_ecy2_dp_m")

use_grid <- spec_grid %>%
  mutate(
    results = pmap(
      list(outcome, shock, interaction, cumulative, iv),
      function(out, shk, int, cum, use_iv) {
        run_analysis(
          data            = df,
          outcome         = out,
          shock           = shk,
          interaction_var = int,
          horizons        = PARAM_HORIZONS,
          cumulative      = cum,
          constant_sample = PARAM_CONSTANT_SAMPLE,
          iv              = use_iv
        )
      },
      .progress = "Running LPs..."
    )
  )

use_grid$results %>%
  map("plot") %>%
  wrap_plots(ncol = 3)
