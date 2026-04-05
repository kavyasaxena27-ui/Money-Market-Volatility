# =============================================================================
# 02c_pairwise_plots.R
# Pairwise scatter plots: cumulative long difference of each outcome
# vs the level of the shock variable at time t
#
# Plot: [y_{t+h} - y_{t-1}] vs x_t
# This is a visual preview of what the LP is estimating at horizon h
#
# SOURCE: After 03b and instrument module — needs df in environment
# =============================================================================

source(here::here("R", "00_setup.R"))
library(patchwork)

# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------

outcomes <- c(
  "log_gdp_sa",
  "log_core_cpi",
  "household_mortgage_rates_sa",
  "gb2yt",
  "gb5yt",
  "gb10yt"
)

outcome_labels <- c(
  log_gdp_sa                  = "Log GDP",
  log_core_cpi                = "Log Core CPI",
  household_mortgage_rates_sa = "Mortgage Rate",
  gb2yt                       = "2Y Gilt Yield",
  gb5yt                       = "5Y Gilt Yield",
  gb10yt                      = "10Y Gilt Yield"
)

shock_vars <- c(
  target       = "MP Shock (target)",
  centered_vol = "Centred Vol"
)

# Horizons to show — pick a few to see how relationship evolves
HORIZONS <- c(1, 6, 12, 24)

# -----------------------------------------------------------------------------
# HELPER: scatter of long difference vs shock at horizon h
# -----------------------------------------------------------------------------

plot_long_diff <- function(df, outcome, shock, h) {

  plot_df <- df %>%
    arrange(date) %>%
    mutate(
      long_diff = lead(.data[[outcome]], h) - lag(.data[[outcome]], 1),
      shock_t   = .data[[shock]]
    ) %>%
    select(date, long_diff, shock_t) %>%
    na.omit()

  ggplot(plot_df, aes(x = shock_t, y = long_diff)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_point(colour = "grey50", size = 0.8, alpha = 0.5) +
    geom_smooth(method = "lm", colour = "firebrick",
                se = TRUE, alpha = 0.15, linewidth = 0.8) +
    labs(
      x        = shock_vars[[shock]],
      y        = glue("y_{{t+{h}}} - y_{{t-1}}"),
      title    = glue("h = {h}"),
      subtitle = outcome_labels[[outcome]]
    ) +
    theme_minimal(base_size = 10)
}

# -----------------------------------------------------------------------------
# GENERATE: one page per outcome, columns = horizons, rows = shock vars
# -----------------------------------------------------------------------------

outcomes %>% walk(function(v) {

  # Row 1: vs target, Row 2: vs centered_vol
  plots <- names(shock_vars) %>%
    map(function(s) {
      HORIZONS %>%
        map(~plot_long_diff(df, v, s, .x))
    }) %>%
    flatten()

  p <- wrap_plots(plots, ncol = length(HORIZONS)) +
    plot_annotation(
      title    = glue("Long differences: {outcome_labels[[v]]}"),
      subtitle = glue(
        "Rows: MP shock (target), Centred vol  |  ",
        "Columns: horizons {paste(HORIZONS, collapse=', ')} months  |  ",
        "Red line = linear fit"
      )
    )

  print(p)
})
