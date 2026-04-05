# =============================================================================
# 02c_pairwise_plots.R
# Pairwise plots: each outcome vs target and vs centered_vol
#
# Two plot types per outcome:
#   1. Time series overlay (dual axis)
#   2. Scatter with smoother
#
# SOURCE: After 03b and instrument module — needs df in environment
# =============================================================================

source(here::here("R", "00_setup.R"))
library(patchwork)
library(scales)

# -----------------------------------------------------------------------------
# CONFIG — update to match your actual column names
# -----------------------------------------------------------------------------

outcomes <- c(
  "log_gdp_sa",
  "log_core_cpi",
  "household_mortgage_rates_sa",
  "gb2yt",
  "gb5yt",
  "gb10yt"
)

shock_vars <- c(
  "target",       # MP shock
  "centered_vol"  # vol measure
)

outcome_labels <- c(
  log_gdp_sa                   = "Log GDP",
  log_core_cpi                 = "Log Core CPI",
  household_mortgage_rates_sa  = "Mortgage Rate",
  gb2yt                        = "2Y Gilt Yield",
  gb5yt                        = "5Y Gilt Yield",
  gb10yt                       = "10Y Gilt Yield"
)

shock_labels <- c(
  target       = "MP Shock (target)",
  centered_vol = "Centred Vol"
)

# -----------------------------------------------------------------------------
# HELPER: dual-axis time series overlay
# Rescales shock to outcome scale for visual comparison
# -----------------------------------------------------------------------------

plot_ts_overlay <- function(df, outcome, shock) {

  plot_df <- df %>%
    select(date, y = all_of(outcome), x = all_of(shock)) %>%
    na.omit()

  # Rescale shock to outcome range for dual-axis
  y_range <- range(plot_df$y, na.rm = TRUE)
  x_range <- range(plot_df$x, na.rm = TRUE)

  scale_factor <- diff(y_range) / diff(x_range)
  offset       <- y_range[1] - x_range[1] * scale_factor

  plot_df <- plot_df %>%
    mutate(x_scaled = x * scale_factor + offset)

  ggplot(plot_df, aes(x = date)) +
    geom_line(aes(y = y), colour = "#0072B2", linewidth = 0.6) +
    geom_line(aes(y = x_scaled), colour = "#E69F00",
              linewidth = 0.5, linetype = "dashed") +
    scale_y_continuous(
      name     = outcome_labels[[outcome]],
      sec.axis = sec_axis(
        ~ (. - offset) / scale_factor,
        name = shock_labels[[shock]]
      )
    ) +
    labs(x = NULL,
         title    = glue("{outcome_labels[[outcome]]} vs {shock_labels[[shock]]}"),
         subtitle = "Blue = outcome (left axis)  |  Orange dashed = shock (right axis)") +
    theme_minimal(base_size = 10) +
    theme(
      axis.title.y       = element_text(colour = "#0072B2"),
      axis.title.y.right = element_text(colour = "#E69F00")
    )
}

# -----------------------------------------------------------------------------
# HELPER: scatter with loess smoother
# -----------------------------------------------------------------------------

plot_scatter <- function(df, outcome, shock) {

  plot_df <- df %>%
    select(y = all_of(outcome), x = all_of(shock)) %>%
    na.omit()

  ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(colour = "grey60", size = 0.8, alpha = 0.6) +
    geom_smooth(method = "loess", colour = "#0072B2",
                fill = "#0072B2", alpha = 0.15, linewidth = 0.8) +
    geom_smooth(method = "lm", colour = "firebrick",
                linetype = "dashed", se = FALSE, linewidth = 0.6) +
    labs(
      x        = shock_labels[[shock]],
      y        = outcome_labels[[outcome]],
      title    = glue("{outcome_labels[[outcome]]} vs {shock_labels[[shock]]}"),
      subtitle = "Blue = loess  |  Red dashed = linear fit"
    ) +
    theme_minimal(base_size = 10)
}

# -----------------------------------------------------------------------------
# GENERATE PLOTS
# -----------------------------------------------------------------------------

# -- vs target ----------------------------------------------------------------
cat("Generating plots vs target...\n")

ts_vs_target <- outcomes %>%
  map(~plot_ts_overlay(df, .x, "target"))

scatter_vs_target <- outcomes %>%
  map(~plot_scatter(df, .x, "target"))

wrap_plots(ts_vs_target, ncol = 2) +
  plot_annotation(title = "Time series: outcomes vs MP shock (target)")

wrap_plots(scatter_vs_target, ncol = 2) +
  plot_annotation(title = "Scatter: outcomes vs MP shock (target)")

# -- vs centered_vol ----------------------------------------------------------
cat("Generating plots vs centered_vol...\n")

ts_vs_vol <- outcomes %>%
  map(~plot_ts_overlay(df, .x, "centered_vol"))

scatter_vs_vol <- outcomes %>%
  map(~plot_scatter(df, .x, "centered_vol"))

wrap_plots(ts_vs_vol, ncol = 2) +
  plot_annotation(title = "Time series: outcomes vs centred vol")

wrap_plots(scatter_vs_vol, ncol = 2) +
  plot_annotation(title = "Scatter: outcomes vs centred vol")

# -- Combined grid: one row per outcome, ts + scatter side by side ------------
cat("Generating combined grid...\n")

combined_plots <- outcomes %>%
  map(function(v) {
    ts_t  <- plot_ts_overlay(df, v, "target")
    sc_t  <- plot_scatter(df, v, "target")
    ts_vol <- plot_ts_overlay(df, v, "centered_vol")
    sc_vol <- plot_scatter(df, v, "centered_vol")
    (ts_t | sc_t | ts_vol | sc_vol) +
      plot_annotation(title = outcome_labels[[v]])
  })

combined_plots %>% walk(print)

cat("Done.\n")
