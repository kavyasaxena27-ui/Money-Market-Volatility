# =============================================================
# instrument_diagnostics_plot.R
# Visual check of Rigobon instruments
# Run after make_instrument() so df has iv1_dummy_r1:r4 etc
# =============================================================

library(tidyverse)
library(patchwork)

df_plot <- as.data.frame(df) %>% arrange(date)

# -------------------------------------------------------------
# 1. All iv1 instruments over time (for sigma)
# Each should be most active in its own regime
# -------------------------------------------------------------

p1 <- df_plot %>%
  dplyr::select(date, regime,
                iv1_dummy_r1, iv1_dummy_r2,
                iv1_dummy_r3, iv1_dummy_r4) %>%
  pivot_longer(-c(date, regime),
               names_to = "instrument", values_to = "value") %>%
  ggplot(aes(x = date, y = value, colour = instrument)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  facet_wrap(~instrument, ncol = 2, scales = "free_y") +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "iv1 instruments (for sigma)",
       subtitle = "Each should be most active in its own regime, near zero elsewhere",
       x = NULL, y = "Instrument value") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

# -------------------------------------------------------------
# 2. iv1_dummy_r1 vs centred vol — should track closely in regime 1
# -------------------------------------------------------------

# Rescale vol for overlay
vol_range  <- range(df_plot[[vol_col]], na.rm = TRUE)
iv1_range  <- range(df_plot$iv1_dummy_r1, na.rm = TRUE)
scale_f    <- diff(iv1_range) / diff(vol_range)
offset     <- iv1_range[1] - vol_range[1] * scale_f

p2 <- df_plot %>%
  mutate(vol_scaled = .data[[vol_col]] * scale_f + offset) %>%
  ggplot(aes(x = date)) +
  geom_rect(
    data = df_plot %>%
      group_by(regime) %>%
      summarise(xmin = min(date), xmax = max(date), .groups = "drop"),
    aes(xmin = xmin, xmax = xmax, fill = regime),
    ymin = -Inf, ymax = Inf, alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(aes(y = iv1_dummy_r1), colour = "firebrick",
            linewidth = 0.6) +
  geom_line(aes(y = vol_scaled), colour = "steelblue",
            linewidth = 0.4, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(
    name     = "iv1_dummy_r1",
    sec.axis = sec_axis(~ (. - offset) / scale_f,
                        name = vol_col)
  ) +
  labs(
    title    = "iv1_dummy_r1 vs centred vol",
    subtitle = "Red = instrument  |  Blue dashed = vol  |  Should co-move in regime 1",
    x = NULL, fill = "Regime"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.title.y       = element_text(colour = "firebrick"),
        axis.title.y.right = element_text(colour = "steelblue"))

# -------------------------------------------------------------
# 3. Instrument summary stats per regime
# Want: large variance in own regime, small elsewhere
# -------------------------------------------------------------

p3 <- df_plot %>%
  dplyr::select(date, regime,
                iv1_dummy_r1, iv1_dummy_r2,
                iv1_dummy_r3, iv1_dummy_r4) %>%
  pivot_longer(-c(date, regime),
               names_to = "instrument", values_to = "value") %>%
  group_by(regime, instrument) %>%
  summarise(
    variance = var(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = regime, y = variance, fill = instrument)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title    = "Instrument variance by regime",
    subtitle = "Each instrument should have highest variance in its own regime",
    x = NULL, y = "Variance", fill = "Instrument"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# -------------------------------------------------------------
# 4. fs1_resid over time — should be mean zero, heteroskedastic
# -------------------------------------------------------------

p4 <- df_plot %>%
  ggplot(aes(x = date, y = fs1_resid)) +
  geom_rect(
    data = df_plot %>%
      group_by(regime) %>%
      summarise(xmin = min(date), xmax = max(date), .groups = "drop"),
    aes(xmin = xmin, xmax = xmax, fill = regime),
    ymin = -Inf, ymax = Inf, alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(colour = "grey30", linewidth = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "First stage residuals (fs1_resid)",
    subtitle = "Should be mean zero and heteroskedastic across regimes",
    x = NULL, y = "Residual", fill = "Regime"
  ) +
  theme_minimal(base_size = 10)

# -------------------------------------------------------------
# Print all
# -------------------------------------------------------------

print(p1)
print(p2)
print(p3)
print(p4)

# Combined
(p2 / p4) | (p1 / p3)
