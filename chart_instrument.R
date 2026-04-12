# After make_instrument() runs
diag_df <- as.data.frame(df) %>%
  arrange(date) %>%
  mutate(
    vol_x_shock  = centered_vol * target,
    fs1_fitted   = vol_x_shock  - attr(df, "fs1")$residuals,  
    # simpler:
    fs1_fitted   = fitted(attr(df, "fs1")),
    fs2_fitted   = fitted(attr(df, "fs2"))
  )

# Plot 1: sigma over time with fitted values
p1 <- ggplot(diag_df, aes(x = date)) +
  geom_rect(data = diag_df %>% group_by(regime) %>%
              summarise(xmin = min(date), xmax = max(date), .groups = "drop"),
            aes(xmin = xmin, xmax = xmax, fill = regime),
            ymin = -Inf, ymax = Inf, alpha = 0.15, inherit.aes = FALSE) +
  geom_line(aes(y = centered_vol), colour = "grey40", linewidth = 0.5) +
  geom_line(aes(y = fs1_fitted),   colour = "firebrick",
            linewidth = 0.8, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Eq (2): sigma ~ X  |  Actual vs fitted",
       x = NULL, y = "centred vol", fill = "Regime") +
  theme_minimal()

# Plot 2: sigma*S over time with fitted values  
p2 <- ggplot(diag_df, aes(x = date)) +
  geom_rect(data = diag_df %>% group_by(regime) %>%
              summarise(xmin = min(date), xmax = max(date), .groups = "drop"),
            aes(xmin = xmin, xmax = xmax, fill = regime),
            ymin = -Inf, ymax = Inf, alpha = 0.15, inherit.aes = FALSE) +
  geom_line(aes(y = vol_x_shock), colour = "grey40", linewidth = 0.5) +
  geom_line(aes(y = fs2_fitted),  colour = "steelblue",
            linewidth = 0.8, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Eq (3): sigma*S ~ X  |  Actual vs fitted",
       x = NULL, y = "vol x shock", fill = "Regime") +
  theme_minimal()

p1 / p2
