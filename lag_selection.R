# Lag selection — simple version
# Run after df is built

library(vars)

# 1. VAR lag selection
lp_vars <- as.data.frame(df) %>%
  dplyr::select(all_of(c(outcomes, "target", "centered_vol"))) %>%
  na.omit()

VARselect(lp_vars, lag.max = 12, type = "const")$selection

# 2. ACF plots — just run this for each outcome manually
# Change outcome and n_lags each time
outcome <- "log_gdp_sa"
n_lags  <- 2

work_data <- as.data.frame(df) %>%
  arrange(date) %>%
  mutate(
    inter_term = centered_vol * target,
    lhs_col    = lead(.data[[outcome]], 0) - lag(.data[[outcome]], 1)
  )

controls <- c("target", "centered_vol", "inter_term",
              glue("diff_{outcome}_lag{1:n_lags}"),
              common_controls)
controls <- controls[controls %in% names(work_data)]

m <- lm(as.formula(paste("lhs_col ~", 
         paste(controls, collapse = " + "))),
        data = work_data, na.action = na.omit)

acf(residuals(m), main = glue("{outcome} — {n_lags} lags"))
