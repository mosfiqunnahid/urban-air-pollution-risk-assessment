############################################################
# Urban Air Pollution Analysis & Public Health Risk Assessment
# Dataset: UCI Air Quality Dataset
# Focus: EDA, Statistical Inference, Probability Analysis,
#        Correlation Analysis, and Regression Modeling
############################################################

# 1. Packages -------------------------------------------------------------

packages <- c(
  "tidyverse", "lubridate", "janitor", "skimr",
  "corrplot", "rstatix", "broom", "scales"
)

install.packages(setdiff(packages, rownames(installed.packages())))

library(tidyverse)
library(lubridate)
library(janitor)
library(skimr)
library(corrplot)
library(rstatix)
library(broom)
library(scales)

# 2. Create Output Folders ------------------------------------------------

dir.create("outputs", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# 3. Import Dataset -------------------------------------------------------

air_raw <- read.csv(
  "data/AirQualityUCI.csv",
  sep = ";",
  dec = ",",
  header = TRUE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

# 4. Data Cleaning --------------------------------------------------------

air_data <- air_raw %>%
  clean_names() %>%
  select(1:15) %>%
  mutate(across(where(is.numeric), ~ na_if(., -200))) %>%
  filter(!is.na(date), !is.na(time)) %>%
  mutate(
    date = dmy(date),
    time_clean = str_replace_all(time, "\\.", ":"),
    datetime = as.POSIXct(
      paste(date, time_clean),
      format = "%Y-%m-%d %H:%M:%S"
    ),
    year = year(date),
    month = month(date),
    hour = hour(datetime),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Autumn",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(datetime), !is.na(season))

# Main analysis dataset
air_clean <- air_data %>%
  drop_na(no2_gt, co_gt, c6h6_gt, t, rh, ah)

# 5. Data Overview --------------------------------------------------------

glimpse(air_clean)
skim(air_clean)

missing_summary <- air_data %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "missing_values"
  ) %>%
  arrange(desc(missing_values))

write_csv(missing_summary, "outputs/missing_value_summary.csv")

# 6. Descriptive Statistics ----------------------------------------------

descriptive_stats <- air_clean %>%
  summarise(
    observations = n(),
    mean_no2 = mean(no2_gt, na.rm = TRUE),
    median_no2 = median(no2_gt, na.rm = TRUE),
    sd_no2 = sd(no2_gt, na.rm = TRUE),
    iqr_no2 = IQR(no2_gt, na.rm = TRUE),
    mean_temperature = mean(t, na.rm = TRUE),
    mean_relative_humidity = mean(rh, na.rm = TRUE),
    mean_absolute_humidity = mean(ah, na.rm = TRUE)
  )

write_csv(descriptive_stats, "outputs/descriptive_statistics.csv")

# 7. Seasonal NO2 Analysis ------------------------------------------------

seasonal_no2 <- air_clean %>%
  group_by(season) %>%
  summarise(
    observations = n(),
    mean_no2 = mean(no2_gt, na.rm = TRUE),
    median_no2 = median(no2_gt, na.rm = TRUE),
    sd_no2 = sd(no2_gt, na.rm = TRUE),
    probability_no2_above_200 = mean(no2_gt > 200, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_no2))

write_csv(seasonal_no2, "outputs/seasonal_no2_summary.csv")

p_seasonal_no2 <- ggplot(seasonal_no2, aes(x = reorder(season, -mean_no2), y = mean_no2)) +
  geom_col() +
  labs(
    title = "Mean NO2 Concentration by Season",
    x = "Season",
    y = "Mean NO2 concentration"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/seasonal_no2_concentration.png", p_seasonal_no2, width = 8, height = 5)

# 8. Probability and Risk Estimation -------------------------------------

risk_summary <- air_clean %>%
  summarise(
    total_observations = n(),
    high_no2_events = sum(no2_gt > 200, na.rm = TRUE),
    probability_high_no2 = mean(no2_gt > 200, na.rm = TRUE)
  )

write_csv(risk_summary, "outputs/high_no2_risk_summary.csv")

risk_by_season <- air_clean %>%
  group_by(season) %>%
  summarise(
    observations = n(),
    high_no2_events = sum(no2_gt > 200, na.rm = TRUE),
    probability_high_no2 = mean(no2_gt > 200, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(probability_high_no2))

write_csv(risk_by_season, "outputs/high_no2_risk_by_season.csv")

# 9. Statistical Inference -----------------------------------------------

# Compare NO2 concentration between winter and summer
winter_summer <- air_clean %>%
  filter(season %in% c("Winter", "Summer"))

wilcox_result <- wilcox_test(
  winter_summer,
  no2_gt ~ season
)

effect_size <- wilcox_effsize(
  winter_summer,
  no2_gt ~ season
)

write_csv(wilcox_result, "outputs/wilcoxon_test_winter_summer.csv")
write_csv(effect_size, "outputs/wilcoxon_effect_size.csv")

# 10. Correlation Analysis ------------------------------------------------

correlation_data <- air_clean %>%
  select(no2_gt, co_gt, c6h6_gt, t, rh, ah) %>%
  drop_na()

cor_matrix <- cor(correlation_data, method = "spearman")

write_csv(
  as.data.frame(cor_matrix) %>% rownames_to_column("variable"),
  "outputs/spearman_correlation_matrix.csv"
)

png("figures/spearman_correlation_heatmap.png", width = 900, height = 700)
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8)
dev.off()

# 11. Regression Modeling -------------------------------------------------

regression_data <- air_clean %>%
  select(no2_gt, co_gt, c6h6_gt, t, rh, ah) %>%
  drop_na()

no2_model <- lm(
  no2_gt ~ co_gt + c6h6_gt + t + rh + ah,
  data = regression_data
)

model_coefficients <- tidy(no2_model)
model_fit <- glance(no2_model)

write_csv(model_coefficients, "outputs/regression_coefficients.csv")
write_csv(model_fit, "outputs/regression_model_fit.csv")

p_regression <- ggplot(regression_data, aes(x = c6h6_gt, y = no2_gt)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Relationship Between Benzene and NO2 Concentration",
    x = "C6H6 concentration",
    y = "NO2 concentration"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/regression_relationship_no2_c6h6.png", p_regression, width = 8, height = 5)

# 12. Export Clean Dataset ------------------------------------------------

write_csv(air_clean, "outputs/clean_air_quality_data.csv")

############################################################
# End of Project Script
############################################################
