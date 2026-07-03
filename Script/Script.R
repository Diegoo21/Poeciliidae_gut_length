############################################################
# Gut length reflects trophic ecology but shows no tradeoff
# with reproductive allotment in livebearing fishes
############################################################

# 1. Load packages -------------------------------------------------------

library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(emmeans)
library(psych)

# 2. Load and prepare data ----------------------------------------------

df <- read_csv("Poeciliid_specimen_data.csv") %>%
  filter(Include == "Yes")

df_mature <- filter(df, Stage == "Female")

# 3. Measurement repeatability ------------------------------------------
# Repeatability of duplicate gut length measurements
icc_results <- ICC(
  df[, c("GutLength1", "GutLength2")]
)

icc_results

# 4. Gut length model ----------------------------------------------------

gut_slope <- lmer(log10(AvgGutLength) ~ log10(SL) * Species + (1|Site), data =df)

anova(gut_slope)
summary(gut_slope)
icc(gut_slope)

# Estimated marginal means at mean standard length
gut_emm <- emmeans(
  gut_slope,
  ~ Species,
  at = list(`log10(SL)` = mean(log10(df$SL), na.rm = TRUE))
)

gut_emm_df <- as.data.frame(gut_emm) %>%
  mutate(
    GutLength_mm = 10^emmean
  )

gut_emm_df

# Pairwise comparisons among species
gut_pairs <- pairs(
  gut_emm,
  adjust = "tukey"
)

summary(gut_pairs)

# 5. Tradeoff model ------------------------------------------------------

tradeoff_model <- lmer(
  log10(AvgGutLength) ~ log10(SL) + log10(EmbryoDryWeight) +
    Species + (1 | Site),
  data = df_mature
)

anova(tradeoff_model)
summary(tradeoff_model)
icc(tradeoff_model)

# 6. Figure 2 ------------------------------------------------------------

figure2 <- ggplot(
  df_plot,
  aes(
    x = SL,
    y = AvgGutLength,
    shape = Species_full,
    linetype = Species_full
  )
) +
  geom_point(
    color = "black",
    size = 2.3,
    alpha = 0.8
  ) +
  geom_smooth(
    aes(group = Species_full),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    fill = "grey75",
    alpha = 0.20,
    linewidth = 0.9,
    fullrange = FALSE
  ) +
  scale_x_log10(
    breaks = c(20, 30, 40, 50, 60, 70)
  ) +
  scale_y_log10(
    breaks = c(10, 20, 30, 50, 100, 200, 300)
  ) +
  labs(
    x = "Standard length (mm)",
    y = "Gut length (mm)",
    shape = "Species",
    linetype = "Species"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )

figure2

ggsave(
  filename = "Figure2_gut_length_allometry.png",
  plot = figure2,
  width = 7,
  height = 5,
  dpi = 300
)

# 7. Export supplementary table -----------------------------------------

table_s1 <- as.data.frame(summary(gut_pairs)) %>%
  mutate(
    Gut_length_ratio = round(10^estimate, 2)
  ) %>%
  select(
    Comparison = contrast,
    Gut_length_ratio,
    p_value = p.value
  )

table_s1

write_csv(table_s1, "TableS1_pairwise_gut_length.csv")