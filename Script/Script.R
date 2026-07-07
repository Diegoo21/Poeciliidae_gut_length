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

df_plot <- df %>%
  filter(
    !is.na(SL),
    !is.na(AvgGutLength),
    !is.na(Species),
    SL > 0,
    AvgGutLength > 0
  ) %>%
  mutate(
    Species_full = recode(
      Species,
      "Acul" = "Alfaro cultratus",
      "Brha" = "Brachyrhaphis rhabdophora",
      "Pann" = "Priapichthys annectens",
      "Pgil" = "Poecilia gillii",
      "Xumb" = "Xenophallus umbratilis"
    )
  )

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
  log10(AvgGutLength) ~ log10(SL) + log10(EmbryoDryWeight) + as.numeric(EmbryoStage) +
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
    size = 2.4,
    alpha = 0.85,
    stroke = 0.8
  ) +
  geom_smooth(
    aes(group = Species_full),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    fill = "grey80",
    alpha = 0.25,
    linewidth = 0.9,
    fullrange = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Alfaro cultratus" = 16,             # filled circle
      "Brachyrhaphis rhabdophora" = 17,   # filled triangle
      "Poecilia gillii" = 15,             # filled square
      "Priapichthys annectens" = 5,       # open diamond
      "Xenophallus umbratilis" = 1        # open circle
    ),
    labels = c(
      expression(italic("Alfaro cultratus")),
      expression(italic("Brachyrhaphis rhabdophora")),
      expression(italic("Poecilia gillii")),
      expression(italic("Priapichthys annectens")),
      expression(italic("Xenophallus umbratilis"))
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Alfaro cultratus" = "solid",
      "Brachyrhaphis rhabdophora" = "dashed",
      "Poecilia gillii" = "longdash",
      "Priapichthys annectens" = "dotdash",
      "Xenophallus umbratilis" = "dotted"
    ),
    labels = c(
      expression(italic("Alfaro cultratus")),
      expression(italic("Brachyrhaphis rhabdophora")),
      expression(italic("Poecilia gillii")),
      expression(italic("Priapichthys annectens")),
      expression(italic("Xenophallus umbratilis"))
    )
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
  guides(
    shape = guide_legend(
      override.aes = list(
        linetype = c("solid", "dashed", "longdash", "dotdash", "dotted"),
        color = "black",
        fill = NA
      )
    ),
    linetype = "none"
  ) +
  theme_bw(base_size = 14, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = c(0.80, 0.18),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11)
  )

figure2

ggsave(
  filename = "Figure2_gut_length_allometry.pdf",
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
