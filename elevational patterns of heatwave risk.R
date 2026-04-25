# =========================================================
# Elevational patterns of heatwave risk across taxa
# =========================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# =========================================================
# 1. Load data
# =========================================================

MAMMALS <- read.csv("F:/huijing/biological/MAMMALS_exposure_sensitivity_suitability_risk.csv")
AMPHIBIANS <- read.csv("F:/huijing/biological/AMPHIBIANS_exposure_sensitivity_suitability_risk.csv")
REPTILES <- read.csv("F:/huijing/biological/REPTILES_exposure_sensitivity_suitability_risk.csv")
BIRD <- read.csv("F:/huijing/biological/Bird_exposure_sensitivity_suitability_risk.csv")

# =========================================================
# 2. Harmonize elevation variable
# =========================================================

MAMMALS$elev    <- MAMMALS$dem_mean
AMPHIBIANS$elev <- AMPHIBIANS$dem_mean
REPTILES$elev   <- REPTILES$elev_mean
BIRDS$elev      <- BIRDS$dem_mean

# =========================================================
# 3. Add taxon labels and merge
# =========================================================

MAMMALS$Taxon    <- "Mammals"
AMPHIBIANS$Taxon <- "Amphibians"
REPTILES$Taxon   <- "Reptiles"
BIRDS$Taxon      <- "Birds"

ALL <- bind_rows(MAMMALS, AMPHIBIANS, REPTILES, BIRDS) %>%
  filter(!is.na(elev))

# =========================================================
# 4. Elevation grouping (quartiles)
# =========================================================

ALL$elev_group <- cut(
  ALL$elev,
  breaks = quantile(ALL$elev, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("Low","Mid-low","Mid-high","High")
)

# =========================================================
# 5. Summary statistics
# =========================================================

elev_summary <- ALL %>%
  group_by(elev_group) %>%
  summarise(
    Exposure_mean    = mean(Exposure, na.rm = TRUE),
    Sensitivity_mean = mean(Sensitivity, na.rm = TRUE),
    Suitability_mean = mean(Suitability, na.rm = TRUE),
    Risk_mean        = mean(Risk_index, na.rm = TRUE),
    Risk_sd          = sd(Risk_index, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

elev_taxon_summary <- ALL %>%
  group_by(Taxon, elev_group) %>%
  summarise(
    Risk_mean = mean(Risk_index, na.rm = TRUE),
    Risk_sd   = sd(Risk_index, na.rm = TRUE),
    .groups = "drop"
  )

# =========================================================
# 6. Plotting
# =========================================================

# Factor labels
labels <- c("Low elevation","Lower-mid","Upper-mid","High elevation")

elev_summary$elev_group <- factor(elev_summary$elev_group, labels = labels)
elev_taxon_summary$elev_group <- factor(elev_taxon_summary$elev_group, labels = labels)

# Theme
base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Panel A: by taxon
p1 <- ggplot(elev_taxon_summary,
             aes(x = elev_group, y = Risk_mean,
                 color = Taxon, group = Taxon)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Risk_mean - Risk_sd,
                    ymax = Risk_mean + Risk_sd),
                width = 0.1) +
  labs(x = "Elevation group", y = "Mean heatwave risk") +
  base_theme

# Panel B: overall
p2 <- ggplot(elev_summary,
             aes(x = elev_group, y = Risk_mean, group = 1)) +
  geom_line(color = "black") +
  geom_point(color = "black") +
  geom_errorbar(aes(ymin = Risk_mean - Risk_sd,
                    ymax = Risk_mean + Risk_sd),
                width = 0.1) +
  labs(x = "Elevation group", y = "Mean heatwave risk") +
  base_theme +
  theme(legend.position = "none")

# Combine
final_plot <- p1 + p2 + plot_layout(ncol = 2)

# =========================================================
# 7. Save outputs
# =========================================================

write.csv(elev_summary, "F:/huijing/biological/elevation_summary.csv", row.names = FALSE)
write.csv(elev_taxon_summary, "F:/huijing/biological/elevation_by_taxon.csv", row.names = FALSE)

ggsave("F:/huijing/biological/elevation_risk_plot.png",
       final_plot, width = 12, height = 5, dpi = 600)
