# =========================================================
# Species-level heatwave risk assessment (ESA framework): taking AMPHIBIANS as a example
# =========================================================

library(data.table)
library(dplyr)
library(janitor)

# =========================================================
# 1. Load data
# =========================================================

risk_part1 <- fread("F:/huijing/biological/AMPHIBIANS_PART1.csv")
risk_part2 <- fread("F:/huijing/biological/AMPHIBIANS_PART2.csv")

risk_data <- bind_rows(risk_part1, risk_part2)

traits <- fread("F:/huijing/biological/Amphibian_traits_avg.csv")

# Clean column names
risk_data <- clean_names(risk_data)
traits    <- clean_names(traits)

# =========================================================
# 2. Merge traits with risk data
# =========================================================

data <- risk_data %>%
  left_join(traits, by = c("species" = "species_full")) %>%
  filter(!is.na(svl))   # keep species with trait data

# =========================================================
# 3. Normalization function (0–1 scaling)
# =========================================================

norm01 <- function(x) {
  if (all(is.na(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

# =========================================================
# 4. Exposure (E)
# =========================================================

exposure_vars <- c(
  "mhw_freq","mhw_total_days","mhw_max_int",
  "mhw_mean_int","mhw_cum_int",
  "slope_intensity","slope_freq","slope_cum",
  "slope_max","slope_total"
)

data[exposure_vars] <- lapply(data[exposure_vars], function(x) norm01(as.numeric(x)))

data$Exposure <- rowMeans(data[exposure_vars], na.rm = TRUE)

# =========================================================
# 5. Sensitivity (S)
# =========================================================

# Remove extinct species
data <- data %>% filter(iucn_category != "EX")

# IUCN scoring
iucn_map <- c("LC"=1, "NT"=2, "VU"=3, "EN"=4, "CR"=5)
data$iucn_num <- norm01(as.numeric(iucn_map[data$iucn_category]))

# Trait-based sensitivity
sens_traits <- c("svl","hll","fll","hl","hw","tl")
sens_traits <- sens_traits[sens_traits %in% names(data)]

data[sens_traits] <- lapply(data[sens_traits], function(x) norm01(as.numeric(x)))

data$Trait_sensitivity <- rowMeans(data[sens_traits], na.rm = TRUE)

# Final sensitivity
data$Sensitivity <- rowMeans(data[, c("iucn_num","Trait_sensitivity")], na.rm = TRUE)

# =========================================================
# 6. Adaptive capacity / Suitability (A)
# =========================================================

adapt_vars <- c("pa_count","pa_quality_mean","area_global","dem_mean","dem_range")
adapt_vars <- adapt_vars[adapt_vars %in% names(data)]

data[adapt_vars] <- lapply(data[adapt_vars], function(x) norm01(as.numeric(x)))

data$Suitability <- rowMeans(data[adapt_vars], na.rm = TRUE)

# =========================================================
# 7. Risk index
# =========================================================

data$Risk_index <- with(data, Exposure * Sensitivity * (1 - Suitability))

# =========================================================
# 8. Risk classification
# =========================================================

breaks <- unique(quantile(data$Risk_index, probs = seq(0,1,0.2), na.rm = TRUE))

if (length(breaks) < 2) {
  data$Risk_level <- "Low"
} else {
  data$Risk_level <- cut(
    data$Risk_index,
    breaks = breaks,
    labels = c("Very low","Low","Medium","High","Very high")[1:(length(breaks)-1)],
    include.lowest = TRUE
  )
}

# =========================================================
# 9. Output
# =========================================================

write.csv(
  data,
  "F:/huijing/biological/output_species_heatwave_risk.csv",
  row.names = FALSE
)

# Summary
summary(data[, c("Exposure","Sensitivity","Suitability","Risk_index")])
table(data$Risk_level)
