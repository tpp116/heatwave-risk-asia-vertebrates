library(terra)
library(dplyr)
library(magrittr)

# ========== Heatwave Indicator Function ==========
Heatwave_Duration <- function(input_seq) {
  threshold <- input_seq[1]; t_seq <- input_seq[-1]
  heat_dates <- which(t_seq > threshold)
  heat_lags <- diff(heat_dates)
  heatwave_events <- split(heat_dates, cumsum(c(1, heat_lags != 1))) %>% 
    Filter(function(x) length(x) >= 3, .)
  heatwave_duration <- lapply(heatwave_events, length) %>% 
    unlist() %>% mean()
  return(heatwave_duration)
}

Heatwave_Total_Days <- function(input_seq) {
  threshold <- input_seq[1]
  t_seq <- input_seq[-1]
  heat_dates <- which(t_seq > threshold)
  heat_lags <- diff(heat_dates)
  events <- split(heat_dates, cumsum(c(1, heat_lags != 1))) %>% 
    Filter(function(x) length(x) >= 3, .)
  if (length(events) == 0) return(0)
  total_days <- sum(sapply(events, length))
  return(total_days)
}

Heatwave_intensity <- function(input_seq) {
  threshold <- input_seq[1]; t_seq <- input_seq[-1]
  heat_dates <- which(t_seq > threshold)
  heat_lags <- diff(heat_dates)
  heatwave_events <- split(heat_dates, cumsum(c(1, heat_lags != 1))) %>% 
    Filter(function(x) length(x) >= 3, .)
  heatwave_tensity <- lapply(heatwave_events, function(x) (mean(t_seq[x] - threshold)) * length(x)) %>%
    unlist() %>% sum(na.rm = TRUE)
  return(heatwave_tensity)
}

Heatwave_MeanIntensity <- function(input_seq) {
  threshold <- input_seq[1]
  t_seq <- input_seq[-1]
  heat_dates <- which(t_seq > threshold)
  heat_lags <- diff(heat_dates)
  events <- split(heat_dates, cumsum(c(1, heat_lags != 1))) %>% 
    Filter(function(x) length(x) >= 3, .)
  if (length(events) == 0) return(NA)
  intensities <- sapply(events, function(x) mean(t_seq[x] - threshold))
  return(mean(intensities, na.rm = TRUE))
}

Heatwave_Tmax <- function(input_seq) {
  threshold <- input_seq[1]
  t_seq <- input_seq[-1]
  heat_dates <- which(t_seq > threshold)
  heat_lags <- diff(heat_dates)
  events <- split(heat_dates, cumsum(c(1, heat_lags != 1))) %>% 
    Filter(function(x) length(x) >= 3, .)
  if (length(events) == 0) return(NA)
  all_days <- unlist(events)
  return(max(t_seq[all_days], na.rm = TRUE))
}

Slight_Heatwave_Frequency <- function(input_seq) {
  threshold <- input_seq[1]; t_seq <- input_seq[-1]
  dates <- which(t_seq > threshold); lags <- diff(dates)
  events <- split(dates, cumsum(c(1, lags != 1))) %>% 
    Filter(function(x) length(x) >= 3 & length(x) <= 4, .)
  return(length(events))
}

Moderate_Heatwave_Frequency <- function(input_seq) {
  threshold <- input_seq[1]; t_seq <- input_seq[-1]
  dates <- which(t_seq > threshold); lags <- diff(dates)
  events <- split(dates, cumsum(c(1, lags != 1))) %>% 
    Filter(function(x) length(x) >= 5 & length(x) <= 7, .)
  return(length(events))
}

Severe_Heatwave_Frequency <- function(input_seq) {
  threshold <- input_seq[1]; t_seq <- input_seq[-1]
  dates <- which(t_seq > threshold); lags <- diff(dates)
  events <- split(dates, cumsum(c(1, lags != 1))) %>% 
    Filter(function(x) length(x) >= 8, .)
  return(length(events))
}

Heatwave_Begin <- function(input_seq) {
  threshold <- input_seq[1]; t_seq <- input_seq[-1]
  dates <- which(t_seq > threshold); lags <- diff(dates)
  events <- split(dates, cumsum(c(1, lags != 1))) %>% 
    Filter(function(x) length(x) >= 3, .)
  return(ifelse(length(events) == 0, NA, min(events[[1]])))
}

Heatwave_End <- function(input_seq) {
  threshold <- input_seq[1]; t_seq <- input_seq[-1]
  dates <- which(t_seq > threshold); lags <- diff(dates)
  events <- split(dates, cumsum(c(1, lags != 1))) %>% 
    Filter(function(x) length(x) >= 3, .)
  return(ifelse(length(events) == 0, NA, max(events[[length(events)]])))
}


# path
library(terra)
library(dplyr)
library(magrittr)


threshold_ras <- rast("D:/heatwave_thresholds/global_land_90p_1979_2014.tif")
years <- 1979:2014
masked_dir <- "D:/heatwave_thresholds/masked"
out_base_dir <- "D:/heatwave_results"

for (sub in c("duration","intensity","slight","moderate","severe","begin","end","total_days","mean_intensity","tmax")) {
  dir.create(file.path(out_base_dir, sub), recursive = TRUE, showWarnings = FALSE)
}

for (year in years) {
  message("processing：", year)
  
  masked_ras <- rast(file.path(masked_dir, paste0("masked_", year, ".tif")))
  masked_aligned <- resample(masked_ras, threshold_ras, method = "bilinear")
  input_ras <- c(threshold_ras, masked_aligned) %>% sds()
  
  # calculate
  duration <- lapp(input_ras, Heatwave_Duration)
  intensity <- lapp(input_ras, Heatwave_intensity)
  slight <- lapp(input_ras, Slight_Heatwave_Frequency)
  moderate <- lapp(input_ras, Moderate_Heatwave_Frequency)
  severe <- lapp(input_ras, Severe_Heatwave_Frequency)
  begin <- lapp(input_ras, Heatwave_Begin)
  end <- lapp(input_ras, Heatwave_End)
  total_days <- lapp(input_ras, Heatwave_Total_Days)
  mean_intensity <- lapp(input_ras, Heatwave_MeanIntensity)
  tmax_hw <- lapp(input_ras, Heatwave_Tmax)
  
  # Save results
  writeRaster(duration, file.path(out_base_dir,"duration", paste0(year, ".tif")), overwrite=TRUE)
  writeRaster(intensity, file.path(out_base_dir,"intensity", paste0(year, ".tif")), overwrite=TRUE)
  writeRaster(slight, file.path(out_base_dir,"slight", paste0(year, ".tif")), overwrite=TRUE)
  writeRaster(moderate, file.path(out_base_dir,"moderate", paste0(year, ".tif")), overwrite=TRUE)
  writeRaster(severe, file.path(out_base_dir,"severe", paste0(year, ".tif")), overwrite=TRUE)
  writeRaster(begin, file.path(out_base_dir,"begin", paste0(year, ".tif")), overwrite=TRUE)
  writeRaster(end, file.path(out_base_dir,"end", paste0(year, ".tif")), overwrite=TRUE)
  writeRaster(total_days, file.path(out_base_dir,"total_days", paste0(year, ".tif")), overwrite=TRUE)
  writeRaster(mean_intensity, file.path(out_base_dir,"mean_intensity", paste0(year, ".tif")), overwrite=TRUE)
  writeRaster(tmax_hw, file.path(out_base_dir,"tmax", paste0(year, ".tif")), overwrite=TRUE)
  
  gc()
}
