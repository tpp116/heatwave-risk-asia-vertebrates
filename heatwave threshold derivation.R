# =========================================================
# Heatwave threshold calculation (90th percentile)
# Author: [Your Name]
# =========================================================

library(terra)
library(sf)
library(dplyr)
library(purrr)

# -----------------------------
# Global settings
# -----------------------------
terraOptions(tempdir = "path_to_tempdir")

base_dir <- "D:/tmax/"   # e.g., yearly folders with .nc files
shp_path <- "D:/tpp/WorldCountries/WorldCountries.shp"           # global land boundary
output_dir <- "D:/heatwave_thresholds/"

masked_dir <- file.path(output_dir, "masked")
tiles_dir  <- file.path(output_dir, "tiles_90p")
final_output <- file.path(output_dir, "global_land_90p.tif")

dir.create(masked_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tiles_dir, showWarnings = FALSE, recursive = TRUE)

years <- 2015:2024

# -----------------------------
# Load template raster
# -----------------------------
sample_file <- list.files(file.path(base_dir, years[1]),
                          pattern = "\\.nc$", full.names = TRUE)[1]
r_template <- rast(sample_file)

# -----------------------------
# Create land mask
# -----------------------------
land <- vect(shp_path)
land <- project(land, crs(r_template))

land_mask <- rasterize(land, r_template, field = 1)
writeRaster(land_mask, file.path(masked_dir, "land_mask.tif"), overwrite = TRUE)

# -----------------------------
# Apply land mask to each year
# -----------------------------
for (year in years) {
  message("Processing year: ", year)
  
  file_paths <- list.files(file.path(base_dir, year),
                           pattern = "\\.nc$", full.names = TRUE)
  file_paths <- sort(file_paths)
  
  r_year <- rast(file_paths)
  r_land <- mask(r_year, land_mask)
  
  writeRaster(r_land,
              file.path(masked_dir, paste0("masked_", year, ".tif")),
              overwrite = TRUE)
  gc()
}

# -----------------------------
# Create 10° × 10° tiles
# -----------------------------
create_global_tiles <- function(dx = 10, dy = 10) {
  lon_seq <- seq(-180, 180 - dx, by = dx)
  lat_seq <- seq(-60, 90 - dy, by = dy)
  
  tiles <- purrr::map_dfr(lon_seq, function(lon) {
    purrr::map_dfr(lat_seq, function(lat) {
      poly <- as.polygons(ext(lon, lon + dx, lat, lat + dy))
      st_as_sf(poly)
    })
  })
  
  vect(tiles)
}

tiles <- create_global_tiles()

# -----------------------------
# Tile-based threshold calculation
# -----------------------------
tile_results <- list()

for (i in seq_len(nrow(tiles))) {
  message("Processing tile ", i, "/", nrow(tiles))
  
  tile <- tiles[i, ]
  stack_list <- list()
  
  for (year in years) {
    r_path <- file.path(masked_dir, paste0("masked_", year, ".tif"))
    if (!file.exists(r_path)) next
    
    r <- rast(r_path)
    r_crop <- crop(r, tile)
    stack_list[[length(stack_list) + 1]] <- r_crop
  }
  
  if (length(stack_list) == 0) next
  
  r_stack <- do.call(c, stack_list)
  
  r_thresh <- app(r_stack, function(x) {
    quantile(x, probs = 0.9, na.rm = TRUE)
  })
  
  out_path <- file.path(tiles_dir, paste0("tile_", i, "_90p.tif"))
  writeRaster(r_thresh, out_path, overwrite = TRUE)
  
  tile_results[[i]] <- rast(out_path)
  gc()
}

# -----------------------------
# Merge tiles
# -----------------------------
message("Merging tiles...")
global_thresh <- do.call(merge, tile_results)

writeRaster(global_thresh, final_output, overwrite = TRUE)

plot(global_thresh, main = "Global 90th Percentile Threshold")