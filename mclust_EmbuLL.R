
# ============================================================
# Script: mclust_EmbuLL.R
# Project: Stratification_Embu
# Purpose: Environmental clustering for weighted stratified
#          sampling design of agroforestry landscape in
#          Embu Living Lab, Kenya
# Author: Sophia Gödde
# Date: April 2026
# ============================================================
# ── PACKAGES ─────────────────────────────────────────────
library(terra)       
library(mclust)     
library(tidyverse)  
library(sf)         
library(corrplot)   
library(RColorBrewer)
library(spsurvey)
library(osmdata) # for masking buidlings and roads when allocating sample points 



# NOTE for self: THIS IS THE CORRECT VERSION connected to GitHub repo (GitHub_Stratification_Embu_LL)


######################## 1. LOAD + Crop RASTERS to LL Boundary ###################
###LL 
ll_boundary <- st_read("data/shapefiles/One_LL_Kenya_V1.shp")
# Check CRS —  match rasters (EPSG:4326)
print(st_crs(ll_boundary))
if (st_crs(ll_boundary)$epsg != 4326) {
  ll_boundary <- st_transform(ll_boundary, 4326)
  cat("Reprojected to EPSG:4326\n")}

# Convert to terra SpatVector 
ll_vect <- vect(ll_boundary)
plot(ll_boundary$geometry, main = "Embu Living Lab boundary")

####### Load RASTERS 
data_path <- "data/rasters/"

# I now downloaded already in LL cropped so this part is not needed 
# clip them before so size and processing smaller 
# load_clip <- function(filename)
#   {r <- rast(paste0(data_path, filename))
#   r <- crop(r, ll_vect)
#   r <- mask(r, ll_vect)
#   return(r)}
# 
# r_temp      <- load_clip("Temperature_MODIS_2000_2025.tif")
# r_rainfall  <- load_clip("Rainfall_CHIRPS_2000_2025.tif")
# r_elevation <- load_clip("Elevation_SRTM.tif")
# r_slope     <- load_clip("Slope_SRTM.tif")
# r_aspect    <- load_clip("Aspect_SRTM.tif")
# r_cropland  <- load_clip("ldsf_kenya_cropland_embuLL.tif")
# r_treecover <- load_clip("ldsf_kenya_tc_embuLL.tif")
# cat("All rasters loaded and clipped to Living Lab boundary\n")
# Dropped: TreeCover_Hansen_pct_2024.tif ; CanopyHeight_ETH_2020_m.tif ; 
#  checking
print(r_rainfall)
print(r_temp)
print(r_treecover)
print(r_cropland)
#print(r_distwater)

# Set all values <= 0 to NA for both LDSF layers
# info  Tor : missing/invalid values are coded as <= 0
r_treecover <- classify(r_treecover, cbind(-Inf, 0, NA))
r_cropland  <- classify(r_cropland,  cbind(-Inf, 0, NA))

# Verify — minimum should now be a small positive number, not -9999 or 0
cat("Tree cover range after NA fix:\n")
print(minmax(r_treecover))

cat("Cropland range after NA fix:\n")
print(minmax(r_cropland))


######################## 2. RESAMPLE ALL TO 100m GRID ######################## 
# 1km = MODIS and 5km CHIRPS need to be downscaled 
# this is acceptable as we are interpolating within known gradients,
# not extrapolating beyond the data range
#   - 100m gives ~10x more pixels
# 0.0009 degrees ≈ 100m at equator (Embu is near 0° latitude)
# Bilinear for continuous variables (smooth interpolation)

ref_grid <- rast(
  ext(ll_vect),
  resolution = 0.0009,
  crs        = "EPSG:4326")
# Continuous variables: bilinear interpolation
r_temp_r      <- resample(r_temp,      ref_grid, method = "bilinear")
r_rainfall_r  <- resample(r_rainfall,  ref_grid, method = "bilinear")
r_elevation_r <- resample(r_elevation, ref_grid, method = "bilinear")
r_slope_r     <- resample(r_slope,     ref_grid, method = "bilinear")
r_aspect_r    <- resample(r_aspect,    ref_grid, method = "bilinear")
r_treecover_r <- resample(r_treecover, ref_grid, method = "bilinear")
r_cropland_r  <- resample(r_cropland,  ref_grid, method = "bilinear")

cat("All layers resampled to 100m reference grid\n")

######################## 3. STACK continous variable LAYERS ######################## 
# Land cover excluded from mclust input (categorical) bcs mclust assumed gaussian distribution 
# landcover used later for post-hoc interpretation/validation  - is that approach methodologically correct??
env_stack <- c(
  r_temp_r,
  r_rainfall_r,
  r_elevation_r,
  r_slope_r,
  r_aspect_r,
  r_treecover_r,
  r_cropland_r)

names(env_stack) <- c(
  "Temp_mean",
  "Rainfall_MAP",
  "Elevation",
  "Slope",
  "Aspect",
  "Tree_cover",
  "Cropland")

env_stack <- mask(env_stack, ll_vect)

print(env_stack)
plot(env_stack) 
cat("Stack bands:", names(env_stack), "\n")




######################## 4. EXTRACT PIXEL VALUES TO DATA FRAME ######################## 
env_df <- as.data.frame(env_stack, xy = TRUE, na.rm = TRUE)
# xy = TRUE keeps the coordinates to map clusters back to raster 
nrow(env_df)   # how many valid pixels: 119217 
cat("Total valid pixels in Living Lab:", nrow(env_df), "\n") #119217




######################## 5. CHECK CORRELATIONS ######################## 
# highly correlated vars (>0.7) will be removed 
vars_to_use <- c("Temp_mean", "Rainfall_MAP", "Elevation", "Slope", "Aspect", "Tree_cover", "Cropland")
cor_matrix <- cor(env_df[, vars_to_use],
                  use    = "complete.obs",
                  method = "pearson")
corrplot(cor_matrix,
         method      = "color",
         type        = "upper",
         tl.col      = "black",
         addCoef.col = "black",
         number.cex  = 0.8,
         title       = "Pearson correlation — Embu LL environmental variables (7)",
         mar         = c(0, 0, 2, 0))
# highest correlation are Temp and Rain ; Temp + Treecover ; Temp + Slope ; cropland +Treecover

# is this code chunk correct? 
high_cor <- which(abs(cor_matrix) > 0.7 & cor_matrix != 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("\n Highly correlated pairs (|r| > 0.7):\n")
  print(data.frame(
    Var1 = rownames(cor_matrix)[high_cor[, 1]],
    Var2 = colnames(cor_matrix)[high_cor[, 2]],
    r    = round(cor_matrix[high_cor], 3)))
  }else {cat("No pairs with |r| > 0.7 — all variables retained\n")}



######################## 6. SCALE VARIABLES ######################## 
# mclust is sensitive to variable scale  -->  mean=0 and SD=1
## different units (°C, mm, m, %) contribute equally to clustering
env_scaled <- scale(env_df[, vars_to_use])
summary(env_scaled)




######################## 7. BIC - OPTIMAL NUMBER OF CLUSTERS ######################## 
# mclust uses BIC (Bayesian Information Criterion) to find best fitting model and select clusters
# LDSF requirement: 3 plots per cluster is minimum 
set.seed(55)
n_pixels <- nrow(env_scaled)
bic_result <- mclustBIC(env_scaled, G = 4:7)
summary(bic_result)
plot(bic_result, main = "BIC — G=4 to 7 for Embu LL (100m)")

#  for 5:10 --> VV10 with BIC -1854266
# for 4:7 --> VVV7 with BIC -1922277
# for sampling 15 days 7 clusters is upper limit since it gives only 7 plots to distribute according to a weighted design


cat("Continuing with G=7 for now.\n")




######################## 8. FIT FINAL MCLUST MODEL ######################## 
optimal_G <- 7

  model_name <- dimnames(bic_result)[[2]][
    which.max(apply(bic_result, 2, max, na.rm = TRUE))]
  
cat("Best covariance model:", model_name, "\n")

# fit on full pixel dataset for BIC 
mc_model <- Mclust(env_scaled, optimal_G, modelNames = model_name)
summary(mc_model)


#cluster ID and uncertainty (prob. for not belonging to same cluster) for each pixel
env_df$cluster <- as.factor(mc_model$classification)
env_df$uncertainty <- mc_model$uncertainty
# ~0 = very confident assignment (good)
# ~0.5 = pixel sits equally between two clusters (transition zone)
cat("Pixels per cluster:\n")
print(table(env_df$cluster))





######################## 9. VISUALIZE + INTERPRETE CLUSTERS ######################## 
cluster_means <- aggregate(
  env_df[, vars_to_use],
  by  = list(Cluster = env_df$cluster),
  FUN = mean)
cat("\nCluster means:\n")
print(round(cluster_means[, -1], 2))
# Use  table to name  strata

# Cluster area summary — kept as separate object from allocation table
cluster_area <- data.frame(
  Cluster  = levels(env_df$cluster),
  N_pixels = as.integer(table(env_df$cluster)))
cluster_area$Pct_area <- round(
  cluster_area$N_pixels / sum(cluster_area$N_pixels) * 100, 1)
cat("\nCluster area breakdown:\n")
print(cluster_area)


# Heatmap of cluster means
cluster_means_long <- cluster_means %>%
  pivot_longer(-Cluster, names_to = "Variable", values_to = "Mean")

ggplot(cluster_means_long,
       aes(x = Variable, y = as.factor(Cluster), fill = Mean)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(Mean, 1)), size = 3.5) +
  scale_fill_distiller(palette = "RdYlGn", direction = 1) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  labs(
    title    = "Mean environmental conditions per cluster",
    subtitle = paste0("Embu Living Lab — ", model_name, ",",
                      optimal_G, " mclust solution"),
    x = "", y = "Cluster"
  )

ggsave("outputs/cluster_means_heatmap.png", width = 9, height = 5, dpi = 300)









######################## 10. MAP CLUSTERS BACK TO RASTER ######################## 
# Convert env_df points to sf object using their coordinates
ref_grid <- rast(
  ext(ll_vect),
  resolution = 0.0009,
  crs        = "EPSG:4326")

cluster_points <- st_as_sf(
  env_df[, c("x", "y", "cluster", "uncertainty")],
  coords = c("x", "y"),
  crs    = 4326)  # WGS84 — same as GEE export and shapefile

# not sure if that chunk is correct?! 

cluster_rast <- rasterize(
  vect(cluster_points),
  ref_grid,
  field  = "cluster",
  fun    = "last")
names(cluster_rast) <- "Cluster"

# Rasterize uncertainty
uncert_rast <- rasterize(
  vect(cluster_points),
  ref_grid,
  field = "uncertainty",
  fun   = "last")
names(uncert_rast) <- "Uncertainty"

# Mask both to LL boundary
cluster_rast <- mask(cluster_rast, ll_vect)
uncert_rast  <- mask(uncert_rast,  ll_vect)

# Checking coverage — should show 3 unique values only
# # Verify full coverage — Non-NA should match nrow(env_df)
cat("Unique cluster values:", unique(values(cluster_rast), na.rm = TRUE), "\n")
cat("Non-NA pixels:", sum(!is.na(values(cluster_rast))),
    "— expected:", nrow(env_df), "\n")
# If these do not match, check CRS alignment between shapefile and rasters




######### Plotting final cluster map #####################
display.brewer.all()
cls_colors <- brewer.pal(max(optimal_G, 8), "Set2")[1:optimal_G]

#cls_colors_6 <- c("#2166ac", "#f4a582", "#1a9641", "darkolivegreen4", "darkorange2", "khaki3")

plot(cluster_rast,
     col    = cls_colors,
     main   = paste0(optimal_G, " Assigned Clusters according to 7 environmental 
     and landscape variables for the GALILEO Living Lab, Embu County, Kenya"),
     legend = TRUE,
     axes   = TRUE)
plot(ll_vect, add = TRUE, border = "black", lwd = 2)

# Uncertainty map — white = confident, red = transition zone
plot(uncert_rast,
     col  = rev(heat.colors(50)),
     main = "Cluster Uncertainty ")
plot(ll_vect, add = TRUE, border = "black", lwd = 2)
# White/pale areas = confident assignment (most pixels)
# Orange/red areas = ambiguous pixels between clusters
# These are NOT sampling errors — they are real landscape transitions





######### Export cluster&uncertainty map in .tif for QGIS/ODK ########################### 
# Unsure about this follwing part 
dir.create("outputs/", showWarnings = FALSE)

writeRaster(cluster_rast,
            "outputs/LL_mclust_7_FINAL.tif",
            overwrite = TRUE,
            datatype  = "INT1U")

writeRaster(uncert_rast,
            "outputs/LL_uncertainty_mclust_7_FINAL.tif",
            overwrite = TRUE)

# maybe useful to export cluster assignments with coordinates as CSV
write_csv(env_df[, c("x", "y", "cluster", "uncertainty")],
          "outputs/Embu_cluster_7_FINAL.csv")

cat("Cluster rasters exported\n")


############# 11. EXCLUSION OF HIGH UNCERTAINTIES ########################### 
cat("\nUncertainty distribution:\n")
print(summary(env_df$uncertainty))
# Max=0.73 

cat("\nPixels with uncertainty >= 0.20 (transition zones):",
    sum(env_df$uncertainty > 0.20), "\n") #24406
cat("Percentage of total pixels:",
    round(sum(env_df$uncertainty > 0.20) / nrow(env_df) * 100, 1), "%\n")
# So 16.3% of total pixcels are in high uncertainty areas and therefore excluded for sampling. 
# median is only 0.021 --> most pixels are assigned with very high confidence! good:)

# Exclude uncertainty > 0.20 from sampling 
#(mainly pixels between clusters and therefore it may not be representable for any of the strata)

#### Confident sampling frame ##
sampling_frame <- env_df %>%
  filter(uncertainty <= 0.20) %>%
  dplyr::select(x, y, cluster)
cat("\nPixels in confident sampling frame:", nrow(sampling_frame), "\n")
#Pixels in confident sampling frame: 125681 
cat("Pixels excluded:", nrow(env_df) - nrow(sampling_frame), "\n")
#Pixels excluded: 24406 
write_csv(sampling_frame, "outputs/EmbuLL_sampling_frame_confident.csv")



# applying uncertainty mask to cluster raster
# This is the raster used for plot placement — uncertain pixels are NA
cluster_rast_confident <- ifel(uncert_rast > 0.20, NA, cluster_rast)
cluster_rast_confident <- mask(cluster_rast_confident, ll_vect)

writeRaster(cluster_rast_confident,
            "outputs/EmbuLL_clusters_VVV7_confident.tif",
            overwrite = TRUE,
            datatype  = "INT1U")




############# 12.  EXCLUSION MASK FROM OSM ############# 
# Excludes roads, buildings and waterbodies from sampling frame
# Gardens/homesteads around buildings are ok
# so building buffer is kept small (10m around structure only) bcs gardens can be included
# download from R didnt work (timeout) so downloaded locally Kenya OSM shapefiles 
# https://download.geofabrik.de/africa/kenya-latest.shp.zip
osm_path <- "data/shapefiles/osm/"

# Load and clip each layer to LL boundary
roads_raw <- st_read(paste0(osm_path, "gis_osm_roads_free_1.shp")) %>%
  st_crop(st_bbox(ll_boundary)) %>%
  filter(fclass %in% c("primary", "secondary", "tertiary",
                       "residential", "unclassified", "track"))

buildings_raw <- st_read(paste0(osm_path,
                                "gis_osm_buildings_a_free_1.shp")) %>%
  st_crop(st_bbox(ll_boundary))


# Example approach
library(sf)
buildings_raw <- st_read(
  paste0(osm_path, "gis_osm_buildings_a_free_1.shp"),
  bbox = st_bbox(ll_boundary) # Limits reading to this area
)


water_poly_raw <- st_read(paste0(osm_path,
                                 "gis_osm_water_a_free_1.shp")) %>%
  st_crop(st_bbox(ll_boundary))

waterways_raw <- st_read(paste0(osm_path,
                                "gis_osm_waterways_free_1.shp")) %>%
  st_crop(st_bbox(ll_boundary))

protectedareas_raw <- st_read(paste0(osm_path,
                                     "gis_osm_protected_areas_a_free_1.shp")) %>%
  st_crop(st_bbox(ll_boundary))

cat("Roads loaded:", nrow(roads_raw), "features\n")
cat("Buildings loaded:", nrow(buildings_raw), "features\n")
cat("Water polygons loaded:", nrow(water_poly_raw), "features\n")
cat("Waterways loaded:", nrow(waterways_raw), "features\n")
cat("Protected areas:", nrow(protectedareas_raw), "features\n")


# Buffer each layer — all in UTM metres then back to WGS84
roads_buf <- roads_raw %>%
  st_transform(32737) %>%
  st_buffer(dist = 30) %>%    # 30m either side of road
  st_union() %>%
  st_transform(4326)

buildings_buf <- buildings_raw %>%
  st_transform(32737) %>%
  st_buffer(dist = 10) %>%    # 10m around structure only
  st_union() %>%              # gardens/homesteads are OK
  st_transform(4326)

water_buf <- water_poly_raw %>%
  st_transform(32737) %>%
  st_buffer(dist = 5) %>%
  st_union() %>%
  st_transform(4326)

waterways_buf <- waterways_raw %>%
  st_transform(32737) %>%
  st_buffer(dist = 10) %>%    # 10m either side of river
  st_union() %>%
  st_transform(4326)

protectedareas_buf <- protectedareas_raw %>%
  st_transform(32737) %>%
  st_buffer(dist = 50) %>%    # 50m from PA boundaries
  st_union() %>%
  st_transform(4326)

# Combine all exclusion zones into one polygon
exclusion_sf <- st_union(roads_buf, buildings_buf) %>%
  st_union(water_buf) %>%
  st_union(waterways_buf) %>%
  st_union(protected_excl) %>%   
  st_make_valid()

cat("Exclusion mask built \n")

# Continue with rasterize step as before
excl_rast <- rasterize(
  vect(exclusion_sf),
  ref_grid,
  field      = 1,
  background = 0)
excl_rast <- mask(excl_rast, ll_vect)

writeRaster(excl_rast,
            "outputs/EmbuLL_exclusion_mask.tif",
            overwrite = TRUE,
            datatype  = "INT1U")

cat("Exclusion mask exported — check in QGIS before continuing\n")





# Combine all exclusion zones
cat("Building exclusion mask...\n")
exclusion_sf <- st_union(roads_buf, buildings_buf) %>%
  st_union(water_exclude) %>%
  st_make_valid()

# Rasterize exclusion mask
excl_rast <- rasterize(
  vect(exclusion_sf),
  ref_grid,
  field      = 1,
  background = 0)
excl_rast <- mask(excl_rast, ll_vect)

# Apply exclusion to confident cluster raster
cluster_rast_final <- ifel(excl_rast == 1, NA, cluster_rast_confident)
cluster_rast_final <- mask(cluster_rast_final, ll_vect)

writeRaster(cluster_rast_final,
            "outputs/EmbuLL_clusters_VVV7_final.tif",
            overwrite = TRUE,
            datatype  = "INT1U")

writeRaster(excl_rast,
            "outputs/EmbuLL_exclusion_mask.tif",
            overwrite = TRUE,
            datatype  = "INT1U")

cat("Exclusion mask exported\n")






#############  13. SAMPLING ALLOCATION ############# 
# 15 field days × 8 subplots/day = 120 subplots
# 4 subplots per plot → 28 plots total
# LDSF minimum: 3 plots per cluster
# Weighted proportionally by available (post-mask) pixel count

masked_df <- as.data.frame(cluster_rast_final, xy = TRUE, na.rm = TRUE)
names(masked_df)[3] <- "cluster"
masked_df$cluster <- as.factor(masked_df$cluster)

cat("\nPixels available per cluster after all masks:\n")
print(table(masked_df$cluster))

total_plots <- 28

cluster_alloc <- masked_df %>%
  count(cluster) %>%
  rename(N_pixels_available = n) %>%
  mutate(
    prop        = N_pixels_available / sum(N_pixels_available),
    alloc_prop  = round(prop * total_plots),
    alloc_final = pmax(alloc_prop, 3),   # LDSF minimum 3 per cluster
    subplots    = alloc_final * 4
  )

# Balance total to exactly 28 plots
while (sum(cluster_alloc$alloc_final) > total_plots) {
  idx <- which.max(cluster_alloc$alloc_final)
  cluster_alloc$alloc_final[idx] <- cluster_alloc$alloc_final[idx] - 1
}
while (sum(cluster_alloc$alloc_final) < total_plots) {
  idx <- which.max(cluster_alloc$prop)
  cluster_alloc$alloc_final[idx] <- cluster_alloc$alloc_final[idx] + 1
}
cluster_alloc$subplots <- cluster_alloc$alloc_final * 4

cat("\nFinal sampling allocation (VVV,7 — 28 plots):\n")
print(cluster_alloc[, c("cluster", "N_pixels_available",
                        "prop", "alloc_final", "subplots")])
cat("Total plots:   ", sum(cluster_alloc$alloc_final), "\n")
cat("Total subplots:", sum(cluster_alloc$subplots), "\n")

write_csv(cluster_alloc, "outputs/EmbuLL_allocation_VVV7.csv")







# ── RANDOM PLOT PLACEMENT WITH 250m MINIMUM SPACING ──────────
# spatSample() draws random points within each cluster mask
# Greedy distance filter enforces 250m minimum between plot centres
# Plot centres are the coordinates pre-loaded into QField GPS

set.seed(55)
all_plots <- list()

for (cl in levels(masked_df$cluster)) {
  
  n_plots <- cluster_alloc$alloc_final[cluster_alloc$cluster == cl]
  if (is.na(n_plots) || n_plots == 0) next
  
  # Isolate this cluster in the final masked raster
  cl_mask <- ifel(cluster_rast_final == as.integer(cl), 1, NA)
  
  # Oversample then filter by spacing
  pts <- spatSample(
    cl_mask,
    size   = n_plots * 15,   # oversample generously
    method = "random",
    na.rm  = TRUE,
    xy     = TRUE
  )
  
  if (nrow(pts) < 1) {
    cat("Warning: no valid pixels for cluster", cl, "\n")
    next
  }
  
  # Project to UTM zone 37S for distance in metres
  pts_sf <- st_as_sf(pts, coords = c("x", "y"), crs = 4326) %>%
    st_transform(32737)
  
  # Greedy selection: keep point only if >250m from all selected
  selected <- pts_sf[1, ]
  for (i in seq_len(nrow(pts_sf))) {
    if (nrow(selected) >= n_plots) break
    dists <- as.numeric(st_distance(pts_sf[i, ], selected))
    if (all(dists > 250)) {
      selected <- rbind(selected, pts_sf[i, ])
    }
  }
  
  if (nrow(selected) < n_plots) {
    cat("Warning: cluster", cl, "— placed", nrow(selected),
        "of", n_plots, "plots (area may be too small for spacing)\n")
  }
  
  selected$cluster <- cl
  selected$plot_id <- paste0("C", cl, "_P",
                             formatC(seq_len(nrow(selected)),
                                     width = 2, flag = "0"))
  all_plots[[cl]] <- selected
  cat("Cluster", cl, "→", nrow(selected), "plots placed\n")
}

plot_centres <- do.call(rbind, all_plots) %>%
  st_transform(4326)

cat("\nTotal plots placed:", nrow(plot_centres), "\n")




#############  14. GENERATE LDSF SUBPLOT COORDINATES ############# 
# LDSF layout: 1 centre subplot + 3 at 11.28m at 0°/120°/240°
# 11.28m = 2 × subplot radius (5.64m) so subplots do not overlap

generate_subplots <- function(centre_sf, plot_id, cluster_id) {
  c_utm  <- st_transform(centre_sf, 32737)
  coords <- st_coordinates(c_utm)
  
  angles <- c(0, 120, 240)
  dist   <- 11.28
  
  all_coords <- rbind(
    coords,
    cbind(
      coords[1] + dist * sin(angles * pi / 180),
      coords[2] + dist * cos(angles * pi / 180)
    )
  )
  
  st_as_sf(
    data.frame(
      plot_id    = plot_id,
      subplot_id = paste0(plot_id, "_S", 1:4),
      subplot_no = 1:4,
      cluster    = cluster_id,
      note       = c("Centre", "North (0°)", "SE (120°)", "SW (240°)"),
      radius_m   = 5.64,
      area_m2    = 100
    ),
    geometry = st_sfc(
      lapply(1:4, function(i) st_point(all_coords[i, ])),
      crs = 32737
    )
  ) %>% st_transform(4326)
}

all_subplots <- mapply(
  generate_subplots,
  centre_sf  = lapply(seq_len(nrow(plot_centres)),
                      function(i) plot_centres[i, ]),
  plot_id    = plot_centres$plot_id,
  cluster_id = plot_centres$cluster,
  SIMPLIFY   = FALSE
)

subplots_all <- do.call(rbind, all_subplots)
cat("Total subplots generated:", nrow(subplots_all), "\n")

#############  15. EXPORT FOR QGIS -> ODK LDSF ############# 
st_write(plot_centres,
         "outputs/EmbuLL_plot_centres_VVV7.gpkg",
         delete_dsn = TRUE)

st_write(subplots_all,
         "outputs/EmbuLL_subplots_VVV7.gpkg",
         delete_dsn = TRUE)

# CSV versions with explicit lat/lon for GPS devices
coords_p <- st_coordinates(plot_centres)
write_csv(
  plot_centres %>%
    st_drop_geometry() %>%
    mutate(longitude = coords_p[, 1], latitude = coords_p[, 2]),
  "outputs/EmbuLL_plot_centres_VVV7.csv"
)

coords_s <- st_coordinates(subplots_all)
write_csv(
  subplots_all %>%
    st_drop_geometry() %>%
    mutate(longitude = coords_s[, 1], latitude = coords_s[, 2]),
  "outputs/EmbuLL_subplots_VVV7.csv"
)

cat("\nAll outputs exported to outputs/\n")
cat("Files for QGIS/QField:\n")
cat("  EmbuLL_clusters_VVV7_final.tif   — cluster map\n")
cat("  EmbuLL_uncertainty_VVV7.tif      — uncertainty map\n")
cat("  EmbuLL_exclusion_mask.tif        — roads/buildings/water\n")
cat("  EmbuLL_plot_centres_VVV7.gpkg    — 28 plot centres\n")
cat("  EmbuLL_subplots_VVV7.gpkg        — 112 subplot points\n")

























######### 13. EXCLUDING HIGH-UNCERTAINTY PIXELS ######## 
# Only sample from pixels where model is confident (uncertainty < 0.2) - reducng missclassification risk
# is that threshold okay? 
sampling_frame <- env_df %>%
  filter(uncertainty < 0.20) %>%
  dplyr::select(x, y, cluster)

cat("Pixels in confident sampling frame:", nrow(sampling_frame), "\n")
# 1131 IS A LOT I guess?! 
cat("Pixels excluded (high uncertainty):", 
    nrow(env_df) - nrow(sampling_frame), "\n")
# 125

write_csv(sampling_frame, "outputs/Embu_sampling_frame_confident.csv")






#########  14. NEXT STEPs ? ################## 
# maybe using Stratified random sample in R or i do it in qgis 
# next steps when all stratification is correct, model well fitted and cluster chosen 
# 1. EmbuLL_clusters_mclust.tif in QGIS 
# 2. Overlaying with roads, villages to check accessability 
# 3. running the sampling design in R - allocate plots and subplots using weighted stratified sampling 
#4. then save as a .gpkg to send for ODK upload or QFIELD upload! 


# maybe in R
#sampleStratified(x, size, exp=10, na.rm=TRUE, xy=FALSE, ext=NULL, sp=FALSE, ...)
# The function may not work well when the size (number of cells) of some strata is relatively small.
