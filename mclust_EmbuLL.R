
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

# NOTE for self: THIS IS THE CORRECT VERSION connected to GitHub


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
#remember: less negative BIC is better model but  when differences between top models are < 2 BIC units,
#     prefer the simpler model (fewer clusters
set.seed(55)
n_pixels <- nrow(env_scaled)
cat("Total pixels for clustering:", n_pixels, "\n")
bic_result_55 <- mclustBIC(env_scaled, G = 5:10)
summary(bic_result_55) # VVV10 with BIC -1854266
plot(bic_result_55, main = "BIC — clusters for Embu LL (100m)")
# ALSO: The BIC values seem to be quite high compared to the literature? So not a good fit. Adding more variables could be better? 
# 
# 
# #### Trying other values 
# # Trying with 6 since 9 might be too much if i sample 3-4 weeks 
# set.seed(66)
# n_pixels <- nrow(env_scaled)
# bic_result_66 <- mclustBIC(env_scaled, G = 1:6)
# summary(bic_result_66)
# plot(bic_result, main = "BIC — clusters for Embu LL")
# # best result: VVV6 VEV 6 and VVV5 with a bic around -12029.64
# 
# # WHICH ONE IS A BETTER CHOICE? 
# # From the old run: I guess the smaller  bcs with only 109 pixels and 6 clusters i would only get 18pixel per cluster. 
# # Some of the smaller clusters maybe even have e.g. 8 pixel and then i cannot allocate "enough" field-clusters to this stratum bcs of 1km2 resolution 
# # then that areas would be undersampled i guess ?! 
# 
# 
# # testing if it will always increase or drop at some point
# set.seed(111)
# n_pixels <- nrow(env_scaled)
# bic_result_111 <- mclustBIC(env_scaled, G = 1:20)
# summary(bic_result_111)
# plot(bic_result, main = "BIC — clusters for Embu LL")
# # seems like 11 or even 20 clusters are fitting well with BIC -11093.19 and -11143.10046 respectively 

cat("Continuing with G=10 for now but VERY unsure if the model is working good enough!\n")





######################## 8. FIT FINAL MCLUST MODEL ######################## 
optimal_G <- 10

  model_name <- dimnames(bic_result_55)[[2]][
    which.max(apply(bic_result_55, 2, max, na.rm = TRUE))]
cat("Best covariance model:", model_name, "\n")
# fit on full pixel dataset not subsample used for BIC 
mc_model_55 <- Mclust(env_scaled, optimal_G, modelNames = "VVV")
summary(mc_model_55)


#cluster ID and uncertainty (prob. for not belonging to same cluster) for each pixel
env_df$cluster <- as.factor(mc_model_55$classification)
env_df$uncertainty <- mc_model_55$uncertainty
# Uncertainty = probability of NOT belonging to assigned cluster
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

#plot(mc_model_55, what = "classification")
#plot(mc_model_55, what = "uncertainty")








######################## 10. MAP CLUSTERS BACK TO RASTER ######################## 
# Convert env_df points to sf object using their coordinates
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




######### PLOTTING CLUSTER MAP #####################
display.brewer.all()
cls_colors <- brewer.pal(max(optimal_G, 10), "Set2")[1:optimal_G]

#cls_colors_6 <- c("#2166ac", "#f4a582", "#1a9641", "darkolivegreen4", "darkorange2", "khaki3")

plot(cluster_rast,
     col    = cls_colors,
     main   = paste0(optimal_G, " Environmental Clusters — Embu Living Lab"),
     legend = TRUE,
     axes   = TRUE)
plot(ll_vect, add = TRUE, border = "black", lwd = 2)

# Uncertainty map — white = confident, red = transition zone
plot(uncert_rast,
     col  = rev(heat.colors(50)),
     main = "Cluster Uncertainty — high values = transition zones")
plot(ll_vect, add = TRUE, border = "black", lwd = 2)
# White/pale areas = confident assignment (most pixels)
# Orange/red areas = ambiguous pixels between clusters
# These are NOT sampling errors — they are real landscape transitions










######### EXPORT FOR QGIS or ODK ########################### 
# Unsure about this follwing part 
dir.create("outputs/", showWarnings = FALSE)

writeRaster(cluster_rast,
            "outputs/LL_mclust_10.tif",
            overwrite = TRUE,
            datatype  = "INT1U")

writeRaster(uncert_rast,
            "outputs/LL_uncertainty_mclust_10.tif",
            overwrite = TRUE)

# maybe useful to export cluster assignments with coordinates as CSV
write_csv(env_df[, c("x", "y", "cluster", "uncertainty")],
          "outputs/Embu_cluster_pixels_9.csv")


######### 12. SAMPLING ALLOCATION PER CLUSTER ######### 
# Proportional allocation = weighted stratified sampling? 
# minimum of 2 plots per stratum? 
# ASSUMPTION
# 3 weeks sampling are 15 field days --> 8 subplots per day on average = 120 subplots in total = 30 plots in total 

cluster_alloc <- cluster_area %>%
  mutate(
    prop           = N_pixels / sum(N_pixels), # should i work with pixel as weighting factor? or size%?
    total_clusters = 30,
    alloc_prop     = round(prop * total_clusters),
    alloc_final    = pmax(alloc_prop, 2))   # minimum 2 per stratum? 
cat("\nSampling allocation per cluster:\n")
print(cluster_alloc)
# alloc_final
#1           5
#2           3
#3           5
#4           7
#5           7
#6           3
cat("Total plots:", sum(cluster_alloc$alloc_final), "\n")
cat("Total sub-plots (x4):", sum(cluster_alloc$alloc_final) * 4, "\n")

write_csv(cluster_alloc, "outputs/EmbuLL_cluster_allocation.csv")




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
