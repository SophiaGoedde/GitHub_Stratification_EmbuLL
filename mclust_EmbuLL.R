
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
library(factoextra) 
library(RColorBrewer)
install.packages("spsurvey")
library(spsurvey)


######################## 1. LOAD + Crop RASTERS to LL Boundary ###################
###LL 
ll_boundary <- st_read("data/shapefiles/One_LL_Kenya_V1.shp")
# Check CRS —  match rasters (EPSG:4326)
print(st_crs(ll_boundary))
# Convert to terra SpatVector for masking
ll_vect <- vect(ll_boundary)
# plot 
plot(ll_boundary$geometry, main = "Embu Living Lab boundary")

##' RASTERS 
data_path <- "data/rasters/"

# i forgot to clip the rasters before to the LL so using a function and then clipping 
load_clip <- function(filename)
  {r <- rast(paste0(data_path, filename))
  r <- crop(r, ll_vect)
  r <- mask(r, ll_vect)
  return(r)}

r_temp      <- load_clip("Temperature_MODIS_2000_2025.tif")
r_rainfall  <- load_clip("Rainfall_CHIRPS_2000_2025.tif")
r_elevation <- load_clip("Elevation_SRTM.tif")
r_slope     <- load_clip("Slope_SRTM.tif")
r_aspect    <- load_clip("Aspect_SRTM.tif")
r_landcover <- load_clip("Landcover_ESA.tif") # IMPORTANT: categroical
# Variables tried but excluded: 
# r_distwater  <- load_clip("Distance_to_Water_JRC.tif")
#   → dropped: zero variance across LL at 1km (all pixels = 0)
# r_canopy     <- load_clip("CanopyHeight_ETH_2020_m.tif")
#   → dropped: r=0.75 with Rainfall_MAP (collinear); kept rainfall instead
# r_treecover  <- load_clip("TreeCover_Hansen_pct_2024.tif")
#   → dropped: Hansen dataset did not load correctly for this area;
#     to be replaced with a more suitable AFS tree cover dataset?!?!? ESpecially Tree cover dataset 

cat("All rasters loaded and clipped to Living Lab boundary\n")

#  checking
print(r_rainfall)
print(r_temp)
#print(r_distwater)

######################## 2. RESAMPLE ALL TO 1km GRID ######################## 
# 1km = MODIS resolution is taken as reference 
# note that CHIRPS was upsampled and the rest downsampled! Is that methodologically okay? 

#  temperature as the reference grid since it is 1km

ref_grid <- r_temp   # reference raster

# Continuous variables: bilinear interpolation
r_rainfall_r  <- resample(r_rainfall,  ref_grid, method = "bilinear")
r_elevation_r <- resample(r_elevation, ref_grid, method = "bilinear")
r_slope_r     <- resample(r_slope,     ref_grid, method = "bilinear")
r_aspect_r    <- resample(r_aspect,    ref_grid, method = "bilinear")
r_landcover_r <- resample(r_landcover, ref_grid, method = "near") # using n nearest neighbor 
#r_distwater_r <- resample(r_distwater, ref_grid, method = "bilinear")
#r_canopy_r    <- resample(r_canopy,    ref_grid, method = "bilinear")

cat("All layers resampled to 1km reference grid\n")

######################## 3. STACK ALL LAYERS (continous ones only)######################## 
# Land cover excluded from mclust input (categorical) bcs mclust assumed gaussian distribution 
# landcover used later for post-hoc interpretation/validation  - is that approach methodologically correct??

env_stack <- c(
  r_temp,
  r_rainfall_r,
  r_elevation_r,
  r_slope_r,
  r_aspect_r)
  #r_distwater_r,
  #r_canopy_r)
# or should i encode landcover as something else and tyr to include it within the clustering? is that possible?

# Renaming bands more simple
names(env_stack) <- c(
  "Temp_mean",
  "Rainfall_MAP",
  "Elevation",
  "Slope",
  "Aspect")
  #"Dist_water",
  #"Canopy_height"
print(env_stack)
plot(env_stack) 

cat("Stack bands:", names(env_stack), "\n")



######################## 4. EXTRACT PIXEL VALUES TO DATA FRAME ######################## 
env_df <- as.data.frame(env_stack, xy = TRUE, na.rm = TRUE)
# i am not entrely sure if this is correct?! 
# xy = TRUE keeps the coordinates 
nrow(env_df)   # how many valid pixels: 1256. --> enough? 
cat("Total valid pixels in Living Lab:", nrow(env_df), "\n") # 1256
# is that number okay? 




######################## 5. CHECK CORRELATIONS ######################## 
# highly correlated vars will be removed - doublecheck again if i add more variables!! 
vars_to_use <- c("Temp_mean", "Rainfall_MAP", "Elevation", "Slope", "Aspect")
cor_matrix <- cor(env_df[, vars_to_use],
                  use    = "complete.obs",
                  method = "pearson")
corrplot(cor_matrix,
         method      = "color",
         type        = "upper",
         tl.col      = "black",
         addCoef.col = "black",
         number.cex  = 0.8,
         title       = "Correlationmatrix for the Embu LL concerning the five environmental variables",
         mar         = c(0, 0, 2, 0))
# RULE: if |r| > 0.7 between two variables, consider dropping one
# for now they all look fine - none has to be removed but when adding new ones probably yes

# is this code chunk correct? 
high_cor <- which(abs(cor_matrix) > 0.7 & cor_matrix != 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("\n Highly correlated pairs (|r| > 0.7):\n")
  print(data.frame(
    Var1 = rownames(cor_matrix)[high_cor[, 1]],
    Var2 = colnames(cor_matrix)[high_cor[, 2]],
    r    = round(cor_matrix[high_cor], 3)))
  }else {cat("No pairs with |r| > 0.7 — all variables retained\n")}
# in the previous run canopy height and rainfall were correlated with r=0.75 
# in this case - which variable to drop?? 




######################## 6. SCALE VARIABLES ######################## 
# mclust is sensitive to variable scale  -->  mean=0 and SD=1
env_scaled <- scale(env_df[, vars_to_use])
summary(env_scaled)




######################## 7. DETERMINE OPTIMAL NUMBER OF CLUSTERS ######################## 
# mclust uses BIC (Bayesian Information Criterion) to find best fitting model and select clusters
#remember: less negative BIC is better model 
set.seed(99)
n_pixels <- nrow(env_scaled)
bic_result_99 <- mclustBIC(env_scaled, G = 1:9)
summary(bic_result_99)
plot(bic_result, main = "BIC — clusters for Embu LL")
# Best result second run: VVV9 or VVV8 or VEV 8 with BIC around -11230.31. 
# 9 clusters seems to be a lot?!?
# Best BIC values in first run were the following (with 109 pixels):  EVE,3 or VVI,6 or  EVI,5

# ALSO: The BIC values seem to be quite high compared to the literature? So not a good fit. Adding more variables could be better? 

# Trying with 6 since 9 might be too much if i sample 3-4 weeks 
set.seed(66)
n_pixels <- nrow(env_scaled)
bic_result_66 <- mclustBIC(env_scaled, G = 1:6)
summary(bic_result_66)
plot(bic_result, main = "BIC — clusters for Embu LL")
# best result: VVV6 VEV 6 and VVV5 with a bic around -12029.64

# WHICH ONE IS A BETTER CHOICE? 
# From the old run: I guess the smaller  bcs with only 109 pixels and 6 clusters i would only get 18pixel per cluster. 
# Some of the smaller clusters maybe even have e.g. 8 pixel and then i cannot allocate "enough" field-clusters to this stratum bcs of 1km2 resolution 
# then that areas would be undersampled i guess ?! 

# trying BIC value with five clusters 
set.seed(55)
n_pixels <- nrow(env_scaled)
bic_result_55 <- mclustBIC(env_scaled, G = 1:5)
summary(bic_result_55)
plot(bic_result, main = "BIC — clusters for Embu LL")
# VVV5 with BIC -12412.38 
# with three it was VVV3 with BIC -12906.36 - not a big difference 

# testing if it will always increase or drop at some point
set.seed(111)
n_pixels <- nrow(env_scaled)
bic_result_111 <- mclustBIC(env_scaled, G = 1:20)
summary(bic_result_111)
plot(bic_result, main = "BIC — clusters for Embu LL")
# seems like 11 or even 20 clusters are fitting well with BIC -11093.19 and -11143.10046 respectively 

cat("Continuing with G=6 for now but VERY unsure!\n")





######################## 8. FIT FINAL MCLUST MODEL ######################## 
mc_model_99 <- Mclust(env_scaled, 9, modelNames = "VVV")
summary(mc_model_99)

#VVV with G=6
mc_model_66 <- Mclust(env_scaled, 6, modelNames = "VVV")
summary(mc_model_66)

#cluster ID and uncertainty (prob. for not belonging to same cluster) for each pixel
env_df$cluster <- as.factor(mc_model_66$classification)
env_df$uncertainty <- mc_model_66$uncertainty
cat("Pixels per cluster:\n")
print(table(env_df$cluster))





######################## 9. VISUALIZE + INTERPRETE CLUSTERS ######################## 
# Cluster means to interpret strata
cluster_means <- aggregate(
  env_df[, vars_to_use],
  by  = list(Cluster = env_df$cluster),
  FUN = mean)
print(round(cluster_means[, -1], 2))

# Area per cluster
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
    subtitle = "Embu Living Lab — EVE,3 mclust solution",
    x = "", y = "Cluster")
ggsave("outputs/cluster_means_heatmap.png", width = 8, height = 4, dpi = 300)

# Classification plot in variable space
plot(mc_model_66, what = "classification")

# Uncertainty plot
plot(mc_model_66, what = "uncertainty")
# it seems like the uncertainty is VERY high 


######################## 10. MAP CLUSTERS BACK TO RASTER ######################## 
# Convert env_df points to sf object using their coordinates
cluster_points <- st_as_sf(
  env_df[, c("x", "y", "cluster", "uncertainty")],
  coords = c("x", "y"),
  crs    = 4326)
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
cat("Unique cluster values:", unique(values(cluster_rast), na.rm = TRUE), "\n")
cat("Non-NA pixels:", sum(!is.na(values(cluster_rast))),
    "— expected:", nrow(env_df), "\n")




######### PLOTTING CLUSTER MAP #####################
display.brewer.all()

cls_colors <- c("#2166ac", "#f4a582", "#1a9641", "darkolivegreen4", "darkorange2", "khaki3")

plot(cluster_rast,
     col     = cls_colors,
     main    = "Six Clusters based on Environmental Variables
     — Embu Living Lab",
     legend  = TRUE,
     axes    = TRUE)
plot(ll_vect, add = TRUE, border = "black", lwd = 2)
# what are the white areas? NA where there. was no data availeble?? 


#plot uncertainty map 
plot(uncert_rast,
     col  = rev(heat.colors(50)),
     main = "Uncertainty for clustering 
     (six clusters based on environmental variables only)")
plot(ll_vect, add = TRUE, border = "black", lwd = 2)






######### EXPORT FOR QGIS or ODK ########################### 
# Unsure about this follwing part 
# exporting as GeoTIFF i assume? 
dir.create("outputs/", showWarnings = FALSE)

writeRaster(cluster_rast,
            "outputs/EmbuLL_clusters_mclust.tif",
            overwrite = TRUE,
            datatype  = "INT1U")

writeRaster(uncert_rast,
            "outputs/EmbuLL_uncertainty_mclust.tif",
            overwrite = TRUE)

# maybe useful to export cluster assignments with coordinates as CSV
write_csv(env_df[, c("x", "y", "cluster", "uncertainty")],
          "outputs/Embu_cluster_pixels.csv")


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
