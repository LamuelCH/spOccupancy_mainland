library(terra)
setwd("RStudio/spOccupancy_NSW_20241219/")
options(mc.core = 10)
###########################################################################
# Create species range mask ---------------------------------------------------
library(terra)
library(dplyr)

# Define study area
aoi <- ext(-100000, 1225969, 5500000, 7400000)  # Define the extent

# Read the species IUCN range shapefile
range = vect("spatial_data/Range data/Dasyurus Maculatus/data_0.shp")
range = project(range, "EPSG:32755")   # reproject to Easting and Northing
range = crop(range,aoi)

# Create a raster mask with 1 km resolution
# Generate an empty raster with the defined extent and resolution
raster = rast(extent = aoi, resolution = 1000, crs = "EPSG:32755")
mask = rasterize(range, raster, field = 1)  # Field = 1 assigns a value of 1 to the rasterized area

# Save the raster mask (optional)
# writeRaster(mask, "RStudio/spOccupancy_NSW_20241219/input/mask_1km.tif", overwrite = TRUE)

# Plot the results (optional)
# plot(mask, main = "1 km Resolution Raster Mask")



# Recall the beta covariates and stacked into single raster
bio = rast("RStudio/spOccupancy_NSW_20241219/input/CHELSA_bio_2011-2040_gfdl-esm4_ssp126_V.2.1_EPSG32755.tif")
bio = bio[[c("bio5", "bio12")]] #we will only use bio5 and bio 12 based on pearson correlation analysis in S0
roads_density = rast("RStudio/spOccupancy_NSW_20241219/input/env_roadDensity_EPSG32755.tif")
pd = rast("RStudio/spOccupancy_NSW_20241219/input/aus_pd_2010-2020_1km_UNadj_EPSG32755.tif")
der = rast("RStudio/spOccupancy_NSW_20241219/input/env_der_EPSG32755.tif")

env_stack <- c(bio,roads_density, pd, der)
env_stack = env_stack * mask
env_stack = as.data.frame(env_stack, xy = TRUE, na.rm = TRUE) #convert to dataframe


bio5.pred = (env_stack$bio5 - mean(data.sfMsPGOcc$occ.covs[, 1])) / sd(data.sfMsPGOcc$occ.covs[, 1])
bio12.pred = (env_stack$bio12 - mean(data.sfMsPGOcc$occ.covs[, 2])) / sd(data.sfMsPGOcc$occ.covs[, 2])
roads.pred = (env_stack$roadLength - mean(data.sfMsPGOcc$occ.covs[, 3])) / sd(data.sfMsPGOcc$occ.covs[, 3])
pd.pred = (env_stack$pd_mean - mean(data.sfMsPGOcc$occ.covs[, 4])) / sd(data.sfMsPGOcc$occ.covs[, 4])
der.pred = (env_stack$der - mean(data.sfMsPGOcc$occ.covs[, 5])) / sd(data.sfMsPGOcc$occ.covs[, 5])

# Order: intercept, predictor (linear), predictor (quadratic)
X.0 <- cbind(1, 
             bio5.pred, bio5.pred^2, 
             bio12.pred, bio12.pred^2,
             roads.pred,
             pd.pred,
             der.pred, der.pred^2)

# Spatial coordinates
coords.0 = as.matrix(env_stack[, c("x", "y")])


# Approximate run time: 30 sec
out.pred <- predict(out.sfMsPGOcc, X.0, coords.0, type = 'occupancy')
str(out.pred)