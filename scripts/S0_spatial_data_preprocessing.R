library(raster)
library(magrittr)
library(terra)
library(slga)
library(ggplot2)
library(sf)
library(geodata)
library(tidyr)

###########################################################################
# Create species range mask ---------------------------------------------------
library(terra)
library(dplyr)

options(mc.core = parallel::detectCores())

# Define study area
aoi_3577 = ext(684147.3, 2225320.5, -4597511.6, -2556437.4)
aoi_32755 <- ext(-100000, 1225969, 5500000, 7400000)
# aoi_polygon <- st_as_sfc(st_bbox(c(xmin = -100000, xmax = 1225969, ymin = 5500000, ymax = 7400000), crs = 32755))

# Step 3: Transform the polygon to EPSG:3577
# aoi_polygon_3577 <- st_transform(aoi_polygon, crs = 3577)

# Step 4: Extract the transformed extent
# aoi_3577 <- st_bbox(aoi_polygon_3577)
# Print the transformed extent
# print(aoi_3577)


# Read the species IUCN range shapefile
# range = vect("spatial_data/Range data/Dasyurus Maculatus/data_0.shp")
# range = project(range, "EPSG:32755")   # reproject to Easting and Northing
# range = crop(range,aoi)

# Create a raster mask with 1 km resolution
# Generate an empty raster with the defined extent and resolution
raster_3577 = rast(extent = aoi_3577, resolution = 1000, crs = "EPSG:3577")
raster_32755 = rast(extent = aoi_32755, resolution = 1000, crs = "EPSG:32755")

mask = rast("input/env_mask_EPSG3577.tif")


# Save the raster mask (optional)
# writeRaster(mask, "RStudio/spOccupancy_NSW_20241219/input/mask_1km.tif", overwrite = TRUE)

# Plot the results (optional)
# plot(mask, main = "1 km Resolution Raster Mask")


###########################################################################
# Raster processing -------------------------------------------------------
# CHELSA bioclimatic data -------------------------------------------------
# Define models and SSPs
models <- c("gfdl-esm4", "mpi-esm1-2-hr", "ukesm1-0-ll", "mri-esm2-0", "ipsl-cm6a-lr")
ssps <- c("ssp126", "ssp370", "ssp585")

# Define directories
base_input_dir <- "data/beta/CHELSA/climatologies/2011-2040/"
output_dir <- "input/raster/"


# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each model and SSP
for (model in models) {
  for (ssp in ssps) {
    
    # Construct input directory path
    input_dir <- file.path(base_input_dir, model, ssp, "bio")
    
    # List all raster files in the directory
    bio_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)
    
    # Proceed only if there are raster files
    if (length(bio_files) > 0) {
      
      # Extract layer names from file names
      layer_names <- sub(".*CHELSA_([^_]+)_.*", "\\1", basename(bio_files))
      
      # Load rasters and stack them
      bio <- rast(bio_files)
      
      # Assign names to layers
      names(bio) <- layer_names
      
      # Reproject the raster using bilinear interpolation
      bio <- terra::project(bio, raster_3577, method = "bilinear") *mask  # Ensure `raster_3577` is defined
      
      # Save the processed raster stack
      output_file <- file.path(output_dir, paste0("CHELSA_bio_2011-2040_", model, "_", ssp, "_V.2.1_EPSG3577.tif"))
      writeRaster(bio, output_file, overwrite = TRUE)
      
      print(paste("Processed and saved:", output_file))
      
    } else {
      warning(paste("No raster files found for", model, ssp))
    }
  }
}


# 1981-2010 CHELSA data
## Define the directory containing the CHELSA rasters
input_dir <- "data/beta/CHELSA/climatologies/1981-2010/bio"
output_dir <- "input/raster"

# List all raster files in the directory
bio <- list.files(input_dir, pattern = "^CHELSA_bio", full.names = TRUE)

# Extract layer names from file names
layer_names <- sub(".*CHELSA_([^_]+)_.*", "\\1", basename(bio))

bio = rast(bio)

# Assign names to the layers in the stack
names(bio) <- layer_names

# Reproject the raster using bilinear method
bio <- terra::project(bio, raster_3577, method = "bilinear")*mask

# Save the stacked raster to a file
output_file = file.path(output_dir, "CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
writeRaster(bio,output_file, overwrite = TRUE)

###########################################################################
# Static layers -----------------------------------------------------------
# DIGITAL ELEVATION MODEL
dem = rast("data/beta/3secSRTM_DEM/DEM_ESRI_GRID_16bit_Integer/dem3s_int/hdr.adf")
dem = terra::project(dem, raster_3577, method = "bilinear")
names(dem) = "dem"

slope = terrain(dem, "slope", unit = "degrees")
aspect = terrain(dem, "aspect", unit="degrees")

# Window size determines the scale of landforms identified
tpi <- terra::focal(dem, w=matrix(1,9,9), 
                    fun=function(x) x[5] - mean(x, na.rm=TRUE))
names(tpi) = "tpi"

# Calculate TRI
tri <- terrain(dem, "TPI")
names(tri) = "tri"

# writeRaster(terrain, 
            # "input/env_terrain_EPSG3577.tif",
            # overwrite=TRUE)

# DEPTH OF REGOLITH ####
# der = rast("data/beta/Soil and Landscape Grid National Soil Attribute Maps/regolith/data/DER_000_999_EV_N_P_AU_NAT_C_20150601.tif")
# der <- terra::project(der, raster_3577, method = "bilinear")
# names(der) = "der"
# writeRaster(der, filename = "input/env_der_EPSG3577.tif", overwrite = TRUE)


# ROCK OUTCROP OCCURRENCE
rock = rast("data/beta/Soil and Landscape Grid National Soil Attribute Maps/rockoutcrop/data/DES_rockoutcrop_N_P_AU_NAT_C_20190901.tif")
rock <- terra::project(rock, raster_3577, method = "bilinear")

names(rock) = "rock"

terrain = c(dem, slope, aspect, tpi, tri, rock)*mask

writeRaster(terrain, filename = "input/env_terrain_EPSG3577.tif", overwrite = TRUE)



# Pearson Correlation Analysis --------------------------------------------
###########################################################################
library(ggplot2)
library(corrplot)

# Bioclimatic variables
bio = rast("input/raster/1981-2010/CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
corr_bio = layerCor(bio, "pearson", na.rm = TRUE)

bio = bio[c("bio5", "bio6")]

# terrain variables
terrain = rast("input/env_terrain_EPSG3577.tif")
corr_terrain = layerCor(terrain, "pearson", na.rm = TRUE)

# plot resutls 
corrplot(corr_terrain$correlation, 
         method = "circle", 
         type = "lower", 
         insig = "blank", 
         addCoef.col = "black", 
         number.cex = 0.6)

corrplot(corr_bio$correlation, 
         method = "circle", 
         type = "lower", 
         insig = "blank", 
         addCoef.col = "black", 
         number.cex = 0.6)

# Assuming your correlation matrix is stored as a matrix or data.frame
# For example: correlation_matrix <- layerCor(your_raster, "pearson", na.rm = TRUE)$pearson

# Convert correlation matrix to a dissimilarity matrix (1 - |correlation|)
# dissimilarity <- as.dist(1 - abs(correlation_matrix$correlation))

# Perform hierarchical clustering
# hc <- hclust(dissimilarity, method = "average")  # "average" or "complete" linkage works well

# Cut the dendrogram to group variables with |r| > 0.7
# groups <- cutree(hc, h = 0.3)  # h = 1 - 0.7 = 0.3 threshold

# View the groups
# correlation_groups <- split(names(groups), groups)
# print(correlation_groups)

# corrplot(correlation_matrix$correlation, method = "square", order = "hclust", hclust.method = "average", 
         # addrect = length(correlation_groups), tl.cex = 0.7)

env_stack = c(subset(bio, c("bio5", "bio6", "bio12", "bio15")), terrain)
correlation_matrix <- layerCor(env_stack, "pearson", na.rm = TRUE)
corrplot(correlation_matrix$correlation, method = "circle", type = "lower", insig = "blank", addCoef.col = "black", number.cex = 0.6)



# Add all your raster layers here

# Repeat this for all your raster layers
###########################################################################


















###########################################################################
# FOREST COVER##########################################################################
#In this section we create the forest cover layer. Because there is obvious pixel bandings and quality issue for 2005, 2010 and 2015 layer,
#we will only use 2000 forest cover 
# WE COMBINE ALL FOREST COVER FROM 2000-2015 AND TAKE THE MEAN VALUE TO EACH PIXEL
# Get a list of all subdirectories ending with "TC_2010"
subdirs <- list.dirs("data/RawData.RAW/Forest Cover/", 
                     recursive = TRUE, 
                     full.names = TRUE)
subdirs <- subdirs[grep("TC_2000$", subdirs)]

# Initialize an empty list to store raster files
GFCC_list <- list()

# Loop through each subdirectory
for (subdir in subdirs) {
  # Get a list of raster files ending with "TC_2010.tif" in this subdirectory
  files <- list.files(path = subdir, pattern = "\\.tif$", full.names = TRUE)
  files <- files[!grepl("_err\\.tif$|_idx\\.tif$|_idx\\.txt$", files)]
  
  # Read the rasters into a list of SpatRaster objects
  GFCC_list <- c(GFCC_list, lapply(files, terra::rast))
}

# REPROJET TO SAME CRS AND RESAMPLE TO 0.01 DEGREE
reprojected_GFCC <- lapply(GFCC_list, function(r) project(r, 
                                                          raster_1km, 
                                                          method = "bilinear")) %>% 
  sprc() %>% 
  mosaic() #%>% 

mask(., env_auMask_1km)

# Filter out FOREST COVER PERCENTAGE > 100%
forestCover_1km <- ifel(reprojected_GFCC > 100,
                        NA, 
                        reprojected_GFCC)

names(forestCover_1km) <- "forestCover"

plot(forestCover_1km)

writeRaster(x = forestCover_1km, filename = "data/env_forestCover_1km.tif", overwrite = TRUE)

###########################################################################
# FOLIAGE PROJECTIVE COVER ####
foliage <- rast("data/raw/raster/Foliage Projective Cover/Woody vegetation cover - Landsat, JRSRP, Australian coverage, 2000-2010/lztmre_aus_y20002011_dm7a2_d20050630_r500m.tif")

#resample to 100m resolution
foliage_1km <- terra::project(foliage, 
                              env_tasMask_1km, 
                              method = "bilinear") * env_tasMask_1km
foliage_1km = clamp((foliage_1km - 100) / (200 - 100) * 100, 0, 100)

names(foliage_1km) = "foliageCover"

plot(foliage_1km)

writeRaster(x = foliage_1km, filename = "data/env_foliage_1km.tif", overwrite = TRUE)

###########################################################################
# Disturbance Layers ------------------------------------------------------
# ROAD DENSITY####
#Road Density
roads <- vect("data/beta/road/2024_12_NationalRoads/NationalRoads2024_12.shp") 
roads = project(roads, y = "EPSG:3577")

# Convert the roads to raster with 1km resolution
roads <- terra::rasterize(roads,raster_3577, field = "Shape_Leng", fun = sum)
# roads = crop(roads, raster) 

# Calculate the total length of road within a 10 km moving window for each cell
roads_density <- terra::focal(roads, 
                              w = 9, #9 cells roughly equal to 10 km
                              fun = sum,
                              na.rm = TRUE)

# Replace NA values in roads_10km_1km with 0
roads_density[is.na(roads_density)] <- 0


# Mask roads_10km_1km with the mask raster
roads_density <- roads_density

plot(roads_density)

names(roads_density) = "roadLength"

writeRaster(x = roads_density,filename = "input/env_roadDensity_EPSG3577.tif", overwrite = TRUE)
###########################################################################
# Population Density ------------------------------------------------------
input_dir <- "data/beta/pd/2010-2020/"
output_dir <- "input/"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define years to process
years <- 2010:2020

# Loop through each year and process the corresponding file
for (year in years) {
  # Construct file name
  input_file <- file.path(input_dir, paste0("aus_pd_", year, "_1km_UNadj.tif"))
  
  # Check if the file exists
  if (file.exists(input_file)) {
    # Load the raster
    pd <- rast(input_file)
    
    # Reproject the raster using bilinear method
    pd_reprojected <- terra::project(pd, raster_3577, method = "bilinear")
    
    # Save the reprojected raster to the output directory
    output_file <- file.path(output_dir, paste0("aus_pd_", year, "_1km_UNadj_EPSG3577.tif"))
    writeRaster(pd_reprojected, output_file, overwrite = TRUE)
    
    # Log the progress
    cat("Processed year:", year, "-> Saved to:", output_file, "\n")
  } else {
    # Log if file is missing
    cat("File missing for year:", year, "\n")
  }
}

cat("Processing complete!")

# Create disturbance index

# DATA READY FOR ANALYSIS


## Define the directory containing the CHELSA rasters
## Define the directory containing the CHELSA rasters
input_dir <- "data/beta/pd/2010-2020/"
output_dir <- "input/"

# List all raster files in the directory
pd <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

pd = rast(pd)
pd = mean(pd)


# Reproject the raster using bilinear method
pd <- terra::project(pd, raster_3577, method = "bilinear")

# Calculate the total density of population within a 10 km moving window for each cell
pd_smooth <- terra::focal(pd, w = 9, fun = sum, na.rm = TRUE)


# Assign names to the layers in the stack
names(pd_smooth) <- "pd_mean"


# Save the stacked raster to a file
output_file = file.path(output_dir, "aus_pd_2010-2020_1km_UNadj_EPSG3577.tif")
writeRaster(pd_smooth,output_file, overwrite = TRUE)



input_dir <- "data/beta/pd/2000-2009/"
output_dir <- "input/raster"

# List all raster files in the directory
pd <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

pd = rast(pd)
pd = mean(pd)


# Reproject the raster using bilinear method
pd <- terra::project(pd, raster_3577, method = "bilinear")

# Calculate the total density of population within a 10 km moving window for each cell
pd_smooth <- terra::focal(pd, w = 9, fun = sum, na.rm = TRUE)


# Assign names to the layers in the stack
names(pd_smooth) <- "pd_mean"


# Save the stacked raster to a file
output_file = file.path(output_dir, "aus_pd_1980-2000_1km_UNadj_EPSG3577.tif")
writeRaster(pd_smooth,output_file, overwrite = TRUE)



