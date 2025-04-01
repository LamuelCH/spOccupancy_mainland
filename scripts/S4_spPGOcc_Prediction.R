# library(terra)
# library(dplyr)
library(spOccupancy)

# Recall the beta covariates and stacked into single raster
bio = rast("input/raster/CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
bio = bio[[c("bio5","bio6", "bio12", "bio15")]] #we will only use bio5 and bio 12 based on pearson correlation analysis in S0

terrain = rast("input/env_terrain_EPSG3577.tif")
env_stack <- c(bio,terrain)
env_stack = as.data.frame(env_stack, xy = TRUE, na.rm = TRUE) #convert to dataframe

load("models/22222/spPGOcc/model_spPGOcc_1981-2010_nthin500_nbatch1e+05_nchain4_nburn1500000.RData")
load("input/list_sfMsPGOcc_1980-2010.RData")


bio5.pred = (env_stack$bio5 - mean(data.sfMsPGOcc$occ.covs[, 1])) / sd(data.sfMsPGOcc$occ.covs[, 1])
bio6.pred = (env_stack$bio6 - mean(data.sfMsPGOcc$occ.covs[, 2])) / sd(data.sfMsPGOcc$occ.covs[, 2])
bio12.pred = (env_stack$bio12 - mean(data.sfMsPGOcc$occ.covs[, 3])) / sd(data.sfMsPGOcc$occ.covs[, 3])
bio15.pred = (env_stack$bio15 - mean(data.sfMsPGOcc$occ.covs[, 4])) / sd(data.sfMsPGOcc$occ.covs[, 4])
dem.pred = (env_stack$dem - mean(data.sfMsPGOcc$occ.covs[, 5])) / sd(data.sfMsPGOcc$occ.covs[, 5])
slope.pred = (env_stack$slope - mean(data.sfMsPGOcc$occ.covs[, 6])) / sd(data.sfMsPGOcc$occ.covs[, 6])
aspect.pred = (env_stack$aspect - mean(data.sfMsPGOcc$occ.covs[, 7])) / sd(data.sfMsPGOcc$occ.covs[, 7])
tpi.pred = (env_stack$tpi - mean(data.sfMsPGOcc$occ.covs[, 8])) / sd(data.sfMsPGOcc$occ.covs[, 8])
tri.pred = (env_stack$tri - mean(data.sfMsPGOcc$occ.covs[, 9])) / sd(data.sfMsPGOcc$occ.covs[, 9])
rock.pred = (env_stack$rock - mean(data.sfMsPGOcc$occ.covs[, 10])) / sd(data.sfMsPGOcc$occ.covs[, 10])

# Order: intercept, elevation (linear), elevation (quadratic)
X.0 <- cbind(1, bio5.pred, bio5.pred^2,
             bio6.pred, bio6.pred^2,
             bio12.pred, bio12.pred^2,
             bio15.pred, bio15.pred^2,
             dem.pred, dem.pred^2,
             slope.pred, slope.pred^2,
             aspect.pred,
             tpi.pred, tpi.pred^2,
             tri.pred, tri.pred^2,
             rock.pred) 

# Spatial coordinates
coords.0 <- as.matrix(env_stack[, c('x', 'y')])

# type = 'occupancy' specified prediction of occupancy (or occurrence). 
# This is also the default.
# Approximate run time: 30 sec
out.pred <- predict(out.spPGOcc, 
                    X.0,  
                    coords.0,
                    n.omp.threads = 192,
                    verbose = TRUE,
                    n.report = 500,
                    type = 'occupancy')

save(out.pred, 
     file = paste0("output/spPGOcc/model_spPGOcc_22222_1981-2010.RData"))


str(out.pred)


# Species richness samples
# rich.pred <- apply(out.pred$z.0.samples, c(1, 3), sum)
# plot.dat <- data.frame(x = hbefElev$Easting, 
                       # y = hbefElev$Northing, 
                       # rich.mean = apply(rich.pred, 2, mean), 
                       # rich.sd = apply(rich.pred, 2, sd))
# Plot species richness of the foliage-gleaning bird community
# across the Hubbard Brook Experimental Forest
# dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
# ggplot() + 
  # geom_stars(data = dat.stars, aes(x = x, y = y, fill = rich.mean)) +
  # scale_fill_viridis_c(na.value = 'transparent') +
  # labs(x = 'Easting', y = 'Northing', fill = '', 
       # title = 'Mean Species Richness') +
  # theme_bw()




# Order: intercept, predictor (linear), predictor (quadratic)
# X.0 <- cbind(1, 
             # bio5.pred, bio5.pred^2, 
             # bio12.pred, bio12.pred^2,
             # roads.pred,
             # pd.pred,
             # der.pred, der.pred^2)

# Spatial coordinates
# coords.0 = as.matrix(env_stack[, c("x", "y")])


# Approximate run time: 30 sec
# out.pred <- predict(out.sfMsPGOcc, X.0, coords.0, type = 'occupancy')
# str(out.pred)