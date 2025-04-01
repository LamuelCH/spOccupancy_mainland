library(terra)
library(dplyr)
library(spOccupancy)

# Create output directories
dir.create("output/22222/sfMsPGOcc/chunks", recursive = TRUE, showWarnings = FALSE)

# Define chunk size for processing
chunk_size <- 1000

# Load the model and data
cat("Loading model and data...\n")
load("models/22222/sfMsPGOcc/model_sfMsPGOcc_1981-2010_nthin50_nbatch10000_nchain4_nburn150000.RData")
load("input/list_sfMsPGOcc_1980-2010.RData")

# Determine number of species
if (!is.null(out.sfMsPGOcc$y)) {
  n_species <- dim(out.sfMsPGOcc$y)[1]
  cat("Detected", n_species, "species in model\n")
  cat("Species names:", rownames(out.sfMsPGOcc$y), "\n")
} else {
  cat("Could not determine species names from model object\n")
  n_species <- 15  # Set manually if not available from model
  cat("Assuming", n_species, "species\n")
}

# Load and prepare environmental data
cat("Loading environmental data...\n")
bio <- rast("input/raster/CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
bio <- bio[[c("bio5","bio6", "bio12", "bio15")]]
terrain <- rast("input/env_terrain_EPSG3577.tif")
env_stack <- c(bio, terrain)

# Convert to data frame, removing NA cells
cat("Converting raster to data frame...\n")
env_stack_df <- as.data.frame(env_stack, xy = TRUE, na.rm = TRUE)
cat("Total prediction locations:", nrow(env_stack_df), "\n")

# Extract means and SDs from original data for standardization
orig_means <- sapply(1:ncol(data.sfMsPGOcc$occ.covs), function(i) mean(data.sfMsPGOcc$occ.covs[, i]))
orig_sds <- sapply(1:ncol(data.sfMsPGOcc$occ.covs), function(i) sd(data.sfMsPGOcc$occ.covs[, i]))

# Standardize variables
cat("Standardizing variables...\n")
env_vars <- c("bio5", "bio6", "bio12", "bio15", "dem", "slope", "aspect", "tpi", "tri", "rock")
std_vars <- paste0(env_vars, "_std")

for (i in 1:length(env_vars)) {
  env_stack_df[[std_vars[i]]] <- (env_stack_df[[env_vars[i]]] - orig_means[i]) / orig_sds[i]
}

# Check for and handle NaN and Inf values
for (var in std_vars) {
  env_stack_df[[var]][is.nan(env_stack_df[[var]]) | is.infinite(env_stack_df[[var]])] <- NA
}

# Create design matrix
cat("Creating design matrix...\n")
X.0 <- cbind(
  1,  # Intercept
  env_stack_df$bio5_std, env_stack_df$bio5_std^2,
  env_stack_df$bio6_std, env_stack_df$bio6_std^2,
  env_stack_df$bio12_std, env_stack_df$bio12_std^2,
  env_stack_df$bio15_std, env_stack_df$bio15_std^2,
  env_stack_df$dem_std, env_stack_df$dem_std^2,
  env_stack_df$slope_std, env_stack_df$slope_std^2,
  env_stack_df$aspect_std,
  env_stack_df$tpi_std, env_stack_df$tpi_std^2,
  env_stack_df$tri_std, env_stack_df$tri_std^2,
  env_stack_df$rock_std
)

# Prepare coordinates
coords.0 <- as.matrix(env_stack_df[, c('x', 'y')])

# Remove rows with NA values
cat("Checking for and removing NA values...\n")
na_rows <- rowSums(is.na(X.0)) > 0
cat("Found", sum(na_rows), "rows with NA values\n")

if (sum(na_rows) > 0) {
  X.0 <- X.0[!na_rows, ]
  coords.0 <- coords.0[!na_rows, ]
  cat("After removing NA rows:", nrow(X.0), "locations remain\n")
}

# Determine number of chunks
n_locations <- nrow(X.0)
n_chunks <- ceiling(n_locations / chunk_size)
cat("Will process data in", n_chunks, "chunks\n")

# Initialize results data frame (empty tidy format)
results_summary <- data.frame(
  x = numeric(0),
  y = numeric(0),
  species = integer(0),
  psi_mean = numeric(0),
  psi_sd = numeric(0),
  z_prob = numeric(0),
  w_mean = numeric(0),
  w_sd = numeric(0)
)

# Process data in chunks
for (i in 1:n_chunks) {
  cat("\n====== Processing chunk", i, "of", n_chunks, "======\n")
  
  # Define indices for this chunk
  start_idx <- (i-1) * chunk_size + 1
  end_idx <- min(i * chunk_size, n_locations)
  
  # Extract chunk data
  X.0.chunk <- X.0[start_idx:end_idx, , drop = FALSE]
  coords.0.chunk <- coords.0[start_idx:end_idx, , drop = FALSE]
  
  cat("Current chunk size:", nrow(X.0.chunk), "locations\n")
  
  # Run prediction
  tryCatch({
    pred.chunk <- predict(out.sfMsPGOcc, 
                          X.0.chunk,  
                          coords.0.chunk,
                          n.omp.threads = 192,
                          verbose = TRUE,
                          n.report = 500,
                          type = 'occupancy')
    
    # Extract and save results for all species and locations
    chunk_summary <- data.frame()
    
    # Get number of locations from the prediction results
    n_loc <- dim(pred.chunk$psi.0.samples)[3]
    
    # Get dimensions of spatial random effect - check if it has species dimension
    w_dims <- dim(pred.chunk$w.0.samples)
    has_species_w <- length(w_dims) >= 3 && w_dims[2] > 1
    
    # Process each location and species
    for (loc in 1:n_loc) {
      # Extract spatial random effect for this location
      if (has_species_w) {
        # Multi-species spatial random effect (uncommon)
        w_means <- apply(pred.chunk$w.0.samples[, , loc], 2, mean)
        w_sds <- apply(pred.chunk$w.0.samples[, , loc], 2, sd)
      } else {
        # Single shared spatial random effect (common)
        w_mean <- mean(pred.chunk$w.0.samples[, 1, loc])
        w_sd <- sd(pred.chunk$w.0.samples[, 1, loc])
      }
      
      for (sp in 1:n_species) {
        # Calculate summary statistics for this species-location combination
        psi_mean <- mean(pred.chunk$psi.0.samples[, sp, loc])
        psi_sd <- sd(pred.chunk$psi.0.samples[, sp, loc])
        z_prob <- mean(pred.chunk$z.0.samples[, sp, loc])  # Posterior occurrence probability
        
        # Get spatial random effect for this species
        if (has_species_w) {
          w_mean_sp <- w_means[sp]
          w_sd_sp <- w_sds[sp]
        } else {
          # All species share the same spatial random effect
          w_mean_sp <- w_mean
          w_sd_sp <- w_sd
        }
        
        # Create a row for this species-location combination
        row_data <- data.frame(
          x = coords.0.chunk[loc, 1],
          y = coords.0.chunk[loc, 2],
          species = sp,
          psi_mean = psi_mean,
          psi_sd = psi_sd,
          z_prob = z_prob,
          w_mean = w_mean_sp,
          w_sd = w_sd_sp
        )
        
        # Add to chunk results
        chunk_summary <- rbind(chunk_summary, row_data)
      }
    }
    
    # Append to results
    results_summary <- rbind(results_summary, chunk_summary)
    
    # Save chunk results
    save(chunk_summary, 
         file = paste0("output/22222/sfMsPGOcc/chunks/chunk_", i, "_of_", n_chunks, ".RData"))
    
    # Clear memory
    rm(pred.chunk, chunk_summary)
    gc()
    
  }, error = function(e) {
    cat("ERROR in chunk", i, ":", conditionMessage(e), "\n")
    
    # If memory error, try with smaller sub-chunks
    if (grepl("cannot allocate", conditionMessage(e))) {
      smaller_chunk_size <- floor(chunk_size / 5)
      cat("Retrying with smaller chunks of size", smaller_chunk_size, "\n")
      
      n_sub_chunks <- ceiling((end_idx - start_idx + 1) / smaller_chunk_size)
      
      for (j in 1:n_sub_chunks) {
        sub_start <- start_idx + (j-1) * smaller_chunk_size
        sub_end <- min(start_idx + j * smaller_chunk_size - 1, end_idx)
        
        cat("Processing sub-chunk", j, "of", n_sub_chunks, "\n")
        
        X.0.sub <- X.0[sub_start:sub_end, , drop = FALSE]
        coords.0.sub <- coords.0[sub_start:sub_end, , drop = FALSE]
        
        tryCatch({
          pred.sub <- predict(out.sfMsPGOcc, 
                              X.0.sub,  
                              coords.0.sub,
                              n.omp.threads = 96,  # Reduced threads for smaller chunks
                              verbose = TRUE,
                              n.report = 100,
                              type = 'occupancy')
          
          # Create summary for sub-chunk
          sub_summary <- data.frame()
          
          # Get number of locations from the prediction results
          n_loc_sub <- dim(pred.sub$psi.0.samples)[3]
          
          # Get dimensions of spatial random effect - check if it has species dimension
          w_dims <- dim(pred.sub$w.0.samples)
          has_species_w <- length(w_dims) >= 3 && w_dims[2] > 1
          
          # Process each location and species
          for (loc in 1:n_loc_sub) {
            # Extract spatial random effect for this location
            if (has_species_w) {
              # Multi-species spatial random effect (uncommon)
              w_means <- apply(pred.sub$w.0.samples[, , loc], 2, mean)
              w_sds <- apply(pred.sub$w.0.samples[, , loc], 2, sd)
            } else {
              # Single shared spatial random effect (common)
              w_mean <- mean(pred.sub$w.0.samples[, 1, loc])
              w_sd <- sd(pred.sub$w.0.samples[, 1, loc])
            }
            
            for (sp in 1:n_species) {
              # Calculate summary statistics
              psi_mean <- mean(pred.sub$psi.0.samples[, sp, loc])
              psi_sd <- sd(pred.sub$psi.0.samples[, sp, loc])
              z_prob <- mean(pred.sub$z.0.samples[, sp, loc])
              
              # Get spatial random effect for this species
              if (has_species_w) {
                w_mean_sp <- w_means[sp]
                w_sd_sp <- w_sds[sp]
              } else {
                # All species share the same spatial random effect
                w_mean_sp <- w_mean
                w_sd_sp <- w_sd
              }
              
              # Create a row for this species-location combination
              row_data <- data.frame(
                x = coords.0.sub[loc, 1],
                y = coords.0.sub[loc, 2],
                species = sp,
                psi_mean = psi_mean,
                psi_sd = psi_sd,
                z_prob = z_prob,
                w_mean = w_mean_sp,
                w_sd = w_sd_sp
              )
              
              # Add to sub-chunk results
              sub_summary <- rbind(sub_summary, row_data)
            }
          }
          
          # Append to overall results
          results_summary <- rbind(results_summary, sub_summary)
          
          # Save sub-chunk results
          save(sub_summary, 
               file = paste0("output/22222/sfMsPGOcc/chunks/chunk_", 
                             i, "_subchunk_", j, "_of_", n_sub_chunks, ".RData"))
          
          rm(pred.sub, sub_summary)
          gc()
          
        }, error = function(e2) {
          cat("ERROR in sub-chunk", j, ":", conditionMessage(e2), "\n")
        })
      }
    }
  })
  
  # Save intermediate results periodically
  if (i %% 5 == 0 || i == n_chunks) {
    cat("Saving intermediate results\n")
    save(results_summary, file = "output/22222/sfMsPGOcc/intermediate_prediction_results.RData")
  }
}

# Save final combined results
save(results_summary, file = "output/22222/sfMsPGOcc/combined_prediction_results.RData")

# Create rasters for each species
if (nrow(results_summary) > 0) {
  cat("Creating prediction rasters for each species...\n")
  
  # Get unique species IDs
  species_ids <- unique(results_summary$species)
  
  for (sp in species_ids) {
    cat("Processing rasters for species", sp, "\n")
    
    # Filter data for this species
    sp_data <- results_summary[results_summary$species == sp, ]
    
    # Create rasters for this species
    psi_mean_raster <- rast(sp_data[, c("x", "y", "psi_mean")], type = "xyz", crs = "EPSG:3577")
    psi_sd_raster <- rast(sp_data[, c("x", "y", "psi_sd")], type = "xyz", crs = "EPSG:3577")
    z_prob_raster <- rast(sp_data[, c("x", "y", "z_prob")], type = "xyz", crs = "EPSG:3577")
    w_mean_raster <- rast(sp_data[, c("x", "y", "w_mean")], type = "xyz", crs = "EPSG:3577")
    w_sd_raster <- rast(sp_data[, c("x", "y", "w_sd")], type = "xyz", crs = "EPSG:3577")
    
    # Create output directory if needed
    species_dir <- paste0("output/22222/sfMsPGOcc/species_", sp)
    dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save rasters for this species
    writeRaster(psi_mean_raster, paste0(species_dir, "/psi_mean.tif"), overwrite = TRUE)
    writeRaster(psi_sd_raster, paste0(species_dir, "/psi_sd.tif"), overwrite = TRUE)
    writeRaster(z_prob_raster, paste0(species_dir, "/z_prob.tif"), overwrite = TRUE)
    writeRaster(w_mean_raster, paste0(species_dir, "/w_mean.tif"), overwrite = TRUE)
    writeRaster(w_sd_raster, paste0(species_dir, "/w_sd.tif"), overwrite = TRUE)
    
    cat("Saved rasters for species", sp, "\n")
  }
}

cat("\nMulti-species prediction process completed!\n")