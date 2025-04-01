library(terra)
library(dplyr)
library(spOccupancy)

# Create output directories
dir.create("output/spPGOcc/chunks", recursive = TRUE, showWarnings = FALSE)

# Define chunk size for processing
chunk_size <- 1000

# Load the model and data
cat("Loading model and data...\n")
load("models/22222/spPGOcc/model_spPGOcc_1981-2010_nthin500_nbatch1e+05_nchain4_nburn1500000.RData")
load("input/list_sfMsPGOcc_1980-2010.RData")

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

# Initialize results data frame
results_summary <- data.frame(
  x = numeric(0),
  y = numeric(0),
  psi.mean = numeric(0),
  psi.sd = numeric(0)
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
    pred.chunk <- predict(out.spPGOcc, 
                          X.0.chunk,  
                          coords.0.chunk,
                          n.omp.threads = 192,
                          verbose = TRUE,
                          n.report = 500,
                          type = 'occupancy')
    
    # Extract and save results
    chunk_summary <- data.frame(
      x = coords.0.chunk[,1],
      y = coords.0.chunk[,2],
      psi.mean = apply(pred.chunk$psi.0.samples, 2, mean),
      psi.sd = apply(pred.chunk$psi.0.samples, 2, sd)
    )
    
    # Append to results
    results_summary <- rbind(results_summary, chunk_summary)
    
    # Save chunk results
    save(chunk_summary, 
         file = paste0("output/spPGOcc/chunks/chunk_", i, "_of_", n_chunks, ".RData"))
    
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
          pred.sub <- predict(out.spPGOcc, 
                              X.0.sub,  
                              coords.0.sub,
                              n.omp.threads = 192,
                              verbose = TRUE,
                              n.report = 500,
                              type = 'occupancy')
          
          sub_summary <- data.frame(
            x = coords.0.sub[,1],
            y = coords.0.sub[,2],
            psi.mean = apply(pred.chunk$psi.0.samples, 2, mean),
            psi.sd = apply(pred.chunk$psi.0.samples, 2, sd)
          )
          
          results_summary <- rbind(results_summary, sub_summary)
          
          save(sub_summary, 
               file = paste0("output/spPGOcc/chunks/chunk_", 
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
  if (i %% 10 == 0 || i == n_chunks) {
    cat("Saving intermediate results\n")
    save(results_summary, file = "output/spPGOcc/intermediate_prediction_results.RData")
  }
}

# Save final combined results
save(results_summary, file = "output/spPGOcc/combined_prediction_results.RData")

# Create rasters if we have results
if (nrow(results_summary) > 0) {
  cat("Creating prediction rasters...\n")
  psi_mean_raster <- rast(results_summary[, c("x", "y", "psi.mean")], type = "xyz", crs = "EPSG:3577")
  psi_sd_raster <- rast(results_summary[, c("x", "y", "psi.sd")], type = "xyz", crs = "EPSG:3577")
  
  writeRaster(psi_mean_raster, "output/spPGOcc/psi_mean_prediction.tif", overwrite = TRUE)
  writeRaster(psi_sd_raster, "output/spPGOcc/psi_sd_prediction.tif", overwrite = TRUE)
  
  cat("Prediction rasters saved\n")
}

cat("\nPrediction process completed!\n")