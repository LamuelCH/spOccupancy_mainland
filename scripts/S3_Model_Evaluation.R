library(spOccupancy)


# Load model --------------------------------------------------------------
load("models/model_sfMsPGOcc_1981-2010_nthin40_nbatch6800_burn10000_reducdedSP.RData")



# inspect the Rhat and ESS values
summary(out.sfMsPGOcc)

# Beta occurrence covariates ------------------------------------------------------------------
# visualize community occurrence covariates convergence
plot(out.sfMsPGOcc, "beta.comm", density =FALSE)

# [Alternative] convergence of individual species
# plot(out.sfMsPGOcc, "beta", density =FALSE)


# Lamba latent loading ----------------------------------------------------
summary(out.sfMsPGOcc$lambda.samples)

# inspect latent factor loading
out.sfMsPGOcc$ESS$lambda #ESS
# Phalangeridae-1      Dasyuridae-1         Canidae-1 
# 0.000000          6.121878          5.062315 
# Felidae-1       Leporidae-1    Macropodidae-1 
# 8.729760          4.512669          7.849480 
# Muridae-1     Peramelidae-1      Petauridae-1 
# 4.002509          4.040657         28.171606 
# Potoroidae-1 Pseudocheiridae-1      Vombatidae-1 
# 7.396990         33.825440          2.579985 s
# Phalangeridae-2      Dasyuridae-2         Canidae-2 
# 0.000000          0.000000          5.378794
# Felidae-2       Leporidae-2    Macropodidae-2 
# 62.325358          5.515756          4.870279 
# Muridae-2     Peramelidae-2      Petauridae-2 
# 3.975408          6.067882         61.970080 
# Potoroidae-2 Pseudocheiridae-2      Vombatidae-2 
# 5.809854         23.895240          5.342259 
out.sfMsPGOcc$rhat$lambda.lower.tri #Rhat
# [1] 12.192902 13.818766  7.474647 10.476819  8.445917
# [6] 13.091519 17.318131  2.410477 11.987074  2.908469
# [11] 20.174667 11.760799  1.479879 11.517220 10.461248
# [16] 13.738965 16.902488  1.389088 11.269297  4.536947
# [21] 10.581946

# Generally when fitting spatial occupancy models, the detection parameters will 
# mix and converge faster than the occupancy parameters (assuming an adequate 
# amount of replication to separate occupancy from detection, which may not be 
# true with a very small number of replicates). In particular, the spatial decay 
# parameters and the occupancy intercepts may show slow mixing and/or convergence, 
# which I discuss more in depth below as to what to do in this case.

plot(out.sfMsPGOcc, 'lambda', density = FALSE) #Visual
# nothing converged atm





# Theta spatial decay parameter -------------------------------------------
# Spatial decay parameter
plot(out.sfMsPGOcc, 'theta', density = FALSE)



summary(out.sfMsPGOcc)

summary(out.sfMsPGOcc$lambda.samples)

# Takes a few seconds to run. 
ppc.occ.out <- ppcOcc(out.sfMsPGOcc, 'freeman-tukey', group = 2)

# Species-specific intercepts
plot(out.sfMsPGOcc$beta.samples[, 1:7], density = FALSE)

plot(out.sfMsPGOcc, "theta", density = FALSE)
# Calculate Bayesian p-values
summary(ppc.occ.out)





# Identify the problematic chain ------------------------------------------
# Create a filtered version of the spOccupancy object
filter_spOccupancy_chains <- function(model, chain_to_remove) {
  # Create a copy of the original model
  filtered_model <- model
  
  # Number of chains
  n_chains <- model$n.chains
  # Iterations per chain
  iterations_per_chain <- dim(model$beta.comm.samples)[1] / n_chains
  
  # Process each parameter that has MCMC samples
  param_list <- c("beta.comm.samples", "tau.sq.beta.samples", "beta.samples", 
                  "alpha.comm.samples", "tau.sq.alpha.samples", "alpha.samples",
                  "theta.samples", "lambda.samples", "phi.samples", "w.samples", 
                  "z.samples", "psi.samples", "sigma.sq.psi.samples")
  
  # Create indices for each chain
  chain_indices <- list()
  for (i in 1:n_chains) {
    start_idx <- (i-1) * iterations_per_chain + 1
    end_idx <- i * iterations_per_chain
    chain_indices[[i]] <- start_idx:end_idx
  }
  
  # Rows to keep (all chains except the one to remove)
  rows_to_keep <- unlist(chain_indices[-chain_to_remove])
  
  # Process each parameter
  for (param in param_list) {
    if (param %in% names(filtered_model)) {
      param_obj <- filtered_model[[param]]
      
      if (inherits(param_obj, "mcmc")) {
        # For mcmc objects (2D: iterations x parameters)
        filtered_obj <- param_obj[rows_to_keep, , drop = FALSE]
        
        # Restore mcmc attributes
        attr(filtered_obj, "mcpar") <- c(1, length(rows_to_keep), 1)
        class(filtered_obj) <- "mcmc"
        
        # Update the model
        filtered_model[[param]] <- filtered_obj
      } 
      else if (is.array(param_obj) && length(dim(param_obj)) == 3) {
        # For 3D arrays (iterations x species x parameters)
        # This handles cases like beta.samples
        n_species <- dim(param_obj)[2]
        n_params <- dim(param_obj)[3]
        
        # Reshape and filter
        filtered_obj <- array(NA, dim = c(length(rows_to_keep), n_species, n_params))
        for (i in 1:n_species) {
          for (j in 1:n_params) {
            filtered_obj[, i, j] <- param_obj[rows_to_keep, i, j]
          }
        }
        
        # Update the model
        filtered_model[[param]] <- filtered_obj
      }
    }
  }
  
  # Update chain count in the model
  filtered_model$n.chains <- n_chains - 1
  filtered_model$n.samples <- length(rows_to_keep)
  
  return(filtered_model)
}

# Apply the function to your model
chain_to_remove <- 2  # The problematic chain (red line in trace plots)
filtered_sfMsPGOcc <- filter_spOccupancy_chains(out.sfMsPGOcc, chain_to_remove)

# Check if filtering worked correctly
original_samples <- dim(out.sfMsPGOcc$beta.comm.samples)[1]
filtered_samples <- dim(filtered_sfMsPGOcc$beta.comm.samples)[1]
cat("Original samples:", original_samples, "\n")
cat("Filtered samples:", filtered_samples, "\n")
cat("Expected filtered samples:", original_samples * (out.sfMsPGOcc$n.chains - 1) / out.sfMsPGOcc$n.chains, "\n")

# Now you can use your existing code with the filtered model
# For example:
summary(filtered_sfMsPGOcc)
plot(filtered_sfMsPGOcc, "beta.comm", density = FALSE)
plot(filtered_sfMsPGOcc, "phi", density = FALSE)
ppc.occ.out <- ppcOcc(filtered_sfMsPGOcc, 'freeman-tukey', group = 2)
summary(ppc.occ.out)



summary(filtered_sfMsPGOcc$lambda.samples)


# Plot beta parameters for the first few species
plot(filtered_sfMsPGOcc, "beta", density = FALSE, species = 1:3)

# Plot lambda parameters 
plot(filtered_sfMsPGOcc, "lambda", density = FALSE)

# Plot theta parameters
plot(filtered_sfMsPGOcc, "theta", density = FALSE)


# Quick check if filtering actually happened
dim_original <- dim(out.sfMsPGOcc$beta.comm.samples)[1]
dim_filtered <- dim(filtered_sfMsPGOcc$beta.comm.samples)[1]

cat("Original dimensions:", dim_original, "\n")
cat("Filtered dimensions:", dim_filtered, "\n")
