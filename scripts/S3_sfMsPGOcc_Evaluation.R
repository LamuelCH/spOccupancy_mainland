# The following script evaluate the convergence of spatial factor multi-species model
library(spOccupancy)


# Defining lamda initial Value --------------------------------------------------
# spOccupancy that use a factor modeling approach can be fairly sensitive to the 
# initial values of the latent factor loadings. This is primarily an issue when 
# there are a large number of rare species. If you encounter difficulties in model
# convergence when running factor models in spOccupancy across multiple chains,
# we recommend first running a single chain of the model for a moderate number of
# iterations until the traceplots look like they are settling around a value 
# (i.e., convergence is closed to being reached). Then extract the estimated mean 
# values for the factor loadings matrix (Λ
# ) and supply these as initial values to the spOccupancy function when running
# the full model across multiple chains. 

load("models/model_sfMsPGOcc_1981-2010_nthin1_nbatch120_nchain1_nburn2000.RData")

# Visual inspect if the value is settling
plot(out.sfMsPGOcc, 'lambda', density = FALSE) 

# If visually settled, extract the mean value for official run
summary(out.sfMsPGOcc$lambda.samples) 

# Iterations = 1:1000
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 1000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean     SD Naive SE Time-series SE
# Wallabia-1       1.00000 0.0000 0.000000        0.00000
# Rattus-1        -1.64119 0.2302 0.007281        0.05122
# Vulpes-1         1.14054 0.1196 0.003783        0.01250
# Canis-1          0.04383 0.1863 0.005891        0.03670
# Trichosurus-1    0.95648 0.1134 0.003586        0.01532
# Vombatus-1       0.03320 0.1083 0.003423        0.01359
# Macropus-1       3.09904 0.3121 0.009869        0.07460
# Notamacropus-1   2.13176 0.2003 0.006334        0.03247
# Felis-1         -0.34928 0.1770 0.005599        0.03484
# Dasyurus-1      -0.40879 0.1509 0.004772        0.02328
# Antechinus-1    -1.11806 0.7962 0.025178        0.40211
# Perameles-1     -1.38463 0.1737 0.005493        0.02862
# Isoodon-1       -1.00118 0.1613 0.005100        0.02475
# Pseudocheirus-1  0.15091 0.1921 0.006075        0.03531
# Oryctolagus-1    0.85159 0.1400 0.004427        0.02077
# Wallabia-2       0.00000 0.0000 0.000000        0.00000
# Rattus-2         1.00000 0.0000 0.000000        0.00000
# Vulpes-2        -2.19305 0.2555 0.008079        0.02860
# Canis-2          0.95029 0.4166 0.013174        0.08019
# Trichosurus-2    0.92489 0.2216 0.007006        0.02923
# Vombatus-2      -4.35135 0.3440 0.010879        0.04622
# Macropus-2      -1.63711 0.3863 0.012217        0.06693
# Notamacropus-2   0.36829 0.2979 0.009420        0.04617
# Felis-2         -0.66094 0.3618 0.011441        0.05740
# Dasyurus-2       1.24508 0.3694 0.011682        0.05500
# Antechinus-2    -0.05014 0.7124 0.022529        0.14944
# Perameles-2      1.27021 0.2612 0.008260        0.03363
# Isoodon-2        2.66528 0.3645 0.011528        0.10220
# Pseudocheirus-2 -1.47847 0.4656 0.014723        0.08836
# Oryctolagus-2   -0.50368 0.2374 0.007508        0.03280

load("models/model_sfMsPGOcc_1981-2010_nthin1_nbatch1360_nchain1_nburn2000.RData")
plot(out.sfMsPGOcc, 'lambda', density = FALSE) 
summary(out.sfMsPGOcc$lambda.samples) 

# Iterations = 1:32000
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 32000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean     SD  Naive SE Time-series SE
# Wallabia-1       1.00000 0.0000 0.0000000       0.000000
# Rattus-1        -1.61435 0.2053 0.0011474       0.008945
# Vulpes-1         1.01789 0.1462 0.0008170       0.006088
# Canis-1          0.04355 0.1895 0.0010592       0.007187
# Trichosurus-1    1.13166 0.1439 0.0008042       0.006332
# Vombatus-1       0.06906 0.1113 0.0006223       0.003754
# Macropus-1       2.88112 0.3116 0.0017418       0.014968
# Notamacropus-1   1.94999 0.2082 0.0011640       0.008607
# Felis-1         -0.35744 0.1921 0.0010737       0.007881
# Dasyurus-1      -0.40666 0.1662 0.0009289       0.005025
# Antechinus-1    -0.70650 0.7112 0.0039760       0.071260
# Perameles-1     -1.35653 0.1717 0.0009599       0.005886
# Isoodon-1       -0.93346 0.1672 0.0009344       0.007343
# Pseudocheirus-1  0.07471 0.2068 0.0011558       0.008130
# Oryctolagus-1    0.80270 0.1429 0.0007990       0.004363
# Wallabia-2       0.00000 0.0000 0.0000000       0.000000
# Rattus-2         1.00000 0.0000 0.0000000       0.000000
# Vulpes-2        -1.59665 0.2439 0.0013634       0.015022
# Canis-2          0.83760 0.3376 0.0018872       0.014783
# Trichosurus-2    0.98227 0.2242 0.0012535       0.015083
# Vombatus-2      -3.95214 0.4491 0.0025106       0.034217
# Macropus-2      -1.54852 0.3471 0.0019405       0.018868
# Notamacropus-2  -0.20547 0.3084 0.0017242       0.025369
# Felis-2         -0.14734 0.3299 0.0018440       0.016495
# Dasyurus-2       1.02661 0.3334 0.0018639       0.020937
# Antechinus-2    -0.66371 0.9514 0.0053185       0.075774
# Perameles-2      1.09775 0.2639 0.0014752       0.013271
# Isoodon-2        1.84275 0.3004 0.0016795       0.026150
# Pseudocheirus-2 -1.28563 0.3862 0.0021588       0.022445
# Oryctolagus-2   -0.46241 0.2329 0.0013019       0.009255



# Load model --------------------------------------------------------------
load("models/model_sfMsPGOcc_1981-2010_nthin1_nbatch160_ncha")


# Overall model performance -----------------------------------------------
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
out.sfMsPGOcc$rhat$lambda.lower.tri #Rhat

plot(out.sfMsPGOcc, 'lambda', density = FALSE) #Visual

# Generally when fitting spatial occupancy models, the detection parameters will 
# mix and converge faster than the occupancy parameters (assuming an adequate 
# amount of replication to separate occupancy from detection, which may not be 
# true with a very small number of replicates). In particular, the spatial decay 
# parameters and the occupancy intercepts may show slow mixing and/or convergence, 
# which I discuss more in depth below as to what to do in this case.



# Theta spatial decay parameter -------------------------------------------
# Spatial decay parameter
plot(out.sfMsPGOcc, 'theta', density = FALSE)




# Posterior Predictive Checks ---------------------------------------------
# The spOccupancy function ppcOcc() performs a posterior predictive check for all 
# spOccupancy model objects as an assessment of Goodness of Fit (GoF). The key idea
# of GoF testing is that a good model should generate data that closely align with 
# the observed data. If there are large differences in the observed data from the
# model-generated data, our model is likely not very useful (Hooten and Hobbs 2015).
# We can use the ppcOcc() and summary() functions to generate a Bayesian p-value 
# as a quick assessment of model fit. A Bayesian p-value that hovers around 0.5 
# indicates adequate model fit, while values less than 0.1 or greater than 0.9 
# suggest our model does not fit the data well. ppcOcc will return an overall
# Bayesian p-value for the entire community, as well as an individual Bayesian 
# p-value for each species. See the introductory spOccupancy vignette and the 
# ppcOcc() help page for additional details. Below we perform a posterior 
# predictive check with the Freeman-Tukey statistic, grouping the data by 
# individual sites (group = 1).

# Takes a few seconds to run. 
ppc.occ.out <- ppcOcc(out.sfMsPGOcc, 'freeman-tukey', group = 1)

# Calculate Bayesian p-values
summary(ppc.occ.out)


























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
