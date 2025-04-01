# For Slurm Script Array  -------------------------------------------------
# Read command-line argument for array index
args <- commandArgs(trailingOnly = TRUE)
# job_index <- as.integer(args[1])  # SLURM_ARRAY_TASK_ID passed here
job_index = 0

data_factor = c("input/list_22222_1980-2010.RData", 
                "input/list_343_1980-2010.RData", 
                "input/list_55_1980-2010.RData")

load (data_factor[job_index+1])

# Define MCMC parameter ---------------------------------------------------
# Define thinning factors and corresponding parameters
# thinning_factors <- c(1, 10, 50, 100)
# sample_factors = c(1000, 2000, 4000, 8000, 16000, 32000)

# n.thin <- thinning_factors[job_index + 1]  # +1 because array starts at 0
n.thin = 1 # use this arg to test init value

#n.sample = sample_factors[job_index + 1] # use this arg to test the init value with single chain
n.sample = 10000

n.burn = 3000 * n.thin # disregard first 2000 posterior samples (not iterations)  

batch.length = 25 # Keep it default
n.batch <- ceiling((n.burn/n.thin + n.sample) * n.thin / 25) # total sample per chain = n_batch * batch.length

n.chains = 1 


# Original script ---------------------------------------------------------
set.seed(123) #so to get the same results
library(spOccupancy)
setwd(".")



# Load the data -----------------------------------------------------------
str(data)
# Species codes
sp.names = rownames(data$y)


# Parametize --------------------------------------------------------------
# check the raw probabilities of occurence from detection-non-dectection data
apply(data$y, 1, mean, na.rm = TRUE)

# Our modelled community is relatively small (N = 33)
# and the species are all quite similar (mammals)
# We have some really rare species and some pretty common species 
# So we try the model with a small number of latent factors (q in statistical notation)
n.factors <- 2



# Ordering species in detection-non detection array ------------------------
# careful consideration of the ordering of species can lead to:
# (1) increased interpretability of the factors and factor loadings; 
# (2) faster model convergence; and 
# (3) improved mixing.

# Current species (here family) ordering
sp.names

# Reorder species:
# 1. Place a common species first
# 2. for the remaning q-1 factors, place species that believed will show different occurrence patterns than the first species and the other species placed before it 
sp.ordered = c("Macropus",
               "Dasyurus",
               "Wallabia",
               "Rattus",
               "Vulpes",
               "Canis",
               "Trichosurus",  
               "Vombatus",
               "Notamacropus" ,
               "Felis",
               "Antechinus",
               "Perameles", 
               "Isoodon",
               "Pseudocheirus",
               "Oryctolagus")

# Create new detection-nondetection data matrix in the new order
y.new = data$y[sp.ordered, ,]
# Create a new data array
data.ordered = data 
# Change the data to the new ordered data
data.ordered$y = y.new
str(data.ordered)

apply(data.ordered$y, 1, mean, na.rm = TRUE)
# Wallabia        Rattus        Vulpes         Canis   Trichosurus      Vombatus      Macropus 
# 0.58737646    0.34321653    0.25864780    0.03942049    0.48528751    0.32457323    0.21922731 
# Notamacropus         Felis      Dasyurus    Antechinus     Perameles       Isoodon Pseudocheirus 
# 0.21945193    0.16183738    0.02392183    0.14746181    0.13376011    0.10220126    0.07176550 
# Oryctolagus 
# 0.07401168 


# Specify initial value and prior distributions ---------------------------
# If not specify, the model will use default value, but model is sensitive to initial, 
# if convergence is not achieved, run a single chain with moderate number of iterations until 
# the traceplots looks settling, then extract the estimated mean values for the factor loading matrix 
# and supply these as initial values.

# Pair-wise distance between all sites
dist.data <- dist(data.ordered$coords)
# Exponential correlation model
cov.model <- "exponential"
# Specify all other initial values identical to lfMsPGOcc() from before
# Number   of species
N <- nrow(data.ordered$y)

# Initiate all lambda initial values to 0.
lambda.inits <- matrix(0, N, n.factors)
# Set diagonal elements to 1
diag(lambda.inits) <- 1
# Set lower triangular elements to random values from a standard normal dist
lambda.inits[lower.tri(lambda.inits)] <- rnorm(sum(lower.tri(lambda.inits)))
# Check it out
lambda.inits

# Initiate lambda value based on preliminary model iterations = 40000
# Extract the mean values from the preliminary model
# lambda_means <- c(
# 1.0000,   # Wallabia-1
# -1.3095,  # Rattus-1
# 1.6395,   # Vulpes-1
# -0.3115,  # Canis-1
#   0.6372,   # Trichosurus-1
#   0.1173,   # Vombatus-1
#   2.8253,   # Macropus-1
#   1.4969,   # Notamacropus-1
#   -0.3714,  # Felis-1
#   -0.7528,  # Dasyurus-1
#   -0.5066,  # Antechinus-1
#   -1.6724,  # Perameles-1
#   -1.5947,  # Isoodon-1
#   0.2345,   # Pseudocheirus-1
#   0.7909    # Oryctolagus-1
# )

# Create the lambda.inits matrix
# lambda.inits <- matrix(lambda_means, nrow = N, ncol = n.factors)


# Create list of initial values.
inits <- list(alpha.comm = 0, #Community level detection coefficients
              beta.comm = 0, #community level occurrence coefficients
              beta = 0, #species-level occurrence coefficients
              alpha = 0, #specie-level detection coefficients 
              tau.sq.beta = 1, #community-level occurrence variance parameters
              tau.sq.alpha = 1, #community-level detection variance parameters
              lambda = lambda.inits, #species-specific factor loading
              phi = 3/30000, #Start with the approximate home range of STQ 30km
              z = apply(data.ordered$y, c(1, 2), max, na.rm = TRUE)) #latent occurrence variables for all species


# Adaptive MCMC Sampler ---------------------------------------------------
accept.rate = 0.43 # leave as default


# Setting Priors ----------------------------------------------------------
# Assume unifor priors for the spatial decay parameter phi 
# We recommend determining the
# bounds of the uniform distribution by computing the smallest distance between sites and the largest distance
# between sites in the observed data set
min.dist <- min(dist.data)
max.dist <- 100000 #100km

priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
               alpha.comm.normal = list(mean = 0, var = 2.72),
               tau.sq.beta.ig = list(a = 0.1, b = 0.1),
               tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
               phi.unif = list(3 / max.dist, 3 / min.dist))

# Tunning parameters
# we also need to specify an initial value for the tuning parameters for the spatial decay and
# smoothness parameters (if applicable). These values are supplied as input in the form of a list with tags
# phi and nu. The initial tuning value can be any value greater than 0, but we recommend starting the value
# out around 0.5. After some initial runs of the model, if you notice the final acceptance rate of a parameter
# is much larger or smaller than the target acceptance rate (accept.rate), you can then change the initial
# tuning value to get closer to the target rate. Here we set the initial tuning value for phi to 1 after some
# initial exploratory runs of the model.
tuning <- list(phi = 0.5)
# if by the end the value is well below the target acceptance rate, recommend rerunning the model with a larger initial tuning parameter (higher than the final reported value in the displayed output of model progress). Opposite vice versa.

# parallelize
n.omp.threads <- 32 # within-chain parallelization
verbose <- TRUE
n.report <- 10 # Report progress at every 200th batch.




###########################################################################
# Model fitting -----------------------------------------------------------
# Specify formula
occ.formula = ~ scale(bio5) + I(scale(bio5)^2) +scale(bio6) + I(scale(bio6)^2) + 
  scale(bio12) + I(scale(bio12)^2) + scale(bio15) + I(scale(bio15)^2) +
  scale(dem) + I(scale(dem)^2) +
  scale(slope) + I(scale(slope)^2) +
  scale(aspect) +
  scale(tpi) + I(scale(tpi)^2) +
  scale(tri) + I(scale(tri)^2) +
  scale(rock)

det.formula = ~ scale(effort) + (1 | project) #assuume no variation in detection probability


# Run Model 
start_time <- Sys.time()
print(paste("Start Time:", start_time))

out.sfMsPGOcc <- sfMsPGOcc(occ.formula = occ.formula,
                           det.formula = det.formula,
                           data = data.ordered,
                           inits = inits,
                           n.batch = n.batch,
                           batch.length = batch.length,
                           accept.rate = 0.43,
                           priors = priors,
                           n.factors = n.factors,
                           cov.model = cov.model,
                           tuning = tuning,
                           n.omp.threads = n.omp.threads,
                           verbose = TRUE,
                           NNGP = TRUE,
                           n.neighbors = 10,
                           n.report = n.report,
                           n.burn = n.burn,
                           n.thin = n.thin,
                           n.chains = n.chains)
#k.fold = 2,
#k.fold.threads = n.omp.threads,
#k.fold.seed = 123

# Extract dataset identifier from data_factor path
dataset_id <- sub(".*list_([^_]+)_.*\\.RData", "\\1", data_factor[job_index+1])

# Create directory if it doesn't exist
dir.create(paste0("models/", dataset_id, "/sfMsPGOcc/init"), recursive = TRUE, showWarnings = FALSE)

# Save with dataset-specific name
save(out.sfMsPGOcc, 
     file = paste0("models/", dataset_id, "/sfMsPGOcc/model_", dataset_id, "_sfMsPGOcc_1981-2010_nthin", 
                   n.thin, "_nbatch", n.batch, "_nchain", n.chains, "_nburn", n.burn, ".RData"))

end_time <- Sys.time()
print(paste("End Time:", end_time))
print(paste("Duration:", end_time - start_time))

