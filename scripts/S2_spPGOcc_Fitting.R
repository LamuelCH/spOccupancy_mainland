# For Slurm Script Array  -------------------------------------------------
tryCatch({
  # Set up error logging
  log_file <- file.path("logs", paste0("job_", Sys.getenv("SLURM_ARRAY_TASK_ID"), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
  
  # Read command-line argument for array index
  args <- commandArgs(trailingOnly = TRUE)
  if(length(args) == 0) {
    job_index <- 0  # Default for local testing
    warning("No job index provided, using default value 0")
  } else {
    job_index <- as.integer(args[1])  # SLURM_ARRAY_TASK_ID passed here
    if(is.na(job_index)) {
      stop("Invalid job index: must be an integer")
    }
  }
  
  # Define MCMC parameter ---------------------------------------------------
  # Define thinning factors and corresponding parameters
  thinning_factors <- c(1, 10, 100, 1000)

    if(job_index + 1 > length(thinning_factors)) {
    stop(sprintf("Job index %d out of range (max: %d)", job_index, length(thinning_factors) - 1))
  }
  
  n.thin <- thinning_factors[job_index + 1]
  n.sample <- 1000
  n.burn <- 3000 * n.thin
  batch.length <- 25
  n.batch <- ceiling((n.burn/n.thin + n.sample) * n.thin / batch.length)
  n.chains <- 4
  n.neighbors <- 10
  
  cat(sprintf("Running with parameters: n.thin=%d, n.sample=%d, n.burn=%d, n.batch=%d\n", 
              n.thin, n.sample, n.burn, n.batch))
  
  # Original script ---------------------------------------------------------
  set.seed(123)
  
  if(!require(spOccupancy)) {
    stop("Package 'spOccupancy' not available. Please install it first.")
  }
  
  # Check working directory and set it correctly
  if(!file.exists("input/data_singleSeason.RData")) {
    # Try to find the right directory
    possible_dirs <- c(".", "..", "../..")
    found <- FALSE
    for(dir in possible_dirs) {
      if(file.exists(file.path(dir, "input/data_singleSeason.RData"))) {
        setwd(dir)
        found <- TRUE
        cat("Setting working directory to:", getwd(), "\n")
        break
      }
    }
    if(!found) {
      stop("Input data file not found. Check working directory.")
    }
  }
  
  # Load the data -----------------------------------------------------------
  data_file <- "input/data_singleSeason.RData"
  if(!file.exists(data_file)) {
    stop(paste("Data file not found:", data_file))
  }
  
  # Load data with error handling
  tryCatch({
    load(data_file)
    if(!exists("data")) {
      stop("Data object not found in loaded file")
    }
  }, error = function(e) {
    stop(paste("Error loading data file:", e$message))
  })
  
  cat("Data loaded successfully\n")
  cat("Structure of data object:\n")
  str(data)
  
  # Species codes
  sp.names <- rownames(data$y)
  cat("Species names:", paste(sp.names, collapse=", "), "\n")
  
  # Parametize --------------------------------------------------------------
  # check the raw probabilities of occurrence from detection-non-detection data
  occ_probs <- apply(data$y, 1, mean, na.rm = TRUE)
  cat("Raw occurrence probabilities:\n")
  print(occ_probs)
  
  # Our modelled community is relatively small
  # and the species are all quite similar (mammals)
  # n.factors <- 2
  
  # Ordering species in detection-non detection array ------------------------
  sp.ordered = c("Dasyurus")
  
  # Check that all ordered species exist in the data
  if(!all(sp.ordered %in% sp.names)) {
    missing_sp <- sp.ordered[!sp.ordered %in% sp.names]
    stop(paste("Some species in sp.ordered are not in the data:", paste(missing_sp, collapse=", ")))
  }
  
  # Create new detection-nondetection data matrix in the new order
  y.new <- data$y[sp.ordered, ,]
  
  # Create a new data array and update with reordered data
  data.ordered <- data
  data.ordered$y <- y.new
  
  cat("New data structure with reordered species:\n")
  str(data.ordered)
  
  # Check for required components in the data
  required_components <- c("y", "occ.covs", "det.covs", "coords")
  missing_components <- required_components[!required_components %in% names(data.ordered)]
  if(length(missing_components) > 0) {
    stop(paste("Missing required data components:", paste(missing_components, collapse=", ")))
  }
  
  # Specify initial value and prior distributions ---------------------------
  # If not specify, the model will use default value, but model is sensitive to initial,
  
  # Pair-wise distance between all sites
  dist.data <- try(dist(data.ordered$coords))
  if(inherits(dist.data, "try-error")) {
    stop("Error calculating distances between sites. Check coords data.")
  }
  
  # Exponential correlation model
  cov.model <- "exponential"
  # Number of species
  N <- nrow(data.ordered$y)
  
  # Create list of initial values
  inits <- list(
    beta = 0,
    alpha = 0,
    phi = 3/30000,
    z = apply(data.ordered$y, 1, max, na.rm = TRUE)
  )
  
  # Setting Priors ----------------------------------------------------------
  min.dist <- min(dist.data, na.rm = TRUE)
  max.dist <- 100000 #100km
  
  if(min.dist <= 0 || is.infinite(min.dist) || is.na(min.dist)) {
    warning("Invalid minimum distance detected. Using default value of 1000.")
    min.dist <- 1000
  }
  
  priors <- list(
    beta.normal = list(mean = 0, var = 2.72),
    alpha.normal = list(mean = 0, var = 2.72),
    phi.unif = c(3 / max.dist, 3 / min.dist)
  )
  
  tuning <- list(phi = 0.5)
  
  # System resources for parallelization
  n.cores.available <- 32
  n.omp.threads <- n.cores.available  # Leave one core free
  cat(sprintf("Using %d of %d available cores for parallelization\n", n.omp.threads, n.cores.available))
  
  verbose <- TRUE
  n.report <- 100
  
  # Model fitting -----------------------------------------------------------
  # Specify formula
  occ.formula = ~ scale(bio5) + I(scale(bio5)^2) + 
    scale(bio6) + I(scale(bio6)^2) + 
    scale(bio12) + I(scale(bio12)^2) + 
    scale(bio15) + I(scale(bio15)^2) +
    scale(dem) + I(scale(dem)^2) +
    scale(eastness) + 
    scale(northness) +
    scale(tpi) + I(scale(tpi)^2) +
    scale(tri) + I(scale(tri)^2) +
    scale(rock)
  
  det.formula <- ~ scale(effort) + (1 | project)
  
  # Check if output directory exists, create if not
  output_dir <- "models/spPGOcc"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Run Model with error handling
  start_time <- Sys.time()
  cat(paste("Start Time:", start_time), "\n")
  
  # Create a temporary file to save intermediate results
  temp_file <- tempfile(pattern = "model_intermediate_", fileext = ".RData")
  
  model_result <- tryCatch({
    spPGOcc(
      occ.formula = occ.formula,
      det.formula = det.formula,
      data = data.ordered,
      inits = inits,
      n.batch = n.batch,
      batch.length = batch.length,
      accept.rate = 0.43,
      priors = priors,
      cov.model = cov.model,
      tuning = tuning,
      n.omp.threads = n.omp.threads,
      verbose = verbose,
      NNGP = TRUE,
      n.neighbors = n.neighbors,
      n.report = n.report,
      n.burn = n.burn,
      n.thin = n.thin,
      n.chains = n.chains
    )
  }, error = function(e) {
    cat("Error in model fitting:", e$message, "\n")
    return(NULL)
  }, warning = function(w) {
    cat("Warning in model fitting:", w$message, "\n")
    return(NULL)
  }, finally = {
    # Calculate runtime regardless of success or failure
    end_time <- Sys.time()
    cat(paste("End Time:", end_time), "\n")
    cat(paste("Duration:", end_time - start_time), "\n")
  })
  
  # Save model if successful
  if(!is.null(model_result)) {
    out.spPGOcc <- model_result
    output_file <- paste0(output_dir, "/model_spPGOcc_2012-2021_nthin", 
                          n.thin, "_nburn", n.burn, "_nchain", n.chains, 
                          "_nsample", n.sample,"_neighbors", n.neighbors, ".RData")
                          # "_nfactors", n.factors, 
                          
    
    tryCatch({
      save(out.spPGOcc, file = output_file)
      cat("Model saved successfully to:", output_file, "\n")
    }, error = function(e) {
      cat("Error saving model:", e$message, "\n")
      # Try saving to a different location
      backup_file <- file.path(tempdir(), basename(output_file))
      cat("Attempting to save to backup location:", backup_file, "\n")
      save(out.spPGOcc, file = backup_file)
    })
  } else {
    cat("Model fitting failed. Check error messages above.\n")
  }
  
}, error = function(e) {
  # Global error handler
  cat(paste("CRITICAL ERROR:", e$message), "\n")
  cat("Stack trace:\n")
  print(sys.calls())
  # Write error to a file for debugging
  error_file <- file.path("logs", paste0("error_", Sys.getenv("SLURM_ARRAY_TASK_ID"), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  dir.create(dirname(error_file), showWarnings = FALSE, recursive = TRUE)
  write(paste("Error:", e$message), file = error_file)
  write(paste("Stack trace:", paste(capture.output(print(sys.calls())), collapse = "\n")), file = error_file, append = TRUE)
  # Return non-zero exit code
  quit(status = 1)
})