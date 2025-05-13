# Load required libraries
library(tidyverse)
library(ggmap)
library(spOccupancy)
library(terra)
library(sf)
library(lubridate)
library(dplyr)

# Setting Project Environment ---------------------------------------------
# Set the path to the root 'data' directory
# setwd(".")

cat("=== Starting data preparation for occupancy modeling ===\n")

# VICTORIA VBA DATA WRANGLING ------------------------------------------------------
cat("\nProcessing VBA data...\n")

# Specify folder path 
path_VBA = "data/occ/VBA"

# Create a list of all VBA data inside the folder
list_VBA = list.files(path = path_VBA, 
                      pattern = "^VBA", 
                      recursive = TRUE, 
                      full.names = TRUE) 

# Create VBA df from the list 
df_VBA = lapply(list_VBA, function(file) {
  df <- read.csv(file, 
                 header = TRUE,
                 skip = 16, # skip the first 16 rows are metadata                   
                 stringsAsFactors = FALSE)
  df$Scientific.Name <- as.character(df$Scientific.Name) # Convert SpeciesCode to character
  df$Project.ID = as.character(df$Project.ID)
  return(df)
}) %>% 
  bind_rows(.) # combine all rows together

# column operation
# Create new columns essential columns
# Add year month day columns
df_VBA = df_VBA %>%
  mutate(
    # Try parsing DateLast with formats focusing on year and month
    Date = parse_date_time(Survey.Start.Date, c(
      "dmy_HMS", "dmy", "ymd_HMS", "ymd", "mdy_HMS", "mdy", # With day and optional time
      "ym", "my", # Year and month only
      "y", "m" # Year or month only (less specific, use with caution)
    ), quiet = TRUE),
    # Extract components from the parsed Date column
    year = year(Date),
    month = month(Date),
    day = day(Date),
    # Add a source column and label all values as "VIC"
    Source = "VIC"
  )

# Check for parsing failures
parse_failures_vba <- df_VBA %>% filter(is.na(Date) & !is.na(Survey.Start.Date))
if(nrow(parse_failures_vba) > 0) {
  warning(paste("Failed to parse", nrow(parse_failures_vba), "dates in VBA data"))
  print(head(parse_failures_vba$Survey.Start.Date))
}

# Add genus column
df_VBA$genus <- sub(" .*", "", df_VBA$Scientific.Name)

# Re-arranging df and retain essential information
df_VBA = df_VBA %>% 
  dplyr::select(Source = Source,
                ProjectID = Project.ID,
                SiteName = Site.Name,
                Year = year,
                Month = month,
                Day = day,
                Genus = genus,
                ScientificName = Scientific.Name,
                Longitude_GDA94 = Longitude.GDA94,
                Latitude_GDA94 = Latitude.GDA94)

# Final filtering
df_VBA = df_VBA %>% 
  distinct() %>%
  filter(!is.na(Longitude_GDA94) & !is.na(Latitude_GDA94)) %>% 
  filter(Year %in% 2012:2021) %>%
  na.omit() # Remove rows with missing values

# Summarize observations by genus
summary_genus_VBA <- df_VBA %>%
  group_by(Genus) %>%
  summarise(Observations = n()) %>%
  arrange(desc(Observations))

# View the summary
cat("\nTop genera in VBA data:\n")
print(summary_genus_VBA, n=20)

# NSW BIONET DATA WRANGLING -----------------------------------------------
cat("\nProcessing BioNet data...\n")

# Specify folder path 
path_BioNet <- "data/occ/BioNet"

# Create a list of all BioNet data inside the folder
list_BioNet <- list.files(path = path_BioNet, 
                          pattern = "^BioNet", 
                          recursive = TRUE, 
                          full.names = TRUE)

# Read all files into a list of DataFrames and ensure SpeciesCode is a character
df_BioNet <- lapply(list_BioNet, function(file) {
  df <- read.delim(file, 
                   header = TRUE, 
                   sep = "\t", 
                   stringsAsFactors = FALSE)
  df$SpeciesCode <- as.character(df$SpeciesCode) # Convert SpeciesCode to character
  return(df)
}) %>% 
  bind_rows(.) %>%
  filter(ObservationType == "Q") # Q = Camera Trap survey

# column operation
# Create new columns essential columns
# Add year month day columns
df_BioNet = df_BioNet %>%
  mutate(
    # Try parsing DateLast with formats focusing on year and month
    Date = parse_date_time(DateLast, c(
      "dmy_HMS", "dmy", "ymd_HMS", "ymd", "mdy_HMS", "mdy", # With day and optional time
      "ym", "my", # Year and month only
      "y", "m" # Year or month only (less specific, use with caution)
    ), quiet = TRUE),
    # Extract components from the parsed Date column
    year = year(Date),
    month = month(Date),
    day = day(Date),
    # Add a source column and label all values as "NSW"
    Source = "NSW"
  )

# Check for parsing failures
parse_failures_bionet <- df_BioNet %>% filter(is.na(Date) & !is.na(DateLast))
if(nrow(parse_failures_bionet) > 0) {
  warning(paste("Failed to parse", nrow(parse_failures_bionet), "dates in BioNet data"))
  print(head(parse_failures_bionet$DateLast))
}

# Add genus column
df_BioNet$genus <- sub(" .*", "", df_BioNet$ScientificName)

# Re-arranging df and retain essential information
df_BioNet = df_BioNet %>% 
  dplyr::select(Source = Source,
                ProjectID = DatasetName,
                SiteName = Description,
                Year = year,
                Month = month,
                Day = day,
                Genus = genus,
                ScientificName = ScientificName,
                Longitude_GDA94 = Longitude_GDA94,
                Latitude_GDA94 = Latitude_GDA94) %>% 
  na.omit()

# Final filtering
df_BioNet = df_BioNet %>% distinct() %>% 
  filter(!is.na(Longitude_GDA94) & !is.na(Latitude_GDA94)) %>% 
  filter(Year %in% 2012:2021) %>%
  na.omit() # Remove rows with missing values

# Summarize observations by genus
summary_genus_BioNet <- df_BioNet %>%
  group_by(Genus) %>%
  summarise(Observations = n()) %>%
  arrange(desc(Observations))

# View the summary
cat("\nTop genera in BioNet data:\n")
print(summary_genus_BioNet, n=20)

# COMBINED OPERATION ------------------------------------------------
cat("\nCombining datasets...\n")
df = bind_rows(df_VBA, df_BioNet) 
cat("Total records:", nrow(df), "\n")

# Convert PCS to EPSG:3577 
df <- df %>%  
  filter(!is.na(Longitude_GDA94) & !is.na(Latitude_GDA94)) %>% 
  st_as_sf(coords = c("Longitude_GDA94", "Latitude_GDA94"), 
           crs = 4283)  %>% 
  st_transform( , crs = 3577) %>% 
  mutate(
    Easting_3577 = st_coordinates(.)[,1],  # X (Easting)
    Northing_3577 = st_coordinates(.)[,2]  # Y (Northing)
  )  %>%
  st_drop_geometry()  # Remove spatial features to keep as a data frame  

# Recall the beta covariates and stacked into single raster
bio = rast("input/CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
bio = bio[[c("bio5","bio6", "bio12", "bio15")]] #we will only use bio5 and bio 12 based on pearson correlation analysis in S0

terrain = rast("input/env_terrain_EPSG3577.tif")
terrain = terrain[[c("dem", "eastness", "northness", "tpi", "tri", "rock")]]
env_stack <- c(bio,terrain)

# SNAP COORDINATES TO CELL CENTERS ----------------------------------------
cat("\nSnapping coordinates to raster cells...\n")

# Snap all unique site coordinates to cell centers
df_processed <- df %>%
  # First, save the original coordinates
  mutate(
    original_Easting_3577 = Easting_3577,
    original_Northing_3577 = Northing_3577
  ) %>%
  # Create a unique identifier for each location
  group_by(Easting_3577, Northing_3577) %>%
  mutate(location_id = cur_group_id()) %>%
  ungroup()

# Get unique locations
unique_locations <- df_processed %>%
  distinct(location_id, Easting_3577, Northing_3577)

# Create spatial points
points <- vect(unique_locations, 
               geom=c("Easting_3577", "Northing_3577"), 
               crs=crs(env_stack))

# Get cell numbers
cell_numbers <- cells(env_stack, points)[,2]

# Get cell center coordinates
cell_centers <- xyFromCell(env_stack, cell_numbers)

# Add to unique locations
unique_locations$cell_center_x <- cell_centers[,1]
unique_locations$cell_center_y <- cell_centers[,2]

# Join back to main dataframe
df_processed <- df_processed %>%
  left_join(unique_locations %>% 
              select(location_id, cell_center_x, cell_center_y), 
            by="location_id") %>%
  # Replace coordinates with cell centers
  mutate(
    Easting_3577 = cell_center_x,
    Northing_3577 = cell_center_y
  ) %>%
  # Clean up
  select(-cell_center_x, -cell_center_y)

# Check snapping distances
df_processed <- df_processed %>%
  mutate(
    snap_distance = sqrt((original_Easting_3577 - Easting_3577)^2 + 
                           (original_Northing_3577 - Northing_3577)^2)
  )

# Check if any points moved too far
max_distance <- max(df_processed$snap_distance, na.rm = TRUE)
mean_distance <- mean(df_processed$snap_distance, na.rm = TRUE)
cat("Maximum snapping distance:", round(max_distance, 2), "meters\n")
cat("Mean snapping distance:", round(mean_distance, 2), "meters\n")

if(max_distance > 1000) {  # If any point moved more than 1km
  warning("Some points moved more than 1km during snapping!")
}


# Visual inspection
# Create a basic plot of site coordinates
ggplot() + 
  geom_point(data = df_processed, 
             aes(x = Easting_3577, y = Northing_3577), 
             color = "blue", 
             alpha = 0.7) +
  theme_minimal() +
  labs(title = "Spatial Distribution of Sampling Sites",
       x = "Easting (EPSG:3577)",
       y = "Northing (EPSG:3577)")

# Scritps for removing outliers
# Upper left would have small X and large Y values
outlier_idx <- which(df_processed$Northing_3577 == max(df_processed$Northing_3577))
df_processed = df_processed[-outlier_idx, ]

# ADD INFORMATION ---------------------------------------------------------
cat("\nPreparing data structure...\n")

# More robust replicate assignment
df_filtered <- df_processed %>%
  mutate(replicate = case_when(
    is.na(Month) ~ NA_integer_,  # Handle missing months
    Year %in% 2012:2013 ~ 1,
    Year %in% 2014:2015 ~ 2,
    Year %in% 2016:2017 ~ 3,
    Year %in% 2018:2019 ~ 4,
    Year %in% 2020:2021 ~ 5,
    TRUE ~ NA_integer_  # Catch any unexpected values
  )) %>%
  filter(!is.na(replicate)) %>%  # Remove records without valid replicates
  group_by(Easting_3577, Northing_3577) %>%  # Group by plotID
  ungroup() %>% 
  arrange(Easting_3577, Northing_3577) %>%  
  mutate(plotID = dense_rank(paste(Easting_3577, Northing_3577, sep = "_"))) %>% 
  mutate(ProjectID = as.integer(factor(ProjectID, levels = unique(ProjectID))))

write.csv(df_filtered, "input/df_singleSeason.csv", row.names = FALSE)

# Count unique projectID values for random effect
num_unique_projects <- length(unique(df_filtered$ProjectID))
cat("Number of unique projectID values:", num_unique_projects, "\n")

# CREATE SPOCCUPANCY OBJECT -----------------------------------------------
cat("\nCreating occupancy arrays...\n")

# Convert df into long format for detection data
y.long <- df_filtered %>%
  group_by(plotID, replicate, Genus) %>%
  summarize(occ = ifelse(n() > 0, 1, 0), .groups = 'drop') %>%
  ungroup()

# Convert df into long format for effort data  
effort.long <- df_filtered %>%
  group_by(plotID, replicate) %>%
  summarize(effort_days = n_distinct(paste(Year, Month, Day, sep = "_")), .groups = 'drop') %>%
  ungroup()

# Species codes, adjust accordingly
genus.codes <- sort(unique(y.long$Genus))
plot.codes <- sort(unique(y.long$plotID))
replicate.codes <- sort(unique(y.long$replicate))

N <- length(genus.codes)
J <- length(plot.codes)
R <- length(replicate.codes)

cat("Data dimensions:\n")
cat("  Species:", N, "\n")
cat("  Sites:", J, "\n")
cat("  Replicates:", R, "\n")

# Create array for detection-nondetection data
y <- array(NA, dim = c(N, J, R))
dimnames(y)[[1]] <- genus.codes
dimnames(y)[[2]] <- as.character(plot.codes)
dimnames(y)[[3]] <- as.character(replicate.codes)

# Create array for effort data
effort.matrix <- array(NA, dim = c(J, R),
                       dimnames = list(as.character(plot.codes), 
                                       as.character(replicate.codes)))

# Fill the y array
for (j in 1:J) {
  if(j %% 100 == 0) cat("  Processing site", j, "of", J, "\n")
    for (r in 1:R) {
      site <- plot.codes[j]
      rep <- replicate.codes[r]
      
      # Check if this site-replicate combination was sampled
      effort_data <- effort.long %>%
        filter(plotID == site, replicate == rep)
      
      if (nrow(effort_data) > 0) {
        # Site was sampled - populate detection data for all species
        for (i in 1:N) {
          genus <- genus.codes[i]
          
          # Check if this species was observed
          species_data <- y.long %>%
            filter(plotID == site, replicate == rep, Genus == genus)
          
          if (nrow(species_data) > 0) {
            y[i, j, r] <- species_data$occ[1]
          } else {
            y[i, j, r] <- 0  # Species not detected
          }
        }
      }
      # If effort_data is empty, leave as NA (default)
    }
}

# Fill the effort matrix
for (j in 1:J) {
    for (r in 1:R) {
      site <- plot.codes[j]
      rep <- replicate.codes[r]
      
      # Get effort for this site and replicate combination
      val <- effort.long %>%
        filter(plotID == site, replicate == rep) %>%
        pull(effort_days)
      
      # If we found a value, populate the array
      if (length(val) > 0) {
        effort.matrix[j, r] <- val[1]
      }
      # Otherwise leave it as NA (default)
    }
}

# Verify NA patterns match
# Check for the first species (all species should have same NA pattern)
na_y <- is.na(y[1, , ])
na_effort <- is.na(effort.matrix)

# Check if patterns match
cat("\nNA pattern match between y and effort:", identical(na_y, na_effort), "\n")

# If not identical, find where they differ
if (!identical(na_y, na_effort)) {
  diff_positions <- which(na_y != na_effort, arr.ind = TRUE)
  cat("Number of mismatched positions:", nrow(diff_positions), "\n")
  if(nrow(diff_positions) > 0) {
    print(head(diff_positions))
  }
}

# Check the total number of observations for each species
species_detections <- apply(y, 1, sum, na.rm = TRUE)
detection_summary <- data.frame(
  species = genus.codes,
  detections = species_detections,
  proportion = species_detections / sum(!is.na(y[1,,]))
)
cat("\nSpecies detection summary:\n")
print(detection_summary)

# DETECTION RANDOM EFFECTS: PROJECT ---------------------------------------
cat("\nProcessing project information...\n")

# Calculate deployment days for each project within each plotID-year-replicate
project_summary <- df_filtered %>%
  group_by(plotID, replicate, ProjectID) %>%
  summarize(
    deployment_days = n_distinct(paste(Year, Month, Day, sep = "_")),
    .groups = 'drop'
  ) %>%
  # For each plotID-year-replicate, keep only the project with max deployment
  group_by(plotID, replicate) %>%
  slice_max(deployment_days, n = 1, with_ties = FALSE) %>%
  ungroup()

# Create the 3D project array
project.matrix <- array(NA, dim = c(J, R),
                        dimnames = list(as.character(plot.codes), 
                                        as.character(replicate.codes)))

# Fill the project matrix with the longest-deployed project
for (i in 1:nrow(project_summary)) {
  row <- project_summary[i, ]
  j <- which(plot.codes == row$plotID)
  r <- which(replicate.codes == row$replicate)
  
  project.matrix[j, r] <- row$ProjectID
}

# Verify the NA pattern matches your detection and effort arrays
na_project <- is.na(project.matrix)
na_y <- is.na(y[1, , ])
cat("NA pattern match between project and y:", identical(na_project, na_y), "\n")

# Check project matrix statistics
n_unique_projects <- length(unique(as.vector(project.matrix[!is.na(project.matrix)])))
cat("Number of unique projects in matrix:", n_unique_projects, "\n")

# Format site coordinates -------------------------------------------------
cat("\nFormatting coordinates...\n")

coords <- df_filtered %>%
  distinct(plotID, .keep_all = TRUE) %>%
  select(Easting_3577, Northing_3577)

# Convert coords into matrix
coords = as.matrix(coords)

# Ensure the matrix is numeric
storage.mode(coords) <- "numeric"

# Assign sequential row names (plotID) as characters
rownames(coords) <- as.character(1:nrow(coords))

# Rename columns
colnames(coords) <- c("X", "Y")

# Verify the structure
cat("Coordinate dimensions:", dim(coords), "\n")

# Format occurrence covariates --------------------------------------------------
cat("\nExtracting environmental covariates...\n")

# Recall the beta covariates
bio = rast("input/CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
bio = bio[[c("bio5","bio6", "bio12", "bio15")]]

terrain = rast("input/env_terrain_EPSG3577.tif")
terrain = terrain[[c("dem", "eastness", "northness", "tpi", "tri", "rock")]]

env_stack <- c(bio,terrain)

beta <- terra::extract(env_stack, coords)

# Check environmental covariate extraction
if(any(is.na(beta))) {
  warning("Some sites have NA values for environmental covariates")
  na_summary <- sapply(beta, function(x) sum(is.na(x)))
  cat("\nNA summary for environmental covariates:\n")
  print(na_summary[na_summary > 0])
}

# Create the occ.covs list
occ.covs <- data.frame(
  bio5 = beta$bio5,
  bio6 = beta$bio6,
  bio12 = beta$bio12,
  bio15 = beta$bio15,
  dem = beta$dem,
  eastness = beta$eastness,
  northness = beta$northness,
  tpi = beta$tpi,
  tri = beta$tri,
  rock = beta$rock
)

# CLEAN UP & PACK INTO OBJECT ---------------------------------------------
cat("\nFinal data cleaning and packaging...\n")

# # Identify rows with NA values in occ.covs
env_vars <- names(occ.covs)

# Check for NAs in the environmental variables
na_check <- sapply(occ.covs[env_vars], function(x) is.na(x))
rows_with_na <- which(rowSums(na_check) > 0)
rows_with_na

if(length(rows_with_na) > 0) {
  cat("Rows with NA values:", paste(rows_with_na, collapse = ", "), "\n")
  cat("Removing", length(rows_with_na), "sites with NA environmental covariates\n")
  
  # Clean other matrices and arrays
  occ.covs <- occ.covs[-rows_with_na, , drop = FALSE]
  effort_clean <- effort.matrix[-rows_with_na, , drop = FALSE]
  project_clean <- project.matrix[-rows_with_na, , drop = FALSE]
  y_clean <- y[, -rows_with_na, , drop = FALSE]
  coords_clean <- coords[-rows_with_na, ,drop = FALSE]
  
} else {
  cat("No NA values found in environmental covariates\n")
  effort_clean <- effort.matrix
  project_clean <- project.matrix
  y_clean <- y
  coords_clean <- coords
}

# Pack all things into a list object
det.covs <- list(
  effort = effort_clean,
  project = project_clean
)

# Create the final data list
data <- list(
  y = y_clean,
  occ.covs = occ.covs,
  det.covs = det.covs,
  coords = coords_clean
)




# Check the structures
cat("\nFinal data dimensions:\n")
cat("Y dimensions:", dim(y_clean), "\n")
cat("Coords dimensions:", dim(coords_clean), "\n")
cat("Effort dimensions:", dim(effort_clean), "\n")
cat("Project dimensions:", dim(project_clean), "\n")

# Verify everything is aligned
cat("\nAlignment check:\n")
cat("Sites in y:", dim(y_clean)[2], "\n")
cat("Sites in coords:", nrow(coords_clean), "\n")
cat("Sites in occ.covs:", length(occ.covs$bio5), "\n")
cat("Sites in det.covs effort:", dim(effort_clean)[1], "\n")
cat("Sites in det.covs project:", dim(project_clean)[1], "\n")

# Comprehensive validation function
validate_occupancy_data <- function(data) {
  cat("\n=== Data Validation ===\n")
  
  # Check dimensions
  y_dims <- dim(data$y)
  cat("Y dimensions:", y_dims, "\n")
  cat("Expected: [species, sites, replicates]\n")
  
  # Check NA patterns
  na_pattern <- apply(data$y, c(2,3), function(x) all(is.na(x)))
  cat("Sites with all NAs:", sum(apply(na_pattern, 1, all)), "\n")
  cat("Replicates with all NAs:", sum(apply(na_pattern, 2, all)), "\n")
  
  # Check detection covariates
  effort_dims <- dim(data$det.covs$effort)
  project_dims <- dim(data$det.covs$project)
  cat("Effort dimensions:", effort_dims, "\n")
  cat("Project dimensions:", project_dims, "\n")
  
  # Check alignment
  cat("\nDimension alignment check:\n")
  cat("Sites match:", y_dims[2] == nrow(data$coords), "\n")
  # cat("Replicates match:", y_dims[2] == n(data$occ.covs$bio5), "\n")
  cat("Effort aligns:", all(effort_dims == y_dims[2:3]), "\n")
  
  # Check for reasonable values
  cat("\nValue range checks:\n")
  cat("Y values:", unique(as.vector(data$y[!is.na(data$y)])), "\n")
  cat("Effort range:", range(data$det.covs$effort, na.rm = TRUE), "\n")
  
  # Check species with zero detections
  species_detections <- apply(data$y, 1, sum, na.rm = TRUE)
  zero_detection_species <- dimnames(data$y)[[1]][species_detections == 0]
  if(length(zero_detection_species) > 0) {
    cat("\nWarning: Species with zero detections:", 
        paste(zero_detection_species, collapse = ", "), "\n")
  }
  
  return(invisible(NULL))
}

# Run validation
validate_occupancy_data(data)

# Save with metadata
metadata <- list(
  date_created = Sys.Date(),
  n_species = N,
  n_sites = nrow(coords_clean),
  n_replicates = R,
  species = genus.codes,
  data_sources = c("VBA", "BioNet"),
  coordinate_system = "EPSG:3577"
)

cat("\nSaving data...\n")
save(data, metadata, file = "input/data_singleSeason.RData")
cat("\nData preparation complete!\n")
cat("Output saved to: input/data_singleSeason.RData\n")


# [OPTIONAL] VISUALIZE SITE -----------------------------------------------
# Convert coordinates to sf object for better mapping
sites_sf <- st_as_sf(as.data.frame(data$coords), 
                     coords = c("X", "Y"), 
                     crs = 3577)

# Plot with spatial context (if you have state boundaries)
ggplot() +
  # Add state boundaries if available
  # geom_sf(data = aus_states, fill = NA, color = "gray70") +
  geom_sf(data = sites_sf, color = "blue", size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Spatial Distribution of Sampling Sites",
       subtitle = paste0("n = ", nrow(data$coords), " sites")) 

# [OPTIONAL] TEST RUN -----------------------------------------------------

# Test run
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


test.sfMsPGOcc <- sfMsPGOcc(
  occ.formula = occ.formula,
  det.formula = det.formula,
  data = data,
  # inits = inits,
  n.batch = 25,
  batch.length = 10,
  accept.rate = 0.43,
  # priors = priors,
  n.factors = 1,
  cov.model = "exponential",
  # tuning = tuning,
  # n.omp.threads = n.omp.thr,
  verbose = 50,
  NNGP = TRUE,
  n.neighbors = 10,
  # n.report = n.repo,
  # n.burn = n.burn,
  n.thin = 1,
  n.chains = 1)

summary(test.sfMsPGOcc, level = 'community')

