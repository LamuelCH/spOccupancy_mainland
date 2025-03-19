library(sf)
library(terra)
library(ggplot2)

# Filtersf# Filter to keep only Dasyuridae
df_dasyuridae <- df %>%
  filter(ScientificName == "Dasyurus maculatus") %>% 
  distinct(x, y, .keep_all = TRUE)


# Save as a new dataframe
df_dasyuridae <- as.data.frame(df_dasyuridae)

# View the new dataframe
head(df_dasyuridae)

# Save the dataframe as a CSV file
write.csv(df_dasyuridae, "data/occ/df_dasyuridae_unique.csv", row.names = FALSE)


bioNet_STQ = read.csv("data/occ/df_dasyuridae_unique.csv")
Atlas_STQ = read.csv("data/occ/Atlas_STQ_2012-2021/Atlas_STQ_2012-2021.csv")


Atlas_STQ = Atlas_STQ %>%
  select(species, decimalLatitude, decimalLongitude, geodeticDatum)

# Convert to sf object
Atlas_STQ_sf <- st_as_sf(Atlas_STQ, 
                         coords = c("decimalLongitude", "decimalLatitude"), 
                         crs = 4326)  # WGS84

bioNet_STQ_sf = st_as_sf(bioNet_STQ,
                         coords = c("x", "y"))

# Transform to EPSG:32755
Atlas_STQ_sf <- st_transform(Atlas_STQ_sf, crs = 32755)

# View the transformed data
head(Atlas_STQ_sf)
head(bioNet_STQ_sf)

# Recall the beta covariates and stacked into single raster
bio = rast("data/beta/CHELSA_bio_2011-2040_gfdl-esm4_ssp126_V.2.1_EPSG32755.tif")
bio = bio[[c("bio5", "bio12")]] #we will only use bio5 and bio 12 based on pearson correlation analysis in S0

roads_density = rast("data/beta/env_roadDensity_EPSG32755.tif")
pd = rast("data/beta/aus_pd_2010-2020_1km_UNadj_EPSG32755.tif")
der = rast("data/beta/env_der_EPSG32755.tif")

env_stack <- c(bio,roads_density, pd, der)

beta_Atlas <- terra::extract(env_stack, Atlas_STQ_sf, df = TRUE)
beta_BioNet = terra::extract(env_stack, bioNet_STQ_sf, df = TRUE)













# Merge both datasets into "Total"
beta_Total <- bind_rows(beta_Atlas, beta_BioNet)


# Standardize values
beta_Total_scaled <- as.data.frame(scale(beta_Total))
beta_BioNet_scaled <- as.data.frame(scale(beta_BioNet))

# beta_Total_scaled <- as.data.frame(beta_Total)
# beta_BioNet_scaled <- as.data.frame(beta_BioNet)


# Add dataset labels for comparison
beta_Total_scaled$Dataset <- "Total (Atlas + BioNet)"
beta_BioNet_scaled$Dataset <- "BioNet"

# Combine into one dataset
beta_combined <- bind_rows(beta_Total_scaled, beta_BioNet_scaled)

# Convert to long format
beta_long <- beta_combined %>%
  pivot_longer(cols = -Dataset, names_to = "Variable", values_to = "Value")


# Density Plot (for distribution comparison)
ggplot(beta_long, aes(x = Value, fill = Dataset)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal() +
  labs(title = "Comparison of Beta BioNet vs Total",
       x = "Value", 
       y = "Density") +
  scale_fill_manual(values = c("Total (Atlas + BioNet)" = "blue", "BioNet" = "red"))



ggplot(beta_long, aes(x = Dataset, y = Value, fill = Dataset)) +
  geom_boxplot(alpha = 0.6) +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal() +
  labs(title = "Box Plot Comparison of Beta BioNet vs Total",
       x = "Dataset", 
       y = "Value") +
  scale_fill_manual(values = c("Total (Atlas + BioNet)" = "blue", "BioNet" = "red"))













# Boxplot comparison
ggplot(beta_long, aes(x = Variable, y = Value, fill = Dataset)) +
  geom_boxplot(outlier.color = "red", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Comparison of Standardized Environmental Variables",
       x = "Environmental Variable",
       y = "Standardized Value (Z-score)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Density plot comparison
ggplot(beta_long, aes(x = Value, color = Dataset, fill = Dataset)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~Variable, scales = "free") +
  theme_minimal() +
  labs(title = "Density Plot Comparison: BioNet vs. Total",
       x = "Standardized Value (Z-score)",
       y = "Density")
