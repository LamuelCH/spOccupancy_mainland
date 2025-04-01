library(dplyr)

# Analysis of replicate counts per site
count_replicates_per_site <- function(df) {
  # Group by site (plotID)
  replicate_counts <- df %>%
    group_by(plotID) %>%
    summarize(
      num_replicates = n_distinct(replicate)
    ) %>%
    ungroup()
  
  # Create a frequency table of sites by number of replicates
  result_table <- table(replicate_counts$num_replicates)
  
  # Calculate percentages
  total_sites <- nrow(replicate_counts)
  percentages <- prop.table(result_table) * 100
  
  # Calculate cumulative counts (sites with at least X replicates)
  cumulative_counts <- sapply(1:5, function(x) {
    sum(replicate_counts$num_replicates >= x)
  })
  
  # Print results
  cat("Distribution of sites by number of replicates:\n")
  for (i in 1:5) {
    count <- sum(replicate_counts$num_replicates == i)
    percentage <- count / total_sites * 100
    cat(sprintf("%d replicate(s): %d sites (%.2f%%)\n", i, count, percentage))
  }
  
  cat("\nCumulative distribution of sites:\n")
  for (i in 1:5) {
    count <- cumulative_counts[i]
    percentage <- count / total_sites * 100
    cat(sprintf("At least %d replicate(s): %d sites (%.2f%%)\n", i, count, percentage))
  }
  
  # Return the detailed count data for further analysis
  invisible(list(
    counts = replicate_counts,
    distribution = result_table,
    percentages = percentages,
    cumulative_counts = cumulative_counts,
    total_sites = total_sites
  ))
}





# Analysis of data distribution across replicates
analyze_data_by_replicate <- function(df) {
  # Count observations per replicate
  replicate_summary <- df %>%
    group_by(replicate) %>%
    summarize(
      num_sites = n_distinct(plotID),
      num_observations = n()
    ) %>%
    ungroup()
  
  # Print results
  cat("Data summary by replicate:\n")
  print(replicate_summary)
  
  # Return the summary
  invisible(replicate_summary)
}




# Create an analysis to examine how different grouping strategies would affect the data
analyze_alternative_groupings <- function(df, years_span = 2012:2021) {
  cat("Analysis of alternative temporal grouping strategies:\n")
  
  # Original grouping: 2 years per replicate (5 replicates)
  cat("\n1. ORIGINAL GROUPING (2 years per replicate, 5 replicates)\n")
  df$replicate_orig <- ((df$Year - min(years_span)) %/% 2) + 1
  
  # Alternative 1: 3 years per replicate (approx. 3 replicates)
  cat("\n2. ALTERNATIVE 1 (3-4 years per replicate, 3 replicates)\n")
  # Create custom year groups: 2012-2014, 2015-2017, 2018-2021
  df$replicate_alt1 <- case_when(
    df$Year %in% 2012:2014 ~ 1,
    df$Year %in% 2015:2017 ~ 2,
    df$Year %in% 2018:2021 ~ 3,
    TRUE ~ NA_integer_
  )
  
  # Alternative 2: 5 years per replicate (2 replicates)  
  cat("\n3. ALTERNATIVE 2 (uneven years per replicate, 4 replicates)\n")
  df$replicate_alt2 <- case_when(
    df$Year %in% 2012:2016 ~ 1,
    df$Year %in% 2017:2021 ~ 2,
    TRUE ~ NA_integer_
  )
  
  # Function to analyze a specific grouping
  analyze_grouping <- function(df, replicate_col) {
    # Count sites by number of replicates
    site_replicates <- df %>%
      group_by(plotID) %>%
      summarize(
        num_replicates = n_distinct(!!sym(replicate_col)),
        .groups = "drop"
      )
    
    # Get the maximum number of replicates in this grouping
    max_replicates <- max(site_replicates$num_replicates, na.rm = TRUE)
    
    # Distribution
    rep_table <- table(site_replicates$num_replicates)
    total_sites <- nrow(site_replicates)
    
    # Print distribution
    cat("Distribution of sites by number of replicates:\n")
    for (i in 1:max_replicates) {
      count <- sum(site_replicates$num_replicates == i, na.rm = TRUE)
      pct <- count / total_sites * 100
      cat(sprintf("%d replicate(s): %d sites (%.2f%%)\n", i, count, pct))
    }
    
    # Print sites with at least X replicates
    cat("\nCumulative distribution:\n")
    for (i in 1:max_replicates) {
      count <- sum(site_replicates$num_replicates >= i, na.rm = TRUE)
      pct <- count / total_sites * 100
      cat(sprintf("At least %d replicate(s): %d sites (%.2f%%)\n", i, count, pct))
    }
    
    # Count observations per replicate
    rep_summary <- df %>%
      group_by(!!sym(replicate_col)) %>%
      summarize(
        num_sites = n_distinct(plotID),
        num_observations = n(),
        .groups = "drop"
      )
    
    cat("\nData summary by replicate:\n")
    print(rep_summary)
    
    invisible(list(site_replicates = site_replicates, rep_summary = rep_summary))
  }
  
  # Analyze each grouping strategy
  cat("\n===== ORIGINAL GROUPING (2 years, 5 replicates) =====\n")
  orig_results <- analyze_grouping(df, "replicate_orig")
  
  cat("\n===== ALTERNATIVE 1 (3-4 years, 3 replicates) =====\n")
  alt1_results <- analyze_grouping(df, "replicate_alt1")
  
  cat("\n===== ALTERNATIVE 2 (5 years, 2 replicates) =====\n")
  alt2_results <- analyze_grouping(df, "replicate_alt2")
  
  # Return results for further analysis
  invisible(list(
    original = orig_results,
    alternative1 = alt1_results,
    alternative2 = alt2_results
  ))
}

# Example usage (you'll run this with your actual data):
# 
# # First analyze the current replicate structure
# replicate_analysis <- count_replicates_per_site(df_filtered)
# data_analysis <- analyze_data_by_replicate(df_filtered)
# 
# # Then explore alternative grouping strategies
# alternative_analysis <- analyze_alternative_groupings(df_filtered)

df_filtered = read.csv("input/df.csv")
# First analyze your current replicate structure
replicate_analysis <- count_replicates_per_site(df_filtered)
data_analysis <- analyze_data_by_replicate(df_filtered)

# Then explore alternative grouping strategies
alternative_analysis <- analyze_alternative_groupings(df_filtered)




# Spatial plot ------------------------------------------------------------
library(sf)
library(ggplot2)
library(patchwork) # For combining plots

# Function to visualize spatial distribution of sites by number of replicates
visualize_spatial_replicates <- function(df) {
  # Ensure we have the spatial coordinates
  if(!all(c("Easting_3577", "Northing_3577", "plotID", "replicate") %in% names(df))) {
    stop("Required columns not found in dataframe: Easting_3577, Northing_3577, plotID, replicate")
  }
  
  # Count the number of distinct replicates per site
  site_replicates <- df %>%
    group_by(plotID, Easting_3577, Northing_3577) %>%
    summarize(
      num_replicates = n_distinct(replicate),
      .groups = "drop"
    )
  
  # Convert to factor for coloring
  site_replicates$num_replicates_factor <- factor(site_replicates$num_replicates, 
                                                  levels = 1:5,
                                                  labels = paste(1:5, "replicate(s)"))
  
  # Create the spatial plot
  p1 <- ggplot(site_replicates, aes(x = Easting_3577, y = Northing_3577, color = num_replicates_factor)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_viridis_d(name = "Number of\nReplicates") +
    theme_minimal() +
    labs(
      title = "Spatial Distribution of Sites by Number of Replicates",
      x = "Easting (m)",
      y = "Northing (m)"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )
  
  # Create a density plot to show the spatial density of sites
  p2 <- ggplot(site_replicates, aes(x = Easting_3577, y = Northing_3577)) +
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_viridis_c(name = "Density") +
    theme_minimal() +
    labs(
      title = "Spatial Density of Sampling Sites",
      x = "Easting (m)",
      y = "Northing (m)"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )
  
  return(list(point_plot = p1, density_plot = p2))
}

# Function to visualize alternate replicate groupings spatially
visualize_alternative_groupings <- function(df, years_span = 2012:2021) {
  # Create alternative groupings
  df_analysis <- df %>%
    mutate(
      # Original grouping (2 years per replicate)
      replicate_orig = ((Year - min(years_span)) %/% 2) + 1,
      
      # Alternative 1: 3 years per replicate (approx. 3 replicates)
      replicate_alt1 = case_when(
        Year %in% 2012:2014 ~ 1,
        Year %in% 2015:2017 ~ 2,
        Year %in% 2018:2021 ~ 3,
        TRUE ~ NA_integer_
      ),
      
      # Alternative 2: 5 years per replicate (2 replicates)
      replicate_alt2 = case_when(
        Year %in% 2012:2016 ~ 1,
        Year %in% 2017:2021 ~ 2,
        TRUE ~ NA_integer_
      )
    )
  
  # Function to create a plot for a specific grouping
  create_replicate_plot <- function(df, replicate_col, title_prefix) {
    # Count the number of distinct replicates per site
    site_replicates <- df %>%
      group_by(plotID, Easting_3577, Northing_3577) %>%
      summarize(
        num_replicates = n_distinct(!!sym(replicate_col)),
        .groups = "drop"
      )
    
    # Convert to factor for coloring
    max_reps <- max(site_replicates$num_replicates, na.rm = TRUE)
    site_replicates$num_replicates_factor <- factor(site_replicates$num_replicates, 
                                                    levels = 1:max_reps,
                                                    labels = paste(1:max_reps, "replicate(s)"))
    
    # Create the spatial plot
    ggplot(site_replicates, aes(x = Easting_3577, y = Northing_3577, color = num_replicates_factor)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_viridis_d(name = "Number of\nReplicates") +
      theme_minimal() +
      labs(
        title = paste(title_prefix, "Grouping"),
        subtitle = paste("Sites with ≥2 replicates:", sum(site_replicates$num_replicates >= 2, na.rm = TRUE)),
        x = "Easting (m)",
        y = "Northing (m)"
      ) +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title = element_text(size = 10)
      )
  }
  
  # Create plots for each grouping
  p_orig <- create_replicate_plot(df_analysis, "replicate_orig", "Original (5 replicates)")
  p_alt1 <- create_replicate_plot(df_analysis, "replicate_alt1", "Alternative 1 (3 replicates)")
  p_alt2 <- create_replicate_plot(df_analysis, "replicate_alt2", "Alternative 2 (2 replicates)")
  
  # Combine plots
  combined_plot <- p_orig + p_alt1 + p_alt2 + plot_layout(ncol = 1)
  
  return(list(original = p_orig, alternative1 = p_alt1, alternative2 = p_alt2, combined = combined_plot))
}

# Visualize the change in spatial coverage over time
visualize_temporal_change <- function(df) {
  # Create separate dataframes for each replicate
  rep_list <- list()
  for (i in unique(df$replicate)) {
    rep_data <- df %>% 
      filter(replicate == i) %>%
      distinct(plotID, Easting_3577, Northing_3577)
    
    rep_data$Replicate <- paste("Replicate", i)
    rep_list[[i]] <- rep_data
  }
  
  # Combine all replicates
  all_reps <- do.call(rbind, rep_list)
  
  # Create plot
  p <- ggplot(all_reps, aes(x = Easting_3577, y = Northing_3577)) +
    geom_point(size = 2, alpha = 0.7) +
    facet_wrap(~ Replicate, ncol = 3) +
    theme_minimal() +
    labs(
      title = "Spatial Distribution of Sampling Sites by Replicate",
      x = "Easting (m)",
      y = "Northing (m)"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  return(p)
}

# Example usage (to run with your actual data):
# plots <- visualize_spatial_replicates(df_filtered)
# print(plots$point_plot)
# print(plots$density_plot)
#
# # Visualize temporal change
# temporal_plot <- visualize_temporal_change(df_filtered)
# print(temporal_plot)
#
# # Visualize alternative groupings
# alt_plots <- visualize_alternative_groupings(df_filtered)
# print(alt_plots$combined)




# CREATE PLOT -------------------------------------------------------------

# Load necessary libraries
library(sf)
library(ggplot2)
library(patchwork)
library(dplyr)

# Basic spatial visualization
plots <- visualize_spatial_replicates(df_filtered)
print(plots$point_plot)
print(plots$density_plot)

# Visualize change in spatial coverage across replicates
temporal_plot <- visualize_temporal_change(df_filtered)
print(temporal_plot)

# Compare alternative replicate grouping strategies
alt_plots <- visualize_alternative_groupings(df_filtered)
print(alt_plots$alternative2)



