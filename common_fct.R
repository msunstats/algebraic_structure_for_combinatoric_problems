# --- common_fct.R: Common Functions and Global Setup for Patient Subgroup Discovery ---

# Load common libraries
suppressMessages(library(dplyr))
suppressMessages(library(rlang))
suppressMessages(library(ggplot2))
suppressMessages(library(fixtheworld)) # Assuming these are from your original environment
suppressMessages(library(rice)) # Assuming these are from your original environment
suppressMessages(library(easyENTIM)) # Assuming these are from your original environment
suppressMessages(library(tidyverse))
suppressMessages(library(aws.s3))
suppressMessages(library(doSNOW))
suppressMessages(library(foreach))
suppressMessages(library(broom))
suppressMessages(library(progress)) # Used for progress bar in GA, but loaded here for consistency

# Set ggplot2 theme
theme_set(theme_bw())

# Custom function to set S3 environment variables
# This function is used to configure AWS S3 environment variables, often for Arvados integration.
setup_environment <- function() {
  # Set S3 environment variables using arvupload package.
  arvupload::set_s3_env_vars()
}

# Custom Hamming Distance Calculation (Moved from bo_method.R to common_fct.R)
# Calculates the Hamming distance between two binary vectors.
calculate_hamming_distance <- function(x, y) {
  if (length(x) != length(y)) {
    stop("Vectors must be of the same length.")
  }
  sum(x != y)
}

# --- Common Analysis Results Generation ---
# This function generates and optionally prints analysis results common to both BO and GA.
# It summarizes metric frequencies, operator distributions for numeric metrics,
# and selected levels for categorical metrics from the identified rules.
generate_analysis_results_common <- function(all_solutions_details_df_sorted, params) {
  # Filter for active rules that have valid metric, operator, and value/levels
  active_rules_for_analysis_df <- all_solutions_details_df_sorted %>%
    filter(!is.na(Metric) & Metric != "" & !is.na(Operator) & Operator != "")
  
  metric_frequency_df <- NULL
  numeric_operator_distribution_df <- NULL
  numeric_cutoff_data_df <- NULL
  categorical_level_frequency_list <- list()
  
  if (nrow(active_rules_for_analysis_df) > 0) {
    # Conditional printing based on verbose_analysis_output parameter
    if (params$verbose_analysis_output) {
      message("\n--- Analysis Results ---")
      message("\n1. Frequency of Each Metric in Defining Rules:")
      metric_frequency_df <- active_rules_for_analysis_df %>%
        group_by(Metric) %>%
        summarise(Count = n(), .groups = 'drop') %>%
        arrange(desc(Count)) %>%
        mutate(Percentage = (Count / sum(Count)) * 100)
      print(metric_frequency_df)
      
      message("\n2. Analysis of Numeric Metrics (Operators and Cutoffs):")
      numeric_cutoff_data_df <- active_rules_for_analysis_df %>%
        filter(Metric %in% params$numeric_metrics) %>%
        mutate(Cutoff = as.numeric(Value_or_Levels))
      
      if (nrow(numeric_cutoff_data_df) > 0) {
        numeric_operator_distribution_df <- numeric_cutoff_data_df %>%
          group_by(Metric, Operator) %>%
          summarise(Count = n(), .groups = 'drop') %>%
          pivot_wider(names_from = Operator, values_from = Count, values_fill = 0) %>%
          mutate(Total_Rules_for_Metric = rowSums(select(., -Metric))) %>%
          arrange(Metric)
        
        message("\n  - Operator Distribution for Numeric Metrics:")
        print(numeric_operator_distribution_df)
        
      } else {
        message("    No numeric metrics were found in defining rules to perform detailed analysis.\n")
      }
      
      message("\n3. Analysis of Categorical Metrics (Selected Levels):")
      categorical_rules_for_analysis_df <- active_rules_for_analysis_df %>%
        filter(Metric %in% params$categorical_metrics)
      
      if (nrow(categorical_rules_for_analysis_df) > 0) {
        all_selected_combinations_by_metric <- list()
        
        for (i in 1:nrow(categorical_rules_for_analysis_df)) {
          row <- categorical_rules_for_analysis_df[i, ]
          metric <- row$Metric
          levels_str <- row$Value_or_Levels
          
          # Parse the levels string (e.g., "Level1, Level2") into an R vector
          # This assumes Value_or_Levels is already formatted as "Level1, Level2"
          # For single levels, it will be just "Level1"
          parsed_levels <- tryCatch({
            # Split by comma and trim whitespace, then unlist to get a character vector
            trimws(unlist(strsplit(levels_str, ",")))
          }, error = function(e) {
            warning(paste("Could not parse levels string for", metric, ":", levels_str, "Error:", e$message))
            return(character(0))
          })
          
          if (length(parsed_levels) > 0) {
            # Standardize combination string for consistent counting (e.g., "Level1, Level2")
            standardized_combo_str <- paste(sort(parsed_levels), collapse = ", ")
            all_selected_combinations_by_metric[[metric]] <- c(all_selected_combinations_by_metric[[metric]], standardized_combo_str)
          }
        }
        
        if (length(all_selected_combinations_by_metric) > 0) {
          for (metric_name in names(all_selected_combinations_by_metric)) {
            combo_counts <- table(all_selected_combinations_by_metric[[metric_name]])
            categorical_level_frequency_list[[metric_name]] <- as.data.frame(combo_counts) %>%
              rename(Combination = Var1, Count = Freq) %>%
              arrange(desc(Count))
            message(sprintf("    Metric: %s", metric_name))
            print(categorical_level_frequency_list[[metric_name]])
            message("\n")
          }
        } else {
          message("    No combinations of categorical levels were successfully parsed from defining rules.")
        }
        
      } else {
        message("    No categorical metrics were used in defining rules to perform detailed analysis.\n")
      }
    } else { # If not verbose, still compute but don't print
      metric_frequency_df <- active_rules_for_analysis_df %>%
        group_by(Metric) %>%
        summarise(Count = n(), .groups = 'drop') %>%
        arrange(desc(Count)) %>%
        mutate(Percentage = (Count / sum(Count)) * 100)
      
      numeric_cutoff_data_df <- active_rules_for_analysis_df %>%
        filter(Metric %in% params$numeric_metrics) %>%
        mutate(Cutoff = as.numeric(Value_or_Levels))
      
      if (nrow(numeric_cutoff_data_df) > 0) {
        numeric_operator_distribution_df <- numeric_cutoff_data_df %>%
          group_by(Metric, Operator) %>%
          summarise(Count = n(), .groups = 'drop') %>%
          pivot_wider(names_from = Operator, values_from = Count, values_fill = 0) %>%
          mutate(Total_Rules_for_Metric = rowSums(select(., -Metric))) %>%
          arrange(Metric)
      }
      
      categorical_rules_for_analysis_df <- active_rules_for_analysis_df %>%
        filter(Metric %in% params$categorical_metrics)
      
      if (nrow(categorical_rules_for_analysis_df) > 0) {
        all_selected_combinations_by_metric <- list()
        for (i in 1:nrow(categorical_rules_for_analysis_df)) {
          row <- categorical_rules_for_analysis_df[i, ]
          metric <- row$Metric
          levels_str <- row$Value_or_Levels
          parsed_levels <- tryCatch({ trimws(unlist(strsplit(levels_str, ","))) }, error = function(e) { character(0) })
          if (length(parsed_levels) > 0) {
            standardized_combo_str <- paste(sort(parsed_levels), collapse = ", ")
            all_selected_combinations_by_metric[[metric]] <- c(all_selected_combinations_by_metric[[metric]], standardized_combo_str)
          }
        }
        if (length(all_selected_combinations_by_metric) > 0) {
          for (metric_name in names(all_selected_combinations_by_metric)) {
            combo_counts <- table(all_selected_combinations_by_metric[[metric_name]])
            categorical_level_frequency_list[[metric_name]] <- as.data.frame(combo_counts) %>%
              rename(Combination = Var1, Count = Freq) %>%
              arrange(desc(Count))
          }
        }
      }
    }
  } else {
    message("\nNo active rules available for detailed analysis.")
  }
  
  return(list(
    metric_frequency = metric_frequency_df,
    numeric_operator_distribution = numeric_operator_distribution_df,
    numeric_cutoff_data = numeric_cutoff_data_df,
    categorical_level_frequencies = categorical_level_frequency_list
  ))
}

# --- Common Plot Fitness vs. Minimum Subgroup Size ---
# This function generates a plot showing the optimal fitness value achieved
# for different minimum subgroup sizes.
plot_fitness_vs_subgroup_size_common <- function(all_results_raw, hv_aval_mean, params) {
  if (!isTRUE(params$benchmark_mode)) {
    message("\n--- Generating Plots ---")
  }
  # Ensure the raw results are treated as a list of lists for consistent processing
  valid_results <- lapply(all_results_raw, function(x) x)
  
  if (length(valid_results) > 0) {
    # Create a dataframe for plotting from the valid results
    plot_data <- data.frame(
      min_subgroup_size = sapply(valid_results, function(x) x$min_subgroup_size),
      optimal_fitness = sapply(valid_results, function(x) x$optimal_fitness)
    )
  } else {
    plot_data <- data.frame(min_subgroup_size = numeric(0), optimal_fitness = numeric(0))
    warning("No valid results were obtained from optimization runs for plotting fitness vs. subgroup size.")
  }
  
  if (nrow(plot_data) > 0) {
    # Calculate first and second derivatives of fitness with respect to subgroup size
    plot_data <- plot_data %>%
      mutate(dFit_dN = c(0, diff(optimal_fitness)) / c(0, diff(min_subgroup_size)),
             d2Fit_dN = c(0, diff(dFit_dN)) / c(0, diff(min_subgroup_size)))
    
    # Construct the plot title
    plot_title <- paste0("Optimal Fitness vs. Minimum Subgroup Size (HV Mean ", params$variable_for_use, ": ", round(hv_aval_mean, 4), ")")
    
    # Create the ggplot object
    p1 <- ggplot(plot_data, aes(x = min_subgroup_size, y = optimal_fitness)) +
      geom_line(color = "blue", size = 1) + # Line connecting the points
      geom_point(color = "red", size = 3) + # Points for each min_subgroup_size
      labs(
        title = plot_title,
        x = "Minimum Subgroup Size",
        y = paste0("Optimal Fitness Value")
      ) +
      theme_minimal() + # Minimal theme for a clean look
      # Set breaks for x-axis to show specific min_subgroup_size values
      scale_x_continuous(breaks = seq(min(params$min_subgroup_sizes_to_test), max(params$min_subgroup_sizes_to_test), by = 5))
    print(p1)
  } else {
    message("Skipping plot generation due to no valid plot data for fitness vs. subgroup size.")
  }
}

# --- Common Generate Numeric Cutoff Plots ---
# This function generates plots showing the distribution of cutoff values for numeric metrics
# used in the identified optimal rules.
generate_numeric_cutoff_plots_common <- function(analysis_results, params) {
  if (params$verbose_analysis_output) { # Only generate plots if verbose output is enabled
    message("\n  - Distribution of Cutoff Values for Numeric Metrics (see plots - Combined by Operator Color, with Transparency):")
  }
  numeric_cutoff_data_df <- analysis_results$numeric_cutoff_data
  
  if (is.null(numeric_cutoff_data_df) || nrow(numeric_cutoff_data_df) == 0) {
    if (params$verbose_analysis_output) {
      message("    No numeric metrics were found in defining rules to plot cutoff distributions.\n")
    }
    return(invisible(NULL)) # Return invisibly if no data to plot
  }
  
  # Aggregate counts for each unique combination of Metric, Cutoff, and Operator
  plot_data_aggregated_counts <- numeric_cutoff_data_df %>%
    group_by(Metric, Cutoff, Operator) %>%
    summarise(TotalFrequency = n(), .groups = 'drop')
  
  # Generate a plot for each unique numeric metric
  for (metric in unique(plot_data_aggregated_counts$Metric)) {
    metric_data_filtered_for_plot <- plot_data_aggregated_counts %>% filter(Metric == metric)
    
    if (nrow(metric_data_filtered_for_plot) == 0) {
      if (params$verbose_analysis_output) {
        message(sprintf("    No data for %s to plot frequency dots.\n", metric))
      }
      next # Skip to the next metric if no data
    }
    
    p <- ggplot(metric_data_filtered_for_plot, aes(x = Cutoff, y = TotalFrequency, color = Operator)) +
      geom_segment(aes(xend = Cutoff, yend = 0), linewidth = 0.5, alpha = 0.7, show.legend = FALSE) + # Vertical lines
      geom_point(size = 3, alpha = 0.7) + # Points at each cutoff frequency
      labs(
        title = paste("Distribution of Cutoff Values for", metric, "by Operator"),
        x = paste0(metric, " Cutoff Value"),
        y = "Frequency",
        color = "Operator"
      ) +
      theme_minimal() +
      # Ensure y-axis breaks are integers for frequency counts
      scale_y_continuous(breaks = function(x) unique(floor(pretty(x))))
    
    if (params$verbose_analysis_output) { # Only print plot if verbose output is enabled
      print(p)
      message("\n") # Add a newline after each plot for better separation
    }
  }
  if (params$verbose_analysis_output) {
    message("    These plots show the exact frequency of each unique cutoff value chosen by the optimization algorithm for each numeric metric, colored by operator, with transparency to reveal stacked points.\n")
  }
}

# --- Common Generate Final Rules Table ---
# This function processes the sorted solutions details into a concise table of optimal rules,
# including statistical significance.
generate_final_rules_table_common <- function(all_solutions_details_df_sorted, params) {
  if (!isTRUE(params$benchmark_mode)) {
    message("\n--- Optimal Rules Table ---")
  }
  
  # Group by MinSubgroupSize, p-value, and OptimalFitness to consolidate rules
  dfRulesNumber <- all_solutions_details_df_sorted %>%
    group_by(MinSubgroupSize, pValVsHv, OptimalFitness) %>%
    summarise(n = n(), # Count of rules for this subgroup
              # Custom rule string generation with conditional formatting for display
              rules = paste(
                Metric,
                Operator,
                # Conditionally format Value_or_Levels for display
                case_when(
                  Operator == "==" ~ gsub("'", "", Value_or_Levels), # Strip quotes for '=='
                  Operator == "%in%" ~ paste0("(", Value_or_Levels, ")"), # Add parentheses for '%in%'
                  TRUE ~ Value_or_Levels # Default for numeric operators
                ),
                collapse = ", "
              ), .groups = 'drop') %>%
    ungroup() %>%
    arrange(n) %>%
    # Apply Holm's method for p-value adjustment (if p-values exist)
    { if(any(!is.na(.$pValVsHv))) mutate(., pValAdjVsHv = p.adjust(pValVsHv, method = "holm")) else mutate(., pValAdjVsHv = NA_real_) } %>%
    # Filter for statistically significant rules (adjusted p-value < 0.05) or if p-value was NA
    filter(pValAdjVsHv < 0.05 | is.na(pValVsHv))
  
  # Further filter to show only rules with 1 or 2 predicates for simplicity (can be adjusted)
  dfRulesFinal <- dfRulesNumber %>%
    filter(n %in% c(1,2))
  
  # Print explanatory messages for the table
  if (!isTRUE(params$benchmark_mode)) {
    message("\nThe table below gives the final rules in selecting the subgroup of patients. They are equally good from data analysis point of view. The ultimate decision should be based on the balance of patient recruitment difficulty and the subgroup size.")
    message("\nHow to read the table:")
    message("- MinSubgroupSize: Patient size of the subgroup of interest")
    message("- n: Number of rules that define this subgroup")
    message(sprintf("- OptimalFitness: Fitness value of the subgroup. In this case defined as fold change of %s of protein %s in the subgroup vs healthy volunteers", params$variable_for_use, params$thisProtein))
    message("- pValVsHv: Raw p value of the subgroup vs healthy volunteers")
    message("- pValAdjVsHv: Adjusted p value of the subgroup vs healthy volunteers")
  }
  
  if (nrow(dfRulesFinal) > 0) {
    # Print the final rules table, rounding numeric values and selecting relevant columns
    if (!isTRUE(params$benchmark_mode)) {
      print(dfRulesFinal %>%
              mutate_if(is.numeric, round, digits = 4) %>%
              select(MinSubgroupSize, n, OptimalFitness, pValVsHv, pValAdjVsHv, rules) %>%
              arrange(rules))
    }
  } else {
    if (!isTRUE(params$benchmark_mode)) {
      message("No significant rules (pValAdjVsHv < 0.05) with 1 or 2 predicates were found for the optimal subgroups.")
    }
  }
  return(dfRulesFinal)
}

# --- Common Plot Protein Expression in Subgroups ---
# This function generates a visualization of protein expression levels within the identified
# subgroups compared to healthy volunteers.
plot_protein_expression_in_subgroups_common <- function(dfRulesFinal, patient_groups_df_sorted, df, params) {
  if (!isTRUE(params$benchmark_mode)) {
    message(sprintf("\n--- Visualization of Differentially Expressed %s in Subgroups ---", params$thisProtein))
    message("Corresponding visualisation of how differentially expressed the subgroup is compared to healthy volunteers, based on the rules chosen above. The y-axis is log2 transformed and the black line is the mean of each subgroup.")
  }
  
  if (nrow(dfRulesFinal) == 0) {
    if (!isTRUE(params$benchmark_mode)) {
      message("Skipping protein expression plot as no final rules were identified.")
    }
    return(invisible(NULL))
  }
  
  # Join patient group assignments with the full clinical data
  dfSubgroupAssigned <- patient_groups_df_sorted %>%
    rename(subgroupAssign = Group) %>%
    distinct() %>% # Ensure unique patient-subgroup assignments
    filter(MinSubgroupSize %in% dfRulesFinal$MinSubgroupSize) %>% # Only include min_subgroup_sizes that have final rules
    inner_join(., df, by = "USUBJID") # Join with original data to get protein values and other metrics
  
  if (nrow(dfSubgroupAssigned) == 0) {
    message("No patients assigned to the final subgroups for plotting.")
    return(invisible(NULL))
  }
  
  # Determine the true unit for the protein variable for plot labeling
  unitThisProtein <- unique(dfSubgroupAssigned$LBSTRESU[dfSubgroupAssigned$LBTESTCD == params$thisProtein])
  trueUnit <- case_when(params$variable_for_use == 'AVAL_BY_LENGTH' ~ paste0(unitThisProtein, '/mm'),
                        params$variable_for_use == 'AVAL_BY_PROTEIN' ~ paste0(unitThisProtein, ''),
                        TRUE ~ NA_character_)
  if (length(trueUnit) == 0 || is.na(trueUnit)) trueUnit <- "Unit" # Fallback unit if not determined
  
  # Generate the plot
  # Ensure variable_for_use is a character string
  if (!is.character(params$variable_for_use)) {
    stop("params$variable_for_use must be a character string, but got: ", class(params$variable_for_use)[1])
  }
  
  # Check if the variable exists in the data
  if (!params$variable_for_use %in% colnames(dfSubgroupAssigned)) {
    stop("Variable '", params$variable_for_use, "' not found in data. Available columns: ", paste(colnames(dfSubgroupAssigned), collapse = ", "))
  }
  
  p2 <- ggplot(dfSubgroupAssigned %>% filter(LBTESTCD == params$thisProtein), # Filter for the specific protein
               aes(x = subgroupAssign, y = .data[[params$variable_for_use]], col = DED_SEVERITY)) +
    geom_jitter(width = 0.2, alpha = 0.5) + # Jittered points to show individual patient values
    scale_y_continuous(transform = 'log2', # Log2 transform the y-axis
                       breaks = scales::breaks_log(n = 10), # Logarithmic breaks
                       labels = scales::label_number()) + # Numeric labels
    facet_wrap(~ MinSubgroupSize, labeller = 'label_both') + # Create separate plots for each min_subgroup_size
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), width = 0.7, size = 0.5, color = "black") + # Add mean lines
    labs(y = paste0(params$variable_for_use, ' (', trueUnit, ')'),
         x = 'Subgroup Assignment',
         title = paste0('Differentially Expressed ', params$thisProtein, ' in Subgroups vs Healthy Volunteers')) +
    theme(plot.title = element_text(hjust = 0.5)) # Center the plot title
  print(p2)
}

# --- Helper Functions for Exhaustive Search ---
# Function to check if exhaustive search is computationally feasible
check_exhaustive_search_feasibility <- function(categorical_metrics, gene_map, max_combinations = 1e6) {
  if (is.null(categorical_metrics) || length(categorical_metrics) == 0) {
    return(list(feasible = FALSE, reason = "No categorical metrics", search_space_size = 0))
  }
  
  search_space_size <- calculate_search_space_size(categorical_metrics, gene_map)
  
  feasible <- search_space_size <= max_combinations
  
  return(list(
    feasible = feasible,
    search_space_size = search_space_size,
    max_allowed = max_combinations,
    reason = ifelse(feasible, "Feasible", "Search space too large")
  ))
}
