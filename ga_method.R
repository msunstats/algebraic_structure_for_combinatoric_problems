# --- ga_method.R: Genetic Algorithm (GA) Method for Patient Subgroup Discovery ---

# Load GA-specific libraries
suppressMessages(library(GA)) # Core GA library

# --- Function: Define Chromosome Encoding and Gene Map for GA ---
# This function defines how clinical tests (numeric and categorical) are represented
# as genes within the GA chromosome, and creates a mapping to decode them later.
define_ga_encoding <- function(de_data, numeric_metrics, categorical_metrics) {
  ga_min <- c() # Minimum value for each gene
  ga_max <- c() # Maximum value for each gene
  ga_type <- c() # Type ("binary" or "real") for each gene
  gene_map <- list() # To map gene indices back to metric names/levels for decoding
  
  current_gene_idx <- 1 # Tracks the current gene index for sequential assignment
  
  # Encode Numeric Metrics: Each numeric metric is represented by 3 genes:
  # 1. Binary: 0 (don't use this metric), 1 (use this metric)
  # 2. Binary: 0 (operator is <=), 1 (operator is >)
  # 3. Real: The threshold value for the metric
  for (metric in numeric_metrics) {
    # Define min/max for the threshold gene based on the data range
    min_val_metric <- min(de_data[[metric]], na.rm = TRUE)
    max_val_metric <- max(de_data[[metric]], na.rm = TRUE)
    
    # Add genes for 'use', 'operator type', and 'threshold value'
    ga_min <- c(ga_min, 0, 0, min_val_metric)
    ga_max <- c(ga_max, 1, 1, max_val_metric)
    ga_type <- c(ga_type, "binary", "binary", "real")
    # Store mapping information for this metric
    gene_map[[metric]] <- list(
      use_idx = current_gene_idx,
      op_idx = current_gene_idx + 1,
      val_idx = current_gene_idx + 2,
      type = "numeric"
    )
    current_gene_idx <- current_gene_idx + 3 # Move to the next available gene index
  }
  
  # Encode Categorical Metrics: Each categorical metric is represented by 1 + K genes:
  # 1. Binary: 0 (don't use this metric), 1 (use this metric)
  # 2. K binary genes: 0 (don't select this level), 1 (select this level) for each unique level
  categorical_levels_map <- list() # Stores unique levels for each categorical metric
  for (metric in categorical_metrics) {
    levels <- unique(de_data[[metric]]) %>% na.omit() %>% as.character() %>% sort()
    if (length(levels) == 0) {
      warning(paste("Skipping categorical metric", metric, "due to no valid levels."))
      next # Skip if no valid levels are found for the metric
    }
    categorical_levels_map[[metric]] <- levels
    levels_count <- length(levels)
    # Add genes for 'use' and for each 'level selection'
    ga_min <- c(ga_min, 0, rep(0, levels_count))
    ga_max <- c(ga_max, 1, rep(1, levels_count))
    ga_type <- c(ga_type, "binary", rep("binary", levels_count))
    # Store mapping information for this metric
    gene_map[[metric]] <- list(
      use_idx = current_gene_idx,
      levels_start_idx = current_gene_idx + 1,
      levels_end_idx = current_gene_idx + levels_count,
      levels = levels,
      type = "categorical"
    )
    current_gene_idx <- current_gene_idx + (1 + levels_count) # Move to the next available gene index
  }
  
  return(list(
    ga_min = ga_min,
    ga_max = ga_max,
    ga_type = ga_type,
    gene_map = gene_map,
    categorical_levels_map = categorical_levels_map # Included for completeness, though gene_map is primary
  ))
}

# --- Function: Fitness Function Generator for GA ---
# This function generates a specific fitness function for the Genetic Algorithm.
# The generated function takes a chromosome (rule combination) and returns its fitness value.
# Fitness is based on protein expression in the defined subgroup relative to healthy volunteers.
generate_ga_fitness_function <- function(current_min_subgroup_size, de_data, hv_aval_mean, numeric_metrics, categorical_metrics, gene_map, protein_variable_name) {
  # The actual fitness function that GA will call
  function(chromosome) {
    current_filter_expressions <- list() # Stores dplyr filter expressions
    
    # Decode Numeric Rules based on the chromosome
    for (metric in numeric_metrics) {
      map_info <- gene_map[[metric]]
      if (is.null(map_info)) next # Skip if metric not in gene_map (e.g. if it had no valid levels)
      
      use_metric <- round(chromosome[map_info$use_idx]) # Round to 0 or 1 for binary gene
      
      if (use_metric == 1) { # If this metric is selected to be used in a rule
        op_type <- round(chromosome[map_info$op_idx]) # Operator type (0 for <=, 1 for >)
        threshold_val <- chromosome[map_info$val_idx] # Real-valued threshold
        
        if (op_type == 0) { # Operator is <=
          current_filter_expressions[[length(current_filter_expressions) + 1]] <- rlang::expr(!!sym(metric) <= !!threshold_val)
        } else { # Operator is >
          current_filter_expressions[[length(current_filter_expressions) + 1]] <- rlang::expr(!!sym(metric) > !!threshold_val)
        }
      }
    }
    
    # Decode Categorical Rules based on the chromosome
    for (metric in categorical_metrics) {
      map_info <- gene_map[[metric]]
      if (is.null(map_info)) next # Skip if metric not in gene_map
      
      use_metric <- round(chromosome[map_info$use_idx])
      
      if (use_metric == 1) { # If this categorical metric is selected to be used
        selected_levels <- c()
        # Check which levels are selected (binary genes for each level)
        for (i in 1:length(map_info$levels)) {
          level_gene_val <- round(chromosome[map_info$levels_start_idx + i - 1])
          if (level_gene_val == 1) {
            selected_levels <- c(selected_levels, map_info$levels[i])
          }
        }
        
        if (length(selected_levels) == 0) {
          # If a categorical metric is "used" (use_metric == 1) but no levels are selected,
          # this rule is invalid, so return 0 fitness.
          return(0)
        }
        
        # Add the filter expression for selected levels
        current_filter_expressions[[length(current_filter_expressions) + 1]] <- rlang::expr(!!sym(metric) %in% !!selected_levels)
      }
    }
    
    # If no rules are active (no metric is selected), return a very small positive fitness.
    # This prevents errors and allows the GA to explore more options.
    if (length(current_filter_expressions) == 0) {
      return(1e-9)
    }
    
    # Apply the generated filters to the disease patient data sequentially
    filtered_de_data <- de_data
    for (filter_expr in current_filter_expressions) {
      filtered_de_data <- filtered_de_data %>% filter(!!filter_expr)
      # If the subgroup size drops below the minimum required at any point, return a low fitness.
      if (nrow(filtered_de_data) < current_min_subgroup_size) {
        return(1e-6)
      }
    }
    
    # Final check on subgroup size after all filters are applied
    if (nrow(filtered_de_data) < current_min_subgroup_size) {
      return(1e-6)
    }
    
    # Calculate the mean protein level for the identified subgroup
    subgroup_aval_mean <- mean(filtered_de_data[[protein_variable_name]], na.rm = TRUE)
    
    # Handle cases where subgroup mean is NaN or non-positive (e.g., all NAs or zero values)
    if (is.nan(subgroup_aval_mean) || subgroup_aval_mean <= 0) {
      return(1e-9)
    }
    
    # Calculate fitness as the ratio of subgroup mean to HV mean
    fitness_val <- subgroup_aval_mean / hv_aval_mean
    
    # Penalize fitness values less than 1 (meaning subgroup mean is lower than HV mean)
    # This encourages finding subgroups with *higher* protein levels than HV.
    if (fitness_val < 1.0) {
      return(fitness_val * 0.5) # Example: halve the fitness for values below 1
    }
    
    return(fitness_val) # Return the calculated fitness value
  }
}


# --- Function: Run Genetic Algorithm for a single min_subgroup_size ---
# This function executes the Genetic Algorithm for a specific minimum subgroup size.
run_genetic_algorithm_for_min_size <- function(current_min_size, ga_params, de_data, hv_aval_mean,
                                               numeric_metrics, categorical_metrics, gene_map,
                                               ga_min, ga_max, ga_type, protein_variable_name) {
  
  # Generate the fitness function specific to the current minimum subgroup size
  current_fitness_function <- generate_ga_fitness_function(
    current_min_size, de_data, hv_aval_mean, numeric_metrics, categorical_metrics, gene_map, protein_variable_name
  )
  
  # Run the Genetic Algorithm using the 'GA' package
  ga_run_results <- GA::ga(
    type = "real-valued", # Chromosomes contain both binary and real values
    fitness = current_fitness_function, # The objective function to maximize
    min = ga_min, # Minimum values for each gene
    max = ga_max, # Maximum values for each gene
    popSize = ga_params$ga_pop_size, # Number of solutions in each generation
    maxiter = ga_params$ga_max_iter, # Maximum number of generations
    run = ga_params$ga_run_limit, # Number of generations without improvement before stopping
    pcrossover = 0.8, # Probability of crossover operation
    pmutation = 0.1, # Probability of mutation operation
    monitor = FALSE, # Do not print progress to console for individual runs (handled by progress bar in orchestrator)
    optim = TRUE, # Attempt to optimize the best solution found at the end
    seed = NULL # Let parallel workers use their own seeds for independent runs
  )
  
  # Return the results of this GA run: min_subgroup_size, best fitness, and the best chromosome
  list(
    min_subgroup_size = current_min_size,
    optimal_fitness = max(ga_run_results@fitness), # The best fitness value found
    best_chromosome = ga_run_results@solution[1, ] # The chromosome that achieved the best fitness (first row if multiple)
  )
}

# --- Function: Orchestrate Parallel GA Runs ---
# Manages the execution of multiple Genetic Algorithm runs in parallel,
# each for a different minimum subgroup size.
orchestrate_parallel_ga <- function(min_subgroup_sizes_to_test, ga_params, de_data, hv_aval_mean,
                                    numeric_metrics, categorical_metrics, gene_map,
                                    ga_min, ga_max, ga_type, protein_variable_name, cl_outer) {
  
  if (!isTRUE(ga_params$benchmark_mode)) {
    message("\n--- Running Genetic Algorithm for different min_subgroup_size values (in parallel) ---\n")
  }
  
  # Register the provided cluster for parallel execution
  registerDoSNOW(cl_outer)
  
  # Export all necessary objects and functions to the workers in the cluster.
  # This ensures each worker has access to the data and functions needed for its task.
  clusterExport(cl_outer, c(
    "run_genetic_algorithm_for_min_size", "generate_ga_fitness_function",
    "de_data", "hv_aval_mean", "numeric_metrics", "categorical_metrics", "gene_map",
    "ga_min", "ga_max", "ga_type", "ga_params", "protein_variable_name"
  ), envir = environment()) # Export from the current environment
  
  # Load required packages within each worker's R session
  clusterEvalQ(cl_outer, {
    library(dplyr)
    library(GA)
    library(rlang)
    # set.seed(Sys.getpid()) # Uncomment if strict reproducibility across parallel runs is not needed
  })
  
  # Setup a progress bar for monitoring the parallel execution
  total_tasks <- length(min_subgroup_sizes_to_test)
  pb <- progress_bar$new(
    format = "[:bar] :percent ETA: :eta (:elapsed)",
    total = total_tasks,
    clear = FALSE,
    width = 60
  )
  progress_update_fun <- function(n) {
    pb$tick()
  }
  
  # Execute the GA runs in parallel using foreach and %dopar%
  all_results_ga_raw <- foreach(
    current_min_size = min_subgroup_sizes_to_test,
    .packages = c("dplyr", "GA", "rlang"), # Packages needed by each worker
    .verbose = FALSE, # Suppress verbose output from foreach
    .errorhandling = 'pass', # If an error occurs in a worker, pass the error object instead of stopping
    .options.snow = list(progress = progress_update_fun) # Pass progress update function to doSNOW
  ) %dopar% {
    # Call the encapsulated GA function for the current min_subgroup_size
    result <- run_genetic_algorithm_for_min_size(
      current_min_size = current_min_size,
      ga_params = ga_params,
      de_data = de_data,
      hv_aval_mean = hv_aval_mean,
      numeric_metrics = numeric_metrics,
      categorical_metrics = categorical_metrics,
      gene_map = gene_map,
      ga_min = ga_min,
      ga_max = ga_max,
      ga_type = ga_type,
      protein_variable_name = protein_variable_name
    )
    return(result)
  }
  # Note: progress_bar objects do not have a close method, so close(pb) is removed.
  if (!isTRUE(ga_params$benchmark_mode)) {
    message("\nAll Genetic Algorithm runs complete.")
  }
  
  # Post-processing: Order and name results for consistency and easier access
  if (is.list(all_results_ga_raw) && length(all_results_ga_raw) > 0) {
    # Check if results are valid and contain 'min_subgroup_size'
    if (all(sapply(all_results_ga_raw, is.list)) && all(sapply(all_results_ga_raw, function(x) !is.null(x$min_subgroup_size)))) {
      result_names <- sapply(all_results_ga_raw, function(x) as.character(x$min_subgroup_size))
      ordered_results <- all_results_ga_raw[order(as.numeric(result_names))]
      names(ordered_results) <- sort(as.numeric(result_names))
      all_results_ga_raw <- ordered_results
    } else {
      warning("Results structure from parallel processing is unexpected. Cannot assign names based on min_subgroup_size.")
    }
  } else {
    warning("all_results_ga_raw is not a list or is empty after parallel processing.")
  }
  
  return(all_results_ga_raw)
}

# --- Function: Process GA Results into Detailed DataFrames ---
# This function decodes the optimal chromosomes from the GA runs
# into human-readable rules and calculates various statistics for each identified subgroup.
# It also assigns patients to "Subgroup", "Rest_DE", and "HV" groups.
process_ga_results <- function(all_results_ga_raw, de_data, hv_data, hv_aval_mean,
                               numeric_metrics, categorical_metrics, gene_map, params) {
  
  # Initialize dataframes to store processed results
  all_solutions_details_df <- data.frame(
    MinSubgroupSize = numeric(),
    OptimalFitness = numeric(),
    SubgroupAVALMean = numeric(),
    DiffFromHV = numeric(),
    PercentDiffFromHV = numeric(),
    NumPatientsInSubgroup = numeric(),
    pValVsHv = numeric(),
    Metric = character(),
    Operator = character(),
    Value_or_Levels = character(),
    stringsAsFactors = FALSE
  )
  
  patient_groups_df <- data.frame(
    MinSubgroupSize = numeric(),
    USUBJID = character(),
    Group = character(), # This column will store "Subgroup", "Rest_DE", "HV"
    stringsAsFactors = FALSE
  )
  
  # Calculate HV USUBJIDs once as they are constant for all solutions
  hv_usubjid_list <- unique(hv_data$USUBJID)
  # This will track which MinSubgroupSizes have already had their HV/Rest_DE added to patient_groups_df
  processed_min_subgroup_sizes_for_hv_rest <- c()
  
  if (!isTRUE(params$benchmark_mode)) {
    message("\n--- Decoding and Analyzing Optimal Rules ---")
  }
  # Ensure the raw results are treated as a list of lists for consistent processing
  valid_results_ga <- lapply(all_results_ga_raw, function(x) x)
  
  if (length(valid_results_ga) > 0) {
    for (i in seq_along(valid_results_ga)) {
      result_entry <- valid_results_ga[[i]]
      current_min_subgroup_size <- result_entry$min_subgroup_size
      current_optimal_fitness <- result_entry$optimal_fitness
      current_best_chromosome_ga <- result_entry$best_chromosome
      
      rules_table_for_current_solution <- data.frame(
        Metric = character(),
        Operator = character(),
        Value_or_Levels = character(),
        stringsAsFactors = FALSE
      )
      final_filter_expressions_apply_i <- list()
      
      if (is.null(current_best_chromosome_ga)) {
        warning(paste("No best chromosome found for min subgroup size:", current_min_subgroup_size))
        next # Skip to the next result if no chromosome was found
      }
      
      # Decode Numeric Rules and add to temporary rules table
      for (metric in numeric_metrics) {
        map_info <- gene_map[[metric]]
        if (is.null(map_info)) next # Skip if metric was not encoded (e.g., no valid levels)
        
        use_metric <- round(current_best_chromosome_ga[map_info$use_idx])
        
        if (use_metric == 1) {
          op_type <- round(current_best_chromosome_ga[map_info$op_idx])
          threshold_val <- current_best_chromosome_ga[map_info$val_idx]
          operator <- ifelse(op_type == 0, "<=", ">")
          
          rules_table_for_current_solution <- rbind(rules_table_for_current_solution, data.frame(
            Metric = metric, Operator = operator, Value_or_Levels = as.character(round(threshold_val, 4)),
            stringsAsFactors = FALSE
          ))
          
          if (op_type == 0) {
            final_filter_expressions_apply_i[[length(final_filter_expressions_apply_i) + 1]] <- rlang::expr(!!sym(metric) <= !!threshold_val)
          } else {
            final_filter_expressions_apply_i[[length(final_filter_expressions_apply_i) + 1]] <- rlang::expr(!!sym(metric) > !!threshold_val)
          }
        }
      }
      
      # Decode Categorical Rules and add to temporary rules table
      for (metric in categorical_metrics) {
        map_info <- gene_map[[metric]]
        if (is.null(map_info)) next # Skip if metric was not encoded
        
        use_metric <- round(current_best_chromosome_ga[map_info$use_idx])
        
        if (use_metric == 1) {
          selected_levels <- c()
          for (j in 1:length(map_info$levels)) {
            level_gene_val <- round(current_best_chromosome_ga[map_info$levels_start_idx + j - 1])
            if (level_gene_val == 1) {
              selected_levels <- c(selected_levels, map_info$levels[j])
            }
          }
          
          if (length(selected_levels) == 0) {
            next # If no levels are selected for a "used" categorical metric, skip this rule
          }
          
          rules_table_for_current_solution <- rbind(rules_table_for_current_solution, data.frame(
            Metric = metric, Operator = "%in%", Value_or_Levels = paste0(selected_levels, collapse = ", "), # Store as comma-separated string
            stringsAsFactors = FALSE
          ))
          
          final_filter_expressions_apply_i[[length(final_filter_expressions_apply_i) + 1]] <- rlang::expr(!!sym(metric) %in% !!selected_levels)
        }
      }
      
      # --- Apply rules and calculate statistics and identify patient groups for this current solution ---
      subgroup_aval_mean_current <- NA
      diff_from_hv_current <- NA
      percent_diff_from_hv_current <- NA
      num_patients_current <- 0
      pVal <- NA
      
      current_subgroup_usubjids <- character(0)
      current_rest_de_usubjids <- character(0)
      
      # Only proceed if there are active rules to apply
      if (nrow(rules_table_for_current_solution) > 0) {
        temp_filtered_de_data <- de_data
        is_valid_subgroup_for_ids <- TRUE
        
        for (filter_expr in final_filter_expressions_apply_i) {
          temp_filtered_de_data <- temp_filtered_de_data %>% filter(!!filter_expr)
          if (nrow(temp_filtered_de_data) < current_min_subgroup_size) {
            is_valid_subgroup_for_ids <- FALSE
            break # Break early if subgroup size constraint is violated
          }
        }
        
        # If the subgroup is valid and meets the minimum size
        if (is_valid_subgroup_for_ids && nrow(temp_filtered_de_data) >= current_min_subgroup_size) {
          subgroup_aval_mean_current <- mean(temp_filtered_de_data[[params$variable_for_use]], na.rm = TRUE)
          
          if (!is.nan(subgroup_aval_mean_current) && subgroup_aval_mean_current > 0) {
            diff_from_hv_current <- subgroup_aval_mean_current - hv_aval_mean
            percent_diff_from_hv_current <- (diff_from_hv_current / hv_aval_mean) * 100
            num_patients_current <- nrow(temp_filtered_de_data)
            
            # Perform Wilcoxon test for statistical significance
            pVal_test <- tryCatch({
              wilcox.test(temp_filtered_de_data[[params$variable_for_use]],
                          hv_data[[params$variable_for_use]],
                          alternative = 'greater', # Test if subgroup mean is GREATER than HV mean
                          paired = FALSE)
            }, error = function(e) {
              warning(paste("Wilcoxon test failed for min_subgroup_size", current_min_subgroup_size, ":", e$message))
              return(NULL) # Return NULL if test fails
            })
            
            if (!is.null(pVal_test)) {
              pVal <- broom::tidy(pVal_test) %>% select(pVal = p.value) %>% pull(pVal)
            }
            
            # Identify USUBJIDs for the current subgroup and the rest of DE patients
            current_subgroup_usubjids <- unique(temp_filtered_de_data$USUBJID)
            all_de_usubjids <- unique(de_data$USUBJID)
            current_rest_de_usubjids <- all_de_usubjids[!all_de_usubjids %in% current_subgroup_usubjids]
            
          } else {
            # If subgroup mean is invalid, clear rules table for this solution
            rules_table_for_current_solution <- rules_table_for_current_solution[0,]
          }
        } else {
          # If subgroup is not valid or too small, clear rules table
          rules_table_for_current_solution <- rules_table_for_current_solution[0,]
        }
      } else {
        # If no rules were generated, clear rules table
        rules_table_for_current_solution <- rules_table_for_current_solution[0,]
      }
      
      # --- Append to all_solutions_details_df (rules table) ---
      # Only add if valid rules were found and subgroup size constraint is met
      if (nrow(rules_table_for_current_solution) > 0 && num_patients_current >= current_min_subgroup_size) {
        common_info <- data.frame(
          MinSubgroupSize = current_min_subgroup_size,
          OptimalFitness = current_optimal_fitness,
          SubgroupAVALMean = subgroup_aval_mean_current,
          DiffFromHV = diff_from_hv_current,
          PercentDiffFromHV = percent_diff_from_hv_current,
          NumPatientsInSubgroup = num_patients_current,
          pValVsHv = pVal,
          stringsAsFactors = FALSE
        )
        
        # Replicate common info for each rule in the current solution and bind
        solution_rows_for_main_df <- cbind(common_info[rep(1, nrow(rules_table_for_current_solution)), ],
                                           rules_table_for_current_solution)
        
        all_solutions_details_df <- rbind(all_solutions_details_df, solution_rows_for_main_df)
      }
      
      # --- Append to patient_groups_df (long format, 3 columns) ---
      # Only add patient groups if the solution successfully formed a valid subgroup AND has actual members
      if(num_patients_current >= current_min_subgroup_size && length(current_subgroup_usubjids) > 0) {
        # Add Subgroup patients
        patient_groups_df <- rbind(patient_groups_df, data.frame(
          MinSubgroupSize = current_min_subgroup_size,
          USUBJID = current_subgroup_usubjids,
          Group = "Subgroup",
          stringsAsFactors = FALSE
        ))
        
        # Add Rest_DE patients (those DE patients not in the subgroup)
        if(length(current_rest_de_usubjids) > 0) {
          patient_groups_df <- rbind(patient_groups_df, data.frame(
            MinSubgroupSize = current_min_subgroup_size,
            USUBJID = current_rest_de_usubjids,
            Group = "Rest_DE",
            stringsAsFactors = FALSE
          ))
        }
        
        # Add HV patients (only once per unique MinSubgroupSize, since they are constant for that size)
        if(!(current_min_subgroup_size %in% processed_min_subgroup_sizes_for_hv_rest) && length(hv_usubjid_list) > 0) {
          patient_groups_df <- rbind(patient_groups_df, data.frame(
            MinSubgroupSize = current_min_subgroup_size,
            USUBJID = hv_usubjid_list,
            Group = "HV",
            stringsAsFactors = FALSE
          ))
          # Mark this MinSubgroupSize as processed for HV/Rest_DE to prevent future duplications
          processed_min_subgroup_sizes_for_hv_rest <- c(processed_min_subgroup_sizes_for_hv_rest, current_min_subgroup_size)
        }
      }
    } # End of for loop through valid_results
  } # End of if(length(valid_results) > 0)
  
  # Sort the final dataframes for readability and consistency
  all_solutions_details_df_sorted <- all_solutions_details_df %>%
    arrange(desc(OptimalFitness), MinSubgroupSize, Metric)
  
  patient_groups_df_sorted <- patient_groups_df %>%
    arrange(MinSubgroupSize, Group, USUBJID) # Sort for readability
  
  return(list(
    all_solutions_details_df_sorted = all_solutions_details_df_sorted,
    patient_groups_df_sorted = patient_groups_df_sorted
  ))
}


# --- Main Orchestrating Function for GA ---
# This is the main function to run the Genetic Algorithm patient subgroup discovery process.
# It orchestrates chromosome encoding, parallel GA runs, results processing, and visualization.
run_genetic_algorithm_analysis <- function(params, df, dfThisProtein, hv_data, de_data, hv_aval_mean, cl_outer = NULL) {
  
  # 1. Define Chromosome Encoding and Gene Map for GA
  # This step determines how clinical tests (numeric and categorical) are represented
  # as genes in the GA chromosome.
  ga_encoding_info <- define_ga_encoding(de_data, params$numeric_metrics, params$categorical_metrics)
  ga_min <- ga_encoding_info$ga_min
  ga_max <- ga_encoding_info$ga_max
  ga_type <- ga_encoding_info$ga_type
  gene_map <- ga_encoding_info$gene_map
  # categorical_levels_map is also returned but not directly used in subsequent functions,
  # its information is embedded in gene_map.
  
  all_results_ga_raw <- list()
  
  # 2. Run or Load Genetic Algorithm Results
  if (params$runSim) {
    # If runSim is TRUE, the GA simulation is executed in parallel for various minimum subgroup sizes.
    all_results_ga_raw <- orchestrate_parallel_ga(
      min_subgroup_sizes_to_test = params$min_subgroup_sizes_to_test,
      ga_params = params, # Pass all GA-related parameters (pop_size, max_iter, run_limit)
      de_data = de_data,
      hv_aval_mean = hv_aval_mean,
      numeric_metrics = params$numeric_metrics,
      categorical_metrics = params$categorical_metrics,
      gene_map = gene_map,
      ga_min = ga_min,
      ga_max = ga_max,
      ga_type = ga_type,
      protein_variable_name = params$variable_for_use,
      cl_outer = cl_outer # Pass the cluster object for parallel execution
    )
    
    if (!isTRUE(params$benchmark_mode)) {
      message("\n--- Simulation results (not saved as requested) ---")
    }
    
  } else {
    # Load pre-computed results
    # If runSim is FALSE, attempt to load previously saved GA results.
    if (!isTRUE(params$benchmark_mode)) {
      message("\n--- Loading pre-computed simulation results ---")
    }
    tryCatch({
      all_results_ga_raw <- aws.s3::s3read_using(FUN = readRDS,
                                                 object = paste0('data/', 'OptSubgroup_Rules_GA_vs_HV_', params$thisProtein, '_', params$variable_for_use, '.rds'),
                                                 bucket = params$.arv_save$collection)
    }, error = function(e) {
      warning(paste("Error loading pre-computed GA data from S3:", e$message, ". Generating dummy results for demonstration."))
      # Fallback: Generate dummy results if loading fails.
      # This requires defining a dummy gene_map and chromosome length if not already available.
      if (!exists("gene_map") || is.null(gene_map)) {
        dummy_encoding <- define_ga_encoding(de_data, params$numeric_metrics, params$categorical_metrics)
        gene_map_dummy <- dummy_encoding$gene_map
        chromosome_length <- length(dummy_encoding$ga_min)
      } else {
        chromosome_length <- length(ga_min)
      }
      # Dummy results for demonstration purposes
      all_results_ga_raw <- list(
        "10" = list(min_subgroup_size = 10, optimal_fitness = 1.5, best_chromosome = sample(0:1, chromosome_length, replace=TRUE)),
        "11" = list(min_subgroup_size = 11, optimal_fitness = 1.6, best_chromosome = sample(0:1, chromosome_length, replace=TRUE)),
        "12" = list(min_subgroup_size = 12, optimal_fitness = 1.7, best_chromosome = sample(0:1, chromosome_length, replace=TRUE))
      )
      # Name the list elements by min_subgroup_size for consistency
      names(all_results_ga_raw) <- sapply(all_results_ga_raw, function(x) as.character(x$min_subgroup_size))
    })
  }
  
  # Stop execution if no GA results are available after simulation or loading.
  if (length(all_results_ga_raw) == 0) {
    stop("No Genetic Algorithm results available for further processing. Check data loading or simulation parameters.")
  }
  
  # 3. Process GA Results into Detailed DataFrames
  # This function decodes the optimal chromosomes into human-readable rules and calculates
  # various statistics for each identified subgroup. It also assigns patients to groups.
  processed_dfs_ga <- process_ga_results(all_results_ga_raw, de_data, hv_data, hv_aval_mean,
                                         params$numeric_metrics, params$categorical_metrics, gene_map, params)
  all_solutions_details_df_sorted_ga <- processed_dfs_ga$all_solutions_details_df_sorted
  patient_groups_df_sorted_ga <- processed_dfs_ga$patient_groups_df_sorted
  
  # 4. Generate Analysis Results (metrics frequencies, etc.)
  # This provides an overview of which clinical tests and their cutoffs/levels were most frequently
  # used in the optimal subgroup definitions across all tested minimum subgroup sizes.
  analysis_results_ga <- generate_analysis_results_common(all_solutions_details_df_sorted_ga, params)
  
  # 5. Plot Fitness vs. Minimum Subgroup Size
  # Visualizes the relationship between the minimum subgroup size constraint and the
  # best fitness achieved by the GA.
  if (!isTRUE(params$benchmark_mode)) {
    plot_fitness_vs_subgroup_size_common(all_results_ga_raw, hv_aval_mean, params)
    
    # 6. Generate Numeric Cutoff Plots (if applicable)
    # Provides visualizations of the distribution of cutoff values for numeric metrics
    # used in the optimal rules.
    if (params$verbose_analysis_output) { # Only plot if verbose output is enabled
      generate_numeric_cutoff_plots_common(analysis_results_ga, params)
    }
  }
  
  # 7. Generate Final Rules Table
  # Presents a concise table of the most relevant optimal rules, including their
  # associated statistics and adjusted p-values.
  dfRulesFinal_ga <- generate_final_rules_table_common(all_solutions_details_df_sorted_ga, params)
  
  # 8. Plot Protein Expression in Subgroups
  # Visualizes the protein expression levels within the identified subgroups
  # compared to healthy volunteers, based on the final selected rules.
  if (!isTRUE(params$benchmark_mode)) {
    plot_protein_expression_in_subgroups_common(dfRulesFinal_ga, patient_groups_df_sorted_ga, df, params)
  }
  
  
  if (!isTRUE(params$benchmark_mode)) {
    message("\n--- Patient Subgroup Discovery (GA) process complete ---")
  }
  
  # Return all generated data frames and analysis results for further programmatic access
  return(list(
    all_solutions_details = all_solutions_details_df_sorted_ga,
    patient_group_details = patient_groups_df_sorted_ga,
    analysis_results = analysis_results_ga,
    optimal_rules_table = dfRulesFinal_ga
  ))
}
