# --- greedy_method.R: Greedy Algorithm for Patient Subgroup Discovery ---

# This script implements a greedy algorithm to find optimal combinations of
# clinical tests to identify patient subgroups.
# It reuses atomic rule generation and fitness evaluation functions from other scripts.

# No specific libraries to load here, as common ones are loaded in common_fct.R
# and bo_method.R which contains generate_all_atomic_rules.

# --- Function: Run Greedy Algorithm for a single min_subgroup_size ---
#' @param current_min_size The minimum number of patients required in the subgroup for this greedy run.
#' @param greedy_params A list of greedy-specific parameters (e.g., max_rules_per_subgroup).
#' @param de_data The dataframe for disease patients.
#' @param hv_aval_mean The mean protein level for healthy volunteers.
#' @param atomic_rules_list A list of all possible atomic rule predicates.
#' @param protein_variable_name The name of the protein variable to use for fitness calculation.
#' @return A list containing the min_subgroup_size, optimal_fitness, and the best_chromosome found.
run_greedy_algorithm_for_min_size <- function(current_min_size, greedy_params, de_data, hv_aval_mean,
                                              atomic_rules_list, protein_variable_name) {
  
  n_atomic_rules_total <- length(atomic_rules_list) # Total number of atomic rules
  
  best_chromosome_overall <- rep(0, n_atomic_rules_total) # Start with no active rules
  optimal_fitness_overall <- evaluate_subgroup_fitness_bo(best_chromosome_overall, de_data, hv_aval_mean, atomic_rules_list, protein_variable_name, current_min_size)
  
  # Ensure initial fitness is not too low if starting with no rules
  if (optimal_fitness_overall == 1e-9) { # If no rules yield 1e-9, set to a baseline
    optimal_fitness_overall <- 0
  }
  
  current_active_rules_indices <- which(best_chromosome_overall == 1)
  
  # Greedy search loop
  # In each iteration, try adding one new rule that maximizes fitness.
  # Stop if no improvement is found or max_rules_per_subgroup is reached.
  for (iter in 1:greedy_params$max_rules_per_subgroup) {
    best_improvement <- -Inf
    best_candidate_rule_idx <- NULL
    
    # Iterate through all atomic rules that are not currently active
    for (rule_idx in 1:n_atomic_rules_total) {
      if (!(rule_idx %in% current_active_rules_indices)) { # Only consider adding non-active rules
        candidate_chromosome <- best_chromosome_overall
        candidate_chromosome[rule_idx] <- 1 # Tentatively add this rule
        
        # Evaluate fitness of the candidate chromosome
        current_fitness <- evaluate_subgroup_fitness_bo(candidate_chromosome, de_data, hv_aval_mean, atomic_rules_list, protein_variable_name, current_min_size)
        
        # Check if this addition improves the fitness
        if (current_fitness > optimal_fitness_overall) {
          if (current_fitness - optimal_fitness_overall > best_improvement) {
            best_improvement <- current_fitness - optimal_fitness_overall
            best_candidate_rule_idx <- rule_idx
          }
        }
      }
    }
    
    # If a beneficial rule was found, add it and update the best solution
    if (!is.null(best_candidate_rule_idx) && best_improvement > 0) {
      best_chromosome_overall[best_candidate_rule_idx] <- 1
      optimal_fitness_overall <- optimal_fitness_overall + best_improvement
      current_active_rules_indices <- which(best_chromosome_overall == 1)
    } else {
      # No further improvement possible by adding a single rule
      break
    }
  }
  
  # Return the results of this greedy run
  list(
    min_subgroup_size = current_min_size,
    optimal_fitness = optimal_fitness_overall,
    best_chromosome = best_chromosome_overall
  )
}

# --- Function: Orchestrate Parallel Greedy Runs ---
#' @param min_subgroup_sizes_to_test A vector of minimum subgroup sizes to test.
#' @param greedy_params A list of greedy-specific parameters.
#' @param de_data The dataframe for disease patients.
#' @param hv_aval_mean The mean protein level for healthy volunteers.
#' @param atomic_rules_list A list of all possible atomic rule predicates.
#' @param protein_variable_name The name of the protein variable to use for fitness calculation.
#' @param cl_outer A pre-existing parallel cluster object (from makeCluster).
#' @return A list of results from each greedy run, indexed by min_subgroup_size.
orchestrate_parallel_greedy <- function(min_subgroup_sizes_to_test, greedy_params, de_data, hv_aval_mean,
                                        atomic_rules_list, protein_variable_name, cl_outer) {
  
  if (!isTRUE(greedy_params$benchmark_mode)) {
    message("\n--- Running Greedy Algorithm for different min_subgroup_size values (in parallel) ---\n")
  }
  
  # Register the provided cluster for parallel execution
  registerDoSNOW(cl_outer)
  
  # Export all necessary objects and functions to the workers in the cluster.
  clusterExport(cl_outer, c(
    "run_greedy_algorithm_for_min_size", "evaluate_subgroup_fitness_bo", # Reusing fitness function from BO
    "greedy_params", "de_data", "hv_aval_mean", "atomic_rules_list", "protein_variable_name"
  ), envir = environment())
  
  # Load required packages within each worker's R session
  clusterEvalQ(cl_outer, {
    library(dplyr)
    library(rlang)
    # No other specific libraries needed for greedy in workers
  })
  
  # Setup a progress bar
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
  
  # Execute the greedy runs in parallel
  all_results_greedy_raw <- foreach(
    current_min_size = min_subgroup_sizes_to_test,
    .packages = c("dplyr", "rlang"),
    .verbose = FALSE,
    .errorhandling = 'pass',
    .options.snow = list(progress = progress_update_fun)
  ) %dopar% {
    result <- run_greedy_algorithm_for_min_size(
      current_min_size = current_min_size,
      greedy_params = greedy_params,
      de_data = de_data,
      hv_aval_mean = hv_aval_mean,
      atomic_rules_list = atomic_rules_list,
      protein_variable_name = protein_variable_name
    )
    return(result)
  }
  if (!isTRUE(greedy_params$benchmark_mode)) {
    message("\nAll Greedy Algorithm runs complete.")
  }
  
  # Post-processing: Order and name results
  if (is.list(all_results_greedy_raw) && length(all_results_greedy_raw) > 0) {
    if (all(sapply(all_results_greedy_raw, is.list)) && all(sapply(all_results_greedy_raw, function(x) !is.null(x$min_subgroup_size)))) {
      result_names <- sapply(all_results_greedy_raw, function(x) as.character(x$min_subgroup_size))
      ordered_results <- all_results_greedy_raw[order(as.numeric(result_names))]
      names(ordered_results) <- sort(as.numeric(result_names))
      all_results_greedy_raw <- ordered_results
    } else {
      warning("Results structure from parallel processing is unexpected. Cannot assign names based on min_subgroup_size.")
    }
  } else {
    warning("all_results_greedy_raw is not a list or is empty after parallel processing.")
  }
  
  return(all_results_greedy_raw)
}

# --- Function: Process Greedy Results into Detailed DataFrames ---
#' @param all_results_greedy_raw Raw results list from parallel Greedy runs.
#' @param de_data The dataframe for disease patients.
#' @param hv_data The dataframe for healthy volunteers.
#' @param hv_aval_mean The mean protein level for healthy volunteers.
#' @param atomic_rules_list A list of all possible atomic rule predicates.
#' @param params A list of parameters for the analysis (includes variable_for_use).
#' @return A list containing two sorted dataframes: all_solutions_details_df_sorted and patient_groups_df_sorted.
process_greedy_results <- function(all_results_greedy_raw, de_data, hv_data, hv_aval_mean, atomic_rules_list, params) {
  
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
    Group = character(),
    stringsAsFactors = FALSE
  )
  
  hv_usubjid_list <- unique(hv_data$USUBJID)
  processed_min_subgroup_sizes_for_hv_rest <- c()
  
  if (!isTRUE(params$benchmark_mode)) {
    message("\n--- Decoding and Analyzing Optimal Rules (Greedy) ---")
  }
  valid_results_greedy <- lapply(all_results_greedy_raw, function(x) x)
  
  if (length(valid_results_greedy) > 0) {
    for (i in seq_along(valid_results_greedy)) {
      result_entry <- valid_results_greedy[[i]]
      current_min_subgroup_size <- result_entry$min_subgroup_size
      current_optimal_fitness <- result_entry$optimal_fitness
      current_best_chromosome_greedy <- result_entry$best_chromosome
      
      rules_table_for_current_solution <- data.frame(
        Metric = character(),
        Operator = character(),
        Value_or_Levels = character(),
        stringsAsFactors = FALSE
      )
      final_filter_expressions_apply_i <- list()
      
      if (is.null(current_best_chromosome_greedy)) {
        warning(paste("No best chromosome found for min subgroup size (Greedy):", current_min_subgroup_size))
        next
      }
      
      # Decode the binary chromosome into human-readable rules
      for (j in 1:length(atomic_rules_list)) {
        if (current_best_chromosome_greedy[j] == 1) {
          pred_str <- atomic_rules_list[[j]]
          
          if (grepl(" %in% ", pred_str, fixed = TRUE)) {
            parts <- strsplit(pred_str, " %in% ", fixed = TRUE)[[1]]
            metric_name <- trimws(parts[1])
            levels_str_raw <- trimws(parts[2])
            levels_parsed <- eval(parse(text = levels_str_raw))
            final_filter_expressions_apply_i[[length(final_filter_expressions_apply_i) + 1]] <- rlang::expr(!!sym(metric_name) %in% !!levels_parsed)
            
            value_for_table <- paste0(levels_parsed, collapse = ", ")
            operator_for_table <- "%in%"
            
          } else if (grepl(" <= ", pred_str, fixed = TRUE)) {
            parts <- strsplit(pred_str, " <= ", fixed = TRUE)[[1]]
            metric_name <- trimws(parts[1])
            val <- as.numeric(trimws(parts[2]))
            final_filter_expressions_apply_i[[length(final_filter_expressions_apply_i) + 1]] <- rlang::expr(!!sym(metric_name) <= !!val)
            
            value_for_table <- as.character(val)
            operator_for_table <- "<="
            
          } else if (grepl(" > ", pred_str, fixed = TRUE)) {
            parts <- strsplit(pred_str, " > ", fixed = TRUE)[[1]]
            metric_name <- trimws(parts[1])
            val <- as.numeric(trimws(parts[2]))
            final_filter_expressions_apply_i[[length(final_filter_expressions_apply_i) + 1]] <- rlang::expr(!!sym(metric_name) > !!val)
            
            value_for_table <- as.character(val)
            operator_for_table <- ">"
            
          } else if (grepl(" == ", pred_str, fixed = TRUE)) {
            parts <- strsplit(pred_str, " == ", fixed = TRUE)[[1]]
            metric_name <- trimws(parts[1])
            level_val <- trimws(parts[2])
            final_filter_expressions_apply_i[[length(final_filter_expressions_apply_i) + 1]] <- rlang::expr(!!sym(metric_name) == !!gsub("'", "", level_val))
            
            value_for_table <- level_val
            operator_for_table <- "=="
            
          } else {
            warning(paste("Unrecognized predicate format:", pred_str))
            next
          }
          
          rules_table_for_current_solution <- rbind(rules_table_for_current_solution, data.frame(
            Metric = metric_name,
            Operator = operator_for_table,
            Value_or_Levels = value_for_table,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      subgroup_aval_mean_current <- NA
      diff_from_hv_current <- NA
      percent_diff_from_hv_current <- NA
      num_patients_current <- 0
      pVal <- NA
      
      current_subgroup_usubjids <- character(0)
      current_rest_de_usubjids <- character(0)
      
      if (nrow(rules_table_for_current_solution) > 0) {
        temp_filtered_de_data <- de_data
        is_valid_subgroup_for_ids <- TRUE
        
        for (filter_expr in final_filter_expressions_apply_i) {
          temp_filtered_de_data <- temp_filtered_de_data %>% filter(!!filter_expr)
          if (nrow(temp_filtered_de_data) < current_min_subgroup_size) {
            is_valid_subgroup_for_ids <- FALSE
            break
          }
        }
        
        if (is_valid_subgroup_for_ids && nrow(temp_filtered_de_data) >= current_min_subgroup_size) {
          subgroup_aval_mean_current <- mean(temp_filtered_de_data[[params$variable_for_use]], na.rm = TRUE)
          
          if (!is.nan(subgroup_aval_mean_current) && subgroup_aval_mean_current > 0) {
            diff_from_hv_current <- subgroup_aval_mean_current - hv_aval_mean
            percent_diff_from_hv_current <- (diff_from_hv_current / hv_aval_mean) * 100
            num_patients_current <- nrow(temp_filtered_de_data)
            
            pVal_test <- tryCatch({
              wilcox.test(temp_filtered_de_data[[params$variable_for_use]],
                          hv_data[[params$variable_for_use]],
                          alternative = 'greater',
                          paired = FALSE)
            }, error = function(e) {
              warning(paste("Wilcoxon test failed for min_subgroup_size", current_min_subgroup_size, ":", e$message))
              return(NULL)
            })
            
            if (!is.null(pVal_test)) {
              pVal <- broom::tidy(pVal_test) %>% select(pVal = p.value) %>% pull(pVal)
            }
            
            current_subgroup_usubjids <- unique(temp_filtered_de_data$USUBJID)
            all_de_usubjids <- unique(de_data$USUBJID)
            current_rest_de_usubjids <- all_de_usubjids[!all_de_usubjids %in% current_subgroup_usubjids]
            
          } else {
            rules_table_for_current_solution <- rules_table_for_current_solution[0,]
          }
        } else {
          rules_table_for_current_solution <- rules_table_for_current_solution[0,]
        }
      } else {
        rules_table_for_current_solution <- rules_table_for_current_solution[0,]
      }
      
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
        
        solution_rows_for_main_df <- cbind(common_info[rep(1, nrow(rules_table_for_current_solution)), ],
                                           rules_table_for_current_solution)
        
        all_solutions_details_df <- rbind(all_solutions_details_df, solution_rows_for_main_df)
      }
      
      if(num_patients_current >= current_min_subgroup_size && length(current_subgroup_usubjids) > 0) {
        patient_groups_df <- rbind(patient_groups_df, data.frame(
          MinSubgroupSize = current_min_subgroup_size,
          USUBJID = current_subgroup_usubjids,
          Group = "Subgroup",
          stringsAsFactors = FALSE
        ))
        
        if(length(current_rest_de_usubjids) > 0) {
          patient_groups_df <- rbind(patient_groups_df, data.frame(
            MinSubgroupSize = current_min_subgroup_size,
            USUBJID = current_rest_de_usubjids,
            Group = "Rest_DE",
            stringsAsFactors = FALSE
          ))
        }
        
        if(!(current_min_subgroup_size %in% processed_min_subgroup_sizes_for_hv_rest) && length(hv_usubjid_list) > 0) {
          patient_groups_df <- rbind(patient_groups_df, data.frame(
            MinSubgroupSize = current_min_subgroup_size,
            USUBJID = hv_usubjid_list,
            Group = "HV",
            stringsAsFactors = FALSE
          ))
          processed_min_subgroup_sizes_for_hv_rest <- c(processed_min_subgroup_sizes_for_hv_rest, current_min_subgroup_size)
        }
      }
    }
  }
  
  all_solutions_details_df_sorted <- all_solutions_details_df %>%
    arrange(desc(OptimalFitness), MinSubgroupSize, Metric)
  
  patient_groups_df_sorted <- patient_groups_df %>%
    arrange(MinSubgroupSize, Group, USUBJID)
  
  return(list(
    all_solutions_details_df_sorted = all_solutions_details_df_sorted,
    patient_groups_df_sorted = patient_groups_df_sorted
  ))
}

# --- Main Orchestrating Function for Greedy ---
#' @param params A list of parameters for the analysis.
#' @param df The full clinical dataframe.
#' @param dfThisProtein The dataframe filtered for the specific protein.
#' @param hv_data The dataframe for healthy volunteers.
#' @param de_data The dataframe for disease patients.
#' @param hv_aval_mean The mean protein level for healthy volunteers.
#' @param atomic_rules_list A list of all possible atomic rule predicates (from generate_all_atomic_rules).
#' @param cl_outer An optional pre-existing parallel cluster object.
#' @return A list containing all generated data frames for further programmatic access.
run_greedy_analysis <- function(params, df, dfThisProtein, hv_data, de_data, hv_aval_mean, atomic_rules_list, cl_outer = NULL) {
  
  all_results_greedy_raw <- list()
  
  if (params$runSim) {
    all_results_greedy_raw <- orchestrate_parallel_greedy(
      min_subgroup_sizes_to_test = params$min_subgroup_sizes_to_test,
      greedy_params = params, # Pass all greedy-related parameters
      de_data = de_data,
      hv_aval_mean = hv_aval_mean,
      atomic_rules_list = atomic_rules_list, # Pass atomic_rules_list
      protein_variable_name = params$variable_for_use,
      cl_outer = cl_outer
    )
    if (!isTRUE(params$benchmark_mode)) {
      message("\n--- Simulation results (not saved as requested) ---")
    }
    
  } else {
    if (!isTRUE(params$benchmark_mode)) {
      message("\n--- Loading pre-computed simulation results for Greedy ---")
    }
    tryCatch({
      all_results_greedy_raw <- aws.s3::s3read_using(FUN = readRDS,
                                                     object = paste0('data/', 'OptSubgroup_Rules_Greedy_vs_HV_', params$thisProtein, '_', params$variable_for_use, '.rds'),
                                                     bucket = params$.arv_save$collection)
    }, error = function(e) {
      warning(paste("Error loading pre-computed Greedy data from S3:", e$message, ". Generating dummy results for demonstration."))
      # Dummy results for demonstration purposes
      # Need a dummy atomic_rules_list length if not generated from real data
      n_atomic_rules_total_dummy <- length(atomic_rules_list) # Use the generated atomic_rules_list length
      
      all_results_greedy_raw <- list(
        "10" = list(min_subgroup_size = 10, optimal_fitness = 1.2, best_chromosome = sample(0:1, n_atomic_rules_total_dummy, replace=TRUE)),
        "11" = list(min_subgroup_size = 11, optimal_fitness = 1.3, best_chromosome = sample(0:1, n_atomic_rules_total_dummy, replace=TRUE)),
        "12" = list(min_subgroup_size = 12, optimal_fitness = 1.4, best_chromosome = sample(0:1, n_atomic_rules_total_dummy, replace=TRUE))
      )
      names(all_results_greedy_raw) <- sapply(all_results_greedy_raw, function(x) as.character(x$min_subgroup_size))
    })
  }
  
  if (length(all_results_greedy_raw) == 0) {
    stop("No Greedy Algorithm results available for further processing. Check data loading or simulation parameters.")
  }
  
  # Process Greedy Results into Detailed DataFrames
  processed_dfs_greedy <- process_greedy_results(all_results_greedy_raw, de_data, hv_data, hv_aval_mean, atomic_rules_list, params)
  all_solutions_details_df_sorted_greedy <- processed_dfs_greedy$all_solutions_details_df_sorted
  patient_groups_df_sorted_greedy <- processed_dfs_greedy$patient_groups_df_sorted
  
  # Generate Analysis Results (metrics frequencies, etc.)
  analysis_results_greedy <- generate_analysis_results_common(all_solutions_details_df_sorted_greedy, params)
  
  # Plot Fitness vs. Minimum Subgroup Size
  if (!isTRUE(params$benchmark_mode)) {
    plot_fitness_vs_subgroup_size_common(all_results_greedy_raw, hv_aval_mean, params)
    
    # Generate Numeric Cutoff Plots (if applicable)
    if (params$verbose_analysis_output) {
      generate_numeric_cutoff_plots_common(analysis_results_greedy, params)
    }
  }
  
  # Generate Final Rules Table
  dfRulesFinal_greedy <- generate_final_rules_table_common(all_solutions_details_df_sorted_greedy, params)
  
  # Plot Protein Expression in Subgroups
  if (!isTRUE(params$benchmark_mode)) {
    plot_protein_expression_in_subgroups_common(dfRulesFinal_greedy, patient_groups_df_sorted_greedy, df, params)
  }
  
  if (!isTRUE(params$benchmark_mode)) {
    message("\n--- Patient Subgroup Discovery (Greedy) process complete ---")
  }
  
  return(list(
    all_solutions_details = all_solutions_details_df_sorted_greedy,
    patient_group_details = patient_groups_df_sorted_greedy,
    analysis_results = analysis_results_greedy,
    optimal_rules_table = dfRulesFinal_greedy
  ))
}

