# --- bo_method_no_orbits.R: Bayesian Optimization (BO) Method for Patient Subgroup Discovery (NO EQUIVALENCE CLASS) ---

# Load BO-specific libraries
suppressMessages(library(dbscan)) # Still needed for consistency, but its use for clustering is removed in this version
suppressMessages(library(stringdist))
suppressMessages(library(DiceKriging))

# --- Global Utility Functions (specific to BO, re-used or slightly adapted) ---

# Custom Hamming Distance Kernel for DiceKriging (re-used)
# Calculates the Hamming distance between two binary vectors.
calculate_hamming_distance <- function(x, y) {
  if (length(x) != length(y)) {
    stop("Vectors must be of the same length.")
  }
  sum(x != y)
}

# Covariance function for DiceKriging based on Hamming distance. (re-used)
covHamming <- function(x1, x2, pars) {
  magnitude <- pars["magnitude"]
  length_scale <- pars["length_scale"]
  
  # Handle potential invalid parameter values
  if (is.null(length_scale) || is.infinite(length_scale) || is.na(length_scale) || length_scale <= 0) {
    length_scale <- 1e-2
  }
  if (is.null(magnitude) || is.infinite(magnitude) || is.na(magnitude) || magnitude <= 0) {
    magnitude <- 1.0
  }
  
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  K <- matrix(0, nrow = n1, ncol = n2)
  
  # Compute the covariance matrix using the exponential kernel with Hamming distance
  for (i in 1:n1) {
    for (j in 1:n2) {
      ham_dist <- calculate_hamming_distance(x1[i, ], x2[j, ])
      K[i, j] <- magnitude * exp(-length_scale * ham_dist)
    }
  }
  return(K)
}

# Expected Improvement Acquisition Function (re-used)
expected_improvement_r <- function(X_candidate, gp_model, X_sample, Y_sample, xi) {
  X_candidate <- as.matrix(X_candidate)
  
  # Predict mean and standard deviation from the Gaussian Process model
  gp_predict <- predict(gp_model, newdata = X_candidate, type = "SK", checkNames = FALSE)
  
  mu <- gp_predict$mean
  sigma <- gp_predict$sd
  
  # Prevent division by zero or very small sigma values
  sigma[sigma < 1e-10] <- 1e-10
  
  f_max <- max(Y_sample) # Current best observed fitness
  
  # Calculate Z-score
  Z <- (mu - f_max - xi) / sigma
  
  # Calculate Expected Improvement
  ei <- sigma * (dnorm(Z) + Z * pnorm(Z))
  
  # Set EI to 0 where sigma is very small (high confidence, little improvement expected)
  ei[sigma <= 1e-10] <- 0
  
  return(ei)
}

# Optimize Acquisition Function (MODIFIED - No Clustering)
# Searches for the next best point (binary rule vector) to sample by maximizing the acquisition function.
# It now primarily relies on random search in the full binary space (exploration).
optimize_acquisition_function_r <- function(acquisition_func, gp_model, X_sample, Y_sample, n_atomic_rules_bo,
                                            n_random_search = 1000, max_active_rules_acquisition = 3, xi_val) {
  
  best_ei <- -Inf
  best_x <- NULL
  
  X_sample <- as.matrix(X_sample)
  
  # Strategy: Random search in the full binary space (Exploration)
  # This ensures diversity and helps discover new promising regions.
  random_vectors_list <- list()
  for (k in 1:n_random_search) {
    # Randomly choose the number of active rules (1s) for the candidate vector
    num_active_rules_rand <- sample(1:max_active_rules_acquisition, 1)
    rand_vec <- rep(0, n_atomic_rules_bo)
    if (num_active_rules_rand > 0) {
      rand_vec[sample(1:n_atomic_rules_bo, num_active_rules_rand)] <- 1
    }
    random_vectors_list[[k]] <- rand_vec
  }
  random_vectors <- do.call(rbind, random_vectors_list)
  random_vectors <- unique(random_vectors) # Remove duplicates
  
  # Filter out candidates already sampled
  random_vectors_to_check <- matrix(nrow = 0, ncol = n_atomic_rules_bo)
  if (nrow(random_vectors) > 0) {
    for (k in 1:nrow(random_vectors)) {
      current_rand_vec <- random_vectors[k, ]
      if (!any(apply(X_sample, 1, function(row) all(row == current_rand_vec)))) { # Check if already sampled
        random_vectors_to_check <- rbind(random_vectors_to_check, current_rand_vec)
      }
    }
  }
  
  # Evaluate EI for new random candidates and update best_x
  if (nrow(random_vectors_to_check) > 0) {
    ei_values_random <- acquisition_func(random_vectors_to_check, gp_model, X_sample, Y_sample, xi = xi_val)
    max_ei_idx_random <- which.max(ei_values_random)
    if (ei_values_random[max_ei_idx_random] > best_ei) {
      best_ei <- ei_values_random[max_ei_idx_random]
      best_x <- random_vectors_to_check[max_ei_idx_random, ]
    }
  }
  
  # Fallback if no new unique point is found (should ideally not happen with enough random search)
  if (is.null(best_x)) {
    best_x <- sample(0:1, size = n_atomic_rules_bo, replace = TRUE)
    attempts <- 0
    while (any(apply(X_sample, 1, function(row) all(row == best_x))) && attempts < 100) {
      best_x <- sample(0:1, size = n_atomic_rules_bo, replace = TRUE)
      attempts <- attempts + 1
    }
    if (attempts == 100) {
      warning("Could not find a unique unsampled point after 100 attempts. Returning a potentially sampled point.")
    }
  }
  
  return(best_x)
}


# --- Function: Generate Atomic Rules (re-used directly from original BO script) ---
generate_all_atomic_rules <- function(de_data, numeric_metrics, categorical_metrics) {
  atomic_rules_list <- list()
  atomic_rules_map <- list()
  current_binary_idx <- 0
  
  # Generate rules for numeric metrics
  for (metric in numeric_metrics) {
    min_val <- floor(min(de_data[[metric]], na.rm = TRUE))
    max_val <- ceiling(max(de_data[[metric]], na.rm = TRUE))
    
    if (is.infinite(min_val) || is.infinite(max_val) || min_val > max_val) {
      warning(paste("Skipping numeric metric", metric, "due to invalid range."))
      next
    }
    
    for (k in min_val:max_val) {
      pred_gt <- paste0(metric, " > ", k)
      atomic_rules_list[[length(atomic_rules_list) + 1]] <- pred_gt
      atomic_rules_map[[pred_gt]] <- current_binary_idx
      current_binary_idx <- current_binary_idx + 1
      
      pred_le <- paste0(metric, " <= ", k)
      atomic_rules_list[[length(atomic_rules_list) + 1]] <- pred_le
      atomic_rules_map[[pred_le]] <- current_binary_idx
      current_binary_idx <- current_binary_idx + 1
    }
  }
  
  # Generate rules for categorical metrics
  for (metric in categorical_metrics) {
    levels <- unique(de_data[[metric]]) %>% na.omit() %>% as.character() %>% sort()
    if (length(levels) == 0) {
      warning(paste("Skipping categorical metric", metric, "due to no valid levels."))
      next
    }
    for (level in levels) {
      pred_eq <- paste0(metric, " == '", level, "'")
      atomic_rules_list[[length(atomic_rules_list) + 1]] <- pred_eq
      atomic_rules_map[[pred_eq]] <- current_binary_idx
      current_binary_idx <- current_binary_idx + 1
    }
  }
  
  n_atomic_rules_total <- length(atomic_rules_list)
  message(sprintf("\nTotal number of atomic rule predicates (N_ATOMIC_RULES_BO): %d\n", n_atomic_rules_total))
  
  return(list(
    atomic_rules_list = atomic_rules_list,
    atomic_rules_map = atomic_rules_map,
    n_atomic_rules_total = n_atomic_rules_total
  ))
}

# --- Function: Objective Function for BO (evaluate_subgroup_fitness_bo - re-used directly) ---
evaluate_subgroup_fitness_bo <- function(binary_rule_vector, de_data, hv_aval_mean, atomic_rules_list, protein_variable_name, min_subgroup_size) {
  binary_rule_vector <- as.logical(binary_rule_vector)
  active_predicate_strings <- atomic_rules_list[binary_rule_vector]
  
  if (length(active_predicate_strings) == 0) {
    return(1e-9)
  }
  
  filter_expressions <- list()
  for (pred_str in active_predicate_strings) {
    if (grepl(" %in% ", pred_str)) {
      parts <- strsplit(pred_str, " %in% ")[[1]]
      metric_name <- trimws(parts[1])
      levels_str_raw <- trimws(parts[2])
      levels_parsed <- eval(parse(text = levels_str_raw))
      filter_expressions[[length(filter_expressions) + 1]] <- rlang::expr(!!sym(metric_name) %in% !!levels_parsed)
    } else if (grepl(" <= ", pred_str)) {
      parts <- strsplit(pred_str, " <= ")[[1]]
      metric_name <- trimws(parts[1])
      val <- as.numeric(trimws(parts[2]))
      filter_expressions[[length(filter_expressions) + 1]] <- rlang::expr(!!sym(metric_name) <= !!val)
    } else if (grepl(" > ", pred_str)) {
      parts <- strsplit(pred_str, " > ")[[1]]
      metric_name <- trimws(parts[1])
      val <- as.numeric(trimws(parts[2]))
      filter_expressions[[length(filter_expressions) + 1]] <- rlang::expr(!!sym(metric_name) > !!val)
    } else if (grepl(" == ", pred_str)) {
      parts <- strsplit(pred_str, " == ")[[1]]
      metric_name <- trimws(parts[1])
      level_val <- gsub("'", "", trimws(parts[2]))
      filter_expressions[[length(filter_expressions) + 1]] <- rlang::expr(!!sym(metric_name) == !!level_val)
    } else {
      warning(paste("Unrecognized predicate format:", pred_str))
      return(1e-9)
    }
  }
  
  filtered_de_data <- de_data
  for (filter_expr in filter_expressions) {
    filtered_de_data <- filtered_de_data %>% filter(!!filter_expr)
    if (nrow(filtered_de_data) < min_subgroup_size) {
      return(1e-6)
    }
  }
  
  if (nrow(filtered_de_data) < min_subgroup_size) {
    return(1e-6)
  }
  
  subgroup_aval_mean <- mean(filtered_de_data[[protein_variable_name]], na.rm = TRUE)
  
  if (is.nan(subgroup_aval_mean) || subgroup_aval_mean <= 0) {
    return(1e-9)
  }
  
  fitness_val <- subgroup_aval_mean / hv_aval_mean
  
  if (fitness_val < 1.0) {
    return(fitness_val * 0.5)
  }
  
  return(fitness_val)
}

# --- Function: Main Bayesian Optimization Loop Function for a single min_subgroup_size (MODIFIED - No Equivalence Class) ---
run_bayesian_optimization_for_min_size <- function(n_atomic_rules_bo, initial_samples, bo_iterations,
                                                   bo_acquisition_random_search, max_initial_active_rules,
                                                   de_data, hv_aval_mean, atomic_rules_list, protein_variable_name,
                                                   min_subgroup_size_fixed, xi) {
  
  X_observed_r <- matrix(nrow = 0, ncol = n_atomic_rules_bo)
  Y_observed_r <- numeric(0)
  
  best_overall_y <- -Inf
  best_overall_x <- NULL
  
  # Phase 1: Initial Sampling and Evaluation
  for (i in 1:initial_samples) {
    num_active_rules <- sample(1:max_initial_active_rules, 1)
    rule_vector <- rep(0, n_atomic_rules_bo)
    if (num_active_rules > 0) {
      rule_vector[sample(1:n_atomic_rules_bo, num_active_rules)] <- 1
    }
    
    # Ensure uniqueness for initial samples (regenerate if duplicate found)
    if (i > 1) {
      current_rule_str <- paste(rule_vector, collapse = "")
      if (current_rule_str %in% apply(X_observed_r, 1, paste, collapse = "")) {
        attempts <- 0
        while (current_rule_str %in% apply(X_observed_r, 1, paste, collapse = "") && attempts < 100) {
          num_active_rules <- sample(1:max_initial_active_rules, 1)
          rule_vector <- rep(0, n_atomic_rules_bo)
          if (num_active_rules > 0) {
            rule_vector[sample(1:n_atomic_rules_bo, num_active_rules)] <- 1
          }
          current_rule_str <- paste(rule_vector, collapse = "")
          attempts <- attempts + 1
        }
        if (attempts == 100) {
          warning(paste("Could not find a unique initial sample after 100 attempts for min_subgroup_size:", min_subgroup_size_fixed))
        }
      }
    }
    
    response <- evaluate_subgroup_fitness_bo(rule_vector, de_data, hv_aval_mean, atomic_rules_list, protein_variable_name, min_subgroup_size_fixed)
    # Add the observed point and its response to the observed data
    X_observed_r <- rbind(X_observed_r, rule_vector)
    Y_observed_r <- c(Y_observed_r, response)
  }
  
  # Handle cases where no initial samples could be generated or evaluated
  if (nrow(X_observed_r) == 0 || length(Y_observed_r) == 0) {
    return(list(best_chromosome = rep(0, n_atomic_rules_bo), optimal_fitness = 0,
                X_observed = X_observed_r, Y_observed = Y_observed_r))
  }
  
  # Remove any remaining duplicates after initial sampling (edge case)
  unique_rows_idx <- !duplicated(X_observed_r)
  X_observed_r <- X_observed_r[unique_rows_idx, , drop = FALSE]
  Y_observed_r <- Y_observed_r[unique_rows_idx]
  
  # Initialize overall best fitness and chromosome
  if (length(Y_observed_r) > 0 && max(Y_observed_r) > -Inf) {
    best_overall_y <- max(Y_observed_r)
    best_overall_x <- X_observed_r[which.max(Y_observed_r), , drop = FALSE]
  } else {
    best_overall_y <- 0
    best_overall_x <- rep(0, n_atomic_rules_bo)
  }
  
  # Bayesian Optimization Loop
  for (i in 1:bo_iterations) {
    # No clustering for equivalence classes in this version
    
    # Train Gaussian Process Surrogate Model
    cov_type_ham <- "custom"
    environment(covHamming) <- environment() # Ensure covHamming can access its environment
    
    km_model <- tryCatch({
      DiceKriging::km(
        formula = ~1, # Constant mean model
        design = X_observed_r, # Observed input points (rule vectors)
        response = Y_observed_r, # Observed output values (fitness)
        covtype = cov_type_ham, # Use the custom Hamming distance kernel
        covpars = c(magnitude = 1.0, length_scale = 0.5), # Initial covariance parameters
        control = list(trace = FALSE, pop.size = 20, max.iter = 10), # Optimization control for hyperparameters
        nugget.estim = TRUE # Estimate nugget effect (noise variance)
      )
    }, error = function(e) {
      warning(paste("DiceKriging km() failed:", e$message, ". Skipping acquisition optimization."))
      return(NULL) # Return NULL if model training fails
    })
    
    # Optimize Acquisition Function to find the next point to sample (No cluster_representatives passed)
    if (is.null(km_model)) {
      # Fallback: if GP model fails, pick a random unsampled point
      next_x_r <- sample(0:1, size = n_atomic_rules_bo, replace = TRUE)
      attempts <- 0
      while (any(apply(X_observed_r, 1, function(row) all(row == next_x_r))) && attempts < 100) {
        next_x_r <- sample(0:1, size = n_atomic_rules_bo, replace = TRUE)
        attempts <- attempts + 1
      }
    } else {
      # Use the acquisition function to find the most promising next point
      next_x_r <- optimize_acquisition_function_r(expected_improvement_r, km_model, X_observed_r, Y_observed_r,
                                                  n_atomic_rules_bo,
                                                  bo_acquisition_random_search, max_initial_active_rules, xi_val = xi)
    }
    
    # Evaluate the new point and update observed data
    new_y_r <- evaluate_subgroup_fitness_bo(next_x_r, de_data, hv_aval_mean, atomic_rules_list, protein_variable_name, min_subgroup_size_fixed)
    
    X_observed_r <- rbind(X_observed_r, next_x_r)
    Y_observed_r <- c(Y_observed_r, new_y_r)
    
    # Update overall best solution
    if (new_y_r > best_overall_y) {
      best_overall_y <- new_y_r
      best_overall_x <- next_x_r
    }
  }
  
  # Return the results for this specific min_subgroup_size
  return(list(
    min_subgroup_size = min_subgroup_size_fixed,
    best_chromosome = best_overall_x,
    optimal_fitness = best_overall_y,
    X_observed = X_observed_r,
    Y_observed = Y_observed_r
  ))
}

# --- Function: Orchestrate Parallel BO Runs (MODIFIED - No Equivalence Class parameters) ---
orchestrate_parallel_bo <- function(min_subgroup_sizes_to_test, params, N_ATOMIC_RULES_BO,
                                    de_data, hv_aval_mean, atomic_rules_list, cl_outer) {
  
  if (!isTRUE(params$benchmark_mode)) {
    message("\n--- Running Bayesian Optimization for different min_subgroup_size values (in parallel - No Equivalence Class) ---\n")
  }
  
  # Register the provided cluster for parallel execution
  registerDoSNOW(cl_outer)
  
  # Export all necessary objects and functions to the workers in the cluster.
  # This ensures each worker has access to the data and functions needed for its task.
  clusterExport(cl_outer, c(
    "run_bayesian_optimization_for_min_size", "evaluate_subgroup_fitness_bo",
    "calculate_hamming_distance", "covHamming", "expected_improvement_r",
    "optimize_acquisition_function_r",
    "de_data", "hv_aval_mean", "atomic_rules_list",
    "N_ATOMIC_RULES_BO", "params"
  ), envir = environment()) # Export from the current environment
  
  # Load required packages within each worker's R session
  clusterEvalQ(cl_outer, {
    library(dplyr)
    library(rlang)
    library(dbscan) # Still load for consistency, even if not used for clustering
    library(DiceKriging)
    library(stringdist)
    library(broom)
    # set.seed(Sys.getpid()) # For distinct random streams per worker, if strict reproducibility not needed across parallel runs
  })
  
  # Setup a text progress bar for monitoring the parallel execution
  total_tasks <- length(min_subgroup_sizes_to_test)
  pb <- utils::txtProgressBar(min = 0, max = total_tasks, style = 3)
  progress_update_fun <- function(n) {
    utils::setTxtProgressBar(pb, n)
  }
  
  # Execute the BO runs in parallel using foreach and %dopar%
  all_results_bo_raw <- foreach(
    current_min_size = min_subgroup_sizes_to_test,
    .packages = c("dplyr", "rlang", "dbscan", "DiceKriging", "stringdist", "broom"), # Packages needed by each worker
    .verbose = FALSE, # Suppress verbose output from foreach
    .errorhandling = 'pass', # If an error occurs in a worker, pass the error object instead of stopping
    .options.snow = list(progress = progress_update_fun) # Pass progress update function to doSNOW
  ) %dopar% {
    # Call the encapsulated BO function for the current min_subgroup_size
    result <- run_bayesian_optimization_for_min_size(
      n_atomic_rules_bo = N_ATOMIC_RULES_BO,
      initial_samples = params$bo_initial_samples,
      bo_iterations = params$bo_iterations,
      bo_acquisition_random_search = params$bo_acquisition_random_search,
      max_initial_active_rules = params$max_initial_active_rules,
      de_data = de_data,
      hv_aval_mean = hv_aval_mean,
      atomic_rules_list = atomic_rules_list,
      protein_variable_name = params$variable_for_use,
      min_subgroup_size_fixed = current_min_size,
      xi = params$bo_xi
    )
    return(result)
  }
  close(pb) # Close the progress bar
  if (!isTRUE(params$benchmark_mode)) {
    message("\nAll Bayesian Optimization runs complete.")
  }
  
  # Post-processing: Order and name results for consistency and easier access
  if (is.list(all_results_bo_raw) && length(all_results_bo_raw) > 0) {
    # Check if results are valid and contain 'min_subgroup_size'
    if (all(sapply(all_results_bo_raw, is.list)) && all(sapply(all_results_bo_raw, function(x) !is.null(x$min_subgroup_size)))) {
      result_names <- sapply(all_results_bo_raw, function(x) as.character(x$min_subgroup_size))
      ordered_results <- all_results_bo_raw[order(as.numeric(result_names))]
      names(ordered_results) <- sort(as.numeric(result_names))
      all_results_bo_raw <- ordered_results
    } else {
      warning("Results structure from parallel processing is unexpected. Cannot assign names based on min_subgroup_size.")
    }
  } else {
    warning("all_results_bo_raw is not a list or is empty after parallel processing.")
  }
  
  return(all_results_bo_raw)
}

# --- Function: Process BO Results into Detailed DataFrames (re-used directly) ---
process_bo_results <- function(all_results_bo_raw, de_data, hv_data, hv_aval_mean, atomic_rules_list, params) {
  
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
    message("\n--- Decoding and Analyzing Optimal Rules ---")
  }
  valid_results_bo <- lapply(all_results_bo_raw, function(x) x)
  
  if (length(valid_results_bo) > 0) {
    for (i in seq_along(valid_results_bo)) {
      result_entry <- valid_results_bo[[i]]
      current_min_subgroup_size <- result_entry$min_subgroup_size
      current_optimal_fitness <- result_entry$optimal_fitness
      current_best_chromosome_bo <- result_entry$best_chromosome
      
      rules_table_for_current_solution <- data.frame(
        Metric = character(),
        Operator = character(),
        Value_or_Levels = character(),
        stringsAsFactors = FALSE
      )
      final_filter_expressions_apply_i <- list()
      
      if (is.null(current_best_chromosome_bo)) {
        warning(paste("No best chromosome found for min subgroup size:", current_min_subgroup_size))
        next
      }
      
      for (j in 1:length(atomic_rules_list)) {
        if (current_best_chromosome_bo[j] == 1) {
          pred_str <- atomic_rules_list[[j]]
          
          if (grepl(" %in% ", pred_str)) {
            parts <- strsplit(pred_str, " %in% ")[[1]]
            metric_name <- trimws(parts[1])
            levels_str_raw <- trimws(parts[2])
            levels_parsed <- eval(parse(text = levels_str_raw))
            final_filter_expressions_apply_i[[length(final_filter_expressions_apply_i) + 1]] <- rlang::expr(!!sym(metric_name) %in% !!levels_parsed)
            
            value_for_table <- paste0(levels_parsed, collapse = ", ")
            operator_for_table <- "%in%"
            
          } else if (grepl(" <= ", pred_str)) {
            parts <- strsplit(pred_str, " <= ")[[1]]
            metric_name <- trimws(parts[1])
            val <- as.numeric(trimws(parts[2]))
            final_filter_expressions_apply_i[[length(final_filter_expressions_apply_i) + 1]] <- rlang::expr(!!sym(metric_name) <= !!val)
            
            value_for_table <- as.character(val)
            operator_for_table <- "<="
            
          } else if (grepl(" > ", pred_str)) {
            parts <- strsplit(pred_str, " > ")[[1]]
            metric_name <- trimws(parts[1])
            val <- as.numeric(trimws(parts[2]))
            final_filter_expressions_apply_i[[length(final_filter_expressions_apply_i) + 1]] <- rlang::expr(!!sym(metric_name) > !!val)
            
            value_for_table <- as.character(val)
            operator_for_table <- ">"
            
          } else if (grepl(" == ", pred_str)) {
            parts <- strsplit(pred_str, " == ")[[1]]
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


# --- Main Orchestrating Function for BO (re-used directly) ---
run_bayesian_optimization_analysis <- function(params, df, dfThisProtein, hv_data, de_data, hv_aval_mean, cl_outer = NULL) {
  
  atomic_rules_info <- generate_all_atomic_rules(de_data, params$numeric_metrics, params$categorical_metrics)
  atomic_rules_list <- atomic_rules_info$atomic_rules_list
  atomic_rules_map <- atomic_rules_info$atomic_rules_map
  N_ATOMIC_RULES_BO <- atomic_rules_info$n_atomic_rules_total
  
  all_results_bo_raw <- list()
  
  if (params$runSim) {
    all_results_bo_raw <- orchestrate_parallel_bo(
      min_subgroup_sizes_to_test = params$min_subgroup_sizes_to_test,
      params = params,
      N_ATOMIC_RULES_BO = N_ATOMIC_RULES_BO,
      de_data = de_data,
      hv_aval_mean = hv_aval_mean,
      atomic_rules_list = atomic_rules_list,
      cl_outer = cl_outer
    )
    if (!isTRUE(params$benchmark_mode)) {
      message("\n--- Simulation results (not saved as requested) ---")
    }
    
  } else {
    if (!isTRUE(params$benchmark_mode)) {
      message("\n--- Loading pre-computed simulation results ---")
    }
    tryCatch({
      all_results_bo_raw <- aws.s3::s3read_using(FUN = readRDS,
                                                 object = paste0('data/', 'OptSubgroup_Rules_BO_vs_HV_', params$thisProtein, '_', params$variable_for_use, '.rds'),
                                                 bucket = params$.arv_save$collection)
    }, error = function(e) {
      warning(paste("Error loading pre-computed data from S3:", e$message, ". Generating dummy results for demonstration."))
      if (!exists("N_ATOMIC_RULES_BO") || is.null(N_ATOMIC_RULES_BO)) {
        N_ATOMIC_RULES_BO_dummy <- 10
      } else {
        N_ATOMIC_RULES_BO_dummy <- N_ATOMIC_RULES_BO
      }
      
      all_results_bo_raw <- list(
        "10" = list(min_subgroup_size = 10, optimal_fitness = 1.8, best_chromosome = sample(0:1, N_ATOMIC_RULES_BO_dummy, replace=TRUE)),
        "11" = list(min_subgroup_size = 11, optimal_fitness = 1.9, best_chromosome = sample(0:1, N_ATOMIC_RULES_BO_dummy, replace=TRUE)),
        "12" = list(min_subgroup_size = 12, optimal_fitness = 2.0, best_chromosome = sample(0:1, N_ATOMIC_RULES_BO_dummy, replace=TRUE))
      )
      names(all_results_bo_raw) <- sapply(all_results_bo_raw, function(x) as.character(x$min_subgroup_size))
    })
  }
  
  if (length(all_results_bo_raw) == 0) {
    stop("No Bayesian Optimization results available for further processing. Check data loading or simulation parameters.")
  }
  
  processed_dfs <- process_bo_results(all_results_bo_raw, de_data, hv_data, hv_aval_mean, atomic_rules_list, params)
  all_solutions_details_df_sorted <- processed_dfs$all_solutions_details_df_sorted
  patient_groups_df_sorted <- processed_dfs$patient_groups_df_sorted
  
  analysis_results <- generate_analysis_results_common(all_solutions_details_df_sorted, params)
  
  if (!isTRUE(params$benchmark_mode)) {
    plot_fitness_vs_subgroup_size_common(all_results_bo_raw, hv_aval_mean, params)
    
    generate_numeric_cutoff_plots_common(analysis_results, params)
  }
  
  dfRulesFinal <- generate_final_rules_table_common(all_solutions_details_df_sorted, params)
  
  if (!isTRUE(params$benchmark_mode)) {
    plot_protein_expression_in_subgroups_common(dfRulesFinal, patient_groups_df_sorted, df, params)
  }
  
  if (!isTRUE(params$benchmark_mode)) {
    message("\n--- Patient Subgroup Discovery (BO) process complete ---")
  }
  
  return(list(
    all_solutions_details = all_solutions_details_df_sorted,
    patient_group_details = patient_groups_df_sorted,
    analysis_results = analysis_results,
    optimal_rules_table = dfRulesFinal
  ))
}
