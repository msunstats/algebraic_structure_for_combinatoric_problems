# Exhaustive Search for Discrete Optimization Cases
# This module provides functions to find the true global optimum by exhausting all possible combinations
# when the search space is discrete (categorical metrics only, no numeric metrics).

library(dplyr)
library(rlang)
library(here)

# --- Function: Calculate Search Space Size ---
# Estimates the total number of combinations for exhaustive search feasibility assessment
calculate_search_space_size <- function(categorical_metrics, gene_map) {
  total_combinations <- 1
  
  cat("Debug: Gene map structure:\n")
  for (metric in categorical_metrics) {
    map_info <- gene_map[[metric]]
    if (is.null(map_info)) {
      cat(sprintf("  %s: NULL (not found in gene map)\n", metric))
      next
    }
    
    cat(sprintf("  %s: use_idx=%s, levels_start_idx=%s, levels=%s\n", 
               metric, 
               ifelse(is.null(map_info$use_idx), "NULL", map_info$use_idx),
               ifelse(is.null(map_info$levels_start_idx), "NULL", map_info$levels_start_idx),
               ifelse(is.null(map_info$levels), "NULL", paste(map_info$levels, collapse=","))))
    
    # For each categorical metric:
    # - 1 choice to not use it (use_metric = 0)
    # - 2^K - 1 choices to use it with different level combinations (excluding empty selection)
    # where K is the number of levels
    K <- length(map_info$levels)
    metric_combinations <- 1 + (2^K - 1)  # 1 for not using + (2^K - 1) for using with various level combinations
    
    total_combinations <- total_combinations * metric_combinations
    
    cat(sprintf("Metric %s: %d levels, %d combinations\n", metric, K, metric_combinations))
  }
  
  cat(sprintf("Total search space size: %s combinations\n", format(total_combinations, big.mark = ",")))
  return(total_combinations)
}

# --- Function: Generate All Valid Chromosomes ---
# Creates all possible valid chromosome combinations for categorical metrics only
generate_all_categorical_chromosomes <- function(categorical_metrics, gene_map) {
  if (length(categorical_metrics) == 0) {
    warning("No categorical metrics provided")
    return(matrix(nrow = 0, ncol = 0))
  }
  
  # Calculate total chromosome length
  total_genes <- 0
  valid_metrics <- c()
  
  cat("Debug: Analyzing gene map for chromosome generation:\n")
  for (metric in categorical_metrics) {
    map_info <- gene_map[[metric]]
    if (!is.null(map_info) && !is.null(map_info$levels_start_idx) && !is.null(map_info$levels)) {
      total_genes <- max(total_genes, map_info$levels_start_idx + length(map_info$levels) - 1)
      valid_metrics <- c(valid_metrics, metric)
      cat(sprintf("  %s: use_idx=%d, levels_start_idx=%d, num_levels=%d\n", 
                 metric, map_info$use_idx, map_info$levels_start_idx, length(map_info$levels)))
    } else {
      cat(sprintf("  %s: INVALID - skipping\n", metric))
    }
  }
  
  if (length(valid_metrics) == 0) {
    warning("No valid metrics found in gene map")
    return(matrix(nrow = 0, ncol = 0))
  }
  
  cat(sprintf("Valid metrics: %s\n", paste(valid_metrics, collapse = ", ")))
  cat(sprintf("Total chromosome length: %d genes\n", total_genes))
  
  # Generate all combinations
  all_chromosomes <- list()
  chromosome_count <- 0
  
  # Recursive function to generate combinations
  generate_combinations <- function(metric_idx, current_chromosome) {
    if (metric_idx > length(valid_metrics)) {
      # Base case: all metrics processed, add chromosome to list
      chromosome_count <<- chromosome_count + 1
      all_chromosomes[[chromosome_count]] <<- current_chromosome
      return()
    }
    
    metric <- valid_metrics[metric_idx]
    map_info <- gene_map[[metric]]
    
    # Option 1: Don't use this metric (use_metric = 0)
    temp_chromosome <- current_chromosome
    temp_chromosome[map_info$use_idx] <- 0
    # Set all level genes to 0 when not using the metric
    for (i in 1:length(map_info$levels)) {
      temp_chromosome[map_info$levels_start_idx + i - 1] <- 0
    }
    generate_combinations(metric_idx + 1, temp_chromosome)
    
    # Option 2: Use this metric (use_metric = 1) with various level combinations
    K <- length(map_info$levels)
    # Generate all possible level combinations (excluding empty set)
    for (combination_id in 1:(2^K - 1)) {
      temp_chromosome <- current_chromosome
      temp_chromosome[map_info$use_idx] <- 1
      
      # Decode which levels are selected in this combination
      for (i in 1:K) {
        bit_value <- bitwAnd(bitwShiftR(combination_id, i - 1), 1)
        temp_chromosome[map_info$levels_start_idx + i - 1] <- bit_value
      }
      
      generate_combinations(metric_idx + 1, temp_chromosome)
    }
  }
  
  # Initialize chromosome with zeros
  initial_chromosome <- rep(0, total_genes)
  generate_combinations(1, initial_chromosome)
  
  # Convert list to matrix
  if (length(all_chromosomes) == 0) {
    return(matrix(nrow = 0, ncol = total_genes))
  }
  
  chromosome_matrix <- matrix(unlist(all_chromosomes), nrow = length(all_chromosomes), byrow = TRUE)
  cat(sprintf("Generated %d valid chromosomes with %d genes each\n", nrow(chromosome_matrix), ncol(chromosome_matrix)))
  
  return(chromosome_matrix)
}

# --- Function: Exhaustive Search for True Global Optimum ---
# Finds the true global optimum by evaluating all possible combinations
exhaustive_search_discrete <- function(current_min_subgroup_size, de_data, hv_aval_mean, 
                                     categorical_metrics, gene_map, protein_variable_name,
                                     max_combinations = 1e6) {
  
  cat("Starting exhaustive search for discrete optimization...\n")
  
  # Check if numeric metrics are present (should be NULL for discrete case)
  if (!is.null(categorical_metrics) && length(categorical_metrics) == 0) {
    stop("No categorical metrics provided for discrete exhaustive search")
  }
  
  # Calculate search space size
  search_space_size <- calculate_search_space_size(categorical_metrics, gene_map)
  
  if (search_space_size > max_combinations) {
    warning(sprintf("Search space size (%s) exceeds maximum allowed (%s). Consider sampling or reducing metrics.", 
                   format(search_space_size, big.mark = ","), 
                   format(max_combinations, big.mark = ",")))
    return(NULL)
  }
  
  # Generate all possible chromosomes
  cat("Generating all possible chromosomes...\n")
  all_chromosomes <- generate_all_categorical_chromosomes(categorical_metrics, gene_map)
  
  if (nrow(all_chromosomes) == 0) {
    stop("No valid chromosomes generated")
  }
  
  # Create fitness function (reuse from GA)
  source(here::here("ga_method.R"))
  fitness_function <- generate_ga_fitness_function(current_min_subgroup_size, de_data, hv_aval_mean, 
                                                  NULL, # numeric_metrics = NULL for discrete case
                                                  categorical_metrics, gene_map, protein_variable_name)
  
  # Evaluate all chromosomes
  cat(sprintf("Evaluating fitness for %d combinations...\n", nrow(all_chromosomes)))
  
  start_time <- Sys.time()
  
  # Initialize results storage
  fitness_values <- numeric(nrow(all_chromosomes))
  
  # Evaluate fitness for each chromosome
  for (i in 1:nrow(all_chromosomes)) {
    if (i %% 1000 == 0) {
      cat(sprintf("Progress: %d/%d (%.1f%%)\n", i, nrow(all_chromosomes), 100*i/nrow(all_chromosomes)))
    }
    
    fitness_values[i] <- fitness_function(all_chromosomes[i, ])
  }
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Find the global optimum
  best_idx <- which.max(fitness_values)
  best_fitness <- fitness_values[best_idx]
  best_chromosome <- all_chromosomes[best_idx, ]
  
  cat(sprintf("\nExhaustive search completed in %.2f seconds\n", total_time))
  cat(sprintf("True global optimum fitness: %.6f\n", best_fitness))
  cat(sprintf("Best chromosome: %s\n", paste(best_chromosome, collapse = " ")))
  
  # Decode the best solution for interpretation
  decoded_solution <- decode_categorical_chromosome(best_chromosome, categorical_metrics, gene_map)
  
  cat("\nBest solution interpretation:\n")
  for (rule in decoded_solution) {
    cat(sprintf("  %s\n", rule))
  }
  
  # Return comprehensive results
  results <- list(
    best_fitness = best_fitness,
    best_chromosome = best_chromosome,
    decoded_solution = decoded_solution,
    all_fitness_values = fitness_values,
    all_chromosomes = all_chromosomes,
    search_space_size = search_space_size,
    evaluation_time = total_time,
    n_evaluations = nrow(all_chromosomes)
  )
  
  return(results)
}

# --- Function: Decode Categorical Chromosome ---
# Converts a chromosome back to human-readable rules
decode_categorical_chromosome <- function(chromosome, categorical_metrics, gene_map) {
  rules <- c()
  
  for (metric in categorical_metrics) {
    map_info <- gene_map[[metric]]
    if (is.null(map_info) || is.null(map_info$use_idx) || is.null(map_info$levels_start_idx) || is.null(map_info$levels)) {
      next
    }
    
    # Check bounds
    if (map_info$use_idx > length(chromosome)) {
      warning(sprintf("use_idx %d out of bounds for chromosome length %d", map_info$use_idx, length(chromosome)))
      next
    }
    
    use_metric <- round(chromosome[map_info$use_idx])
    
    if (use_metric == 1) {
      selected_levels <- c()
      for (i in 1:length(map_info$levels)) {
        level_idx <- map_info$levels_start_idx + i - 1
        if (level_idx > length(chromosome)) {
          warning(sprintf("level_idx %d out of bounds for chromosome length %d", level_idx, length(chromosome)))
          next
        }
        
        level_gene_val <- round(chromosome[level_idx])
        if (level_gene_val == 1) {
          selected_levels <- c(selected_levels, map_info$levels[i])
        }
      }
      
      if (length(selected_levels) > 0) {
        rule_text <- sprintf("%s %s {%s}", metric, 
                           ifelse(length(selected_levels) == 1, "==", "in"), 
                           paste(selected_levels, collapse = ", "))
        rules <- c(rules, rule_text)
      }
    }
  }
  
  if (length(rules) == 0) {
    rules <- "No rules (empty solution)"
  }
  
  return(rules)
}

# --- Function: Compare Methods to True Optimum ---
# Compares optimization method results against the true global optimum
compare_to_true_optimum <- function(method_results, true_optimum_fitness) {
  comparison_results <- list()
  
  for (method_name in names(method_results)) {
    method_data <- method_results[[method_name]]
    
    if (is.null(method_data) || !("best_fitness" %in% names(method_data))) {
      next
    }
    
    best_method_fitness <- method_data$best_fitness
    optimality_gap <- (true_optimum_fitness - best_method_fitness) / true_optimum_fitness * 100
    performance_ratio <- best_method_fitness / true_optimum_fitness
    
    comparison_results[[method_name]] <- list(
      method_fitness = best_method_fitness,
      true_optimum_fitness = true_optimum_fitness,
      optimality_gap_percent = optimality_gap,
      performance_ratio = performance_ratio,
      is_global_optimum = abs(best_method_fitness - true_optimum_fitness) < 1e-10
    )
    
    cat(sprintf("%s: %.6f (%.2f%% of optimum, gap: %.2f%%)\n", 
               method_name, best_method_fitness, performance_ratio * 100, optimality_gap))
  }
  
  return(comparison_results)
}

# Export functions
if (exists("ga_method.R", envir = .GlobalEnv)) {
  cat("Exhaustive search functions loaded successfully\n")
}
