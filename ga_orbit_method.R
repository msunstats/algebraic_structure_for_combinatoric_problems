# --- ga_orbit_method.R: Genetic Algorithm with Orbit Awareness for Patient Subgroup Discovery ---

# This script implements a Genetic Algorithm that incorporates "orbit awareness"
# through niche elite preservation. It aims to maintain diversity in the population
# by ensuring the best solutions from different "orbits" (clusters of similar rules)
# are carried over to the next generation.

# Load necessary libraries
suppressMessages(library(GA)) # Core GA library functions (still needed for some utilities if any, but main loop is custom)
suppressMessages(library(dbscan)) # For clustering to identify orbits
# rlang, dplyr are loaded via common_fct.R and exported to cluster workers

# --- Custom GA Operators (replacing internal GA package functions) ---

# Custom Population Initialization
# Generates a random initial population respecting gene min/max and types.
custom_ga_init <- function(popSize, ga_min, ga_max, ga_type) {
  n_genes <- length(ga_min)
  population <- matrix(nrow = popSize, ncol = n_genes)
  for (i in 1:popSize) {
    for (j in 1:n_genes) {
      if (ga_type[j] == "binary") {
        population[i, j] <- sample(0:1, 1)
      } else if (ga_type[j] == "real") {
        population[i, j] <- runif(1, ga_min[j], ga_max[j])
      } else {
        stop("Unsupported gene type in custom_ga_init")
      }
    }
  }
  return(population)
}

# Custom Roulette Wheel Selection
# Selects individuals for reproduction based on their fitness values.
custom_ga_selection_rw <- function(fitness_values, num_to_select) {
  if (all(fitness_values <= 0 | is.na(fitness_values))) { # Handle all non-positive or NA fitness
    warning("All fitness values are non-positive or NA during selection. Selecting randomly.")
    return(sample(1:length(fitness_values), num_to_select, replace = TRUE))
  }
  
  # Normalize fitness values to probabilities
  # Shift fitness values to be non-negative if there are negative values
  min_fitness <- min(fitness_values, na.rm = TRUE)
  if (min_fitness < 0) {
    shifted_fitness <- fitness_values - min_fitness # Shift so smallest is 0
  } else {
    shifted_fitness <- fitness_values
  }
  
  # Replace NA fitness with 0 for selection purposes
  shifted_fitness[is.na(shifted_fitness)] <- 0
  
  total_fitness <- sum(shifted_fitness)
  if (total_fitness == 0) { # If after shifting, all are zero (e.g., all were same negative value)
    warning("Total fitness is zero after shifting. Selecting randomly.")
    return(sample(1:length(fitness_values), num_to_select, replace = TRUE))
  }
  
  selection_probs <- shifted_fitness / total_fitness
  
  # Sample indices based on probabilities
  selected_indices <- sample(1:length(fitness_values), num_to_select, replace = TRUE, prob = selection_probs)
  return(selected_indices)
}

# Custom Single Point Crossover
# Combines two parent chromosomes to create two offspring.
custom_ga_crossover_sp <- function(parent1, parent2, ga_min, ga_max, ga_type) {
  n_genes <- length(parent1)
  crossover_point <- sample(1:(n_genes - 1), 1) # Choose a random crossover point
  
  child1 <- c(parent1[1:crossover_point], parent2[(crossover_point + 1):n_genes])
  child2 <- c(parent2[1:crossover_point], parent1[(crossover_point + 1):n_genes])
  
  return(rbind(child1, child2))
}

# Custom Mutation
# Randomly alters genes in a chromosome.
custom_ga_mutation <- function(chromosome, ga_min, ga_max, ga_type, pmutation) {
  n_genes <- length(chromosome)
  mutated_chromosome <- chromosome
  for (j in 1:n_genes) {
    if (runif(1) < pmutation) { # Apply mutation with probability pmutation
      if (ga_type[j] == "binary") {
        mutated_chromosome[j] <- 1 - mutated_chromosome[j] # Flip the bit
      } else if (ga_type[j] == "real") {
        # Mutate real value within its min/max range
        # A simple mutation: add a small random perturbation, then clamp
        perturbation <- runif(1, -0.1 * (ga_max[j] - ga_min[j]), 0.1 * (ga_max[j] - ga_min[j]))
        mutated_chromosome[j] <- mutated_chromosome[j] + perturbation
        mutated_chromosome[j] <- max(ga_min[j], min(ga_max[j], mutated_chromosome[j])) # Clamp to bounds
      }
    }
  }
  return(mutated_chromosome)
}


# --- Function: Define Chromosome Encoding and Gene Map for GA (Re-used from ga_method.R) ---
# This function defines how clinical tests (numeric and categorical) are represented
# as genes within the GA chromosome, and creates a mapping to decode them later.
define_ga_encoding_orbit <- function(de_data, numeric_metrics, categorical_metrics) {
  ga_min <- c() # Minimum value for each gene
  ga_max <- c() # Maximum value for each gene
  ga_type <- c() # Type ("binary" or "real") for each gene
  gene_map <- list() # To map gene indices back to metric names/levels for decoding
  
  current_gene_idx <- 1 # Tracks the current gene index for sequential assignment
  
  # Encode Numeric Metrics: Each numeric metric is represented by 3 genes:
  # 1. Binary: 0 (don't use this metric), 1 (use this metric)
  # 2. Binary: 0 (operator is <=), 1 (operator is >)
  # 3. Real: The threshold value for the metric
  if (!is.null(numeric_metrics) && length(numeric_metrics) > 0) { # Only process if numeric_metrics is not NULL and has elements
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
  }
  
  # Encode Categorical Metrics: Each categorical metric is represented by 1 + K genes:
  # 1. Binary: 0 (don't use), 1 (use)
  # 2. K binary genes: 0 (don't select level), 1 (select level) for each unique level
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

# --- Function: Convert GA Chromosome to Atomic Rule Vector ---
# This helper function converts a GA chromosome (mix of binary and real)
# into the binary atomic rule vector representation used by evaluate_subgroup_fitness_bo.
# This is crucial for clustering in the atomic rule space.
convert_chromosome_to_atomic_vector <- function(chromosome, gene_map, atomic_rules_map, atomic_rules_list, n_atomic_rules_total) {
  atomic_vector <- rep(0, n_atomic_rules_total)
  
  # Decode Numeric Rules (this block will be skipped if numeric_metrics is NULL)
  for (metric in names(gene_map)) {
    map_info <- gene_map[[metric]]
    if (map_info$type == "numeric") {
      # Ensure indices are within chromosome bounds
      if (map_info$use_idx > length(chromosome) || map_info$op_idx > length(chromosome) || map_info$val_idx > length(chromosome)) {
        warning(paste("convert_chromosome_to_atomic_vector: Numeric gene indices out of bounds for metric", metric))
        next
      }
      use_metric <- round(chromosome[map_info$use_idx])
      if (use_metric == 1) {
        op_type <- round(chromosome[map_info$op_idx])
        threshold_val <- chromosome[map_info$val_idx]
        
        k_val <- round(threshold_val)
        
        pred_str <- if (op_type == 0) { # <= operator
          paste0(metric, " <= ", k_val)
        } else { # > operator
          paste0(metric, " > ", k_val)
        }
        
        if (pred_str %in% names(atomic_rules_map)) {
          atomic_vector[atomic_rules_map[[pred_str]] + 1] <- 1 # +1 because R indices are 1-based
        }
      }
    } else if (map_info$type == "categorical") {
      # Ensure use_idx is within chromosome bounds
      if (map_info$use_idx > length(chromosome)) {
        warning(paste("convert_chromosome_to_atomic_vector: Categorical use_idx out of bounds for metric", metric))
        next
      }
      use_metric <- round(chromosome[map_info$use_idx])
      if (use_metric == 1) {
        selected_levels <- c()
        for (i in 1:length(map_info$levels)) {
          level_gene_idx <- map_info$levels_start_idx + i - 1
          # Ensure level_gene_idx is within chromosome bounds
          if (level_gene_idx > length(chromosome)) {
            warning(paste("convert_chromosome_to_atomic_vector: Categorical level_gene_idx out of bounds for metric", metric, "level", map_info$levels[i]))
            next
          }
          level_gene_val <- round(chromosome[level_gene_idx])
          if (level_gene_val == 1) {
            selected_levels <- c(selected_levels, map_info$levels[i])
          }
        }
        if (length(selected_levels) > 0) {
          for(level in selected_levels) {
            pred_str <- paste0(metric, " == '", level, "'")
            if (pred_str %in% names(atomic_rules_map)) {
              atomic_vector[atomic_rules_map[[pred_str]] + 1] <- 1
            }
          }
        } else {
          # If a categorical metric is "used" (use_metric == 1) but no levels are selected,
          # this chromosome essentially forms no rule for this metric.
          # We don't add to atomic_vector, and it will result in lower fitness.
        }
      }
    }
  }
  return(atomic_vector)
}


# --- Function: Run Orbit-Aware Genetic Algorithm for a single min_subgroup_size ---
#' @param current_min_size The minimum number of patients required in the subgroup for this GA run.
#' @param ga_params A list of GA-specific parameters (pop_size, max_iter, run_limit, orbit_check_interval, dbscan_eps, dbscan_minpts).
#' @param de_data The dataframe for disease patients.
#' @param hv_aval_mean The mean protein level for healthy volunteers.
#' @param numeric_metrics A character vector of numeric clinical metrics.
#' @param categorical_metrics A character vector of categorical clinical metrics.
#' @param gene_map The mapping from gene indices to metric names/levels.
#' @param ga_min Vector of minimum values for each gene.
#' @param ga_max Vector of maximum values for each gene.
#' @param ga_type Vector of types ("binary", "real") for each gene.
#' @param protein_variable_name The name of the protein variable to use for fitness calculation.
#' @param atomic_rules_map The map from atomic rule strings to their binary indices.
#' @param atomic_rules_list The list of atomic rule strings.
#' @param n_atomic_rules_total Total number of atomic rules.
#' @return A list containing the min_subgroup_size, optimal_fitness, and the best_chromosome found.
run_genetic_algorithm_orbit_for_min_size <- function(current_min_size, ga_params, de_data, hv_aval_mean,
                                                     numeric_metrics, categorical_metrics, gene_map,
                                                     ga_min, ga_max, ga_type, protein_variable_name,
                                                     atomic_rules_map, atomic_rules_list, n_atomic_rules_total) {
  
  # Validate inputs
  tryCatch({
    if (is.null(de_data) || nrow(de_data) == 0) {
      stop(paste("Invalid de_data for min_size", current_min_size))
    }
    if (is.null(atomic_rules_map) || length(atomic_rules_map) == 0) {
      stop(paste("Invalid atomic_rules_map for min_size", current_min_size))
    }
    if (is.null(atomic_rules_list) || length(atomic_rules_list) == 0) {
      stop(paste("Invalid atomic_rules_list for min_size", current_min_size))
    }
    
    # Check if required columns exist in de_data
    missing_numeric <- NULL
    missing_categorical <- NULL
    if (!is.null(numeric_metrics)) {
      missing_numeric <- setdiff(numeric_metrics, names(de_data))
    }
    if (!is.null(categorical_metrics)) {
      missing_categorical <- setdiff(categorical_metrics, names(de_data))
    }
    if (length(missing_numeric) > 0 || length(missing_categorical) > 0) {
      stop(paste("Missing columns in de_data for min_size", current_min_size, 
                ": numeric =", paste(missing_numeric, collapse=","), 
                "; categorical =", paste(missing_categorical, collapse=",")))
    }
  
    # Generate the fitness function specific to the current minimum subgroup size
    current_fitness_function <- generate_ga_fitness_function(
      current_min_size, de_data, hv_aval_mean, numeric_metrics, categorical_metrics, gene_map, protein_variable_name
    )
  
  popSize <- ga_params$ga_pop_size
  maxiter <- ga_params$ga_max_iter # Corrected from ga_params$max_iter
  pcrossover <- ga_params$pcrossover 
  pmutation <- ga_params$pmutation 
  
  # Orbit-specific parameters
  orbit_check_interval <- ga_params$orbit_check_interval
  dbscan_eps <- ga_params$dbscan_eps
  dbscan_minpts <- ga_params$dbscan_minpts
  
  # Initialize population randomly using custom function
  population <- custom_ga_init(popSize, ga_min, ga_max, ga_type)
  
  # Evaluate initial population
  fitness_values <- apply(population, 1, current_fitness_function)
  
  # Initialize best_overall_chromosome to a default valid chromosome (e.g., all zeros)
  # This prevents returning NULL if no "better" solutions are found later, or if initial fitness is problematic.
  best_overall_chromosome_default <- rep(0, length(ga_min))
  best_overall_fitness_default <- evaluate_subgroup_fitness_bo(
    convert_chromosome_to_atomic_vector(best_overall_chromosome_default, gene_map, atomic_rules_map, atomic_rules_list, n_atomic_rules_total),
    de_data, hv_aval_mean, atomic_rules_list, protein_variable_name, current_min_size
  )
  best_overall_chromosome <- best_overall_chromosome_default
  best_overall_fitness <- best_overall_fitness_default
  
  # If any initial population member is better, update best_overall
  if (length(fitness_values) > 0 && max(fitness_values, na.rm = TRUE) > best_overall_fitness) { 
    best_overall_fitness <- max(fitness_values, na.rm = TRUE)
    best_overall_chromosome <- population[which.max(fitness_values), , drop = FALSE]
  }
  
  # Main GA loop
  for (iter in 1:maxiter) {
    # message(sprintf("GA Orbit (MinSize %d): Iteration %d/%d, Current Pop Max Fitness: %.4f, Overall Best: %.4f", current_min_size, iter, maxiter, max(fitness_values, na.rm = TRUE), best_overall_fitness)) 
    
    # 1. Orbit Discovery and Niche Elite Preservation
    niche_elites <- NULL
    if (iter == 1 || (iter %% orbit_check_interval == 0 && iter < maxiter)) { # Check periodically
      # Convert current population to atomic rule vectors for clustering
      atomic_vectors_pop <- t(apply(population, 1, convert_chromosome_to_atomic_vector, 
                                    gene_map, atomic_rules_map, atomic_rules_list, n_atomic_rules_total))
      
      # Filter out rows that are all zeros (no active rules) or all NAs (conversion failed)
      valid_rows_idx <- apply(atomic_vectors_pop, 1, function(row) !all(row == 0) && !any(is.na(row)))
      
      # message(sprintf("GA Orbit (MinSize %d): Iter %d, Valid atomic vectors for DBSCAN: %d / %d (from popSize %d)", current_min_size, iter, sum(valid_rows_idx), popSize, popSize)) 
      
      if (sum(valid_rows_idx) >= dbscan_minpts) { # Ensure enough valid points for clustering
        atomic_vectors_for_dbscan <- atomic_vectors_pop[valid_rows_idx, , drop = FALSE]
        
        # Perform DBSCAN clustering using Hamming distance
        dbscan_result <- tryCatch({
          # Compute Hamming distance matrix for binary vectors
          # Hamming distance = number of positions where vectors differ
          n_vectors <- nrow(atomic_vectors_for_dbscan)
          hamming_dist <- as.dist(apply(combn(n_vectors, 2), 2, function(idx) {
            sum(atomic_vectors_for_dbscan[idx[1], ] != atomic_vectors_for_dbscan[idx[2], ])
          }))
          
          # Apply DBSCAN with pre-computed Hamming distance matrix
          dbscan::dbscan(hamming_dist, eps = dbscan_eps, minPts = dbscan_minpts)
        }, error = function(e) {
          warning(paste("DBSCAN failed in GA orbit detection for MinSize", current_min_size, "Iter", iter, ":", e$message, ". Skipping niching for this iteration."))
          return(NULL)
        })
        
        if (!is.null(dbscan_result)) {
          cluster_labels <- dbscan_result$cluster
          unique_clusters <- unique(cluster_labels)
          
          # message(sprintf("GA Orbit (MinSize %d): Iter %d, Unique clusters found: %d", current_min_size, iter, length(unique_clusters))) 
          
          # Identify niche elites (best individual per cluster)
          for (cluster_id in unique_clusters) {
            if (cluster_id != 0) { # Exclude noise points (cluster 0)
              cluster_indices_in_dbscan_data <- which(cluster_labels == cluster_id)
              # Map back to original population indices
              original_pop_indices <- which(valid_rows_idx)[cluster_indices_in_dbscan_data]
              
              if (length(original_pop_indices) > 0) {
                best_in_cluster_idx <- original_pop_indices[which.max(fitness_values[original_pop_indices])]
                niche_elites <- rbind(niche_elites, population[best_in_cluster_idx, , drop = FALSE])
              }
            }
          }
        }
      } else {
        # message(sprintf("GA Orbit (MinSize %d): Not enough valid points (%d) for DBSCAN at iteration %d. Skipping niching.", current_min_size, sum(valid_rows_idx), iter))
      }
    }
    
    # Initialize next generation population
    next_population <- matrix(nrow = 0, ncol = ncol(population))
    
    # Add niche elites to the next generation first
    if (!is.null(niche_elites)) {
      # Remove duplicates among niche elites if any (e.g., if multiple clusters point to same best solution)
      niche_elites <- unique(niche_elites)
      next_population <- niche_elites
      # message(sprintf("GA Orbit (MinSize %d): Iter %d, Niche elites added: %d", current_min_size, iter, nrow(niche_elites))) 
    }
    
    # Determine how many slots are left to fill
    remaining_slots <- popSize - nrow(next_population)
    
    if (remaining_slots > 0) {
      # 2. Selection, Crossover, Mutation for the rest of the population
      # Use GA package's internal functions
      
      # Selection (e.g., Roulette Wheel or Tournament)
      # Handle cases where fitness_values might be all 0 or very low, leading to issues with selection
      if (sum(fitness_values, na.rm = TRUE) > 0 && !all(is.na(fitness_values))) { # Only proceed with selection if there's some valid fitness
        selected_indices <- custom_ga_selection_rw(fitness_values, remaining_slots) # Use custom selection
        selected_parents <- population[selected_indices, , drop = FALSE]
      } else {
        # If all fitness is zero or NA, select randomly to try and introduce diversity
        # message(sprintf("GA Orbit (MinSize %d): Iter %d, All fitness zero/NA, selecting random parents.", current_min_size, iter)) 
        selected_parents <- population[sample(1:popSize, remaining_slots, replace = TRUE), , drop = FALSE]
      }
      
      # Crossover and Mutation
      offspring <- matrix(nrow = remaining_slots, ncol = ncol(population))
      # Ensure there are at least two parents for crossover
      if (nrow(selected_parents) >= 2) {
        for (i in seq(1, remaining_slots, by = 2)) {
          parent1_idx <- sample(1:nrow(selected_parents), 1)
          parent2_idx <- sample(1:nrow(selected_parents), 1)
          
          parent1 <- selected_parents[parent1_idx, , drop = FALSE]
          parent2 <- selected_parents[parent2_idx, , drop = FALSE]
          
          # Crossover
          if (runif(1) < pcrossover) {
            crossed <- custom_ga_crossover_sp(parent1, parent2, ga_min, ga_max, ga_type) # Use custom crossover
            child1 <- crossed[1, , drop = FALSE]
            child2 <- crossed[2, , drop = FALSE]
          } else {
            child1 <- parent1
            child2 <- parent2
          }
          
          # Mutation
          child1 <- custom_ga_mutation(child1, ga_min, ga_max, ga_type, pmutation) # Use custom mutation
          child2 <- custom_ga_mutation(child2, ga_min, ga_max, ga_type, pmutation) # Use custom mutation
          
          offspring[i, ] <- child1
          if (i + 1 <= remaining_slots) {
            offspring[i+1, ] <- child2
          }
        }
      } else { # Fallback if not enough parents for crossover (e.g., remaining_slots is 1 or selected_parents has < 2 rows)
        # message(sprintf("GA Orbit (MinSize %d): Iter %d, Not enough parents for crossover, mutating remaining slots.", current_min_size, iter)) 
        # Just mutate the available parents to fill remaining slots
        if (nrow(selected_parents) > 0) {
          offspring[1:remaining_slots, ] <- custom_ga_mutation(selected_parents[sample(1:nrow(selected_parents), remaining_slots, replace=TRUE), , drop = FALSE], ga_min, ga_max, ga_type, pmutation)
        } else { # If no parents at all, generate random offspring
          # message(sprintf("GA Orbit (MinSize %d): Iter %d, No parents available, generating random offspring.", current_min_size, iter)) 
          offspring[1:remaining_slots, ] <- custom_ga_init(remaining_slots, ga_min, ga_max, ga_type) # Use custom init
        }
      }
      
      # Add offspring to the next population
      next_population <- rbind(next_population, offspring)
    }
    
    # Update population and fitness values for the next iteration
    population <- next_population
    fitness_values <- apply(population, 1, current_fitness_function)
    
    # Update overall best solution found so far
    current_best_fitness_this_iter <- max(fitness_values, na.rm = TRUE) 
    if (current_best_fitness_this_iter > best_overall_fitness) {
      best_overall_fitness <- current_best_fitness_this_iter
      best_overall_chromosome <- population[which.max(fitness_values), , drop = FALSE]
    }
    
    # Check for convergence (optional, based on ga_run_limit)
    # This part is more complex to implement manually; for simplicity, we run for maxiter.
    # If a `run` limit is desired, it would involve tracking generations without improvement.
  }
  
  # message(sprintf("GA Orbit (Final): MinSize %d, Optimal Fitness: %.4f", current_min_size, best_overall_fitness)) 
  
  # Return the results of this GA run
  return(list(
    min_subgroup_size = current_min_size,
    optimal_fitness = best_overall_fitness,
    best_chromosome = best_overall_chromosome
  ))
  
  }, error = function(e) {
    # Return error information that can be caught by parallel processing
    error_msg <- paste("GA Orbit failed for min_size", current_min_size, ":", e$message)
    warning(error_msg)
    # Return a recognizable error structure instead of stopping
    return(structure(list(message = error_msg, call = sys.call()), class = "try-error"))
  })
}

# --- Function: Orchestrate Parallel Orbit-Aware GA Runs ---
#' @param min_subgroup_sizes_to_test A vector of minimum subgroup sizes to test.
#' @param ga_params A list of GA-specific parameters.
#' @param de_data The dataframe for disease patients.
#' @param hv_aval_mean The mean protein level for healthy volunteers.
#' @param numeric_metrics A character vector of numeric clinical metrics.
#' @param categorical_metrics A character vector of categorical clinical metrics.
#' @param gene_map The mapping from gene indices to metric names/levels.
#' @param ga_min Vector of minimum values for each gene.
#' @param ga_max Vector of maximum values for each gene.
#' @param ga_type Vector of types ("binary", "real") for each gene.
#' @param protein_variable_name The name of the protein variable to use for fitness calculation.
#' @param atomic_rules_map The map from atomic rule strings to their binary indices.
#' @param atomic_rules_list The list of atomic rule strings.
#' @param n_atomic_rules_total Total number of atomic rules.
#' @param cl_outer A pre-existing parallel cluster object (from makeCluster).
#' @return A list of results from each GA run, indexed by min_subgroup_size.
orchestrate_parallel_ga_orbit <- function(min_subgroup_sizes_to_test, ga_params, de_data, hv_aval_mean,
                                          numeric_metrics, categorical_metrics, gene_map,
                                          ga_min, ga_max, ga_type, protein_variable_name,
                                          atomic_rules_map, atomic_rules_list, n_atomic_rules_total, cl_outer) {
  
  if (!isTRUE(ga_params$benchmark_mode)) {
    message("\n--- Running Orbit-Aware Genetic Algorithm for different min_subgroup_size values (in parallel) ---\n")
  }
  
  registerDoSNOW(cl_outer)
  
  # Export all necessary objects and functions to the workers in the cluster.
  clusterExport(cl_outer, c(
    "run_genetic_algorithm_orbit_for_min_size", "define_ga_encoding_orbit", 
    "generate_ga_fitness_function", "convert_chromosome_to_atomic_vector",
    "custom_ga_init", "custom_ga_selection_rw", "custom_ga_crossover_sp", "custom_ga_mutation", 
    "de_data", "hv_aval_mean", "numeric_metrics", "categorical_metrics", "gene_map",
    "ga_min", "ga_max", "ga_type", "ga_params", "protein_variable_name",
    "atomic_rules_map", "atomic_rules_list", "n_atomic_rules_total"
  ), envir = environment())
  
  # Load required packages within each worker's R session
  clusterEvalQ(cl_outer, {
    library(dplyr)
    library(GA) 
    library(rlang)
    library(dbscan) 
    library(stringdist) 
    source(here::here('common_fct.R')) 
  })
  
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
  
  all_results_ga_orbit_raw <- foreach(
    current_min_size = min_subgroup_sizes_to_test,
    .packages = c("dplyr", "GA", "rlang", "dbscan", "stringdist"),
    .verbose = FALSE,
    .errorhandling = 'pass', 
    .options.snow = list(progress = progress_update_fun)
  ) %dopar% {
    result <- run_genetic_algorithm_orbit_for_min_size(
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
      protein_variable_name = protein_variable_name,
      atomic_rules_map = atomic_rules_map,
      atomic_rules_list = atomic_rules_list,
      n_atomic_rules_total = n_atomic_rules_total
    )
    return(result)
  }
  if (!isTRUE(ga_params$benchmark_mode)) {
    message("\nAll Orbit-Aware Genetic Algorithm runs complete.")
  }
  
  if (is.list(all_results_ga_orbit_raw) && length(all_results_ga_orbit_raw) > 0) {
    # Check for errors in parallel processing results
    error_indices <- sapply(all_results_ga_orbit_raw, inherits, "try-error")
    if (any(error_indices)) {
      error_messages <- sapply(all_results_ga_orbit_raw[error_indices], function(x) x$message)
      warning(paste0("GA Orbit parallel processing encountered ", sum(error_indices), " errors:\n", 
                    paste(unique(error_messages), collapse = "\n")))
    }
    
    # Filter out any error objects from parallel processing
    all_results_ga_orbit_raw_filtered <- all_results_ga_orbit_raw[!error_indices]
    
    if (length(all_results_ga_orbit_raw_filtered) > 0 && all(sapply(all_results_ga_orbit_raw_filtered, is.list)) && all(sapply(all_results_ga_orbit_raw_filtered, function(x) !is.null(x$min_subgroup_size)))) {
      result_names <- sapply(all_results_ga_orbit_raw_filtered, function(x) as.character(x$min_subgroup_size))
      ordered_results <- all_results_ga_orbit_raw_filtered[order(as.numeric(result_names))]
      names(ordered_results) <- sort(as.numeric(result_names))
      all_results_ga_orbit_raw <- ordered_results 
    } else {
      warning("Results structure from parallel processing is unexpected or empty after filtering. Cannot assign names based on min_subgroup_size.")
      all_results_ga_orbit_raw <- list() 
    }
  } else {
    warning("all_results_ga_orbit_raw is not a list or is empty after parallel processing.")
    all_results_ga_orbit_raw <- list() 
  }
  
  return(all_results_ga_orbit_raw)
}

# --- Function: Process Orbit-Aware GA Results into Detailed DataFrames (re-used from ga_method.R) ---
# This function decodes the optimal chromosomes from the GA runs
# into human-readable rules and calculates various statistics for each identified subgroup.
# It also assigns patients to "Subgroup", "Rest_DE", and "HV" groups.
process_ga_orbit_results <- function(all_results_ga_orbit_raw, de_data, hv_data, hv_aval_mean,
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
    message("\n--- Decoding and Analyzing Optimal Rules (Orbit-Aware GA) ---")
  }
  valid_results_ga_orbit <- lapply(all_results_ga_orbit_raw, function(x) x)
  
  # message(sprintf("process_ga_orbit_results: Number of raw results entries received: %d", length(valid_results_ga_orbit))) 
  
  if (length(valid_results_ga_orbit) > 0) {
    for (i in seq_along(valid_results_ga_orbit)) {
      result_entry <- valid_results_ga_orbit[[i]]
      current_min_subgroup_size <- result_entry$min_subgroup_size
      current_optimal_fitness <- result_entry$optimal_fitness
      current_best_chromosome_ga_orbit <- result_entry$best_chromosome
      
      # message(sprintf("  Processing MinSize %d (Optimal Fitness: %.4f)", current_min_subgroup_size, current_optimal_fitness)) 
      # message(paste("  Best Chromosome (first 10 genes):", paste(round(current_best_chromosome_ga_orbit[1:min(10, length(current_best_chromosome_ga_orbit))], 2), collapse = " "))) 
      
      
      rules_table_for_current_solution <- data.frame(
        Metric = character(),
        Operator = character(),
        Value_or_Levels = character(),
        stringsAsFactors = FALSE
      )
      final_filter_expressions_apply_i <- list()
      
      # Check if best_chromosome is valid before proceeding with decoding
      if (is.null(current_best_chromosome_ga_orbit) || length(current_best_chromosome_ga_orbit) == 0 || all(is.na(current_best_chromosome_ga_orbit))) {
        warning(paste("No valid best chromosome found for min subgroup size (Orbit-Aware GA):", current_min_subgroup_size, ". Skipping decoding."))
        next 
      }
      
      # Decode Numeric Rules and add to temporary rules table
      for (metric in numeric_metrics) {
        map_info <- gene_map[[metric]]
        if (is.null(map_info)) next 
        
        use_metric <- round(current_best_chromosome_ga_orbit[map_info$use_idx])
        
        if (use_metric == 1) {
          op_type <- round(current_best_chromosome_ga_orbit[map_info$op_idx])
          threshold_val <- current_best_chromosome_ga_orbit[map_info$val_idx]
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
        if (is.null(map_info)) next 
        
        use_metric <- round(current_best_chromosome_ga_orbit[map_info$use_idx])
        
        if (use_metric == 1) {
          selected_levels <- c()
          for (j in 1:length(map_info$levels)) {
            level_gene_val <- round(current_best_chromosome_ga_orbit[map_info$levels_start_idx + j - 1])
            if (level_gene_val == 1) {
              selected_levels <- c(selected_levels, map_info$levels[j])
            }
          }
          
          if (length(selected_levels) == 0) {
            # message(sprintf("    MinSize %d: Categorical metric %s used (gene=1) but no levels selected. Skipping rule.", current_min_subgroup_size, metric)) 
            next 
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
      
      if (nrow(rules_table_for_current_solution) > 0) {
        temp_filtered_de_data <- de_data
        is_valid_subgroup_for_ids <- TRUE
        
        for (filter_expr in final_filter_expressions_apply_i) {
          temp_filtered_de_data <- temp_filtered_de_data %>% filter(!!filter_expr)
          if (nrow(temp_filtered_de_data) < current_min_subgroup_size) {
            is_valid_subgroup_for_ids <- FALSE
            # message(sprintf("    MinSize %d: Subgroup size (%d) fell below min_subgroup_size (%d) during filtering. Invalid subgroup.", current_min_subgroup_size, nrow(temp_filtered_de_data), current_min_subgroup_size)) 
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
            
            # message(sprintf("    MinSize %d: Valid subgroup found. Patients: %d, Fitness: %.4f, pVal: %.4f", current_min_subgroup_size, num_patients_current, current_optimal_fitness, pVal)) 
            
          } else {
            # message(sprintf("    MinSize %d: Subgroup mean is NaN or non-positive. Invalidating rules.", current_min_subgroup_size)) 
            rules_table_for_current_solution <- rules_table_for_current_solution[0,]
          }
        } else {
          # message(sprintf("    MinSize %d: Final subgroup size (%d) below min_subgroup_size (%d). Invalidating rules.", current_min_subgroup_size, nrow(temp_filtered_de_data), current_min_subgroup_size)) 
          rules_table_for_current_solution <- rules_table_for_current_solution[0,]
        }
      } else {
        # message(sprintf("    MinSize %d: No active rules generated from chromosome. Invalidating rules.", current_min_subgroup_size)) 
        rules_table_for_current_solution <- rules_table_for_current_solution[0,]
      }
      
      # message(sprintf("  MinSize %d: Rules table rows after processing: %d", current_min_subgroup_size, nrow(rules_table_for_current_solution))) 
      
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
  # message(sprintf("process_ga_orbit_results: Total solutions added to details_df: %d", nrow(all_solutions_details_df))) 
  
  all_solutions_details_df_sorted <- all_solutions_details_df %>%
    arrange(desc(OptimalFitness), MinSubgroupSize, Metric)
  
  patient_groups_df_sorted <- patient_groups_df %>%
    arrange(MinSubgroupSize, Group, USUBJID)
  
  return(list(
    all_solutions_details_df_sorted = all_solutions_details_df_sorted,
    patient_groups_df_sorted = patient_groups_df_sorted
  ))
}


# --- Main Orchestrating Function for Orbit-Aware GA ---
#' @param params A list of parameters for the analysis.
#' @param df The full clinical dataframe.
#' @param dfThisProtein The dataframe filtered for the specific protein.
#' @param hv_data The dataframe for healthy volunteers.
#' @param de_data The dataframe for disease patients.
#' @param hv_aval_mean The mean protein level for healthy volunteers.
#' @param atomic_rules_map The map from atomic rule strings to their binary indices.
#' @param atomic_rules_list The list of atomic rule strings.
#' @param n_atomic_rules_total Total number of atomic rules.
#' @param cl_outer An optional pre-existing parallel cluster object.
#' @return A list containing all generated data frames for further programmatic access.
run_genetic_algorithm_orbit_analysis <- function(params, df, dfThisProtein, hv_data, de_data, hv_aval_mean,
                                                 atomic_rules_map, atomic_rules_list, n_atomic_rules_total, cl_outer = NULL) {
  
  # 1. Define Chromosome Encoding and Gene Map for GA
  ga_encoding_info <- define_ga_encoding_orbit(de_data, params$numeric_metrics, params$categorical_metrics)
  ga_min <- ga_encoding_info$ga_min
  ga_max <- ga_encoding_info$ga_max
  ga_type <- ga_encoding_info$ga_type
  gene_map <- ga_encoding_info$gene_map
  
  all_results_ga_orbit_raw <- list()
  
  if (params$runSim) {
    # 2. Orchestrate Parallel Orbit-Aware GA Runs
    all_results_ga_orbit_raw <- orchestrate_parallel_ga_orbit(
      min_subgroup_sizes_to_test = params$min_subgroup_sizes_to_test,
      ga_params = params, # Pass all GA-related parameters (pop_size, max_iter, run_limit, orbit_check_interval, etc.)
      de_data = de_data,
      hv_aval_mean = hv_aval_mean,
      numeric_metrics = params$numeric_metrics,
      categorical_metrics = params$categorical_metrics,
      gene_map = gene_map,
      ga_min = ga_min,
      ga_max = ga_max,
      ga_type = ga_type,
      protein_variable_name = params$variable_for_use,
      atomic_rules_map = atomic_rules_map,
      atomic_rules_list = atomic_rules_list,
      n_atomic_rules_total = n_atomic_rules_total,
      cl_outer = cl_outer
    )
    message("\n--- Simulation results (not saved as requested) ---")
    
  } else {
    message("\n--- Loading pre-computed simulation results for Orbit-Aware GA ---")
    tryCatch({
      all_results_ga_orbit_raw <- aws.s3::s3read_using(FUN = readRDS,
                                                       object = paste0('data/', 'OptSubgroup_Rules_GA_Orbit_vs_HV_', params$thisProtein, '_', params$variable_for_use, '.rds'),
                                                       bucket = params$.arv_save$collection)
    }, error = function(e) {
      warning(paste("Error loading pre-computed Orbit-Aware GA data from S3:", e$message, ". Generating dummy results for demonstration."))
      # Fallback: Generate dummy results if loading fails.
      chromosome_length <- length(ga_min) # Use the length from encoding
      all_results_ga_orbit_raw <- list(
        "10" = list(min_subgroup_size = 10, optimal_fitness = 1.6, best_chromosome = sample(0:1, chromosome_length, replace=TRUE)),
        "11" = list(min_subgroup_size = 11, optimal_fitness = 1.7, best_chromosome = sample(0:1, chromosome_length, replace=TRUE)),
        "12" = list(min_subgroup_size = 12, optimal_fitness = 1.8, best_chromosome = sample(0:1, chromosome_length, replace=TRUE))
      )
      names(all_results_ga_orbit_raw) <- sapply(all_results_ga_orbit_raw, function(x) as.character(x$min_subgroup_size))
    })
  }
  
  if (length(all_results_ga_orbit_raw) == 0) {
    stop("No Orbit-Aware Genetic Algorithm results available for further processing. Check data loading or simulation parameters.")
  }
  
  # 3. Process Orbit-Aware GA Results into Detailed DataFrames
  processed_dfs_ga_orbit <- process_ga_orbit_results(all_results_ga_orbit_raw, de_data, hv_data, hv_aval_mean,
                                                     params$numeric_metrics, params$categorical_metrics, gene_map, params)
  all_solutions_details_df_sorted_ga_orbit <- processed_dfs_ga_orbit$all_solutions_details_df_sorted
  patient_groups_df_sorted_ga_orbit <- processed_dfs_ga_orbit$patient_groups_df_sorted
  
  # 4. Generate Analysis Results
  analysis_results_ga_orbit <- generate_analysis_results_common(all_solutions_details_df_sorted_ga_orbit, params)
  
  # 5. Plot Fitness vs. Minimum Subgroup Size
  if (!isTRUE(params$benchmark_mode)) {
    plot_fitness_vs_subgroup_size_common(all_results_ga_orbit_raw, hv_aval_mean, params)
    
    # 6. Generate Numeric Cutoff Plots (if applicable)
    if (params$verbose_analysis_output) {
      generate_numeric_cutoff_plots_common(analysis_results_ga_orbit, params)
    }
  }
  
  # 7. Generate Final Rules Table
  dfRulesFinal_ga_orbit <- generate_final_rules_table_common(all_solutions_details_df_sorted_ga_orbit, params)
  
  # 8. Plot Protein Expression in Subgroups
  if (!isTRUE(params$benchmark_mode)) {
    plot_protein_expression_in_subgroups_common(dfRulesFinal_ga_orbit, patient_groups_df_sorted_ga_orbit, df, params)
  }
  
  if (!isTRUE(params$benchmark_mode)) {
    message("\n--- Patient Subgroup Discovery (Orbit-Aware GA) process complete ---")
  }
  
  return(list(
    all_solutions_details = all_solutions_details_df_sorted_ga_orbit,
    patient_group_details = patient_groups_df_sorted_ga_orbit,
    analysis_results = analysis_results_ga_orbit,
    optimal_rules_table = dfRulesFinal_ga_orbit
  ))
}
