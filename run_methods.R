# --- run_analysis.R: Main Script to Run Patient Subgroup Discovery Methods ---

# Load the 'here' package for robust file path management
suppressMessages(library(here))
suppressMessages(library(tictoc)) # For measuring execution time
suppressMessages(library(ggplot2)) # Ensure ggplot2 is loaded for comparison plot
suppressMessages(library(dplyr)) # Ensure dplyr is loaded for data manipulation
suppressMessages(library(tidyr)) # For pivot_longer
suppressMessages(library(plotly)) # For interactive 3D plotting

# Source common functions and libraries
# This loads shared utility functions and common R packages.
source(here::here('common_fct.R'))

# Source method-specific functions
# This loads functions unique to Bayesian Optimization and Genetic Algorithm.
# The original BO with orbits is now bo_orbit_method.R
source(here::here('bo_orbit_method.R')) # With orbits
source(here::here('bo_method.R'))      # Without orbits)
source(here::here('ga_orbit_method.R')) # Orbit-Aware GA method
source(here::here('ga_method.R'))     # GA without orbits
source(here::here('greedy_method.R')) # Source greedy method
source(here::here('exhaustive_search.R')) # Source exhaustive search for ground truth validation

# --- Define Base Parameters (common across all runs) ---
base_params <- list(
  study = "BP44617",
  # Setting numeric_metrics to NULL to focus on categorical/discrete search space
  # numeric_metrics = c('NEISCORE', 'OCFSCORE', 'TFBUT'), # Can be a subset of c('LENGTH', 'NEISCORE', 'OCFSCORE', 'TFBUT') or NULL
  numeric_metrics = NULL, # Can be a subset of c('LENGTH', 'NEISCORE', 'OCFSCORE', 'TFBUT') or NULL
  categorical_metrics = c('MGLSRT', 'DED_SEVERITY', 'DED_TYPE'),
  variable_for_use = 'AVAL_BY_LENGTH',
  thisProtein = 'ABC',
  runSim = TRUE, # Set to TRUE to run simulation, FALSE to load saved results
  nCore = 4, # Number of cores for parallel processing
  min_subgroup_sizes_to_test = 5:40, # Range of minimum subgroup sizes to test
  # BO-specific parameters that are NOT varied in the systematic comparison, but need to be passed
  # These are for the original BO (now bo_orbit_method.R) and the new bo_method.R
  max_initial_active_rules = 5,
  bo_epsilon_clustering = 1, # Only relevant for bo_orbit_method.R
  bo_min_samples_clustering = 5, # Only relevant for bo_orbit_method.R
  # Greedy-specific parameters that are NOT varied in the systematic comparison, but need to be passed
  max_rules_per_subgroup = 2,
  # Orbit-Aware GA specific parameters
  orbit_check_interval = 10, # How often to check for orbits (every X generations)
  dbscan_eps = 1, # Epsilon for DBSCAN in orbit detection (Hamming distance)
  dbscan_minpts = 5, # MinPts for DBSCAN in orbit detection
  # General display parameters
  show_rules = FALSE,
  verbose_analysis_output = FALSE,
  benchmark_mode = TRUE, # Set to FALSE to show detailed rule tables after each run
  useRealData = FALSE, # Initialize to FALSE, set to TRUE if real data loads successfully
  .arv_save = list(
    collection = 'ABC',
    filename = "reports/ABC", # Generic filename, specific reports will be generated per run
    report_html = "./reports/ABC.html"
  )
)

# --- Define Stability Testing Parameters ---
stability_params <- list(
  n_stability_runs = 20, # Number of times to run each method for stability testing
  stability_min_subgroup_sizes = c(10, 20, 30), # Selected subgroup sizes to test
  stability_seed_base = 1234, # Base seed for reproducible stability testing
  n_param_sets_per_method = 5 # Number of parameter sets to randomly sample per method for stability testing
)

# --- Define Hyperparameter Ranges for Systematic Comparison ---

# Bayesian Optimization Parameter Ranges (these will be varied)
# These parameters will be used for BOTH bo_method.R (no orbits) and bo_orbit_method.R (with orbits)
bo_initial_samples_vals <- c(100, 300, 500) # Expanded range
bo_iterations_vals <- c(100, 500, 1000) # Expanded range
bo_acquisition_random_search_vals <- c(2000) # Expanded range
bo_xi_vals <- c(0.01, 0.1, 1.0, 10.0) # Expanded range for exploration parameter

# Genetic Algorithm Parameter Ranges (these will be varied)
ga_pop_size_vals <- c(50, 100)
ga_max_iter_vals <- c(20, 50, 100, 150)
ga_run_limit_vals <- c(10, 30, 60, 100)
# Add crossover and mutation probabilities to GA params
ga_pcrossover_vals <- c(0.8) # Fixed for now, could be varied
ga_pmutation_vals <- c(0.1) # Fixed for now, could be varied

# Orbit-Aware Genetic Algorithm Parameter Ranges (these will be varied)
ga_orbit_pop_size_vals <- c(50, 100) # Expanded range
ga_orbit_max_iter_vals <- c(20, 50, 100, 150) # Expanded range
ga_orbit_run_limit_vals <- c(10, 30, 60, 100) # Expanded range
ga_orbit_pcrossover_vals <- c(0.8)
ga_orbit_pmutation_vals <- c(0.1)


# Greedy Algorithm Parameter Ranges (these will be varied)
greedy_max_rules_per_subgroup_vals <- seq(2, 30, 2) # Varying the max number of rules


# --- Generate All Hyperparameter Combinations ---

# Generate BO Parameter Sets (used for both BO (No Orbits) and BO (With Orbits))
bo_param_grid <- expand.grid(
  bo_initial_samples = bo_initial_samples_vals,
  bo_iterations = bo_iterations_vals,
  bo_acquisition_random_search = bo_acquisition_random_search_vals,
  bo_xi = bo_xi_vals, # Include xi in the grid
  stringsAsFactors = FALSE
)
bo_param_sets <- lapply(1:nrow(bo_param_grid), function(i) as.list(bo_param_grid[i, ]))
names(bo_param_sets) <- paste0(
  "BO_IS", bo_param_grid$bo_initial_samples,
  "_I", bo_param_grid$bo_iterations,
  "_ARS", bo_param_grid$bo_acquisition_random_search,
  "_XI", bo_param_grid$bo_xi
)

# Generate GA Parameter Sets
ga_param_grid <- expand.grid(
  ga_pop_size = ga_pop_size_vals,
  ga_max_iter = ga_max_iter_vals,
  ga_run_limit = ga_run_limit_vals,
  pcrossover = ga_pcrossover_vals, # Add to GA params
  pmutation = ga_pmutation_vals, # Add to GA params
  stringsAsFactors = FALSE
)
ga_param_sets <- lapply(1:nrow(ga_param_grid), function(i) as.list(ga_param_grid[i, ]))
names(ga_param_sets) <- paste0(
  "GA_PS", ga_param_grid$ga_pop_size,
  "_MI", ga_param_grid$ga_max_iter,
  "_RL", ga_param_grid$ga_run_limit
)

# Generate Orbit-Aware GA Parameter Sets
ga_orbit_param_grid <- expand.grid(
  ga_pop_size = ga_orbit_pop_size_vals,
  ga_max_iter = ga_orbit_max_iter_vals,
  ga_run_limit = ga_orbit_run_limit_vals,
  pcrossover = ga_orbit_pcrossover_vals,
  pmutation = ga_orbit_pmutation_vals,
  stringsAsFactors = FALSE
)
ga_orbit_param_sets <- lapply(1:nrow(ga_orbit_param_grid), function(i) as.list(ga_orbit_param_grid[i, ]))
names(ga_orbit_param_sets) <- paste0(
  "GA_Orbit_PS", ga_orbit_param_grid$ga_pop_size,
  "_MI", ga_orbit_param_grid$ga_max_iter,
  "_RL", ga_orbit_param_grid$ga_run_limit
)


# Generate Greedy Parameter Sets
greedy_param_grid <- expand.grid(
  max_rules_per_subgroup = greedy_max_rules_per_subgroup_vals,
  stringsAsFactors = FALSE
)
greedy_param_sets <- lapply(1:nrow(greedy_param_grid), function(i) as.list(greedy_param_grid[i, ]))
names(greedy_param_sets) <- paste0(
  "Greedy_MRPS", greedy_param_grid$max_rules_per_subgroup
)


# --- Global Data Loading and Preprocessing ---
message("--- Global Data Loading and Preprocessing ---")

setup_environment()

df_global <- NULL
dfThisProtein_global <- NULL
hv_data_global <- NULL
de_data_global <- NULL
hv_aval_mean_global <- NULL

if (base_params$useRealData) {
  # Try to load real data from S3
  tryCatch({
    filename_clin_global <- paste("data/", base_params$study, "_ZOC_LB_Clin", ".rds", sep="")
    filename_dic_global <- paste("data/", base_params$study, "_ZOC_LB_dic", ".rds", sep="")
    
    df_global <- aws.s3::s3read_using(FUN = readRDS, object = filename_clin_global, bucket = base_params$.arv_save$collection)
    dfOeDic_global <- aws.s3::s3read_using(FUN = readRDS, object = filename_dic_global, bucket = base_params$.arv_save$collection)
    message("Global data loaded successfully from S3.")
  }, error = function(e) {
    stop(paste("Error loading real data from S3:", e$message, 
               "\nReal data was requested (useRealData = TRUE) but could not be loaded.",
               "\nSet useRealData = FALSE to use dummy data instead."))
  })
} else {
  # Generate dummy data directly
  message("Using dummy data for analysis (useRealData = FALSE)")
  set.seed(123)
  n_patients <- 200
  df_dummy_global <- data.frame(
    USUBJID = paste0("P", 1:n_patients),
    ARMCD = sample(c("HV", "DE"), n_patients, replace = TRUE, prob = c(0.2, 0.8)),
    LBTESTCD = base_params$thisProtein,
    # Categorical metrics
    MGLSRT = sample(c("0", "1", "2", "3", "4", "5"), n_patients, replace = TRUE),
    DED_SEVERITY = sample(c("Mild", "Moderate", "Severe"), n_patients, replace = TRUE),
    DED_TYPE = sample(c("Aqueous Deficient", "Evaporative", "MIXED TYPE"), n_patients, replace = TRUE),
    LBSTRESU = "ng/mL"
  )
  
  # Add numeric metrics columns if they are specified
  if (!is.null(base_params$numeric_metrics) && length(base_params$numeric_metrics) > 0) {
    for (metric in base_params$numeric_metrics) {
      # Generate realistic dummy values for numeric metrics
      if (metric == "NEISCORE") {
        df_dummy_global[[metric]] <- sample(c(0, 1, 2, 3, 4, 5), n_patients, replace = TRUE, prob = c(0.1, 0.15, 0.2, 0.25, 0.2, 0.1))
      } else if (metric == "OCFSCORE") {
        df_dummy_global[[metric]] <- sample(c(0, 1, 2, 3, 4), n_patients, replace = TRUE, prob = c(0.2, 0.25, 0.3, 0.2, 0.05))
      } else if (metric == "TFBUT") {
        # Tear Film Break-Up Time in seconds - typically between 1-20 seconds
        df_dummy_global[[metric]] <- round(runif(n_patients, min = 1, max = 20), 1)
      } else {
        # Generic numeric values for any other metrics
        df_dummy_global[[metric]] <- round(runif(n_patients, min = 0, max = 10), 2)
      }
    }
  }
  
  # Add the main variable for analysis
  df_dummy_global[[base_params$variable_for_use]] <- ifelse(df_dummy_global$ARMCD == "HV", runif(n_patients, 10, 50), runif(n_patients, 20, 100))
  # Introduce some structure for rules to find based on categorical variables
  df_dummy_global[[base_params$variable_for_use]] <- ifelse(
    df_dummy_global$DED_SEVERITY == "Severe" & df_dummy_global$DED_TYPE == "Aqueous Deficient",
    df_dummy_global[[base_params$variable_for_use]] * 2,
    df_dummy_global[[base_params$variable_for_use]]
  )
  
  # Add some structure to numeric metrics to make rules more discoverable
  if (!is.null(base_params$numeric_metrics) && length(base_params$numeric_metrics) > 0) {
    for (metric in base_params$numeric_metrics) {
      if (metric %in% names(df_dummy_global)) {
        # Create some correlation between numeric metrics and the outcome
        # Higher severity cases tend to have higher numeric scores
        severity_multiplier <- ifelse(df_dummy_global$DED_SEVERITY == "Severe", 1.3,
                                      ifelse(df_dummy_global$DED_SEVERITY == "Moderate", 1.1, 1.0))
        
        if (metric == "TFBUT") {
          # For TFBUT, lower values are worse, so invert the relationship
          df_dummy_global[[metric]] <- df_dummy_global[[metric]] / severity_multiplier
        } else {
          # For other metrics, higher values are worse
          df_dummy_global[[metric]] <- round(df_dummy_global[[metric]] * severity_multiplier, 
                                             ifelse(metric == "TFBUT", 1, 0))
        }
      }
    }
  }
  df_global <- df_dummy_global
}

if (!is.null(df_global)) {
  df_global = df_global %>%
    # Ensure categorical metrics are factors or characters as needed
    mutate(MGLSRT = as.character(MGLSRT),
           DED_SEVERITY = as.character(DED_SEVERITY),
           DED_TYPE = as.character(DED_TYPE))
  dfThisProtein_global <- df_global %>% filter(LBTESTCD == base_params$thisProtein)
  if (nrow(dfThisProtein_global) == 0) {
    stop(paste("No data found for protein:", base_params$thisProtein, "in global data loading."))
  }
  dfThisProtein_global[[base_params$variable_for_use]] <- as.numeric(dfThisProtein_global[[base_params$variable_for_use]])
  hv_data_global <- dfThisProtein_global %>% filter(ARMCD == "HV")
  de_data_global <- dfThisProtein_global %>% filter(ARMCD == "DE")
  if (nrow(hv_data_global) == 0) {
    stop("No Healthy Volunteer (HV) data found in global data loading.")
  }
  if (nrow(de_data_global) == 0) {
    stop("No Disease (DE) data found in global data loading.")
  }
  hv_aval_mean_global <- mean(hv_data_global[[base_params$variable_for_use]], na.rm = TRUE)
  message(sprintf("\nMean %s for Healthy Volunteers (HV) in global data: %.4f\n\n", base_params$variable_for_use, hv_aval_mean_global))
} else {
  stop("Global data frame (df_global) is NULL. Cannot proceed with preprocessing.")
}

# --- Determine Output Directory and Data Type Suffix for File Names ---
# Include flag for whether numeric metrics are used
numeric_flag <- ifelse(is.null(base_params$numeric_metrics) || length(base_params$numeric_metrics) == 0, "no_numeric", "with_numeric")
data_type_suffix <- paste0(ifelse(base_params$useRealData, "real", "dummy"), "_", numeric_flag, "_data")
output_subdir <- data_type_suffix  # This will be used as the subdirectory name
output_dir <- here::here(paste0('output/', output_subdir))
message(sprintf("Using %s for analysis (%s metrics)", 
                ifelse(base_params$useRealData, "real data", "dummy data"),
                ifelse(numeric_flag == "with_numeric", "with numeric", "categorical only")))
message(sprintf("Results will be saved to: output/%s/", output_subdir))

# --- External Cluster Management (only if running simulations) ---
if (base_params$runSim == TRUE) {
  # Determine the maximum number of cores needed across all planned runs.
  # Assuming nCore is constant in base_params, so we just use that.
  max_n_cores_needed <- base_params$nCore
  message(sprintf("Creating parallel cluster with %d cores...", max_n_cores_needed))
  cl <- makeCluster(max_n_cores_needed, type = "SOCK")
  registerDoSNOW(cl)  # Register the cluster for parallel processing
  message("Parallel cluster created and registered.")
} else {
  cl <- NULL  # No cluster needed for loading results
}

# --- Store Results for Comparison ---
all_comparison_fitness_data <- data.frame()
all_comparison_times_data <- data.frame(Method = character(), ParamSet = character(), Time_Seconds = numeric(), stringsAsFactors = FALSE)

# ======================================================================
#                   EXHAUSTIVE SEARCH VALIDATION
# ======================================================================
# This section runs exhaustive search for ground truth validation
# before any optimization methods are evaluated (for discrete cases only)

message("\n==============================================")
message("      EXHAUSTIVE SEARCH VALIDATION           ")
message("==============================================\n")

# Check if this is a discrete case (no numeric metrics)
is_discrete_case <- is.null(base_params$numeric_metrics) || length(base_params$numeric_metrics) == 0

if (is_discrete_case && base_params$runSim) {
  message("Discrete optimization case detected - running exhaustive search for ground truth validation...")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate GA encoding to get gene map
  ga_encoding_info <- define_ga_encoding(de_data_global, base_params$numeric_metrics, base_params$categorical_metrics)
  gene_map_global <- ga_encoding_info$gene_map
  
  # Check feasibility
  feasibility <- check_exhaustive_search_feasibility(base_params$categorical_metrics, gene_map_global, max_combinations = 1e6)
  
  if (feasibility$feasible) {
    message(sprintf("Search space size: %s combinations - proceeding with exhaustive search", 
                    format(feasibility$search_space_size, big.mark = ",")))
    
    message(sprintf("Running exhaustive search in parallel for %d minimum subgroup sizes...", 
                    length(base_params$min_subgroup_sizes_to_test)))
    
    # Use foreach for parallel computation across subgroup sizes
    exhaustive_results_list <- foreach(min_size = base_params$min_subgroup_sizes_to_test, 
                                       .combine = rbind,
                                       .packages = c("dplyr", "rlang", "here"),
                                       .export = c("exhaustive_search_discrete", "calculate_search_space_size", 
                                                   "generate_all_categorical_chromosomes", "decode_categorical_chromosome",
                                                   "generate_ga_fitness_function", "define_ga_encoding", 
                                                   "de_data_global", "hv_aval_mean_global", "gene_map_global")) %dopar% {
                                                     
                                                     # Run exhaustive search for this minimum subgroup size
                                                     tryCatch({
                                                       exhaustive_result <- exhaustive_search_discrete(
                                                         current_min_subgroup_size = min_size,
                                                         de_data = de_data_global,
                                                         hv_aval_mean = hv_aval_mean_global,
                                                         categorical_metrics = base_params$categorical_metrics,
                                                         gene_map = gene_map_global,
                                                         protein_variable_name = base_params$variable_for_use,
                                                         max_combinations = 1e6
                                                       )
                                                       
                                                       if (!is.null(exhaustive_result)) {
                                                         # Create result row
                                                         result_row <- data.frame(
                                                           MinSubgroupSize = min_size,
                                                           TrueGlobalOptimum = exhaustive_result$best_fitness,
                                                           BestSolutionRules = paste(exhaustive_result$decoded_solution, collapse = " AND "),
                                                           SearchSpaceSize = exhaustive_result$search_space_size,
                                                           EvaluationTime = exhaustive_result$evaluation_time,
                                                           NumEvaluations = exhaustive_result$n_evaluations,
                                                           stringsAsFactors = FALSE
                                                         )
                                                         
                                                         return(result_row)
                                                       } else {
                                                         # Return row with NA values if exhaustive search failed
                                                         return(data.frame(
                                                           MinSubgroupSize = min_size,
                                                           TrueGlobalOptimum = NA,
                                                           BestSolutionRules = "Exhaustive search failed",
                                                           SearchSpaceSize = feasibility$search_space_size,
                                                           EvaluationTime = NA,
                                                           NumEvaluations = NA,
                                                           stringsAsFactors = FALSE
                                                         ))
                                                       }
                                                       
                                                     }, error = function(e) {
                                                       warning(sprintf("Exhaustive search failed for min_size %d: %s", min_size, e$message))
                                                       return(data.frame(
                                                         MinSubgroupSize = min_size,
                                                         TrueGlobalOptimum = NA,
                                                         BestSolutionRules = paste("Error:", e$message),
                                                         SearchSpaceSize = feasibility$search_space_size,
                                                         EvaluationTime = NA,
                                                         NumEvaluations = NA,
                                                         stringsAsFactors = FALSE
                                                       ))
                                                     })
                                                   }
    
    # Convert results to data frame if it isn't already
    if (!is.data.frame(exhaustive_results_list)) {
      exhaustive_global_optima <- data.frame(exhaustive_results_list, stringsAsFactors = FALSE)
    } else {
      exhaustive_global_optima <- exhaustive_results_list
    }
    
    # Save exhaustive search results
    exhaustive_file <- file.path(output_dir, 'exhaustive_search_global_optima.rds')
    saveRDS(exhaustive_global_optima, exhaustive_file)
    
    # Display summary
    message("\n=== EXHAUSTIVE SEARCH RESULTS SUMMARY ===")
    successful_runs <- sum(!is.na(exhaustive_global_optima$TrueGlobalOptimum))
    total_time <- sum(exhaustive_global_optima$EvaluationTime, na.rm = TRUE)
    
    message(sprintf("Successful runs: %d/%d", successful_runs, nrow(exhaustive_global_optima)))
    message(sprintf("Total computation time: %.2f seconds", total_time))
    
    if (successful_runs > 0) {
      message(sprintf("Global optimum range: %.6f - %.6f", 
                      min(exhaustive_global_optima$TrueGlobalOptimum, na.rm = TRUE),
                      max(exhaustive_global_optima$TrueGlobalOptimum, na.rm = TRUE)))
      
      # Show some example results
      message("\nExample results:")
      sample_rows <- head(exhaustive_global_optima[!is.na(exhaustive_global_optima$TrueGlobalOptimum), ], 3)
      for (i in 1:nrow(sample_rows)) {
        row <- sample_rows[i, ]
        message(sprintf("  Min size %d: Optimum = %.6f, Rules = %s", 
                        row$MinSubgroupSize, row$TrueGlobalOptimum, 
                        substr(row$BestSolutionRules, 1, 60)))
      }
    }
    
    message(sprintf("Results saved to: %s", exhaustive_file))
    
  } else {
    message(sprintf("Search space too large (%s combinations) - skipping exhaustive search", 
                    format(feasibility$search_space_size, big.mark = ",")))
    message("Consider reducing categorical metrics or using sampling-based validation")
    
    # Save empty results to indicate exhaustive search was not feasible
    exhaustive_global_optima <- data.frame(
      MinSubgroupSize = base_params$min_subgroup_sizes_to_test,
      TrueGlobalOptimum = NA,
      BestSolutionRules = "Search space too large",
      SearchSpaceSize = feasibility$search_space_size,
      EvaluationTime = NA,
      NumEvaluations = NA,
      stringsAsFactors = FALSE
    )
    
    saveRDS(exhaustive_global_optima, file.path(output_dir, 'exhaustive_search_global_optima.rds'))
  }
  
} else if (!is_discrete_case && base_params$runSim) {
  message("Mixed discrete/continuous optimization case detected - exhaustive search not applicable")
  message("Exhaustive search only works for discrete cases (categorical metrics only)")
  
  # Save empty results to indicate exhaustive search was not applicable
  exhaustive_global_optima <- data.frame(
    MinSubgroupSize = base_params$min_subgroup_sizes_to_test,
    TrueGlobalOptimum = NA,
    BestSolutionRules = "Not applicable (numeric metrics present)",
    SearchSpaceSize = NA,
    EvaluationTime = NA,
    NumEvaluations = NA,
    stringsAsFactors = FALSE
  )
  
  saveRDS(exhaustive_global_optima, file.path(output_dir, 'exhaustive_search_global_optima.rds'))
  
} else {
  message("Exhaustive search validation skipped (results loaded from disk or runSim=FALSE)")
  
  # Try to load existing exhaustive search results
  exhaustive_file <- file.path(output_dir, 'exhaustive_search_global_optima.rds')
  if (file.exists(exhaustive_file)) {
    exhaustive_global_optima <- readRDS(exhaustive_file)
    message(sprintf("Loaded existing exhaustive search results from: %s", exhaustive_file))
  } else {
    message("No existing exhaustive search results found")
    exhaustive_global_optima <- NULL
  }
}

message("\n=== EXHAUSTIVE SEARCH VALIDATION COMPLETED ===\n")

# ======================================================================
#               MAIN OPTIMIZATION METHODS EVALUATION
# ======================================================================

# --- Check runSim Parameter ---
if (base_params$runSim == FALSE) {
  message("\n==============================================")
  message("          LOADING SAVED RESULTS MODE         ")
  message("==============================================\n")
  
  # Determine which data type files to load
  numeric_flag <- ifelse(is.null(base_params$numeric_metrics) || length(base_params$numeric_metrics) == 0, "no_numeric", "with_numeric")
  data_type_suffix <- paste0(ifelse(base_params$useRealData, "real", "dummy"), "_", numeric_flag, "_data")
  output_subdir <- data_type_suffix  # This will be used as the subdirectory name
  
  # Check if saved results exist
  fitness_file <- here::here(paste0('output/', output_subdir, '/all_comparison_fitness_data.rds'))
  if (!file.exists(fitness_file)) {
    stop(paste0("Saved results not found for ", ifelse(base_params$useRealData, "real data", "dummy data"), 
                ". Please run with runSim = TRUE first to generate the results."))
  }
  
  # Load the saved datasets
  message(paste0("Loading saved analysis results for ", ifelse(base_params$useRealData, "real data", "dummy data"), "..."))
  all_comparison_fitness_data <- readRDS(here::here(paste0('output/', output_subdir, '/all_comparison_fitness_data.rds')))
  all_comparison_times_data <- readRDS(here::here(paste0('output/', output_subdir, '/all_comparison_times_data.rds')))
  all_comparison_data_for_3d_plot <- readRDS(here::here(paste0('output/', output_subdir, '/all_comparison_data_for_3d_plot.rds')))
  
  # Load exhaustive search results if available
  exhaustive_file <- here::here(paste0('output/', output_subdir, '/exhaustive_search_global_optima.rds'))
  if (file.exists(exhaustive_file)) {
    exhaustive_global_optima <- readRDS(exhaustive_file)
    message(sprintf("  - Exhaustive search results: %d subgroup sizes", nrow(exhaustive_global_optima)))
    
    # Check if enhanced comparison data with gaps exists
    gaps_file <- here::here(paste0('output/', output_subdir, '/all_comparison_fitness_data_with_gaps.rds'))
    if (file.exists(gaps_file)) {
      all_comparison_fitness_data_with_gaps <- readRDS(gaps_file)
      message(sprintf("  - Comparison data with optimality gaps: %d rows", nrow(all_comparison_fitness_data_with_gaps)))
    }
  } else {
    message("  - No exhaustive search results found")
  }
  
  message("Successfully loaded:")
  message(sprintf("  - Fitness data: %d rows", nrow(all_comparison_fitness_data)))
  message(sprintf("  - Times data: %d rows", nrow(all_comparison_times_data))) 
  message(sprintf("  - 3D plot data: %d rows", nrow(all_comparison_data_for_3d_plot)))
  
  # Skip to plotting section
  message("\nSkipping analysis - jumping to visualization...")
  
} else {
  message("\n==============================================")
  message("          RUNNING ANALYSIS MODE              ")
  message("==============================================\n")
  
  # --- Generate Atomic Rules (needed for Greedy as well as GA Orbit) ---
  # This needs to be done once globally as it defines the search space for all methods.
  atomic_rules_info_global <- generate_all_atomic_rules(de_data_global, base_params$numeric_metrics, base_params$categorical_metrics)
  atomic_rules_list_global <- atomic_rules_info_global$atomic_rules_list
  N_ATOMIC_RULES_BO_global <- atomic_rules_info_global$n_atomic_rules_total
  
  # --- Benchmark Progress Tracking ---
  if (base_params$benchmark_mode) {
    total_runs <- length(bo_param_sets) * 2 + length(ga_param_sets) + length(ga_orbit_param_sets) + length(greedy_param_sets)
    message(sprintf("BENCHMARK MODE: Running %d total parameter combinations across 5 methods\n", total_runs))
  }
  
  # --- Run Bayesian Optimization (No Orbits) Analyses for different parameter sets ---
  # This uses the new bo_method.R script
  bo_param_count <- 0
  bo_total_params <- length(bo_param_sets)
  for (param_set_name in names(bo_param_sets)) {
    bo_param_count <- bo_param_count + 1
    
    if (base_params$benchmark_mode) {
      message(sprintf("[%d/%d] BO (No Orbits): %s", bo_param_count, bo_total_params, param_set_name))
    } else {
      message(sprintf("\n=============================================="))
      message(sprintf("  Running Bayesian Optimization (No Orbits): %s", param_set_name))
      message(sprintf("==============================================\n"))
    }
    
    current_bo_params <- c(base_params, bo_param_sets[[param_set_name]])
    
    tic(paste("Bayesian Optimization (No Orbits) Run -", param_set_name))
    results_bo_no_orbits <- run_bayesian_optimization_analysis( # This calls the function from the new bo_method.R
      params = current_bo_params,
      df = df_global,
      dfThisProtein = dfThisProtein_global,
      hv_data = hv_data_global,
      de_data = de_data_global,
      hv_aval_mean = hv_aval_mean_global,
      cl_outer = cl
    )
    bo_no_orbits_time <- toc()
    
    # Collect fitness data
    if (!is.null(results_bo_no_orbits$all_solutions_details) && nrow(results_bo_no_orbits$all_solutions_details) > 0) {
      current_fitness_data <- results_bo_no_orbits$all_solutions_details %>%
        select(MinSubgroupSize, OptimalFitness) %>%
        distinct() %>%
        mutate(Method = "Bayesian Optimization (No Orbits)", ParamSet = param_set_name)
      all_comparison_fitness_data <- bind_rows(all_comparison_fitness_data, current_fitness_data)
    } else {
      warning(sprintf("No valid fitness data for BO (No Orbits) run: %s", param_set_name))
    }
    
    # Collect execution time
    all_comparison_times_data <- bind_rows(all_comparison_times_data, data.frame(
      Method = "Bayesian Optimization (No Orbits)",
      ParamSet = param_set_name,
      Time_Seconds = bo_no_orbits_time$toc - bo_no_orbits_time$tic
    ))
    
    # Print detailed results only if not in benchmark mode
    if (!base_params$benchmark_mode) {
      print(results_bo_no_orbits$optimal_rules_table)
    }
  }
  
  # --- Run Bayesian Optimization (With Orbits) Analyses for different parameter sets ---
  # This uses the original bo_method.R, now renamed bo_orbit_method.R
  bo_orbit_param_count <- 0
  bo_orbit_total_params <- length(bo_param_sets) # Using the same param sets as the no-orbits version
  for (param_set_name in names(bo_param_sets)) { 
    bo_orbit_param_count <- bo_orbit_param_count + 1
    
    if (base_params$benchmark_mode) {
      message(sprintf("[%d/%d] BO (With Orbits): %s", bo_orbit_param_count, bo_orbit_total_params, param_set_name))
    } else {
      message(sprintf("\n=============================================="))
      message(sprintf("  Running Bayesian Optimization (With Orbits): %s", param_set_name))
      message(sprintf("==============================================\n"))
    }
    
    current_bo_params <- c(base_params, bo_param_sets[[param_set_name]])
    
    tic(paste("Bayesian Optimization (With Orbits) Run -", param_set_name))
    results_bo_with_orbits <- run_bayesian_optimization_analysis_orbit( # This calls the function from bo_orbit_method.R
      params = current_bo_params,
      df = df_global,
      dfThisProtein = dfThisProtein_global,
      hv_data = hv_data_global,
      de_data = de_data_global,
      hv_aval_mean = hv_aval_mean_global,
      cl_outer = cl
    )
    bo_with_orbits_time <- toc()
    
    # Collect fitness data
    if (!is.null(results_bo_with_orbits$all_solutions_details) && nrow(results_bo_with_orbits$all_solutions_details) > 0) {
      current_fitness_data <- results_bo_with_orbits$all_solutions_details %>%
        select(MinSubgroupSize, OptimalFitness) %>%
        distinct() %>%
        mutate(Method = "Bayesian Optimization (With Orbits)", ParamSet = param_set_name)
      all_comparison_fitness_data <- bind_rows(all_comparison_fitness_data, current_fitness_data)
    } else {
      warning(sprintf("No valid fitness data for BO (With Orbits) run: %s", param_set_name))
    }
    
    # Collect execution time
    all_comparison_times_data <- bind_rows(all_comparison_times_data, data.frame(
      Method = "Bayesian Optimization (With Orbits)",
      ParamSet = param_set_name,
      Time_Seconds = bo_with_orbits_time$toc - bo_with_orbits_time$tic
    ))
    
    # Print detailed results only if not in benchmark mode
    if (!base_params$benchmark_mode) {
      print(results_bo_with_orbits$optimal_rules_table)
    }
  }
  
  
  # --- Run Genetic Algorithm Analyses for different parameter sets ---
  ga_param_count <- 0
  ga_total_params <- length(ga_param_sets)
  for (param_set_name in names(ga_param_sets)) {
    ga_param_count <- ga_param_count + 1
    
    if (base_params$benchmark_mode) {
      message(sprintf("[%d/%d] GA: %s", ga_param_count, ga_total_params, param_set_name))
    } else {
      message(sprintf("\n=============================================="))
      message(sprintf("  Running Genetic Algorithm: %s", param_set_name))
      message(sprintf("==============================================\n"))
    }
    
    current_ga_params <- c(base_params, ga_param_sets[[param_set_name]])
    
    tic(paste("Genetic Algorithm Run -", param_set_name))
    results_ga <- run_genetic_algorithm_analysis(
      params = current_ga_params,
      df = df_global,
      dfThisProtein = dfThisProtein_global,
      hv_data = hv_data_global,
      de_data = de_data_global,
      hv_aval_mean = hv_aval_mean_global,
      cl_outer = cl
    )
    ga_time <- toc()
    
    # Collect fitness data
    if (!is.null(results_ga$all_solutions_details) && nrow(results_ga$all_solutions_details) > 0) {
      current_fitness_data <- results_ga$all_solutions_details %>%
        select(MinSubgroupSize, OptimalFitness) %>%
        distinct() %>%
        mutate(Method = "Genetic Algorithm", ParamSet = param_set_name)
      all_comparison_fitness_data <- bind_rows(all_comparison_fitness_data, current_fitness_data)
    } else {
      warning(sprintf("No valid fitness data for GA run: %s", param_set_name))
    }
    
    # Collect execution time
    all_comparison_times_data <- bind_rows(all_comparison_times_data, data.frame(
      Method = "Genetic Algorithm",
      ParamSet = param_set_name,
      Time_Seconds = ga_time$toc - ga_time$tic
    ))
    
    # Print detailed results only if not in benchmark mode
    if (!base_params$benchmark_mode) {
      print(results_ga$optimal_rules_table)
    }
  }
  
  # --- Run Orbit-Aware Genetic Algorithm Analyses for different parameter sets ---
  ga_orbit_param_count <- 0
  ga_orbit_total_params <- length(ga_orbit_param_sets)
  for (param_set_name in names(ga_orbit_param_sets)) {
    ga_orbit_param_count <- ga_orbit_param_count + 1
    
    if (base_params$benchmark_mode) {
      message(sprintf("[%d/%d] GA (Orbit): %s", ga_orbit_param_count, ga_orbit_total_params, param_set_name))
    } else {
      message(sprintf("\n=============================================="))
      message(sprintf("  Running Orbit-Aware Genetic Algorithm: %s", param_set_name))
      message(sprintf("==============================================\n"))
    }
    
    current_ga_orbit_params <- c(base_params, ga_orbit_param_sets[[param_set_name]])
    
    tic(paste("Orbit-Aware GA Run -", param_set_name))
    results_ga_orbit <- run_genetic_algorithm_orbit_analysis(
      params = current_ga_orbit_params,
      df = df_global,
      dfThisProtein = dfThisProtein_global,
      hv_data = hv_data_global,
      de_data = de_data_global,
      hv_aval_mean = hv_aval_mean_global,
      atomic_rules_map = atomic_rules_info_global$atomic_rules_map, # Pass map
      atomic_rules_list = atomic_rules_info_global$atomic_rules_list, # Pass list
      n_atomic_rules_total = atomic_rules_info_global$n_atomic_rules_total, # Pass total
      cl_outer = cl
    )
    ga_orbit_time <- toc()
    
    # Collect fitness data
    if (!is.null(results_ga_orbit$all_solutions_details) && nrow(results_ga_orbit$all_solutions_details) > 0) {
      current_fitness_data <- results_ga_orbit$all_solutions_details %>%
        select(MinSubgroupSize, OptimalFitness) %>%
        distinct() %>%
        mutate(Method = "Orbit-Aware GA", ParamSet = param_set_name)
      all_comparison_fitness_data <- bind_rows(all_comparison_fitness_data, current_fitness_data)
    } else {
      warning(sprintf("No valid fitness data for Orbit-Aware GA run: %s", param_set_name))
    }
    
    # Collect execution time
    all_comparison_times_data <- bind_rows(all_comparison_times_data, data.frame(
      Method = "Orbit-Aware GA",
      ParamSet = param_set_name,
      Time_Seconds = ga_orbit_time$toc - ga_orbit_time$tic
    ))
    
    # Print detailed results only if not in benchmark mode
    if (!base_params$benchmark_mode) {
      print(results_ga_orbit$optimal_rules_table)
    }
  }
  
  
  # --- Run Greedy Algorithm Analyses ---
  # Loop through greedy_param_sets (now multiple scenarios)
  greedy_param_count <- 0
  greedy_total_params <- length(greedy_param_sets)
  for (param_set_name in names(greedy_param_sets)) {
    greedy_param_count <- greedy_param_count + 1
    
    if (base_params$benchmark_mode) {
      message(sprintf("[%d/%d] Greedy: %s", greedy_param_count, greedy_total_params, param_set_name))
    } else {
      message(sprintf("\n=============================================="))
      message(sprintf("  Running Greedy Algorithm: %s", param_set_name))
      message(sprintf("==============================================\n"))
    }
    
    # Combine base_params with the specific greedy param set for this run
    current_greedy_params <- c(base_params, greedy_param_sets[[param_set_name]])
    
    tic(paste("Greedy Algorithm Run -", param_set_name))
    results_greedy <- run_greedy_analysis(
      params = current_greedy_params,
      df = df_global,
      dfThisProtein = dfThisProtein_global,
      hv_data = hv_data_global,
      de_data = de_data_global,
      hv_aval_mean = hv_aval_mean_global,
      atomic_rules_list = atomic_rules_list_global, # Pass the globally generated atomic rules
      cl_outer = cl
    )
    greedy_time <- toc()
    
    # Collect fitness data
    if (!is.null(results_greedy$all_solutions_details) && nrow(results_greedy$all_solutions_details) > 0) {
      current_fitness_data <- results_greedy$all_solutions_details %>%
        select(MinSubgroupSize, OptimalFitness) %>%
        distinct() %>%
        mutate(Method = "Greedy Algorithm", ParamSet = param_set_name)
      all_comparison_fitness_data <- bind_rows(all_comparison_fitness_data, current_fitness_data)
    } else {
      warning(sprintf("No valid fitness data for Greedy run: %s", param_set_name))
    }
    
    # Collect execution time
    all_comparison_times_data <- bind_rows(all_comparison_times_data, data.frame(
      Method = "Greedy Algorithm",
      ParamSet = param_set_name,
      Time_Seconds = greedy_time$toc - greedy_time$tic
    ))
    
    # Print detailed results only if not in benchmark mode  
    if (!base_params$benchmark_mode) {
      print(results_greedy$optimal_rules_table)
    }
  }
  
  # --- Prepare Data for Visualization (for runSim = TRUE) ---
  if (base_params$runSim == TRUE) {
    # Join time data to fitness data and ensure proper ordering for lines
    all_comparison_data_for_3d_plot <- all_comparison_fitness_data %>%
      left_join(all_comparison_times_data %>% select(Method, ParamSet, Time_Seconds), by = c("Method", "ParamSet")) %>%
      # Crucial for correct line drawing: order by MinSubgroupSize within each unique series
      arrange(Method, ParamSet, MinSubgroupSize)
    
  }
  
  
  # --- Save Analysis Results to Local Disk (only if runSim = TRUE) ---
  if (base_params$runSim == TRUE) {
    message("\nSaving analysis results to local disk...")
    
    # Ensure output directory exists (should already be created in exhaustive search section)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save the main comparison datasets in subdirectory
    saveRDS(all_comparison_fitness_data, file.path(output_dir, 'all_comparison_fitness_data.rds'))
    saveRDS(all_comparison_times_data, file.path(output_dir, 'all_comparison_times_data.rds'))
    saveRDS(all_comparison_data_for_3d_plot, file.path(output_dir, 'all_comparison_data_for_3d_plot.rds'))
    
    # --- Add Optimality Gap Analysis (if exhaustive search was performed) ---
    exhaustive_file <- file.path(output_dir, 'exhaustive_search_global_optima.rds')
    if (file.exists(exhaustive_file) && nrow(all_comparison_fitness_data) > 0) {
      message("\n=== OPTIMALITY GAP ANALYSIS ===")
      
      # Load exhaustive search results
      exhaustive_global_optima <- readRDS(exhaustive_file)
      
      # Check if we have valid exhaustive search results
      successful_exhaustive <- sum(!is.na(exhaustive_global_optima$TrueGlobalOptimum))
      
      if (successful_exhaustive > 0) {
        # Merge with exhaustive search results to calculate optimality gaps
        comparison_with_gaps <- all_comparison_fitness_data %>%
          left_join(exhaustive_global_optima %>% select(MinSubgroupSize, TrueGlobalOptimum), 
                    by = "MinSubgroupSize") %>%
          mutate(
            OptimalityGap = ifelse(!is.na(TrueGlobalOptimum), 
                                   (TrueGlobalOptimum - OptimalFitness) / TrueGlobalOptimum * 100, 
                                   NA),
            PerformanceRatio = ifelse(!is.na(TrueGlobalOptimum), 
                                      OptimalFitness / TrueGlobalOptimum, 
                                      NA)
          )
        
        # Save enhanced comparison data
        saveRDS(comparison_with_gaps, file.path(output_dir, 'all_comparison_fitness_data_with_gaps.rds'))
        
        # Display optimality gap summary
        gap_summary <- comparison_with_gaps %>%
          filter(!is.na(OptimalityGap)) %>%
          group_by(Method) %>%
          summarise(
            MeanGap = mean(OptimalityGap, na.rm = TRUE),
            MedianGap = median(OptimalityGap, na.rm = TRUE),
            MinGap = min(OptimalityGap, na.rm = TRUE),
            MaxGap = max(OptimalityGap, na.rm = TRUE),
            MeanPerformanceRatio = mean(PerformanceRatio, na.rm = TRUE),
            .groups = 'drop'
          ) %>%
          arrange(MeanGap)
        
        message("Optimality gap summary by method (lower is better):")
        for (i in 1:nrow(gap_summary)) {
          row <- gap_summary[i, ]
          message(sprintf("  %s: Mean gap = %.2f%%, Performance ratio = %.3f", 
                          row$Method, row$MeanGap, row$MeanPerformanceRatio))
        }
        
        message(sprintf("Enhanced comparison data saved to: %s", 
                        file.path(output_dir, 'all_comparison_fitness_data_with_gaps.rds')))
        
      } else {
        message("No valid exhaustive search results available for optimality gap analysis")
      }
    } else {
      message("Exhaustive search results not available - skipping optimality gap analysis")
    }
    
    message(paste0(
      "Saved results for ", ifelse(base_params$useRealData, "real data", "dummy data"), " in subdirectory '", output_subdir, "':\n",
      "  - output/", output_subdir, "/all_comparison_fitness_data.rds\n",
      "  - output/", output_subdir, "/all_comparison_times_data.rds\n", 
      "  - output/", output_subdir, "/all_comparison_data_for_3d_plot.rds"
    ))
    
  } else {
    message("\nSkipping save - results were loaded from disk.")
  }
  
  
} # End of runSim == TRUE condition


# --- Stability Testing (Reproducibility Analysis) ---
# This section is decoupled into two parts:
# 1. Data Generation: Generate stability data if needed
# 2. Analysis & Visualization: Analyze stability data (can be run independently)

# Check if stability data files exist
stability_data_file <- here::here(paste0('output/', output_subdir, '/stability_raw_results.rds'))


standardize_rule_format <- function(rule_string) {
  if (is.na(rule_string) || rule_string == "" || rule_string == "No rules extracted") {
    return(rule_string)
  }
  
  # Convert "VARIABLE == value" to "VARIABLE %in% {value}" 
  standardized <- gsub("([A-Za-z_][A-Za-z0-9_.]*) == ([^,]+)", "\\1 %in% {\\2}", rule_string)
  
  # Clean up spacing around curly brackets in %in% expressions
  # Remove extra spaces before and after opening/closing braces
  standardized <- gsub("%in%\\s*\\{\\s*", "%in% {", standardized)  # standardize opening
  standardized <- gsub("\\s*\\}", "}", standardized)              # remove space before closing
  
  # Clean up spaces around commas within braces - simpler approach
  # First normalize multiple spaces to single space
  standardized <- gsub("\\s+", " ", standardized)
  # Then fix comma spacing within braces: remove space before comma, ensure space after
  standardized <- gsub("\\s*,\\s*", ", ", standardized)
  # Remove any leading/trailing spaces within braces
  standardized <- gsub("\\{\\s+", "{", standardized)
  standardized <- gsub("\\s+\\}", "}", standardized)
  
  # Remove unnecessary quotes around numeric values within braces
  # This handles cases like MGLSRT %in% {'1'} -> MGLSRT %in% {1}
  standardized <- gsub("'([0-9.]+)'", "\\1", standardized)
  standardized <- gsub("\"([0-9.]+)\"", "\\1", standardized)
  
  return(standardized)
}

message("\n==============================================")
message("          STABILITY TESTING PHASE           ")
message("==============================================\n")

if (base_params$runSim == TRUE) {
  message("=== GENERATING STABILITY DATA ===")
  
  # Initialize stability results storage
  all_stability_results <- data.frame()
  
  # Select multiple representative parameter sets from each method for stability testing
  # (randomly sampling parameter sets from each method for more robust testing)
  set.seed(stability_params$stability_seed_base) # Ensure reproducible parameter selection
  n_param_sets_to_test <- min(stability_params$n_param_sets_per_method, length(bo_param_sets), length(ga_param_sets), 
                              length(ga_orbit_param_sets), length(greedy_param_sets))
  
  stability_param_sets <- list(
    "Bayesian Optimization (No Orbits)" = list(
      param_sets = sample(bo_param_sets, n_param_sets_to_test, replace = FALSE),
      run_function = run_bayesian_optimization_analysis
    ),
    "Bayesian Optimization (With Orbits)" = list(
      param_sets = sample(bo_param_sets, n_param_sets_to_test, replace = FALSE),
      run_function = run_bayesian_optimization_analysis_orbit
    ),
    "Genetic Algorithm" = list(
      param_sets = sample(ga_param_sets, n_param_sets_to_test, replace = FALSE),
      run_function = run_genetic_algorithm_analysis
    ),
    "Orbit-Aware GA" = list(
      param_sets = sample(ga_orbit_param_sets, n_param_sets_to_test, replace = FALSE),
      run_function = run_genetic_algorithm_orbit_analysis
    ),
    "Greedy Algorithm" = list(
      param_sets = sample(greedy_param_sets, n_param_sets_to_test, replace = FALSE),
      run_function = run_greedy_analysis
    )
  )
  
  message(sprintf("Testing stability with %d randomly selected parameter sets per method", n_param_sets_to_test))
  
  # Loop through each method for stability testing
  for (method_name in names(stability_param_sets)) {
    message(sprintf("\nTesting stability for: %s", method_name))
    
    method_info <- stability_param_sets[[method_name]]
    
    # Loop through the selected parameter sets for this method
    for (param_set_idx in 1:length(method_info$param_sets)) {
      current_param_set <- method_info$param_sets[[param_set_idx]]
      param_set_name <- names(method_info$param_sets)[param_set_idx]
      
      # Loop through selected minimum subgroup sizes
      for (min_size in stability_params$stability_min_subgroup_sizes) {
        message(sprintf("  %s, size %d: Running %d stability runs...", param_set_name, min_size, stability_params$n_stability_runs), appendLF = FALSE)
        
        # Store results for this method-paramset-size combination
        size_results <- data.frame()
        
        # NOTE: Currently running stability runs sequentially within each parameter set
        # The parallel cluster (cl) is available and used within individual method runs
        # For better parallelization, consider using foreach() here in future versions
        # Run the method multiple times
        for (run_idx in 1:stability_params$n_stability_runs) {
          # Set seed for reproducible randomness tracking
          set.seed(stability_params$stability_seed_base + param_set_idx * 1000 + min_size * 100 + run_idx)
          
          # Create parameters for this specific run (only test the selected min_subgroup_size)
          stability_run_params <- c(base_params, current_param_set)
          stability_run_params$min_subgroup_sizes_to_test <- min_size
          stability_run_params$verbose_analysis_output <- FALSE # Suppress verbose output
          stability_run_params$show_rules <- FALSE # Suppress rule display
          
          # Suppress method-specific verbose output and warnings
          suppressMessages(suppressWarnings({
            tryCatch({
              # Run the method
              if (method_name == "Orbit-Aware GA") {
                # Special handling for Orbit-Aware GA which needs atomic rules
                result <- method_info$run_function(
                  params = stability_run_params,
                  df = df_global,
                  dfThisProtein = dfThisProtein_global,
                  hv_data = hv_data_global,
                  de_data = de_data_global,
                  hv_aval_mean = hv_aval_mean_global,
                  atomic_rules_map = atomic_rules_info_global$atomic_rules_map,
                  atomic_rules_list = atomic_rules_info_global$atomic_rules_list,
                  n_atomic_rules_total = atomic_rules_info_global$n_atomic_rules_total,
                  cl_outer = cl
                )
              } else if (method_name == "Greedy Algorithm") {
                # Special handling for Greedy Algorithm which needs atomic rules
                result <- method_info$run_function(
                  params = stability_run_params,
                  df = df_global,
                  dfThisProtein = dfThisProtein_global,
                  hv_data = hv_data_global,
                  de_data = de_data_global,
                  hv_aval_mean = hv_aval_mean_global,
                  atomic_rules_list = atomic_rules_list_global,
                  cl_outer = cl
                )
              } else {
                # Standard handling for BO methods and GA
                result <- method_info$run_function(
                  params = stability_run_params,
                  df = df_global,
                  dfThisProtein = dfThisProtein_global,
                  hv_data = hv_data_global,
                  de_data = de_data_global,
                  hv_aval_mean = hv_aval_mean_global,
                  cl_outer = cl
                )
              }
              
              # Extract fitness value and solution details for this run
              if (!is.null(result$all_solutions_details) && nrow(result$all_solutions_details) > 0) {
                fitness_value <- result$all_solutions_details$OptimalFitness[1]
                
                # Extract the optimal solution rules with standardized formatting
                optimal_solution_rules <- ""
                if (!is.null(result$optimal_rules_table) && nrow(result$optimal_rules_table) > 0) {
                  # Method 1: Extract from optimal_rules_table (BO and Greedy methods)
                  rules_for_size <- result$optimal_rules_table %>%
                    filter(MinSubgroupSize == min_size)
                  
                  if (nrow(rules_for_size) > 0) {
                    # Check available columns and use appropriate rule column
                    if ("rules" %in% colnames(rules_for_size)) {
                      # Use the 'rules' column (lowercase) - this is what we found in the result structure
                      raw_rules <- rules_for_size$rules[1] # Take the first (and likely only) rule set
                      # Standardize format: convert any %in% to == format for consistency
                      optimal_solution_rules <- standardize_rule_format(raw_rules)
                    } else if ("Rules" %in% colnames(rules_for_size)) {
                      # Fallback to 'Rules' column (uppercase)
                      rules_list <- sort(rules_for_size$Rules)
                      raw_rules <- paste(rules_list, collapse = " | ")
                      optimal_solution_rules <- standardize_rule_format(raw_rules)
                    } else if ("RulesToDisplay" %in% colnames(rules_for_size)) {
                      # Another potential column name
                      rules_list <- sort(rules_for_size$RulesToDisplay)
                      raw_rules <- paste(rules_list, collapse = " | ")
                      optimal_solution_rules <- standardize_rule_format(raw_rules)
                    } else {
                      optimal_solution_rules <- "No recognized rule column found"
                    }
                  } else {
                    optimal_solution_rules <- "No rules found"
                  }
                } else if (!is.null(result$all_solutions_details) && nrow(result$all_solutions_details) > 0) {
                  # Method 2: Extract from all_solutions_details (GA methods)
                  # Filter for the specific subgroup size
                  solutions_for_size <- result$all_solutions_details %>%
                    filter(MinSubgroupSize == min_size)
                  
                  if (nrow(solutions_for_size) > 0) {
                    # Check if we have rule components (Metric, Operator, Value_or_Levels)
                    if (all(c("Metric", "Operator", "Value_or_Levels") %in% colnames(solutions_for_size))) {
                      # Build rules from individual components with standardized formatting
                      rule_components <- solutions_for_size %>%
                        select(Metric, Operator, Value_or_Levels) %>%
                        distinct() %>%
                        mutate(
                          # Create standardized rule string for each component
                          rule_part = case_when(
                            Operator == "==" ~ paste(Metric, "==", Value_or_Levels),
                            Operator == "%in%" ~ {
                              # Standardize %in% format to be cleaner and more consistent
                              values_clean <- gsub("^\\s*|\\s*$", "", Value_or_Levels) # trim whitespace
                              paste(Metric, "%in% {", values_clean, "}")
                            },
                            TRUE ~ paste(Metric, Operator, Value_or_Levels)
                          )
                        ) %>%
                        arrange(Metric) # Sort for consistency
                      
                      # Combine all rule parts
                      if (nrow(rule_components) > 0) {
                        optimal_solution_rules <- paste(rule_components$rule_part, collapse = ", ")
                      }
                    } else {
                      optimal_solution_rules <- "No recognized rule components found"
                    }
                  } else {
                    optimal_solution_rules <- "No solutions found for this size"
                  }
                } else {
                  optimal_solution_rules <- "No solution details available"
                }
                
                # Store the result with both fitness and solution (always store, even if rules are empty)
                size_results <- rbind(size_results, data.frame(
                  Method = method_name,
                  ParamSetName = param_set_name,
                  MinSubgroupSize = min_size,
                  RunIndex = run_idx,
                  OptimalFitness = fitness_value,
                  OptimalRules = optimal_solution_rules,
                  stringsAsFactors = FALSE
                ))
              } else {
                # If no all_solutions_details, still create a row to see what's missing
                if (run_idx == 1) { # Only show this once per parameter set
                  message(sprintf("   No all_solutions_details found for %s, paramset %s", method_name, param_set_name))
                }
              }
              
            }, error = function(e) {
              # Report the actual error for debugging
              if (run_idx == 1) { # Only report on first error to avoid spam
                warning(sprintf("Errors occurred in stability runs for %s (paramset: %s, size %d): %s", 
                                method_name, param_set_name, min_size, e$message), call. = FALSE)
              }
            })
          }))
          
          # Show progress dots every 5 runs
          if (run_idx %% 5 == 0) {
            cat(".")
          }
        }
        
        message(" Done!")
        
        # Combine results
        all_stability_results <- rbind(all_stability_results, size_results)
      }
    }
  }
  
  # Standardize rule formats before saving
  message("\n=== STANDARDIZING RULE FORMATS ===")
  all_stability_results2 = all_stability_results
  
  # Apply standardization to all rules
  if ("OptimalRules" %in% colnames(all_stability_results)) {
    original_count <- sum(!is.na(all_stability_results$OptimalRules) & 
                            all_stability_results$OptimalRules != "" &
                            all_stability_results$OptimalRules != "No rules extracted")
    
    all_stability_results$OptimalRules <- sapply(all_stability_results$OptimalRules, standardize_rule_format)
    
    standardized_count <- sum(!is.na(all_stability_results$OptimalRules) & 
                                all_stability_results$OptimalRules != "" &
                                all_stability_results$OptimalRules != "No rules extracted")
    
    message(sprintf("Rules standardized: %d rules processed (%d non-empty)", nrow(all_stability_results), original_count))
  }
  
  # Save the standardized stability data
  message("\n=== SAVING STABILITY DATA ===")
  output_dir <- here::here(paste0('output/', output_subdir))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  saveRDS(all_stability_results, stability_data_file)
  message(paste0("Stability raw data saved: ", stability_data_file))
  message("All stability analysis and visualization should be done in benchmark_visualization.Rmd")
  
} else {
  message("=== LOADING SAVED STABILITY DATA ===")
  
  # Check if stability data files exist
  if (!file.exists(stability_data_file)) {
    stop(paste0("Saved stability results not found for ", ifelse(base_params$useRealData, "real data", "dummy data"), 
                ". Please run with runSim = TRUE first to generate the stability results."))
  }
  
  message(paste0("Stability data file exists for ", ifelse(base_params$useRealData, "real data", "dummy data")))
  message("All stability analysis and visualization should be done in benchmark_visualization.Rmd")
}

message("\n=== STABILITY DATA GENERATION COMPLETED ===\n")
message("Next steps:")
message("1. Use benchmark_visualization.Rmd for stability analysis and plots")
message("2. The raw stability data is available at: ", stability_data_file)

# Clean up parallel cluster
if (!is.null(cl)) {
  stopCluster(cl)
  message("Parallel cluster stopped.")
}

# --- END OF SCRIPT ---
# All further analysis should be done in benchmark_visualization.Rmd
