# Quotient Space Learning for Patient Subgroup Discovery when label is not fixed

A comprehensive R-based framework for discovering optimal patient subgroups using advanced optimization methods, including novel **quotient space (quotient-space-aware)** approaches. This project benchmarks multiple optimization algorithms for identifying clinically meaningful patient subgroups based on biomarker expression levels.

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Optimization Methods](#optimization-methods)
- [File Descriptions](#file-descriptions)
- [Output](#output)
- [Visualization](#visualization)
- [Advanced Usage](#advanced-usage)
- [Citation](#citation)

---

## Overview

This framework addresses the challenge of identifying patient subgroups that exhibit significantly different biomarker expression levels compared to healthy volunteers. Here the main difficulty is that, the label of whether a patient belongs to the subgroup is not fixed (dependent on the rules for which the subgroup is defined). The key innovation is the identification of monoid structure and the incorporation of **quotient space learning** that leverage equivalence classes of clinically similar rule combinations to improve optimization efficiency and solution quality.

### Problem Statement

Given a dataset of patients with:
- Clinical metrics (categorical and/or numeric)
- Biomarker expression values
- Patient classification (Disease vs. Healthy Volunteer)

The goal is to find rule combinations (e.g., "DED_SEVERITY == 'Severe' AND DED_TYPE == 'Aqueous Deficient'") that define patient subgroups with maximally elevated biomarker levels relative to healthy volunteers.

### Fitness Function

The optimization maximizes:

```
Fitness = (Mean biomarker in subgroup) / (Mean biomarker in healthy volunteers)
```

Subject to a minimum subgroup size constraint.

---

## Key Features

- **Multiple Optimization Methods**: Bayesian Optimization, Genetic Algorithms, and Greedy Search
- **Quotient Space (quotient-Aware) Variants**: Novel implementations that cluster similar solutions to maintain diversity
- **Exhaustive Search Validation**: Ground truth computation for discrete (categorical-only) cases
- **Parallel Processing**: Efficient execution across multiple CPU cores
- **Comprehensive Benchmarking**: Systematic comparison across hyperparameter configurations
- **Interactive Visualization**: 3D performance landscapes and detailed comparison plots
- **Flexible Data Support**: Works with both real clinical data and simulated dummy data. However real clinical data is not available to public as it contains confidential data

---

## Project Structure

```
code_public_github/
├── run_methods.R              # Main entry point - orchestrates all analyses
├── common_fct.R               # Shared utility functions and libraries
├── bo_method.R                # Bayesian Optimization (standard)
├── bo_orbit_method.R          # Bayesian Optimization with quotient space awareness
├── ga_method.R                # Genetic Algorithm (standard)
├── ga_orbit_method.R          # Genetic Algorithm with quotient space awareness
├── greedy_method.R            # Greedy search algorithm
├── exhaustive_search.R        # Exhaustive search for ground truth validation
├── benchmark_visualization.Rmd # R Markdown for visualizing benchmark results
├── output/                    # Generated results (.rds files)
├── plots/                     # Generated plots and figures
├── tables/                    # Generated data tables
└── report/                    # Generated HTML reports
```

---

## Installation

### Prerequisites

Ensure you have R 4.3.0 installed along with the following packages:

```r
# Core packages
install.packages(c(
  "here",           # File path management
  "dplyr",          # Data manipulation
  "tidyr",          # Data tidying
  "ggplot2",        # Plotting
  "plotly",         # Interactive 3D plots
  "rlang",          # Tidy evaluation
  "foreach",        # Parallel loops
  "doSNOW",         # Parallel backend
  "progress",       # Progress bars
  "tictoc",         # Timing
  "knitr",          # Report generation
  "kableExtra",     # Table formatting
  "gridExtra",      # Grid arrangements
  "htmlwidgets",    # Interactive widgets
  "broom"           # Tidy model outputs
))

# Method-specific packages
install.packages(c(
  "GA",             # Genetic Algorithm
  "DiceKriging",    # Gaussian Process for BO
  "dbscan",         # DBSCAN clustering for orbit detection
  "stringdist"      # String distance calculations
))
```

### Optional Packages

For real clinical data integration (if using AWS S3):
```r
install.packages("aws.s3")
# Note: Additional internal packages (fixtheworld, rice, easyENTIM, arvupload) 
# may be required for real data loading
```


## Quick Start

### 1. Basic Execution with Dummy Data

```r
# Set working directory to project root
setwd("/path/to/code_public_github")

# Run the main analysis script
source("run_methods.R")
```

This will:
1. Generate synthetic patient data
2. Run all optimization methods with default parameters
3. Save results to `output/dummy_no_numeric_data/`
4. Display performance summaries

### 2. Generate Benchmark Visualization Report

```r
# Render the benchmark visualization report
rmarkdown::render(
  "benchmark_visualization.Rmd",
  params = list(
    useRealData = FALSE,
    useNumericMetrics = FALSE
  )
)
```

The HTML report will be saved to `report/benchmark_visualisation_dummy_no_numeric_data.html`.

---

## Configuration

### Main Parameters (in `run_methods.R`)

Edit the `base_params` list to customize the analysis:

```r
base_params <- list(
  # Data settings
  useRealData = FALSE,          # TRUE for real data, FALSE for dummy data
  study = "BP44617",            # Study identifier (for real data)
  
  # Metrics to use
  numeric_metrics = NULL,       # Set to c('NEISCORE', 'OCFSCORE', 'TFBUT') for numeric
  categorical_metrics = c('MGLSRT', 'DED_SEVERITY', 'DED_TYPE'),
  
  # Target variable
  variable_for_use = 'AVAL_BY_LENGTH',  # Biomarker variable
  thisProtein = 'ABC',                   # Protein identifier
  
  # Execution settings
  runSim = TRUE,                # TRUE to run simulations, FALSE to load saved results
  nCore = 4,                    # Number of parallel cores
  min_subgroup_sizes_to_test = 5:40,  # Range of subgroup sizes
  
  # Display settings
  benchmark_mode = TRUE,        # Suppress detailed output during benchmarking
  verbose_analysis_output = FALSE
)
```

### Hyperparameter Grids

The script automatically generates parameter combinations for systematic comparison:

```r
# Bayesian Optimization parameters
bo_initial_samples_vals <- c(100)
bo_iterations_vals <- c(100)
bo_xi_vals <- c(0.01, 0.1)

# Genetic Algorithm parameters
ga_pop_size_vals <- c(50, 60)
ga_max_iter_vals <- c(20)
ga_run_limit_vals <- c(10)

# Greedy parameters
greedy_max_rules_per_subgroup_vals <- seq(2, 30, 2)
```

---

## Optimization Methods

### Bayesian Optimization (BO)

**Files**: `bo_method.R`, `bo_orbit_method.R`

Uses Gaussian Process regression with a Hamming distance kernel to model the fitness landscape over the binary rule space.

| Variant | File | Description |
|---------|------|-------------|
| BO (Standard) | `bo_method.R` | Random search in binary space for acquisition optimization |
| BO (With Quotient Space) | `bo_orbit_method.R` | Clusters similar solutions and searches within equivalence classes |

### Genetic Algorithm (GA)

**Files**: `ga_method.R`, `ga_orbit_method.R`

Evolutionary approach with chromosome encoding for both numeric and categorical rules.

| Variant | File | Description |
|---------|------|-------------|
| GA (Standard) | `ga_method.R` | Uses R's `GA` package with standard operators |
| GA (With Quotient Space) | `ga_orbit_method.R` | Custom GA with DBSCAN-based orbit detection and niche elite preservation |

### Greedy Search

**File**: `greedy_method.R`

Simple iterative approach that adds one rule at a time to maximize fitness.

### Exhaustive Search (Ground Truth)

**File**: `exhaustive_search.R`

For discrete (categorical-only) cases, evaluates all possible rule combinations to find the true global optimum.

---

## File Descriptions

### Core Files

| File | Purpose |
|------|---------|
| `run_methods.R` | **Main entry point**. Loads data, configures parameters, runs all optimization methods, and saves results. |
| `common_fct.R` | Shared utilities: library loading, AWS setup, analysis result generation, plotting functions. |

### Method Implementations

| File | Method | Quotient Space |
|------|--------|----------------|
| `bo_method.R` | Bayesian Optimization | No |
| `bo_orbit_method.R` | Bayesian Optimization | Yes (cluster-based search) |
| `ga_method.R` | Genetic Algorithm | No |
| `ga_orbit_method.R` | Genetic Algorithm | Yes (niche elite preservation) |
| `greedy_method.R` | Greedy Search | No |
| `exhaustive_search.R` | Exhaustive Search | N/A (ground truth) |

### Visualization

| File | Purpose |
|------|---------|
| `benchmark_visualization.Rmd` | R Markdown document generating comprehensive benchmark reports with interactive plots |

---

## Output

### Directory Structure

Results are organized by data type and metric configuration:

```
output/
├── dummy_no_numeric_data/       # Dummy data, categorical metrics only
│   ├── all_comparison_fitness_data.rds
│   ├── all_comparison_times_data.rds
│   ├── all_comparison_data_for_3d_plot.rds
│   └── exhaustive_search_global_optima.rds
├── dummy_with_numeric_data/     # Dummy data, with numeric metrics
└── real_no_numeric_data/        # Real data, categorical metrics only
```

### Output Files

| File | Contents |
|------|----------|
| `all_comparison_fitness_data.rds` | Fitness values for all methods, parameter sets, and subgroup sizes |
| `all_comparison_times_data.rds` | Execution times for each method and parameter set |
| `all_comparison_data_for_3d_plot.rds` | Combined data for 3D visualization |
| `exhaustive_search_global_optima.rds` | True global optima from exhaustive search (discrete cases) |

---

## Visualization

### Running the Visualization Report

```r
# For dummy data without numeric metrics
rmarkdown::render(
  "benchmark_visualization.Rmd",
  params = list(useRealData = FALSE, useNumericMetrics = FALSE)
)

# For dummy data with numeric metrics
rmarkdown::render(
  "benchmark_visualization.Rmd",
  params = list(useRealData = FALSE, useNumericMetrics = TRUE)
)
```

### Generated Visualizations

1. **3D Performance Landscape**: Interactive plot showing fitness vs. subgroup size vs. execution time
2. **Method Comparison Plots**: GA vs GA (Quotient Space), BO vs BO (Quotient Space)
3. **Optimality Gap Analysis**: Distance from true global optimum (when exhaustive search is available)
4. **Performance Summary Tables**: Mean fitness, execution time, and consistency metrics

---
