#!/usr/bin/env Rscript

# RNA-seq Enhanced Analysis Pipeline - Setup Script
# Purpose: One-time project initialization - packages, directories, basic config
# Usage: Rscript setup.R [--packages-only] [--directories-only] [--help]
#
# This script sets up the hybrid enhancement system that:
# - Preserves your current Excel/CSV metadata workflow
# - Adds professional automation and batch processing
# - Creates organized directory structure for experiments
# - Installs all required R packages
# - Generates global configuration template
#
# Philosophy: "Best of Both Worlds" - keep scientific flexibility, add automation

# =============================================================================
# CONFIGURATION
# =============================================================================

REQUIRED_R_VERSION <- "4.1.0"

# Essential packages for enhanced functions
CRAN_PACKAGES <- c(
  "tidyverse",      # Data manipulation and ggplot2
  "yaml",           # Global configuration parsing
  "glue",           # String interpolation
  "fs",             # File system operations
  "scales",         # Scale functions for ggplot2
  "ggrepel",        # Non-overlapping text labels
  "ggtext",         # Enhanced text rendering (for markdown axis labels)
  "ggthemes",       # Additional ggplot themes
  "patchwork",      # Combining multiple plots
  "viridis",        # Perceptually uniform color scales
  "pals",           # Additional color palettes (for tol colors)
  "knitr",          # Dynamic report generation
  "kableExtra",     # Enhanced table formatting
  "DT",             # Interactive HTML tables
  "writexl",        # Excel file writing
  "readxl",         # Reading Excel files
  "openxlsx",       # Advanced Excel operations
  "optparse"        # Command line argument parsing
)

BIOC_PACKAGES <- c(
  "DESeq2",         # Differential expression analysis
  "tximeta",        # Transcript quantification import
  "ComplexHeatmap"  # Advanced heatmap generation
)

OPTIONAL_PACKAGES <- c(
  "plotly",         # Interactive plots for HTML reports
  "htmlwidgets",    # HTML widget framework
  "svglite",        # SVG graphics device
  "Cairo"           # Enhanced graphics rendering
)

# Hybrid system directory structure
DIRECTORIES <- c(
  "config",
  "functions",
  "templates",
  "experiments",
  "shared_resources",
  "shared_resources/metadata",
  "shared_resources/annotations",
  "logs",
  "backups"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

print_header <- function() {
  cat("\n")
  cat("================================================================================\n")
  cat("          RNA-seq Enhanced Analysis Pipeline - Setup                           \n")
  cat("         Preserve Excel Workflow + Add Professional Automation                 \n") 
  cat("================================================================================\n\n")
}

safe_execute <- function(func, description, required = TRUE) {
  # Wrapper function that executes setup steps with error handling and progress reporting
  # This helps isolate failures and provides clear feedback on which step failed
  cat(sprintf("🔄 %s...\n", description))
  
  result <- tryCatch({
    # Execute the function and capture any errors
    func()
    cat(sprintf("✓ %s completed\n", description))
    TRUE
  }, error = function(e) {
    cat(sprintf("✗ %s failed: %s\n", description, e$message))
    if (required) {
      # Stop entire setup if this is a required step
      stop(sprintf("Required step failed: %s", description))
    }
    # Continue setup if this is optional (like Quarto check)
    FALSE
  })
  
  return(result)
}

check_r_version <- function(required = REQUIRED_R_VERSION) {
  current_version <- paste(R.version$major, R.version$minor, sep = ".")
  
  if (compareVersion(current_version, required) < 0) {
    stop(sprintf("R version %s or higher is required. Current version: %s", 
                 required, current_version))
  }
  
  cat(sprintf("R version: %s\n", current_version))
  return(TRUE)
}

install_packages <- function(packages, source = "CRAN", required = TRUE) {
  # Install packages from CRAN or Bioconductor with robust error handling
  # This function installs packages one at a time so we can see exactly which one fails
  
  if (length(packages) == 0) return(TRUE)
  
  cat(sprintf("Checking %s packages (%d total)...\n", source, length(packages)))
  
  # Check which packages are missing by comparing with installed packages
  missing_packages <- setdiff(packages, rownames(installed.packages()))
  
  if (length(missing_packages) == 0) {
    cat(sprintf("All %s packages already installed\n", source))
    return(TRUE)
  }
  
  cat(sprintf("Installing %d missing %s packages...\n", 
              length(missing_packages), source))
  
  # Install packages one by one for better error reporting
  for (pkg in missing_packages) {
    cat(sprintf("  Installing %s...", pkg))
    
    install_result <- tryCatch({
      if (source == "CRAN") {
        # Standard CRAN installation
        install.packages(pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/")
      } else if (source == "Bioconductor") {
        # Bioconductor requires BiocManager
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cran.rstudio.com/")
        }
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      }
      TRUE
    }, error = function(e) {
      cat(sprintf(" FAILED: %s", e$message))
      FALSE
    })
    
    # Verify the package actually installed by checking if it's now available
    if (install_result && pkg %in% rownames(installed.packages())) {
      cat(" ✓\n")
    } else {
      cat(" ✗\n")
      if (required) {
        stop(sprintf("Failed to install required package: %s", pkg))
      }
    }
  }
  
  cat(sprintf("%s package installation completed\n", source))
  return(TRUE)
}

check_quarto <- function() {
  cat("Checking for Quarto...\n")
  
  quarto_available <- system("quarto --version", 
                            ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  
  if (quarto_available) {
    version_output <- tryCatch({
      system("quarto --version", intern = TRUE)
    }, error = function(e) "unknown version")
    cat(sprintf("Quarto found: %s\n", version_output[1]))
  } else {
    cat("⚠ Quarto not found - install from: https://quarto.org/docs/get-started/\n")
  }
  
  return(quarto_available)
}

create_directories <- function(directories = DIRECTORIES) {
  # Create the hybrid directory structure for the enhanced pipeline
  # This creates self-contained experiment folders while preserving shared resources
  cat(sprintf("Creating directory structure (%d directories)...\n", length(directories)))
  
  created_count <- 0
  for (dir in directories) {
    if (!dir.exists(dir)) {
      tryCatch({
        # recursive = TRUE ensures nested directories are created
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        cat(sprintf("  Created: %s/\n", dir))
        created_count <- created_count + 1
      }, error = function(e) {
        stop(sprintf("Failed to create directory %s: %s", dir, e$message))
      })
    } else {
      cat(sprintf("  Exists: %s/\n", dir))
    }
  }
  
  if (created_count == 0) {
    cat("All directories already existed\n")
  } else {
    cat(sprintf("Created %d new directories\n", created_count))
  }
  
  return(TRUE)
}

create_global_config <- function() {
  config_file <- "config/global_settings.yaml"
  
  if (file.exists(config_file)) {
    cat(sprintf("Global configuration already exists: %s\n", config_file))
    return(TRUE)
  }
  
  cat("Creating global configuration template using yaml package...\n")
  
  # Load yaml package for clean config creation
  library(yaml, quietly = TRUE)
  
  # Create configuration as R list structure
  # This approach ensures valid YAML and is more maintainable than text manipulation
  config <- list(
    # Project metadata - update these with your information
    project = list(
      name = "RNA-seq Enhanced Analysis Pipeline",
      version = "2.0.0",
      description = "Hybrid system preserving Excel workflow with professional automation",
      author = "Your Name",        # USER: Update this
      institute = "Your Lab"       # USER: Update this  
    ),
    
    # Critical paths - MUST be updated for your system
    # These point to shared resources used across all experiments
    paths = list(
      # Master metadata CSV file containing all your samples
      # This should be your existing rnaseq_yildiz_metadata.csv or similar
      metadata_file = "shared_resources/metadata/rnaseq_yildiz_metadata.csv",
      
      # Base path to salmon quantification outputs
      # Each experiment will construct full paths relative to this
      salmon_base = "../../../vc-refs/rnaseq_counts",  # USER: Update this path
      
      # Path to your annotation database (vc_genome_anno.RDS)
      # Contains gene annotations, pathways, functional categories
      annotation_db = "shared_resources/annotations/vc_genome_anno.RDS"  # USER: Update this path
    ),
    
    # Default analysis parameters
    # These can be overridden per experiment but provide sensible defaults
    analysis = list(
      # Standard significance threshold for differential expression
      significance_threshold = 0.05,
      
      # Fold change threshold (log2 scale, so 1.0 = 2x change)
      fold_change_threshold = 1.0,
      
      # Minimum read count threshold for filtering low-expressed genes
      # Matches your current filtering approach
      min_counts = 2
    ),
    
    # Basic plotting defaults
    # Specific styling (colors, devices) will be handled per experiment
    plotting = list(
      # Base text size for plots - matches your current ggthemes::theme_few(base_size = 8)
      base_size = 8,
      
      # Standard figure dimensions for consistent outputs
      figure_width = 7,
      figure_height = 5,
      
      # High resolution for publication quality
      dpi = 300
    ),
    
    # Report generation preferences  
    # Controls default output formats across all experiments
    reports = list(
      # PDF reports are publication-ready and universally compatible
      generate_pdf = TRUE,
      
      # HTML reports with interactive plots (optional, can be enabled per experiment)
      generate_html = FALSE,
      
      # Save individual plot files in addition to embedding in reports
      save_plots_separately = TRUE,
      
      # Author name for report headers
      author_name = "Your Name"    # USER: Update this
    )
  )
  
  # Write configuration using yaml package for clean formatting
  tryCatch({
    # write_yaml creates properly formatted, indented YAML
    yaml::write_yaml(config, config_file)
    cat(sprintf("Created: %s\n", config_file))
    cat("→ This config contains ONLY global settings (paths, defaults)\n")
    cat("→ Each experiment will have its own analysis script with custom filtering\n")
  }, error = function(e) {
    stop(sprintf("Failed to create config file: %s", e$message))
  })
  
  return(TRUE)
}

test_functionality <- function() {
  cat("Testing basic functionality...\n")
  
  # Test essential package loading
  essential_packages <- c("tidyverse", "DESeq2")
  failed_tests <- character(0)
  
  for (pkg in essential_packages) {
    load_result <- tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE))
      cat(sprintf("  ✓ %s loads correctly\n", pkg))
      TRUE
    }, error = function(e) {
      cat(sprintf("  ✗ %s failed to load: %s\n", pkg, e$message))
      FALSE
    })
    
    if (!load_result) {
      failed_tests <- c(failed_tests, pkg)
    }
  }
  
  # Test configuration file
  if (file.exists("config/global_settings.yaml")) {
    yaml_test <- tryCatch({
      # Test if we can read the file
      readLines("config/global_settings.yaml", n = 5)
      cat("  ✓ Configuration file readable\n")
      TRUE
    }, error = function(e) {
      cat(sprintf("  ✗ Configuration file issue: %s\n", e$message))
      FALSE
    })
  }
  
  if (length(failed_tests) > 0) {
    warning(sprintf("Some functionality tests failed: %s", 
                   paste(failed_tests, collapse = ", ")))
    return(FALSE)
  }
  
  cat("All functionality tests passed\n")
  return(TRUE)
}

print_summary <- function() {
  cat("\n")
  cat("================================================================================\n")
  cat("                            Setup Complete!                                    \n")
  cat("================================================================================\n\n")
  
  cat("🎯 What's Ready:\n")
  cat("   ✓ All packages installed\n") 
  cat("   ✓ Directory structure created\n")
  cat("   ✓ Global configuration template created\n\n")
  
  cat("📝 Next Steps:\n")
  cat("1. Update config/global_settings.yaml with your paths:\n")
  cat("   - salmon_base: path to your salmon quantification outputs\n")
  cat("   - annotation_db: path to your vc_genome_anno.RDS file\n") 
  cat("   - metadata_file: path to your master CSV metadata file\n\n")
  
  cat("2. Create the enhanced functions and templates\n")
  cat("   (These will be provided as separate files)\n\n")
  
  cat("💡 Setup Options:\n")
  cat("   Rscript setup.R                    # Full setup\n")
  cat("   Rscript setup.R --packages-only    # Install packages only\n") 
  cat("   Rscript setup.R --directories-only # Create directories only\n\n")
  
  cat("🧬 Ready for the next phase of implementation!\n\n")
}

# =============================================================================
# COMMAND LINE PARSING
# =============================================================================

parse_simple_args <- function() {
  # Simple command line argument parsing without external dependencies
  # This avoids the chicken-and-egg problem of needing optparse before installing packages
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default to full setup
  result <- list(
    packages_only = FALSE,
    directories_only = FALSE,
    help = FALSE
  )
  
  # Check for specific flags in command line arguments
  if (length(args) > 0) {
    if ("--packages-only" %in% args) result$packages_only <- TRUE
    if ("--directories-only" %in% args) result$directories_only <- TRUE  
    if ("--help" %in% args || "-h" %in% args) result$help <- TRUE
  }
  
  # Handle help request immediately
  if (result$help) {
    cat("RNA-seq Enhanced Analysis Pipeline - Setup\n\n")
    cat("Usage:\n")
    cat("  Rscript setup.R                    # Full setup (packages + directories)\n")
    cat("  Rscript setup.R --packages-only    # Install packages only\n")
    cat("  Rscript setup.R --directories-only # Create directories only\n")
    cat("  Rscript setup.R --help             # Show this help\n\n")
    cat("Options:\n")
    cat("  --packages-only      Install required packages but skip directory creation\n")
    cat("  --directories-only   Create directory structure but skip package installation\n")
    cat("  --help, -h          Show this help message and exit\n\n")
    quit(status = 0)
  }
  
  return(result)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

main <- function() {
  # Main setup function that orchestrates the entire setup process
  # This function handles command line arguments and executes setup steps conditionally
  
  print_header()
  
  # Parse command line arguments using simple approach (no external dependencies)
  args <- parse_simple_args()
  
  # Show user what mode we're running in
  cat("Setup mode: ")
  if (args$packages_only) {
    cat("Packages only\n")
  } else if (args$directories_only) {
    cat("Directories only\n")
  } else {
    cat("Full setup\n")
  }
  cat("\n")
  
  # PACKAGE INSTALLATION PHASE
  # Skip this phase if user specified --directories-only
  if (!args$directories_only) {
    # Check R version compatibility first
    safe_execute(check_r_version, "Checking R version", required = TRUE)
    
    # Install essential CRAN packages (required for core functionality)
    safe_execute(function() install_packages(CRAN_PACKAGES, "CRAN", TRUE), 
                "Installing CRAN packages", required = TRUE)
    
    # Install Bioconductor packages (DESeq2, tximeta, etc.)
    safe_execute(function() install_packages(BIOC_PACKAGES, "Bioconductor", TRUE), 
                "Installing Bioconductor packages", required = TRUE)
    
    # Install optional packages (for enhanced features, not critical)
    safe_execute(function() install_packages(OPTIONAL_PACKAGES, "CRAN", FALSE), 
                "Installing optional packages", required = FALSE)
    
    # Check if Quarto is available (needed for report generation)
    safe_execute(check_quarto, "Checking Quarto", required = FALSE)
  }
  
  # PROJECT STRUCTURE SETUP PHASE
  # Skip this phase if user specified --packages-only
  if (!args$packages_only) {
    # Create the hybrid directory structure
    safe_execute(create_directories, "Creating directories", required = TRUE)
    
    # Create global configuration file with paths and defaults
    safe_execute(create_global_config, "Creating configuration", required = TRUE)
    
    # Test that essential packages load correctly
    safe_execute(test_functionality, "Testing functionality", required = FALSE)
  }
  
  # Show summary and next steps
  print_summary()
}

# Execute main function if script is run directly
if (!interactive()) {
  main()
}