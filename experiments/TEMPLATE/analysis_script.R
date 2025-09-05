#!/usr/bin/env Rscript

# RNA-seq Analysis: experiment_template
# Creation date: CREATION_DATE

# =============================================================================
# SETUP & CONFIGURATION
# =============================================================================

library(tidyverse)
library(DESeq2)
library(tximeta)
library(yaml)
library(glue)
library(here)

# Experiment name - CHANGE THIS!!
experiment_name <- "experiment_template"

# Establish project root
here::i_am(paste0("experiments/", experiment_name, "/analysis_script.R"))

# Source enhanced functions
source(here("functions", "rnaseq_analysis.R"))
source(here("functions", "rnaseq_plotting.R"))
source(here("functions", "rnaseq_output_management.R"))

# =============================================================================
# EXPERIMENT CONFIGURATION
# =============================================================================

output_label <- experiment_name
cat(glue("Running experiment: {experiment_name}\n"))

# Load global config using here()
config <- read_yaml(here("config", "global_settings.yaml"))

# Experiment-specific overrides (customize as needed)
# config$analysis$min_counts <- 1
# config$analysis$significance_threshold <- 0.01
# config$plotting$base_size <- 10

# Determine what processing steps to run
pipeline_step <- Sys.getenv("PIPELINE_STEP", "full")  # full, dds-only, analysis-only, reports-only
cat(glue("Running {experiment_name} - Step: {pipeline_step}\n"))

# =============================================================================
# COLDATA FILTERING LOGIC
# =============================================================================
# Customize this filtering logic for your experiment
# This section can be run interactively for testing

# Load master metadata
coldata.in <- read.csv(config$paths$metadata_file)

coldata <- coldata.in %>%
  filter(experiment_group %in% c(
    "iron"
    # "aki_rvvA",
    # "lb_rvvA"
  )) %>%
  filter(!strain %in% c(
    "dryhB"
    # "dfur",
    # "drvvA"
  )) %>%
  # filter(!sample %in% c("ID_dvxrB_3", "IR_dvxrB_2")) %>%
  mutate(
    strain = fct_relevel(strain, "WT"),
    condition = fct_relevel(condition, "iron_replete"),
    group = fct_inorder(sampleLabel),
    names = sample,
    files = file.path(
      config$paths$salmon_base,
      experiment_folder,
      sample_FolderName,
      "salmon_quants",
      "quant.sf"
    )
  ) %>%
  select(sample, strain, condition, group, files, temperature, names)

# =============================================================================
# CONTRAST DEFINITIONS
# =============================================================================
# Define your comparisons - customize for your experiment

contrast_list <- list(
  c("group", "IR_drvvA", "IR_WT"),
  c("group", "ID_drvvA", "ID_WT"),
  c("group", "IR_dfur", "IR_WT"),
  c("group", "ID_dfur", "ID_WT"),
  c("group", "IR_dfur", "IR_WT"),
  c("group", "ID_dfur", "ID_WT")
  # Add more contrasts as needed
)

# =============================================================================
# STEP 1: DDS GENERATION
# =============================================================================

if (pipeline_step %in% c("full", "dds-only")) {
  cat("Generating DDS object...\n")
  
  # Validate files exist
  missing_files <- coldata$files[!file.exists(coldata$files)]
  if (length(missing_files) > 0) {
    cat("ERROR: Missing salmon files:\n")
    cat(paste(missing_files, collapse = "\n"))
    stop("Missing salmon quantification files. Check paths.")
  }
  
  cat(glue("Found {nrow(coldata)} samples for analysis\n"))
  
  # Create DDS object using enhanced functions
  dds <- create_dds_from_salmon(
    coldata = coldata, 
    design_formula = ~group,
    min_counts = config$analysis$min_counts,
    experiment_name = experiment_name
  )
  
  # Save DDS and metadata
  save_dds_with_metadata(dds, experiment_name, coldata, config)
  
  cat(glue("DDS saved for {experiment_name}\n"))
}

# =============================================================================
# STEP 2: DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

if (pipeline_step %in% c("full", "analysis-only")) {
  cat("Running differential expression analysis...\n")
  
  # Load DDS if not already in memory
  if (!exists("dds")) {
    dds <- load_dds_with_validation(experiment_name)
  }
  
  # Load annotation database
  annotations <- readRDS(config$paths$annotation_db)
  sub_anots <- prepare_annotations(annotations)
  
  # Run enhanced differential expression analysis
  # This replaces your current deseq_makeOutput function
  results <- enhanced_deseq_makeOutput(
    contrast_list = contrast_list,
    dds = dds,
    sub_anots = sub_anots,
    merge_anno = TRUE,
    experiment_name = experiment_name,
    config = config
  )
  
  # Save organized results
  save_experiment_results(results, experiment_name, config)
  
  cat(glue("Analysis results saved for {experiment_name}\n"))
}

# =============================================================================
# STEP 3: REPORT GENERATION  
# =============================================================================

if (pipeline_step %in% c("full", "analysis-only", "reports-only")) {
  cat("Generating reports...\n")
  
  # Load results if not in memory
  if (!exists("results")) {
    results <- load_experiment_results(experiment_name)
  }
  if (!exists("dds")) {
    dds <- load_dds_with_validation(experiment_name)
  }
  
  # Generate all comparison reports
  generate_comparison_reports(
    results = results,
    dds = dds, 
    experiment_name = experiment_name,
    config = config
  )
  
  # Generate experiment summary
  generate_experiment_summary_report(
    results = results,
    dds = dds,
    experiment_name = experiment_name, 
    config = config
  )
  
  # Generate experiment README
  generate_experiment_readme(experiment_name, results, config)
  
  cat(glue("Reports generated for {experiment_name}\n"))
}

# =============================================================================
# COMPLETION
# =============================================================================

cat(glue("{experiment_name} analysis complete!\n"))
cat(glue("Results in: experiments/{experiment_name}/outputs/\n"))

# Print summary stats if available
if (exists("results") && exists("dds")) {
  cat("\nEXPERIMENT SUMMARY:\n")
  cat(glue("   • Samples analyzed: {ncol(dds)}\n"))
  cat(glue("   • Transcripts: {nrow(dds)}\n"))
  cat(glue("   • Comparisons: {length(results$xlsx) - 1}\n"))  # -1 for sigHits_compiled
  cat(glue("   • Significant genes (any comparison): {nrow(results$xlsx$sigHits_compiled)}\n"))
}