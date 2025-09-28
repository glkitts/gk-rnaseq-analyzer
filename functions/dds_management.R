#!/usr/bin/env Rscript

# DDS generation and management - Streamlined version
# Creation date: 2025-09-05

library(tidyverse)
library(DESeq2)
library(tximeta)
library(glue)
library(yaml)
library(cli)
library(here)

here::i_am("functions/dds_management.R")

# ============================================================================= #
# CORE DDS FUNCTIONS ----
# ============================================================================= #

#' Create DDS from salmon quantification with enhanced error handling
#'
#' This function creates a DESeq2 dataset from salmon quantification files,
#' runs the DESeq2 analysis, and filters low-count transcripts.
#'
#' @param coldata Data frame with sample metadata and file paths. Must contain
#'   'files' and 'names' columns.
#' @param design_formula Formula for DESeq2 design (default: ~group)
#' @param min_counts Minimum count threshold for filtering low-expressed transcripts (default: 2)
#' @param experiment_name Optional experiment name for logging and validation
#'
#' @return DESeqDataSet object with filtering stats and original coldata attached as attributes
#'
#' @examples
#' \dontrun{
#' dds <- create_dds_from_salmon(coldata, design_formula = ~condition, min_counts = 2)
#' }
create_dds_from_salmon <- function(coldata,
                                   design_formula = ~ group,
                                   min_counts = 2,
                                   experiment_name = NULL) {
  
  # Basic validation
  required_cols <- c("files", "names")
  missing_cols <- setdiff(required_cols, colnames(coldata))
  if (length(missing_cols) > 0) {
    cli_abort("Missing required columns in coldata: {.field {missing_cols}}")
  }
  if (nrow(coldata) == 0) {
    cli_abort("coldata cannot be empty")
  }
  
  cli_inform("Creating DDS for {experiment_name} with {nrow(coldata)} samples")
  
  # Create tximeta object
  gse <- tximeta(coldata, skipMeta = TRUE)
  
  # Create DESeq2 object and run analysis
  dds <- DESeqDataSet(gse, design = design_formula)
  dds <- DESeq(dds)
  
  # Filter low counts and attach metadata
  filter_result <- filter_low_counts(dds, min_counts = min_counts)
  
  # Attach metadata to DDS
  attr(filter_result$dds, "filtering_stats") <- filter_result$filtering_stats
  attr(filter_result$dds, "original_coldata") <- coldata
  
  return(filter_result$dds)
}

#' Filter low count transcripts with reporting
#'
#' Removes transcripts with low read counts and reports filtering statistics.
#'
#' @param dds DESeqDataSet object
#' @param min_counts Minimum total count threshold across all samples (default: 2)
#'
#' @return List containing:
#'   \itemize{
#'     \item dds: Filtered DESeqDataSet object
#'     \item filtering_stats: List with filtering statistics
#'   }
filter_low_counts <- function(dds, min_counts = 2) {
  nrow_before <- nrow(dds)
  keep <- rowSums(counts(dds)) > min_counts
  dds_filtered <- dds[keep, ]
  nrow_after <- nrow(dds_filtered)
  nrow_removed <- nrow_before - nrow_after
  
  cli_inform("Filtering: {nrow_removed} of {nrow_before} transcripts removed (≤{min_counts} counts)")
  cli_inform("{nrow_after} transcripts retained for analysis")
  
  filtering_stats <- list(
    transcripts_before_filtering = nrow_before,
    transcripts_after_filtering = nrow_after,
    transcripts_removed = nrow_removed,
    min_counts_threshold = min_counts
  )
  
  list(dds = dds_filtered, filtering_stats = filtering_stats)
}

#' Load DDS object from file
#'
#' @param experiment_name String identifier for the experiment
#'
#' @return DESeqDataSet object
load_dds <- function(experiment_name) {
  dds_path <- here("experiments", experiment_name, "outputs", "technical", "R", paste0(experiment_name, ".dds.RDS"))
  
  if (!file.exists(dds_path)) {
    stop(glue("DDS not found. Run: Rscript run_experiment.R {experiment_name} --dds-only"))
  }
  
  return(read_rds(dds_path))
}

#' Save DDS object with comprehensive metadata
#'
#' Saves the DESeq2 dataset object along with metadata about the analysis,
#' quality metrics, and sample information. Also saves coldata as a separate CSV file.
#'
#' @param dds DESeqDataSet object (must be processed with DESeq())
#' @param experiment_name String identifier for the experiment (used for file naming)
#' @param coldata Optional data frame with original sample metadata. If NULL,
#'   will use coldata attached to the DDS object.
#' @param config Optional list with analysis configuration parameters
#'
#' @return List containing paths to saved files
save_dds_with_metadata <- function(dds,
                                   experiment_name,
                                   coldata = NULL,
                                   config = NULL) {
  # Basic validation
  if (!is(dds, "DESeqDataSet")) {
    cli_abort("dds must be a DESeqDataSet object")
  }
  if (is.null(experiment_name) || experiment_name == "") {
    cli_abort("experiment_name is required and cannot be empty")
  }
  
  # Create output directories and file paths
  outputs_dir <- ensure_experiment_outputs(experiment_name)
  dds_path <- file.path(outputs_dir, "technical", "R", paste0(experiment_name, ".dds.RDS"))
  metadata_path <- file.path(outputs_dir, "technical/metadata", paste0(experiment_name, ".metadata.yaml"))
  coldata_path <- file.path(outputs_dir, "technical/metadata", paste0(experiment_name, ".coldata.csv"))
  
  # Extract coldata if not provided
  if (is.null(coldata)) {
    coldata <- attr(dds, "original_coldata")
  }
  
  # Build metadata
  metadata <- build_metadata(dds, experiment_name, coldata, config, 
                           list(coldata_path = coldata_path))
  
  # Save files
  if (!is.null(coldata)) {
    cli_inform("Saving coldata: {basename(coldata_path)}")
    write_csv(coldata, coldata_path)
  }
  
  cli_inform("Saving DDS: {basename(dds_path)}")
  saveRDS(dds, dds_path)
  
  cli_inform("Saving metadata: {basename(metadata_path)}")
  yaml::write_yaml(metadata, metadata_path)
  
  # Report success
  file_size_mb <- round(file.size(dds_path) / 1024^2, 1)
  cli_inform("DDS saved for {experiment_name} ({file_size_mb} MB) | Samples: {ncol(dds)} | Transcripts: {nrow(dds)}")
  
  # Return paths
  list(
    dds_path = dds_path,
    metadata_path = metadata_path,
    coldata_path = coldata_path
  )
}

# ============================================================================= #
# METADATA AND UTILITY FUNCTIONS ----
# ============================================================================= #

#' Extract sample group counts from coldata and design formula
#'
#' @param coldata Sample metadata
#' @param design_formula DESeq2 design formula
#'
#' @return List with sample counts per group, or NULL if extraction fails
extract_sample_group_counts <- function(coldata, design_formula) {
  if (is.null(coldata))
    return(NULL)
  
  design_factors <- all.vars(design_formula)
  if (length(design_factors) == 0)
    return(NULL)
  
  primary_factor <- design_factors[length(design_factors)]
  if (!primary_factor %in% colnames(coldata))
    return(NULL)
  
  coldata %>%
    group_by(!!sym(primary_factor)) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    {
      setNames(.$n_samples, .[[primary_factor]])
    } %>%
    as.list()
}

#' Build quality metrics from DDS object
#'
#' @param dds DESeqDataSet object
#'
#' @return List with quality metrics, or NULL if none available
build_quality_metrics <- function(dds) {
  metrics <- list()
  
  if (!is.null(sizeFactors(dds))) {
    metrics$size_factors_range <- round(range(sizeFactors(dds)), 3)
    metrics$size_factors_median <- round(median(sizeFactors(dds)), 3)
    metrics$library_sizes_range <- range(colSums(counts(dds)))
    metrics$library_sizes_median <- median(colSums(counts(dds)))
  }
  
  if (!is.null(dispersions(dds))) {
    metrics$dispersion_median <- round(median(dispersions(dds), na.rm = TRUE), 4)
    metrics$dispersion_range <- round(range(dispersions(dds), na.rm = TRUE), 4)
  }
  
  if (length(metrics) == 0)
    NULL
  else
    metrics
}

#' Build complete metadata object for experiment
#'
#' @param dds DESeqDataSet object
#' @param experiment_name String identifier
#' @param coldata Sample metadata
#' @param config Optional configuration
#' @param paths File paths list
#'
#' @return List with comprehensive experiment metadata
build_metadata <- function(dds, experiment_name, coldata, config, paths) {
  filtering_stats <- attr(dds, "filtering_stats")
  
  list(
    experiment = list(
      name = experiment_name,
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      deseq2_version = as.character(packageVersion("DESeq2")),
      tximeta_version = as.character(packageVersion("tximeta"))
    ),
    
    dataset = list(
      n_samples = ncol(dds),
      n_transcripts = nrow(dds),
      design_formula = as.character(design(dds))[2],
      coldata_file = if (!is.null(paths$coldata_path)) basename(paths$coldata_path) else NULL,
      sample_group_counts = extract_sample_group_counts(coldata, design(dds))
    ),
    
    filtering = filtering_stats,
    
    analysis_params = list(
      min_counts = if (!is.null(config)) config$analysis$min_counts else 2,
      significance_threshold = if (!is.null(config)) config$analysis$significance_threshold else 0.05,
      fold_change_threshold = if (!is.null(config)) config$analysis$fold_change_threshold else 1.0
    ),
    
    source_paths = if (!is.null(config)) {
      list(
        metadata_file = config$paths$metadata_file,
        salmon_base = config$paths$salmon_base,
        annotation_db = config$paths$annotation_db
      )
    } else NULL,
    
    quality_metrics = build_quality_metrics(dds)
  )
}

#' Ensure experiment output directories exist
#'
#' Creates the required directory structure for an experiment if it doesn't exist.
#' Sources the function from rnaseq_output_management.R to avoid code duplication.
#'
#' @param experiment_name String identifier for the experiment
#'
#' @return String path to the outputs directory
ensure_experiment_outputs <- function(experiment_name) {
  source(here::here("functions", "output_management.R"))
  ensure_experiment_outputs(experiment_name)
}