#!/usr/bin/env Rscript

# DDS generation and management
# Creation date: 2025-09-04

library(tidyverse)
library(DESeq2)
library(tximeta)
library(glue)
library(yaml)
library(cli)
library(here)

here::i_am("functions/dds_management.R")

# ============================================================================= #
# VALIDATION HELPERS ----
# ============================================================================= #

#' Validate coldata for DDS creation
#' @param coldata Data frame with sample metadata
#' @noRd
validate_coldata <- function(coldata) {
  required_cols <- c("files", "names")
  missing_cols <- setdiff(required_cols, colnames(coldata))
  
  if (length(missing_cols) > 0) {
    cli_abort("Missing required columns in coldata: {.field {missing_cols}}")
  }
  
  if (nrow(coldata) == 0) {
    cli_abort("coldata cannot be empty")
  }
  
  invisible(coldata)
}

#' Validate DDS object for saving
#' @param dds DESeqDataSet object
#' @noRd
validate_dds_for_saving <- function(dds) {
  if (!is(dds, "DESeqDataSet")) {
    cli_abort("dds must be a DESeqDataSet object")
  }
  invisible(dds)
}

#' Validate experiment name
#' @param experiment_name String identifier
#' @noRd
validate_experiment_name <- function(experiment_name) {
  if (is.null(experiment_name) || experiment_name == "") {
    cli_abort("experiment_name is required and cannot be empty")
  }
  invisible(experiment_name)
}

# ============================================================================= #
# DDS CREATION HELPERS ----
# ============================================================================= #

#' Create tximeta object from coldata
#' @param coldata Validated coldata
#' @noRd
create_tximeta_object <- function(coldata) {
  # cli_inform("📊 Creating tximeta object from {nrow(coldata)} samples")
  tximeta(coldata, skipMeta = TRUE)
}

#' Create DESeqDataSet from tximeta object
#' @param gse tximeta object
#' @param design_formula Design formula
#' @noRd
create_deseq_object <- function(gse, design_formula) {
  # cli_inform("🧬 Creating DESeqDataSet object")
  DESeqDataSet(gse, design = design_formula)
}

#' Run DESeq2 analysis
#' @param dds DESeqDataSet object
#' @noRd
run_deseq_analysis <- function(dds) {
  # cli_inform("🔬 Running DESeq2 analysis")
  DESeq(dds)
}

#' Attach analysis metadata to DDS object
#' @param dds Processed DDS object
#' @param filtering_stats Filtering statistics
#' @param coldata Original coldata
#' @noRd
attach_dds_metadata <- function(dds, filtering_stats, coldata) {
  attr(dds, "filtering_stats") <- filtering_stats
  attr(dds, "original_coldata") <- coldata
  dds
}

# ============================================================================= #
# DDS CREATION AND MANAGEMENT ----
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
  # Validate inputs
  validate_coldata(coldata)
  
  cli_inform("🔬 Creating DDS for {nrow(coldata)} samples")
  
  # Create pipeline: tximeta -> DESeq2 -> filter -> attach metadata
  gse <- create_tximeta_object(coldata)
  dds <- create_deseq_object(gse, design_formula)
  dds <- run_deseq_analysis(dds)
  
  # Filter and get stats
  filter_result <- filter_low_counts(dds, min_counts = min_counts)
  
  # Attach metadata and return
  attach_dds_metadata(filter_result$dds, filter_result$filtering_stats, coldata)
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
  
  cli_inform(
    "🗂 Filtering: {nrow_removed} of {nrow_before} transcripts removed (≤{min_counts} counts)"
  )
  cli_inform("✅ {nrow_after} transcripts retained for analysis")
  
  filtering_stats <- list(
    transcripts_before_filtering = nrow_before,
    transcripts_after_filtering = nrow_after,
    transcripts_removed = nrow_removed,
    min_counts_threshold = min_counts
  )
  
  list(dds = dds_filtered, filtering_stats = filtering_stats)
}

# ============================================================================= #
# SAVE HELPERS ----
# ============================================================================= #

#' Create output file paths for experiment
#' @param experiment_name String identifier
#' @noRd
create_output_paths <- function(experiment_name) {
  outputs_dir <- ensure_experiment_outputs(experiment_name)
  
  list(
    outputs_dir = outputs_dir,
    dds_path = file.path(outputs_dir, paste0(experiment_name, ".dds.RDS")),
    metadata_path = file.path(outputs_dir, paste0(experiment_name, ".metadata.yaml")),
    coldata_path = file.path(outputs_dir, paste0(experiment_name, ".coldata.csv"))
  )
}

#' Extract coldata from DDS or provided argument
#' @param dds DESeqDataSet object
#' @param coldata Optional provided coldata
#' @noRd
extract_coldata <- function(dds, coldata = NULL) {
  if (is.null(coldata)) {
    attr(dds, "original_coldata")
  } else {
    coldata
  }
}

#' Extract sample group counts from coldata and design
#' @param coldata Sample metadata
#' @param design_formula DESeq2 design formula
#' @noRd
extract_sample_group_counts <- function(coldata, design_formula) {
  if (is.null(coldata))
    return(NULL)
  
  design_factors <- all.vars(design_formula)
  if (length(design_factors) == 0)
    return(NULL)
  
  primary_factor <- design_factors[length(design_factors)]
  if (!primary_factor %in% colnames(coldata))
    return(NULL)
  
  # Use original working logic
  coldata %>%
    group_by(!!sym(primary_factor)) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    {
      setNames(.$n_samples, .[[primary_factor]])
    } %>%
    as.list()
}

#' Build quality metrics from DDS
#' @param dds DESeqDataSet object
#' @noRd
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

#' Build complete metadata object
#' @param dds DESeqDataSet object
#' @param experiment_name String identifier
#' @param coldata Sample metadata
#' @param config Optional configuration
#' @param paths File paths
#' @noRd
build_metadata <- function(dds,
                           experiment_name,
                           coldata,
                           config,
                           paths) {
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
      coldata_file = basename(paths$coldata_path),
      sample_group_counts = extract_sample_group_counts(coldata, design(dds))
    ),
    
    filtering = filtering_stats,
    
    # analysis_params = list(
    #   min_counts = if (!is.null(config))
    #     config$analysis$min_counts
    #   else
    #     2,
    #   significance_threshold = if (!is.null(config))
    #     config$analysis$significance_threshold
    #   else
    #     0.05,
    #   fold_change_threshold = if (!is.null(config))
    #     config$analysis$fold_change_threshold
    #   else
    #     1.0
    # ),
    analysis_params = list(
      min_counts = if_else(!is.null(config),
        config$analysis$min_counts,
        2),
      significance_threshold = if_else(!is.null(config),
        config$analysis$significance_threshold,
        0.05),
      fold_change_threshold = if_else(!is.null(config),
        config$analysis$fold_change_threshold,
        1.0)
    ),
    
    source_paths = if (!is.null(config)) {
      list(
        metadata_file = config$paths$metadata_file,
        salmon_base = config$paths$salmon_base,
        annotation_db = config$paths$annotation_db
      )
    } else
      NULL,
    
    quality_metrics = build_quality_metrics(dds)
  )
}

#' Save DDS object to file
#' @param dds DESeqDataSet object
#' @param dds_path File path
#' @noRd
save_dds_object <- function(dds, dds_path) {
  cli_inform("💾 Saving DDS: {basename(dds_path)}")
  saveRDS(dds, dds_path)
}

#' Save coldata to CSV file
#' @param coldata Sample metadata
#' @param coldata_path File path
#' @noRd
save_coldata_file <- function(coldata, coldata_path) {
  if (!is.null(coldata)) {
    cli_inform("📊 Saving coldata: {basename(coldata_path)}")
    write_csv(coldata, coldata_path)
  }
}

#' Save metadata to YAML file
#' @param metadata Metadata object
#' @param metadata_path File path
#' @noRd
save_metadata_file <- function(metadata, metadata_path) {
  cli_inform("📋 Saving metadata: {basename(metadata_path)}")
  yaml::write_yaml(metadata, metadata_path)
}

#' Report successful save operation
#' @param experiment_name String identifier
#' @param dds DESeqDataSet object
#' @param dds_path Path to saved DDS file
#' @noRd
report_save_success <- function(experiment_name, dds, dds_path) {
  file_size_mb <- round(file.size(dds_path) / 1024^2, 1)
  cli_inform(
    "DDS saved for {experiment_name} ({file_size_mb} MB) | Samples: {ncol(dds)} | Transcripts: {nrow(dds)}"
  )
}

# ============================================================================= #
# SAVE AND LOAD DDS ----
# ============================================================================= #

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
#' @return List containing paths to saved files:
#'   \itemize{
#'     \item dds_path: Path to saved DDS RDS file
#'     \item metadata_path: Path to saved metadata YAML file
#'     \item coldata_path: Path to saved coldata CSV file
#'   }
#'
#' @details
#' Creates three files in the experiment's outputs directory:
#' \itemize{
#'   \item {experiment_name}_dds.RDS: The DESeq2 object
#'   \item {experiment_name}_metadata.yaml: Analysis metadata and quality metrics
#'   \item {experiment_name}_coldata.csv: Sample metadata
#' }
#'
#' @examples
#' \dontrun{
#' result <- save_dds_with_metadata(dds, "my_experiment", coldata, config)
#' }
save_dds_with_metadata <- function(dds,
                                   experiment_name,
                                   coldata = NULL,
                                   config = NULL) {
  # Validate inputs
  validate_dds_for_saving(dds)
  validate_experiment_name(experiment_name)
  
  # Setup paths and extract data
  paths <- create_output_paths(experiment_name)
  coldata <- extract_coldata(dds, coldata)
  
  # Build metadata
  metadata <- build_metadata(dds, experiment_name, coldata, config, paths)
  
  # Save all files
  save_coldata_file(coldata, paths$coldata_path)
  save_dds_object(dds, paths$dds_path)
  save_metadata_file(metadata, paths$metadata_path)
  
  # Report success
  report_save_success(experiment_name, dds, paths$dds_path)
  
  # Return paths (excluding outputs_dir for cleaner interface)
  list(
    dds_path = paths$dds_path,
    metadata_path = paths$metadata_path,
    coldata_path = paths$coldata_path
  )
}

load_dds <- function(experiment_name) {
  dds_path <- here("experiments", experiment_name, "outputs", paste0(experiment_name, ".dds.RDS"))
  
  if (!file.exists(dds_path)) {
    stop(glue("DDS not found. Run: Rscript run_experiment.R {experiment_name} --dds-only"))
  }
  
  return(read_rds(dds_path))
}


# ============================================================================= #
# UTILITY FUNCTIONS ----
# ============================================================================= #

#' Ensure experiment output directories exist
#'
#' Creates the required directory structure for an experiment if it doesn't exist.
#' Sources the function from rnaseq_output_management.R to avoid code duplication.
#'
#' @param experiment_name String identifier for the experiment
#'
#' @return String path to the outputs directory
ensure_experiment_outputs <- function(experiment_name) {
  source(here::here("functions", "rnaseq_output_management.R"))
  ensure_experiment_outputs(experiment_name)
}