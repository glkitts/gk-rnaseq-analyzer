# RNAseq differential analysis functions
# Creation date: 2025-09-04

library(tidyverse)
library(DESeq2)
library(tximeta)
library(glue)
library(yaml)
library(here)

here::i_am("functions/input_validation.R")

# In functions/input_validation.R or at top of differential_expression.R

validate_inputs <- function(..., .context = "analysis") {
  #' Validate any combination of pipeline inputs
  #' 
  #' @param ... Named arguments to validate (dds, contrast, annotations, etc.)
  #' @param .context String describing where validation is called from (for better error messages)
  #' 
  #' @examples
  #' validate_inputs(dds = my_dds, contrast = my_contrast)
  #' validate_inputs(annotations = annot_df, gene_sets = gs_list, .context = "compilation")
  
  inputs <- list(...)
  if (length(inputs) == 0) return(TRUE)
  
  # Validate each provided input
  for (input_name in names(inputs)) {
    input_value <- inputs[[input_name]]
    
    # Skip NULL inputs (optional parameters)
    if (is.null(input_value)) next
    
    # Dispatch to appropriate validator
    tryCatch({
      switch(input_name,
             "dds" = validate_dds(input_value, .context),
             "contrast" = validate_contrast_vector(input_value, .context),
             "contrast_list" = validate_contrast_list(input_value, .context),
             "annotations" = validate_annotations(input_value, .context),
             "gene_sets" = validate_gene_sets(input_value, .context),
             "config" = validate_config(input_value, .context),
             "file_path" = validate_file_path(input_value, .context),
             "output_dir" = validate_output_dir(input_value, .context),
             # Add more as needed
             cli_warn("Unknown input type for validation: {input_name}")
      )
    }, error = function(e) {
      # Re-throw with context about which input failed
      cli_abort("Validation failed for {input_name} in {.context}: {e$message}")
    })
  }
  
  # Special cross-input validations (when multiple inputs interact)
  if (all(c("dds", "contrast") %in% names(inputs))) {
    validate_contrast_against_dds(inputs$contrast, inputs$dds, .context)
  }
  
  if (all(c("dds", "contrast_list") %in% names(inputs))) {
    walk(inputs$contrast_list, validate_contrast_against_dds, 
         dds = inputs$dds, .context = .context)
  }
  
  TRUE
}

# Individual validation functions
validate_dds <- function(dds, context) {
  if (!is(dds, "DESeqDataSet")) {
    cli_abort("dds must be a DESeqDataSet object. Got {class(dds)[1]}")
  }
  
  if (nrow(dds) == 0) {
    cli_abort("DDS object is empty (0 transcripts)")
  }
  
  if (ncol(dds) == 0) {
    cli_abort("DDS object has no samples")
  }
  
  TRUE
}

validate_contrast_vector <- function(contrast, context) {
  if (!is.character(contrast) || length(contrast) != 3) {
    cli_abort("contrast must be character vector of length 3: c(factor, treatment, control)")
  }
  
  if (any(contrast == "")) {
    cli_abort("contrast elements cannot be empty strings")
  }
  
  TRUE
}

validate_contrast_list <- function(contrast_list, context) {
  if (!is.list(contrast_list) || length(contrast_list) == 0) {
    cli_abort("contrast_list must be a non-empty list")
  }
  
  # Validate each contrast in the list
  iwalk(contrast_list, function(contrast, idx) {
    tryCatch({
      validate_contrast_vector(contrast, context)
    }, error = function(e) {
      cli_abort("Invalid contrast at position {idx}: {e$message}")
    })
  })
  
  TRUE
}

validate_annotations <- function(annotations, context) {
  if (!is.data.frame(annotations)) {
    cli_abort("annotations must be a data frame. Got {class(annotations)[1]}")
  }
  
  if (!"Label" %in% colnames(annotations)) {
    available_cols <- colnames(annotations)[1:min(5, ncol(annotations))]
    cli_abort(c(
      "annotations must contain 'Label' column for joining",
      "i" = "Available columns: {paste(available_cols, collapse = ', ')}{if(ncol(annotations) > 5) '...' else ''}"
    ))
  }
  
  if (nrow(annotations) == 0) {
    cli_abort("annotations data frame is empty")
  }
  
  TRUE
}

validate_gene_sets <- function(gene_sets, context) {
  if (!is.list(gene_sets)) {
    cli_abort("gene_sets must be a list. Got {class(gene_sets)[1]}")
  }
  
  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    cli_abort("gene_sets must be a named list (all elements need names)")
  }
  
  # Check that gene sets contain character vectors
  invalid_sets <- map_lgl(gene_sets, function(gs) {
    !is.character(gs) || length(gs) == 0
  })
  
  if (any(invalid_sets)) {
    invalid_names <- names(gene_sets)[invalid_sets]
    cli_abort("Invalid gene sets (must be non-empty character vectors): {paste(invalid_names, collapse = ', ')}")
  }
  
  TRUE
}

validate_config <- function(config, context) {
  if (!is.list(config)) {
    cli_abort("config must be a list. Got {class(config)[1]}")
  }
  
  # Check for required sections
  required_sections <- c("analysis", "paths")
  missing_sections <- setdiff(required_sections, names(config))
  
  if (length(missing_sections) > 0) {
    cli_abort("config missing required sections: {paste(missing_sections, collapse = ', ')}")
  }
  
  # Check for required analysis parameters
  required_analysis <- c("significance_threshold", "fold_change_threshold")
  missing_analysis <- setdiff(required_analysis, names(config$analysis))
  
  if (length(missing_analysis) > 0) {
    cli_abort("config$analysis missing required parameters: {paste(missing_analysis, collapse = ', ')}")
  }
  
  TRUE
}

validate_file_path <- function(file_path, context) {
  if (!is.character(file_path) || length(file_path) != 1) {
    cli_abort("file_path must be a single character string")
  }
  
  if (!file.exists(file_path)) {
    cli_abort("File does not exist: {file_path}")
  }
  
  TRUE
}

validate_output_dir <- function(output_dir, context) {
  if (!is.character(output_dir) || length(output_dir) != 1) {
    cli_abort("output_dir must be a single character string")
  }
  
  # Try to create if it doesn't exist
  if (!dir.exists(output_dir)) {
    tryCatch({
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }, error = function(e) {
      cli_abort("Cannot create output directory: {output_dir}\n{e$message}")
    })
  }
  
  TRUE
}

# Cross-input validation
validate_contrast_against_dds <- function(contrast, dds, context) {
  factor_name <- contrast[1]
  treatment <- contrast[2] 
  control <- contrast[3]
  
  if (!factor_name %in% names(colData(dds))) {
    available_factors <- names(colData(dds))
    cli_abort(c(
      "Factor '{factor_name}' not found in DDS colData",
      "i" = "Available factors: {paste(available_factors, collapse = ', ')}"
    ))
  }
  
  factor_levels <- levels(colData(dds)[[factor_name]])
  missing_levels <- setdiff(c(treatment, control), factor_levels)
  
  if (length(missing_levels) > 0) {
    cli_abort(c(
      "Contrast levels not found: {paste(missing_levels, collapse = ', ')}",
      "i" = "Available levels for {factor_name}: {paste(factor_levels, collapse = ', ')}"
    ))
  }
  
  TRUE
}