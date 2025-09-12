#!/usr/bin/env Rscript

# RNAseq analysis functions
# Creation date: 2025-09-04

library(tidyverse)
library(DESeq2)
library(tximeta)
library(glue)
library(openxlsx2)
library(here)
library(yaml)

here::i_am("functions/output_management.R")

ensure_experiment_outputs <- function(experiment_name) {
  #' Create complete experiment output directory structure if missing
  #' 
  #' Called by any function that needs to write files.
  #' 
  #' @param experiment_name String: experiment identifier
  #' @return String: path to outputs directory
  
  # Base experiment directory
  experiment_dir <- here("experiments", experiment_name)
  if (!dir.exists(experiment_dir)) {
    stop(sprintf(
      "Experiment directory not found: %s\n💡 Create it with: Rscript new_experiment_setup.R %s",
      experiment_name, experiment_name
    ))
  }
  
  # Base outputs directory
  outputs_dir <- file.path(experiment_dir, "outputs")
  
  # Complete directory structure (matches new_experiment_setup.R)
  required_dirs <- c(
    "outputs",
    "outputs/results", 
    "outputs/results/individual_comparisons",
    "outputs/plots",
    "outputs/R", 
    "outputs/metadata",
    "outputs/logs"
  )
  
  dirs_created <- 0
  for (rel_dir in required_dirs) {
    full_dir <- file.path(experiment_dir, rel_dir)
    if (!dir.exists(full_dir)) {
      dir.create(full_dir, recursive = TRUE, showWarnings = FALSE)
      dirs_created <- dirs_created + 1
    }
  }
  
  # Only announce if we had to create missing directories
  if (dirs_created > 0) {
    cat(sprintf("📁 Created %d missing output directories\n", dirs_created))
  }
  
  return(outputs_dir)
}


save_individual_csvs <- function(res.l,
                                 experiment_name,
                                 comparison_name,
                                 output_dir = NULL,
                                 types = c("all", "DE")) {
  ensure_experiment_outputs(experiment_name)
  
  if (is.null(output_dir)) {
    output_dir <- here("experiments",
                       experiment_name,
                       "outputs/results/individual_comparisons/")
  }
  
  filtered_res.l <- res.l[names(res.l) %in% types]
  
  # Return file paths instead of data
  file_paths <- filtered_res.l %>%
    imap_chr(function(res, name) {
      output_path <- paste0(output_dir, comparison_name, "_", name, ".csv")
      
      res %>%
        readr::write_excel_csv(file = output_path, na = "")
      
      return(output_path)
    })
  
  return(file_paths)
}

save_master_table <- function(allHits_compiled,
                              experiment_name,
                              output_dir_rds = NULL,
                              output_dir_csv = NULL,
                              save_rds = T,
                              save_csv = F) {
  ensure_experiment_outputs(experiment_name)
  
  if (save_rds == T) {
    if (is.null(output_dir_rds)) {
      output_dir <- here("experiments", experiment_name, "outputs/R/")
    }
    output_name <- paste0(output_dir, experiment_name, ".masterTable.RDS")
    allHits_compiled %>%
      write_rds(file = output_name)
  }
  
  
  
  if (save_csv == T) {
    if (is.null(output_dir_csv)) {
      output_dir <- here("experiments", experiment_name, "outputs/results/")
    }
    output_name <- paste0(output_dir, experiment_name, ".masterTable.csv")
    allHits_compiled %>%
      readr::write_excel_csv(file = output_name, na = "")
  }
  
  return(allHits_compiled)
  
}


# ============================================================================= #
# PREPARE FOR EXCEL OUTPUT ----
# ============================================================================= #

#' Save Excel file with enhanced formatting using openxlsx2
#'
#' Creates Excel files with intelligent column sizing, number formatting,
#' and professional styling. Optimized for differential expression results.
#'
#' @param data_list Named list of data frames (from create_de_compilation)
#' @param file_path String, output file path (.xlsx extension)
#' @param experiment_name String, for metadata and sheet naming
#' @param freeze_panes Logical, freeze first row and first two columns (default: TRUE)
#' @param add_filters Logical, add autofilters to data sheets (default: TRUE)
#'
#' @return String file path of created Excel file
save_excel_enhanced <- function(data_list,
                                file_path,
                                experiment_name,
                                freeze_panes = TRUE,
                                add_filters = TRUE) {
  
  # Validate inputs
  if (!is.list(data_list) || is.null(names(data_list))) {
    cli_abort("data_list must be a named list of data frames")
  }
  if (length(data_list) == 0) {
    cli_abort("data_list cannot be empty")
  }
  if (!str_ends(file_path, "\\.xlsx$")) {
    cli_abort("file_path must end with .xlsx")
  }
  
  # Filter out empty data frames
  valid_sheets <- data_list %>%
    keep(~ is.data.frame(.x) && nrow(.x) > 0)
  
  if (length(valid_sheets) == 0) {
    cli_abort("No valid data frames found in data_list")
  }
  
  # Create workbook
  wb <- wb_workbook()
  
  # Add worksheets and data
  iwalk(valid_sheets, ~ {
    sheet_name <- .y
    sheet_data <- .x
    
    # Add worksheet and data
    wb <<- wb$add_worksheet(sheet = sheet_name)$
      add_data(sheet = sheet_name, x = sheet_data, na.strings = "")
  })
  
  # Apply formatting to all sheets
  iwalk(valid_sheets, ~ {
    wb <<- apply_enhanced_formatting(wb, .y, .x, freeze_panes, add_filters)
  })
  
  # Save workbook
  tryCatch({
    wb$save(file = file_path, overwrite = TRUE)
  }, error = function(e) {
    cli_abort("Failed to save Excel file: {e$message}")
  })
  
  return(file_path)
}

#' Apply enhanced formatting to a worksheet
#'
#' @param wb Workbook object
#' @param sheet_name String, worksheet name
#' @param sheet_data Data frame being written
#' @param freeze_panes Logical
#' @param add_filters Logical
#'
#' @return Modified workbook object
apply_enhanced_formatting <- function(wb, sheet_name, sheet_data, 
                                      freeze_panes, add_filters) {
  
  n_rows <- nrow(sheet_data) + 1  # +1 for header
  n_cols <- ncol(sheet_data)
  
  # 1. HEADER FORMATTING - bold headers with light background
  tryCatch({
    header_dims <- wb_dims(rows = 1, cols = 1:n_cols)
    wb <- wb$add_font(sheet = sheet_name, dims = header_dims, bold = TRUE)$
      add_fill(sheet = sheet_name, dims = header_dims, color = wb_color("lightgray"))
  }, error = function(e) {
    # Fallback to just bold if fill fails
    tryCatch({
      wb <<- wb$add_font(sheet = sheet_name, dims = header_dims, bold = TRUE)
    }, error = function(e2) {
      # Skip header formatting entirely if it fails
    })
  })
  
  # 2. COLUMN WIDTHS with special handling for README sheet
  if (str_detect(sheet_name, "README")) {
    # Auto-fit README sheet columns
    tryCatch({
      wb <- wb$set_col_widths(sheet = sheet_name, cols = 1:n_cols, widths = "auto")
    }, error = function(e) {
      # Skip if auto-fit fails
    })
  } else {
    # Use optimized widths for data sheets
    column_widths <- calculate_optimized_column_widths(sheet_data)
    
    if (length(column_widths) > 0) {
      tryCatch({
        wb <- wb$set_col_widths(sheet = sheet_name, cols = seq_along(column_widths), widths = column_widths)
      }, error = function(e) {
        # Skip column widths if it fails
      })
    }
  }
  
  # 3. NUMBER FORMATTING for log2FC and padj columns
  format_columns <- detect_format_columns(sheet_data)
  wb <- apply_number_formatting(wb, sheet_name, format_columns, n_rows)
  
  # 4. FREEZE PANES (first row and first two columns)
  if (freeze_panes && n_rows > 1 && n_cols > 2) {
    tryCatch({
      wb <- wb$freeze_pane(sheet = sheet_name, first_active_row = 2, first_active_col = 3)
    }, error = function(e) {
      # Skip freeze panes if it fails
    })
  }
  
  # 5. AUTO FILTERS
  if (add_filters && n_rows > 1) {
    tryCatch({
      filter_dims <- wb_dims(rows = 1, cols = 1:n_cols)
      wb <- wb$add_filter(sheet = sheet_name, dims = filter_dims)
    }, error = function(e) {
      # Skip filters if it fails
    })
  }
  
  return(wb)
}

#' Determine if column is a log2FC type column using fixed patterns
#'
#' @param col_name String, column name
#' @return Logical
is_log2fc_column <- function(col_name) {
  str_detect(col_name, "log2FC|log2FoldChange|logFC") || 
    str_ends(col_name, fixed(".log2FC"))
}

#' Determine if column is a padj type column using fixed patterns
#'
#' @param col_name String, column name
#' @return Logical
is_padj_column <- function(col_name) {
  str_detect(col_name, "padj") || 
    str_ends(col_name, fixed(".padj")) ||
    str_detect(col_name, "p\\.adj|pvalue|pval")
}

#' Calculate optimized column widths with FIXED widths for specific types
#'
#' @param data Data frame
#' @return Numeric vector of column widths
calculate_optimized_column_widths <- function(data) {
  
  col_names <- colnames(data)
  n_cols <- length(col_names)
  
  # Calculate width for each column individually
  widths <- map_dbl(seq_len(n_cols), ~ {
    col_name <- col_names[.x]
    col_data <- data[[.x]]
    width.log2FC_padj <- 20
    
    # Base width on header length
    header_width <- nchar(col_name) * 1.2
    
    # FIXED WIDTHS for specific column types (don't use max with header_width)
    if (is_log2fc_column(col_name)) {
      # log2FC columns: FIXED width of 8 (ignore header length)
      return(width.log2FC_padj)
    } else if (is_padj_column(col_name)) {
      # P-value columns: FIXED width of 10 (ignore header length)
      return(width.log2FC_padj)
    } else if (str_detect(col_name, "GO terms|GO_")) {
      # GO terms: long content (unchanged)
      return(40)
    } else if (str_detect(col_name, "Product|Description|Function")) {
      # Product descriptions: 20% wider than before (35 -> 42)
      return(42)
    } else if (str_detect(col_name, "yildiz_|pathway|Pathway")) {
      # Pathway annotations: 2/3 of original (25 -> 17)
      return(17)
    } else if (str_detect(col_name, "Label|Gene|gene")) {
      # Gene identifiers (unchanged)
      return(max(header_width, 15))
    } else if (str_detect(col_name, "baseMean|Mean")) {
      # Mean expression columns (unchanged)
      return(max(header_width, 11))
    } else if (str_detect(col_name, "nc_")) {
      # Normalized count columns (unchanged)
      return(max(header_width, 10))
    } else if (str_detect(col_name, "Type|Begin|End|Length")) {
      # Metadata columns (unchanged)
      return(max(header_width, 10))
    } else {
      # Default: adaptive based on content sample
      return(calculate_adaptive_width_simple(col_data, header_width))
    }
  })
  
  # Ensure reasonable bounds
  widths <- pmax(pmin(widths, 50), 6)  # Min 6, Max 50
  
  return(widths)
}

#' Calculate adaptive width for a single column
#'
#' @param col_data Vector of column data
#' @param header_width Numeric, width based on header
#' @return Numeric width
calculate_adaptive_width_simple <- function(col_data, header_width) {
  
  if (!is.character(col_data) || length(col_data) == 0) {
    return(max(header_width, 11))
  }
  
  # Sample non-NA content for performance
  sample_data <- col_data[!is.na(col_data)]
  if (length(sample_data) == 0) {
    return(max(header_width, 12))
  }
  
  # Sample for efficiency on large datasets
  sample_size <- min(10, length(sample_data))
  content_sample <- sample(sample_data, sample_size)
  avg_content_width <- mean(nchar(content_sample), na.rm = TRUE)
  content_width <- min(avg_content_width * 1.1, 30)  # Cap at 30
  
  return(max(header_width, content_width, 8))  # Minimum 8
}

#' Detect columns needing special number formatting
#'
#' Uses the same logic as column width calculation for consistency
#'
#' @param data Data frame
#' @return List with log2fc and padj column indices
detect_format_columns <- function(data) {
  col_names <- colnames(data)
  
  list(
    # Use the same detection functions for consistency
    log2fc = which(map_lgl(col_names, is_log2fc_column)),
    padj = which(map_lgl(col_names, is_padj_column))
  )
}

#' Apply number formatting to specific column types
#'
#' @param wb Workbook object
#' @param sheet_name String, worksheet name
#' @param format_columns List with column indices to format
#' @param n_rows Integer, total rows including header
#' @return Modified workbook object
apply_number_formatting <- function(wb, sheet_name, format_columns, n_rows) {
  
  # Only apply formatting if we have data rows
  if (n_rows <= 2) {
    return(wb)
  }
  
  # Format log2FC columns (2 decimal places)
  if (length(format_columns$log2fc) > 0) {
    walk(format_columns$log2fc, ~ {
      tryCatch({
        data_dims <- wb_dims(rows = 2:n_rows, cols = .x)
        wb <<- wb$add_numfmt(sheet = sheet_name, dims = data_dims, numfmt = "0.00")
      }, error = function(e) {
        # Skip this column formatting if it fails
      })
    })
  }
  
  # Format padj columns (scientific notation)
  if (length(format_columns$padj) > 0) {
    walk(format_columns$padj, ~ {
      tryCatch({
        data_dims <- wb_dims(rows = 2:n_rows, cols = .x)
        wb <<- wb$add_numfmt(sheet = sheet_name, dims = data_dims, numfmt = "0.00E+00")
      }, error = function(e) {
        # Skip this column formatting if it fails
      })
    })
  }
  
  return(wb)
}

#' Wrapper function that integrates with the existing pipeline
#'
#' Updated to save directly to outputs folder as primary result
#'
#' @param de_compilation Named list from create_de_compilation()
#' @param experiment_name String, experiment identifier
#' @param output_dir String, directory for output file (optional, uses ensure_experiment_outputs if NULL)
#' @param filename String, optional custom filename (without extension)
#'
#' @return String path to created Excel file
save_de_compilation_excel <- function(de_compilation,
                                      experiment_name,
                                      output_dir = NULL,
                                      filename = NULL) {
  
  # Create output directory - save directly to outputs as primary result
  if (is.null(output_dir)) {
    output_dir <- ensure_experiment_outputs(experiment_name)
    # Save directly to outputs folder instead of outputs/results
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate filename
  if (is.null(filename)) {
    filename <- paste0(experiment_name, "_DEcompilation.xlsx")
  } else if (!str_ends(filename, "\\.xlsx$")) {
    filename <- paste0(filename, ".xlsx")
  }
  
  file_path <- file.path(output_dir, filename)
  
  # Save with enhanced formatting
  save_excel_enhanced(
    data_list = de_compilation,
    file_path = file_path,
    experiment_name = experiment_name
  )
  
  return(file_path)
}


