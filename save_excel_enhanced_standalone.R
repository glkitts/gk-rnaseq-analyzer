#!/usr/bin/env Rscript


# Enhanced Excel Export with Smart Formatting
# ============================================================================= -
#
# Exports named lists of data frames to Excel with intelligent column sizing
# and formatting. Optimized for RNA-seq data but works with any tabular data.
#
# DEPENDENCIES:
#   library(openxlsx2)
#   library(tidyverse)
#
# USAGE:
#   source("save_excel_enhanced_standalone.R")
#
#   my_data <- list(
#     "Results" = tibble(Label = c("VC0001", "VC0002"),
#                        log2FoldChange = c(2.5, -1.3),
#                        padj = c(0.001, 0.05)),
#     "Summary" = tibble(comparison = c("A_v_B"), n_genes = c(150))
#   )
#
#   save_excel_enhanced(
#     data_list = my_data,
#     file_path = "output.xlsx",
#     freeze_panes = TRUE,
#     add_filters = TRUE
#   )
#
# ============================================================================= -

library(openxlsx2)
library(tidyverse)


#' Save Excel file with enhanced formatting
#'
#' Creates Excel files with intelligent column sizing and number formatting.
#' Applies RNA-seq specific formatting to recognized column types (log2FC, padj, etc.)
#' and adaptive formatting to other columns.
#'
#' @param data_list Named list of data frames
#' @param file_path String, output file path (.xlsx extension)
#' @param freeze_panes Logical, freeze first row and first two columns (default: TRUE)
#' @param add_filters Logical, add autofilters to data sheets (default: TRUE)
#'
#' @return String file path of created Excel file
save_excel_enhanced <- function(data_list,
                                file_path,
                                freeze_panes = TRUE,
                                add_filters = TRUE) {

  # Filter out empty data frames
  valid_sheets <- data_list %>%
    keep(~ is.data.frame(.x) && nrow(.x) > 0)

  # Create workbook
  wb <- wb_workbook()

  # Add worksheets and data
  iwalk(valid_sheets, ~ {
    wb <<- wb$add_worksheet(sheet = .y)$
      add_data(sheet = .y, x = .x, na.strings = NULL)
  })

  # Apply formatting to all sheets
  iwalk(valid_sheets, ~ {
    wb <<- apply_enhanced_formatting(wb, .y, .x, freeze_panes, add_filters)
  })

  # Save workbook
  wb$save(file = file_path, overwrite = TRUE)

  return(file_path)
}


#' Apply enhanced formatting to a worksheet
#'
#' Consolidates all formatting logic: column widths, number formats,
#' header styling, freeze panes, and filters.
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
  col_names <- colnames(sheet_data)

  # ======================================================================== #
  # HEADER FORMATTING
  # ======================================================================== #
  header_dims <- wb_dims(rows = 1, cols = 1:n_cols)
  wb <- wb$add_font(sheet = sheet_name, dims = header_dims, bold = TRUE)$
    add_fill(sheet = sheet_name, dims = header_dims, color = wb_color("lightgray"))

  # ======================================================================== #
  # COLUMN WIDTHS - Smart sizing based on column type
  # ======================================================================== #
  if (str_detect(sheet_name, "README")) {
    # Auto-fit README sheets
    wb <- wb$set_col_widths(sheet = sheet_name, cols = 1:n_cols, widths = "auto")
  } else {
    # Calculate optimized widths for data sheets
    column_widths <- map_dbl(seq_len(n_cols), ~ {
      col_name <- col_names[.x]
      col_data <- sheet_data[[.x]]
      header_width <- nchar(col_name) * 1.2

      # RNA-seq specific column widths
      if (str_detect(col_name, "log2FC|log2FoldChange|logFC") ||
          str_ends(col_name, fixed(".log2FC"))) {
        return(20)  # log2FC columns
      } else if (str_detect(col_name, "padj|p\\.adj|pvalue|pval") ||
                 str_ends(col_name, fixed(".padj"))) {
        return(20)  # P-value columns
      } else if (str_detect(col_name, "GO terms|GO_")) {
        return(40)  # GO terms
      } else if (str_detect(col_name, "Product|Description|Function")) {
        return(42)  # Product descriptions
      } else if (str_detect(col_name, "yildiz_|pathway|Pathway")) {
        return(17)  # Pathway annotations
      } else if (str_detect(col_name, "Label|Gene|gene")) {
        return(max(header_width, 15))  # Gene identifiers
      } else if (str_detect(col_name, "baseMean|Mean")) {
        return(max(header_width, 11))  # Mean expression
      } else if (str_detect(col_name, "nc_")) {
        return(max(header_width, 10))  # Normalized counts
      } else if (str_detect(col_name, "vst_")) {
        return(max(header_width, 10))  # VST counts
      } else if (str_detect(col_name, "Type|Begin|End|Length")) {
        return(max(header_width, 10))  # Metadata columns
      } else {
        # Adaptive width for unknown columns
        if (is.character(col_data) && length(col_data) > 0) {
          sample_data <- col_data[!is.na(col_data)]
          if (length(sample_data) > 0) {
            sample_size <- min(10, length(sample_data))
            content_sample <- sample(sample_data, sample_size)
            avg_width <- mean(nchar(content_sample), na.rm = TRUE) * 1.1
            return(max(header_width, min(avg_width, 30), 8))
          }
        }
        return(max(header_width, 11))
      }
    })

    # Apply bounds (min 6, max 50)
    column_widths <- pmax(pmin(column_widths, 50), 6)
    wb <- wb$set_col_widths(sheet = sheet_name,
                           cols = seq_along(column_widths),
                           widths = column_widths)
  }

  # ======================================================================== #
  # NUMBER FORMATTING
  # ======================================================================== #
  if (n_rows > 2) {
    # Detect and format log2FC columns (2 decimal places)
    log2fc_cols <- which(map_lgl(col_names, ~ {
      str_detect(.x, "log2FC|log2FoldChange|logFC") || str_ends(.x, fixed(".log2FC"))
    }))

    walk(log2fc_cols, ~ {
      data_dims <- wb_dims(rows = 2:n_rows, cols = .x)
      wb <<- wb$add_numfmt(sheet = sheet_name, dims = data_dims, numfmt = "0.00")
    })

    # Detect and format padj columns (scientific notation)
    padj_cols <- which(map_lgl(col_names, ~ {
      str_detect(.x, "padj|p\\.adj|pvalue|pval") || str_ends(.x, fixed(".padj"))
    }))

    walk(padj_cols, ~ {
      data_dims <- wb_dims(rows = 2:n_rows, cols = .x)
      wb <<- wb$add_numfmt(sheet = sheet_name, dims = data_dims, numfmt = "0.00E+00")
    })
  }

  # ======================================================================== #
  # FREEZE PANES & FILTERS
  # ======================================================================== #
  if (freeze_panes && n_rows > 1 && n_cols > 2) {
    wb <- wb$freeze_pane(sheet = sheet_name, first_active_row = 2, first_active_col = 3)
  }

  if (add_filters && n_rows > 1) {
    wb <- wb$add_filter(sheet = sheet_name, rows = 1, cols = 1:n_cols)
  }

  return(wb)
}
