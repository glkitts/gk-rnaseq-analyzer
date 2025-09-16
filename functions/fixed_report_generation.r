#!/usr/bin/env Rscript

# Fixed report generation functions
# Saves temporary files to outputs/R directory as requested

library(tidyverse)
library(DESeq2)
library(here)
library(yaml)
library(cli)

here::i_am("functions/fixed_report_generation.R")

# ============================================================================= #
# REPORT DATA PREPARATION ----
# ============================================================================= #

#' Prepare comprehensive data package for report generation
prepare_report_data <- function(res.l.all, 
                                dds, 
                                coldata, 
                                contrast_list,
                                config = NULL,
                                annotations = NULL,
                                gsea_results = NULL) {
  
  # Basic experiment metadata
  experiment_metadata <- list(
    n_samples = ncol(dds),
    n_transcripts = nrow(dds),
    n_comparisons = length(res.l.all),
    comparison_names = names(res.l.all),
    design_formula = as.character(design(dds))[2]
  )
  
  # Prepare PCA data with fixed column naming
  pca_data <- prepare_pca_data_fixed(dds)
  
  # Prepare heatmap data (top variable genes)
  heatmap_data <- prepare_heatmap_data_fixed(dds, top_n = 20)
  
  # Summary statistics per comparison
  comparison_summary <- res.l.all %>%
    imap_dfr(function(comp_data, comp_name) {
      tibble(
        comparison = comp_name,
        total_genes = nrow(comp_data$all),
        sig_genes = nrow(comp_data$sig),
        de_genes = nrow(comp_data$DE),
        max_log2fc = max(abs(comp_data$all$log2FoldChange), na.rm = TRUE),
        min_padj = if(nrow(comp_data$sig) > 0) min(comp_data$sig$padj, na.rm = TRUE) else NA
      )
    })
  
  # Prepare individual comparison data for detailed reports
  comparison_data <- res.l.all %>%
    imap(function(comp_data, comp_name) {
      # Find corresponding contrast
      contrast_idx <- which(sapply(contrast_list, function(c) {
        format_comparison_name(c) == comp_name
      }))
      
      contrast_info <- if (length(contrast_idx) > 0) {
        contrast_list[[contrast_idx[1]]]
      } else {
        c("group", str_split(comp_name, "_v_")[[1]])
      }
      
      list(
        name = comp_name,
        contrast = contrast_info,
        results = comp_data,
        volcano_data = prepare_volcano_data_fixed(comp_data$all),
        gsea_data = if (!is.null(gsea_results)) gsea_results[[comp_name]] else NULL
      )
    })
  
  # Return comprehensive data package
  list(
    experiment = list(
      metadata = experiment_metadata,
      coldata = coldata,
      config = config,
      analysis_date = Sys.Date()
    ),
    plotting = list(
      pca = pca_data,
      heatmap = heatmap_data,
      comparison_summary = comparison_summary
    ),
    comparisons = comparison_data,
    annotations = annotations
  )
}

#' Fixed PCA data preparation avoiding column name conflicts
prepare_pca_data_fixed <- function(dds) {
  tryCatch({
    # Use variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    
    # Calculate PCA
    pca_result <- prcomp(t(assay(vsd)))
    
    # Extract variance explained
    variance_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
    
    # Create plotting data with unique column names
    pca_coords <- pca_result$x %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
    
    # Get sample metadata with unique naming
    sample_metadata <- as.data.frame(colData(dds)) %>%
      rownames_to_column("sample_id") %>%
      # Rename any conflicting columns
      rename_with(~ paste0("meta_", .x), -sample_id)
    
    # Join PCA coordinates with metadata
    pca_data <- pca_coords %>%
      left_join(sample_metadata, by = "sample_id")
    
    list(
      data = pca_data,
      variance_explained = variance_explained[1:4],  # PC1-4
      loadings = pca_result$rotation[, 1:4]
    )
  }, error = function(e) {
    cli_warn("PCA preparation failed: {e$message}")
    return(NULL)
  })
}

#' Fixed heatmap data preparation
prepare_heatmap_data_fixed <- function(dds, top_n = 20) {
  tryCatch({
    # Get normalized counts
    norm_counts <- counts(dds, normalized = TRUE)
    
    # Calculate row variance and select top genes
    row_vars <- apply(norm_counts, 1, var, na.rm = TRUE)
    top_genes <- names(sort(row_vars, decreasing = TRUE))[1:min(top_n, length(row_vars))]
    
    # Extract and format data
    heatmap_matrix <- norm_counts[top_genes, , drop = FALSE] %>%
      log2() %>%
      t() %>%  # Samples as rows
      scale() %>%  # Z-score normalization
      t()  # Back to genes as rows
    
    # Clean gene names (extract VC#### from rownames)
    gene_labels <- rownames(heatmap_matrix) %>%
      str_extract("VC\\d+") %>%
      coalesce(rownames(heatmap_matrix))
    
    rownames(heatmap_matrix) <- gene_labels
    
    list(
      matrix = heatmap_matrix,
      sample_annotation = as.data.frame(colData(dds)),
      gene_labels = gene_labels
    )
  }, error = function(e) {
    cli_warn("Heatmap data preparation failed: {e$message}")
    return(NULL)
  })
}

#' Fixed volcano plot data preparation
prepare_volcano_data_fixed <- function(results_df) {
  tryCatch({
    # Clean data for plotting
    volcano_df <- results_df %>%
      filter(!is.na(log2FoldChange), !is.na(padj)) %>%
      mutate(
        # Handle zero p-values for visualization
        padj_plot = pmax(padj, 1e-300),
        # Significance categories
        sig_category = case_when(
          padj <= 0.05 & abs(log2FoldChange) >= 1 ~ "DE",
          padj <= 0.05 ~ "Significant",
          TRUE ~ "Not Significant"
        )
      )
    
    # Calculate dynamic axis limits
    x_range <- range(volcano_df$log2FoldChange, na.rm = TRUE)
    x_buffer <- diff(x_range) * 0.1
    x_limits <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
    
    y_max <- max(-log10(volcano_df$padj_plot), na.rm = TRUE)
    y_limits <- c(0, y_max * 1.1)
    
    list(
      data = volcano_df,
      x_limits = x_limits,
      y_limits = y_limits,
      summary = list(
        total_genes = nrow(volcano_df),
        de_genes = sum(volcano_df$sig_category == "DE"),
        sig_genes = sum(volcano_df$sig_category %in% c("DE", "Significant"))
      )
    )
  }, error = function(e) {
    cli_warn("Volcano data preparation failed: {e$message}")
    return(NULL)
  })
}

# ============================================================================= #
# REPORT GENERATION WITH OUTPUTS/R DIRECTORY ----
# ============================================================================= #

#' Generate reports using outputs/R directory for temporary files
generate_reports_fixed <- function(report_data, 
                                   experiment_name, 
                                   config = NULL,
                                   formats = c("html")) {  # Start with HTML only
  
  # Create output directories
  reports_dir <- here("experiments", experiment_name, "outputs", "reports")
  temp_dir <- here("experiments", experiment_name, "outputs", "R")
  
  if (!dir.exists(reports_dir)) {
    dir.create(reports_dir, recursive = TRUE)
  }
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }
  
  generated_files <- list()
  
  # Check if simple template exists
  simple_template <- here("templates", "simple_experimental_summary.qmd")
  
  if (file.exists(simple_template)) {
    cli_inform("Generating experimental summary report...")
    
    # Save report data to R directory
    temp_data_file <- file.path(temp_dir, paste0(experiment_name, "_report_data.RDS"))
    saveRDS(report_data, temp_data_file)
    
    # Generate summary report
    summary_files <- generate_summary_fixed(
      report_data, experiment_name, reports_dir, temp_data_file, formats
    )
    generated_files <- c(generated_files, summary_files)
    
  } else {
    cli_warn("Quarto template not found, using minimal HTML reports...")
    
    # Fall back to minimal reports
    source(here("functions", "fixed_minimal_reports.R"))
    minimal_files <- generate_minimal_reports(
      experiment_name, 
      # Extract results from report_data
      map(report_data$comparisons, "results"),
      # We don't have DDS here, create dummy
      list(ncol = report_data$experiment$metadata$n_samples, 
           nrow = report_data$experiment$metadata$n_transcripts),
      report_data$experiment$coldata,
      map(report_data$comparisons, "contrast")
    )
    generated_files <- c(generated_files, minimal_files)
  }
  
  cli_inform("Generated {length(generated_files)} report files")
  return(generated_files)
}

#' Generate experimental summary with better error handling
generate_summary_fixed <- function(report_data, experiment_name, reports_dir, temp_data_file, formats) {
  
  template_path <- here("templates", "simple_experimental_summary.qmd")
  generated_files <- list()
  
  for (format in formats) {
    output_file <- file.path(
      reports_dir, 
      paste0(experiment_name, "_summary.", if (format == "html") "html" else "pdf")
    )
    
    # Simple quarto command
    cmd <- paste(
      "quarto render",
      shQuote(template_path),
      "--output", shQuote(output_file),
      "--to", format,
      "-P", paste0("experiment_name:", experiment_name),
      "-P", paste0("temp_data_file:", temp_data_file),
      "-P", paste0("output_format:", format)
    )
    
    tryCatch({
      cli_inform("Running: {cmd}")
      result <- system(cmd, intern = TRUE)
      
      if (file.exists(output_file)) {
        cli_inform("Generated: {basename(output_file)}")
        generated_files[[paste0("summary_", format)]] <- output_file
      } else {
        cli_warn("Output file not created: {basename(output_file)}")
        cli_warn("Quarto output: {paste(result, collapse = '; ')}")
      }
      
    }, error = function(e) {
      cli_warn("Failed to generate {basename(output_file)}: {e$message}")
    })
  }
  
  # Clean up temp file
  if (file.exists(temp_data_file)) {
    unlink(temp_data_file)
  }
  
  return(generated_files)
}

# ============================================================================= #
# INTEGRATION FUNCTION ----
# ============================================================================= #

#' Main function to call from analysis_script.R
generate_experiment_reports <- function(experiment_name, res.l.all, dds, coldata, contrast_list, config = NULL, gsea_results = NULL) {
  
  cli_inform("Preparing report data for {experiment_name}...")
  
  tryCatch({
    # Prepare comprehensive data package
    report_data <- prepare_report_data(
      res.l.all = res.l.all,
      dds = dds,
      coldata = coldata,
      contrast_list = contrast_list,
      config = config,
      gsea_results = gsea_results
    )
    
    # Generate reports
    generated_files <- generate_reports_fixed(
      report_data = report_data,
      experiment_name = experiment_name,
      config = config,
      formats = c("html")  # Start with HTML only
    )
    
    cli_inform("Report generation completed!")
    return(generated_files)
    
  }, error = function(e) {
    cli_abort("Report generation failed: {e$message}")
  })
}