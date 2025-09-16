#!/usr/bin/env Rscript

# Report generation functions for RNA-seq analysis pipeline
# Creation date: 2025-09-15

library(tidyverse)
library(DESeq2)
library(here)
library(yaml)
library(glue)

here::i_am("functions/report_generation.R")

# ============================================================================= #
# REPORT DATA PREPARATION ----
# ============================================================================= #

#' Prepare comprehensive data package for report generation
#'
#' Creates a structured list containing all data needed for both experimental
#' summary and individual comparison reports. This approach allows reports to
#' be extended later without changing the pipeline.
#'
#' @param res.l.all Named list of comparison results from main workflow
#' @param dds DESeqDataSet object 
#' @param coldata Sample metadata
#' @param contrast_list List of contrasts used in analysis
#' @param config Global configuration
#' @param annotations Gene annotations (optional)
#' @param gsea_results GSEA results if available (optional)
#'
#' @return List with organized data for report templates
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
  
  # Prepare PCA data
  pca_data <- prepare_pca_data(dds)
  
  # Prepare heatmap data (top variable genes)
  heatmap_data <- prepare_heatmap_data(dds, top_n = 20)
  
  # Summary statistics per comparison
  comparison_summary <- res.l.all %>%
    imap_dfr(function(comp_data, comp_name) {
      tibble(
        comparison = comp_name,
        total_genes = nrow(comp_data$all),
        sig_genes = nrow(comp_data$sig),
        de_genes = nrow(comp_data$DE),
        max_log2fc = max(abs(comp_data$all$log2FoldChange), na.rm = TRUE),
        min_padj = min(comp_data$sig$padj, na.rm = TRUE)
      )
    })
  
  # Prepare individual comparison data for detailed reports
  comparison_data <- res.l.all %>%
    imap(function(comp_data, comp_name) {
      list(
        name = comp_name,
        contrast = contrast_list[[which(names(res.l.all) == comp_name)]],
        results = comp_data,
        # Volcano plot data with dynamic limits
        volcano_data = prepare_volcano_data(comp_data$all),
        # GSEA data if available
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

#' Prepare PCA data for plotting
prepare_pca_data <- function(dds) {
  # Use variance stabilizing transformation
  vsd <- vst(dds, blind = FALSE)
  
  # Calculate PCA
  pca_result <- prcomp(t(assay(vsd)))
  
  # Extract variance explained
  variance_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
  
  # Create plotting data
  pca_data <- pca_result$x %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    bind_cols(as.data.frame(colData(dds)))
  
  list(
    data = pca_data,
    variance_explained = variance_explained[1:4],  # PC1-4
    loadings = pca_result$rotation[, 1:4]
  )
}

#' Prepare heatmap data for top variable genes
prepare_heatmap_data <- function(dds, top_n = 20) {
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Calculate row variance and select top genes
  row_vars <- apply(norm_counts, 1, var)
  top_genes <- names(sort(row_vars, decreasing = TRUE))[1:top_n]
  
  # Extract and format data
  heatmap_matrix <- norm_counts[top_genes, ] %>%
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
}

#' Prepare volcano plot data with dynamic limits
prepare_volcano_data <- function(results_df) {
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
}

# ============================================================================= #
# REPORT GENERATION ----
# ============================================================================= #

#' Generate all reports for an experiment
#'
#' Creates both experimental summary and individual comparison reports
#' in HTML and PDF formats. HTML reports include cross-linking.
#'
#' @param report_data Prepared data from prepare_report_data()
#' @param experiment_name String identifier
#' @param config Global configuration
#' @param formats Vector of output formats (default: c("html", "pdf"))
#'
#' @return List of generated file paths
generate_reports <- function(report_data, 
                            experiment_name, 
                            config = NULL,
                            formats = c("html", "pdf")) {
  
  # Create reports directory
  reports_dir <- ensure_reports_directory(experiment_name)
  
  generated_files <- list()
  
  # Generate experimental summary report
  cli_inform("Generating experimental summary report...")
  summary_files <- generate_experimental_summary(
    report_data, experiment_name, reports_dir, formats
  )
  generated_files <- c(generated_files, summary_files)
  
  # Generate individual comparison reports
  cli_inform("Generating {length(report_data$comparisons)} comparison reports...")
  comparison_files <- report_data$comparisons %>%
    map(function(comp_data) {
      generate_comparison_report(
        comp_data, report_data$experiment, reports_dir, formats
      )
    }) %>%
    flatten()
  
  generated_files <- c(generated_files, comparison_files)
  
  # Generate index page for HTML reports (if HTML format requested)
  if ("html" %in% formats) {
    index_file <- generate_report_index(report_data, experiment_name, reports_dir)
    generated_files <- c(generated_files, list(index = index_file))
  }
  
  cli_inform("Generated {length(generated_files)} report files")
  return(generated_files)
}

#' Generate experimental summary report
generate_experimental_summary <- function(report_data, experiment_name, reports_dir, formats) {
  
  template_path <- here("templates", "experimental_summary.qmd")
  
  if (!file.exists(template_path)) {
    cli_warn("Template not found: {template_path}")
    return(list())
  }
  
  generated_files <- list()
  
  for (format in formats) {
    output_file <- file.path(
      reports_dir, 
      paste0(experiment_name, "_summary.", if (format == "html") "html" else "pdf")
    )
    
    # Render with quarto
    quarto_render_safe(
      input = template_path,
      output_file = output_file,
      output_format = format,
      execute_params = list(
        experiment_name = experiment_name,
        report_data = report_data,
        output_format = format
      )
    )
    
    generated_files[[paste0("summary_", format)]] <- output_file
  }
  
  return(generated_files)
}

#' Generate individual comparison report
generate_comparison_report <- function(comp_data, experiment_info, reports_dir, formats) {
  
  template_path <- here("templates", "comparison_report.qmd")
  
  if (!file.exists(template_path)) {
    cli_warn("Template not found: {template_path}")
    return(list())
  }
  
  generated_files <- list()
  comparison_name <- comp_data$name
  
  for (format in formats) {
    output_file <- file.path(
      reports_dir, 
      paste0(comparison_name, "_report.", if (format == "html") "html" else "pdf")
    )
    
    # Render with quarto
    quarto_render_safe(
      input = template_path,
      output_file = output_file,
      output_format = format,
      execute_params = list(
        comparison_name = comparison_name,
        comparison_data = comp_data,
        experiment_info = experiment_info,
        output_format = format
      )
    )
    
    generated_files[[paste0(comparison_name, "_", format)]] <- output_file
  }
  
  return(generated_files)
}

# ============================================================================= #
# UTILITY FUNCTIONS ----
# ============================================================================= #

#' Ensure reports directory exists
ensure_reports_directory <- function(experiment_name) {
  reports_dir <- here("experiments", experiment_name, "outputs", "reports")
  
  if (!dir.exists(reports_dir)) {
    dir.create(reports_dir, recursive = TRUE)
  }
  
  return(reports_dir)
}

#' Safe quarto rendering with error handling
quarto_render_safe <- function(input, output_file, output_format, execute_params) {
  
  tryCatch({
    # Use quarto::quarto_render if available, otherwise system call
    if (requireNamespace("quarto", quietly = TRUE)) {
      quarto::quarto_render(
        input = input,
        output_file = output_file,
        output_format = output_format,
        execute_params = execute_params,
        quiet = TRUE
      )
    } else {
      # Fallback to system call
      param_yaml <- tempfile(fileext = ".yml")
      yaml::write_yaml(execute_params, param_yaml)
      
      system(glue(
        "quarto render {input} ",
        "--output {output_file} ",
        "--to {output_format} ",
        "--execute-params {param_yaml}"
      ))
      
      unlink(param_yaml)
    }
    
    cli_inform("Generated: {basename(output_file)}")
    
  }, error = function(e) {
    cli_warn("Failed to generate {basename(output_file)}: {e$message}")
  })
}

#' Generate HTML index page linking all reports
generate_report_index <- function(report_data, experiment_name, reports_dir) {
  
  index_content <- glue(
    "# {experiment_name} - Analysis Reports\n\n",
    "Generated: {Sys.Date()}\n\n",
    "## Experimental Summary\n",
    "- [Summary Report]({experiment_name}_summary.html)\n\n",
    "## Individual Comparisons\n",
    "{paste(map_chr(report_data$comparisons, ~ glue('- [{.x$name}]({.x$name}_report.html)')), collapse = '\n')}\n\n",
    "## Analysis Overview\n",
    "- **Samples:** {report_data$experiment$metadata$n_samples}\n",
    "- **Transcripts:** {report_data$experiment$metadata$n_transcripts}\n", 
    "- **Comparisons:** {report_data$experiment$metadata$n_comparisons}\n"
  )
  
  index_file <- file.path(reports_dir, "index.html")
  
  # Convert markdown to HTML (simple approach)
  if (requireNamespace("markdown", quietly = TRUE)) {
    html_content <- markdown::markdownToHTML(text = index_content, fragment.only = TRUE)
    writeLines(html_content, index_file)
  } else {
    # Fallback: write as markdown
    index_file <- file.path(reports_dir, "index.md")
    writeLines(index_content, index_file)
  }
  
  return(index_file)
}