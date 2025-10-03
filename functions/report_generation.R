#!/usr/bin/env Rscript

# Report generation functions for RNA-seq analysis pipeline
# Creation date: 2025-09-27
# Following gk_functions.R tidyverse style and report_generation_design.md specification

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(here)
  library(yaml)
  library(glue)
  library(cli)
})

here::i_am("functions/report_generation.R")

# Source output management functions for ensure_experiment_outputs
source(here("functions/output_management.R"))

# ============================================================================= #
# CORE REPORT GENERATION FUNCTIONS ----
# ============================================================================= #

#' Compile comprehensive report data from analysis results
#'
#' Creates single RDS file containing all data needed for report generation
#' Following the single-file approach from design specification
#'
#' @param dds DESeqDataSet object with normalized counts
#' @param res.l.all Named list of comparison results
#' @param contrast_list List of contrasts used in analysis
#' @param coldata Sample metadata
#' @param config Global configuration from analysis
#' @param gsea_results GSEA results if available (optional)
#'
#' @return List with organized data for report templates
compile_report_data <- function(dds, res.l.all, contrast_list, coldata, config, gsea_results = NULL) {

  cli_inform("Compiling report data for {length(res.l.all)} comparisons...")

  # Analysis metadata
  analysis_metadata <- list(
    experiment_name = config$experiment_name,
    n_samples = ncol(dds),
    n_genes = nrow(dds),
    n_comparisons = length(res.l.all),
    comparison_names = names(res.l.all),
    design_formula = as.character(design(dds))[2],
    analysis_date = Sys.Date(),
    pipeline_version = "2025-09-27"
  )

  # Comparison summary for overview using imap
  comparison_summary <- res.l.all %>%
    imap_dfr(~ tibble(
      comparison = .y,
      total_genes = nrow(.x$all),
      sig_genes = nrow(.x$sig),
      de_genes = nrow(.x$DE),
      max_log2fc = round(max(abs(.x$all$log2FoldChange), na.rm = TRUE), 2),
      min_padj = if(nrow(.x$sig) > 0) min(.x$sig$padj, na.rm = TRUE) else NA_real_
    ))

  # Compile comprehensive report data package
  list(
    # Core analysis results
    dds = dds,
    res.l.all = res.l.all,
    contrast_list = contrast_list,
    coldata = coldata,
    config_snapshot = config,

    # Summary information for overview
    metadata = analysis_metadata,
    comparison_summary = comparison_summary,

    # GSEA results if available
    gsea_results = gsea_results,

    # Generation timestamp
    generated_at = Sys.time()
  ) %>%
    {cli_inform("Report data compiled successfully"); .}
}

#' Save report data to technical directory
save_report_data <- function(experiment_name, report_data) {

  # Ensure proper directory structure exists
  ensure_experiment_outputs(experiment_name)

  # Determine output path following new directory structure
  output_dir <- here("experiments", experiment_name, "outputs", "technical", "R")
  output_file <- file.path(output_dir, paste0(experiment_name, "_data.RDS"))

  report_data %>%
    write_rds(output_file)

  cli_inform("Report data saved to: {output_file}")
  return(output_file)
}

#' Load report data from technical directory
load_report_data <- function(experiment_name) {

  input_file <- here("experiments", experiment_name, "outputs", "technical", "R",
                    paste0(experiment_name, "_data.RDS"))

  if (!file.exists(input_file)) {
    cli_abort("Report data not found: {input_file}")
  }

  cli_inform("Loading report data from: {input_file}")
  read_rds(input_file)
}

# ============================================================================= #
# INDIVIDUAL REPORT GENERATORS ----
# ============================================================================= #

#' Generate overview landing page
#'
#' Creates main overview.html file at outputs root for collaborator access
#'
#' @param experiment_name Experiment identifier
#' @param report_data Compiled report data from compile_report_data()
#' @param output_formats Vector of formats (default: "html")
generate_overview_page <- function(experiment_name, report_data, output_formats = "html") {

  cli_inform("Generating overview page...")

  template_path <- here("templates", "reports", "overview.qmd")
  if (!file.exists(template_path)) {
    cli_warn("Overview template not found: {template_path}")
    return(NULL)
  }

  # Output directly to outputs/ directory (not reports/)
  output_dir <- here("experiments", experiment_name, "outputs")

  output_formats %>%
    keep(~ .x == "html") %>%  # Only HTML for overview
    map(~ {
      output_file <- file.path(output_dir, "overview.html")

      success <- render_quarto_template(
        template_path = template_path,
        output_file = output_file,
        output_format = "html",
        params = list(
          experiment_name = experiment_name,
          data_dir = normalizePath(file.path(output_dir, "technical", "R"))
        )
      )

      if (success) {
        cli_inform("Generated overview: {output_file}")
        list(overview = output_file)
      } else {
        list()
      }
    }) %>%
    flatten()
}

#' Generate experimental summary report
#'
#' Creates comprehensive experimental overview with QC plots
#'
#' @param experiment_name Experiment identifier
#' @param report_data Compiled report data
#' @param output_formats Vector of formats (default: c("html", "pdf"))
generate_experiment_summary_report <- function(experiment_name, report_data, output_formats = c("html", "pdf")) {

  cli_inform("Generating experimental summary report...")

  template_path <- here("templates", "reports", "experimental_summary.qmd")
  if (!file.exists(template_path)) {
    cli_warn("Experimental summary template not found: {template_path}")
    return(list())
  }

  # Output to reports/html/ and reports/pdf/
  output_dir <- here("experiments", experiment_name, "outputs")
  data_dir <- file.path(output_dir, "technical", "R")

  output_formats %>%
    set_names(output_formats) %>%
    imap(~ {
      format_dir <- file.path(output_dir, "reports", .x)
      if (!dir.exists(format_dir)) {
        dir.create(format_dir, recursive = TRUE)
      }

      output_file <- file.path(format_dir, paste0(experiment_name, "_experimental_summary.", .x))

      success <- render_quarto_template(
        template_path = template_path,
        output_file = output_file,
        output_format = .x,
        params = list(
          experiment_name = experiment_name,
          data_dir = normalizePath(data_dir)
        )
      )

      if (success) {
        cli_inform("Generated experimental summary ({.x}): {output_file}")
        list(output_file)
      } else {
        list()
      }
    }) %>%
    set_names(paste0("summary_", names(.))) %>%
    flatten()
}

#' Generate individual comparison report
#'
#' Creates detailed DE analysis report for single comparison
#'
#' @param experiment_name Experiment identifier
#' @param comparison_name Name of specific comparison
#' @param report_data Compiled report data
#' @param output_formats Vector of formats (default: c("html", "pdf"))
generate_comparison_report <- function(experiment_name, comparison_name, report_data, output_formats = c("html", "pdf")) {

  cli_inform("Generating comparison report: {comparison_name}")

  template_path <- here("templates", "reports", "comparison_report.qmd")
  if (!file.exists(template_path)) {
    cli_warn("Comparison report template not found: {template_path}")
    return(list())
  }

  # Verify comparison exists
  if (!comparison_name %in% names(report_data$res.l.all)) {
    cli_warn("Comparison '{comparison_name}' not found in results")
    return(list())
  }

  output_dir <- here("experiments", experiment_name, "outputs")
  data_dir <- file.path(output_dir, "technical", "R")

  output_formats %>%
    set_names(output_formats) %>%
    imap(~ {
      format_dir <- file.path(output_dir, "reports", .x)
      if (!dir.exists(format_dir)) {
        dir.create(format_dir, recursive = TRUE)
      }

      output_file <- file.path(format_dir, paste0(comparison_name, "_report.", .x))

      success <- render_quarto_template(
        template_path = template_path,
        output_file = output_file,
        output_format = .x,
        params = list(
          experiment_name = experiment_name,
          comparison_name = comparison_name,
          data_dir = normalizePath(data_dir)
        )
      )

      if (success) {
        cli_inform("Generated {comparison_name} report ({.x}): {output_file}")
        list(output_file)
      } else {
        list()
      }
    }) %>%
    set_names(paste0(comparison_name, "_", names(.))) %>%
    flatten()
}

#' Generate interactive heatmap page
#'
#' Creates separate interactive heatmap for all genes
#'
#' @param experiment_name Experiment identifier
#' @param report_data Compiled report data
#' @param output_formats Vector of formats (default: "html")
generate_interactive_heatmap <- function(experiment_name, report_data, output_formats = "html") {

  cli_inform("Generating interactive heatmap...")

  template_path <- here("templates", "reports", "interactive_heatmap.qmd")
  if (!file.exists(template_path)) {
    cli_warn("Interactive heatmap template not found: {template_path}")
    return(list())
  }

  output_dir <- here("experiments", experiment_name, "outputs")
  data_dir <- file.path(output_dir, "technical", "R")

  output_formats %>%
    keep(~ .x == "html") %>%  # Only HTML for interactive content
    map(~ {
      format_dir <- file.path(output_dir, "reports", "html")
      if (!dir.exists(format_dir)) {
        dir.create(format_dir, recursive = TRUE)
      }

      output_file <- file.path(format_dir, paste0(experiment_name, "_interactive_heatmap.html"))

      success <- render_quarto_template(
        template_path = template_path,
        output_file = output_file,
        output_format = "html",
        params = list(
          experiment_name = experiment_name,
          data_dir = normalizePath(data_dir)
        )
      )

      if (success) {
        cli_inform("Generated interactive heatmap: {output_file}")
        list(interactive_heatmap = output_file)
      } else {
        list()
      }
    }) %>%
    flatten()
}

#' Generate all reports for an experiment
#'
#' Master function that creates complete report suite using tidyverse patterns
#'
#' @param experiment_name Experiment identifier
#' @param report_data Compiled report data from compile_report_data()
#' @param output_formats Vector of formats (default: c("html", "pdf"))
generate_all_reports <- function(experiment_name, report_data, output_formats = c("html", "pdf")) {

  cli_inform("Generating complete report suite for {experiment_name}")
  cli_inform("Output formats: {paste(output_formats, collapse = ', ')}")

  # Generate all report types using tidyverse patterns
  all_generated_files <- list()

  # 1. Overview page (HTML only)
  if ("html" %in% output_formats) {
    cli_inform("Generating overview reports...")
    overview_files <- generate_overview_page(experiment_name, report_data, "html")
    all_generated_files <- c(all_generated_files, overview_files)
  }

  # 2. Experimental summary report
  cli_inform("Generating summary reports...")
  summary_files <- generate_experiment_summary_report(experiment_name, report_data, output_formats)
  all_generated_files <- c(all_generated_files, summary_files)

  # 3. Interactive heatmap (HTML only)
  if ("html" %in% output_formats) {
    cli_inform("Generating heatmap reports...")
    heatmap_files <- generate_interactive_heatmap(experiment_name, report_data, "html")
    all_generated_files <- c(all_generated_files, heatmap_files)
  }

  # 4. Individual comparison reports
  cli_inform("Generating comparison reports...")
  comparison_files <- names(report_data$res.l.all) %>%
    set_names(names(report_data$res.l.all)) %>%
    map(~ generate_comparison_report(experiment_name, .x, report_data, output_formats)) %>%
    flatten()
  all_generated_files <- c(all_generated_files, comparison_files)

  cli_inform("Report generation complete: {length(all_generated_files)} files generated")

  # Log generated files for debugging using walk
  all_generated_files %>%
    iwalk(~ cli_inform("  {.y}: {.x}"))

  all_generated_files
}

# ============================================================================= #
# UTILITY FUNCTIONS ----
# ============================================================================= #

#' Safe Quarto template rendering with error handling
#'
#' Renders Quarto template with parameters and proper error handling
#'
#' @param template_path Path to .qmd template file
#' @param output_file Full path for output file
#' @param output_format Format ("html" or "pdf")
#' @param params List of parameters to pass to template
#'
#' @return TRUE if successful, FALSE if failed
render_quarto_template <- function(template_path, output_file, output_format, params) {

  tryCatch({

    # Ensure proper directory structure exists (if experiment_name is in params)
    if ("experiment_name" %in% names(params)) {
      ensure_experiment_outputs(params$experiment_name)
    }

    # Ensure output directory exists
    output_file %>%
      dirname() %>%
      {if (!dir.exists(.)) dir.create(., recursive = TRUE)}

    # Use system call for Quarto rendering (R package has CLI path issues)
    # The R quarto package causes "--output option cannot specify a relative or absolute path" errors
    if (FALSE && requireNamespace("quarto", quietly = TRUE)) {

      # Use absolute paths for better reliability
      abs_template_path <- normalizePath(template_path)
      abs_output_file <- normalizePath(dirname(output_file))
      output_filename <- basename(output_file)

      quarto::quarto_render(
        input = abs_template_path,
        output_file = file.path(abs_output_file, output_filename),
        output_format = output_format,
        execute_params = params,
        quiet = FALSE  # Enable output for debugging
      )

    } else {

      # Fallback to system call with temporary parameter file
      param_file <- tempfile(fileext = ".yml")
      params %>% write_yaml(param_file)

      # Use absolute paths and render in template directory, then move
      abs_template_path <- normalizePath(template_path)
      template_dir <- dirname(abs_template_path)
      template_name <- basename(abs_template_path)
      abs_param_file <- normalizePath(param_file)

      # Get expected output filename (Quarto will create this in template directory)
      base_name <- tools::file_path_sans_ext(template_name)
      expected_output <- file.path(template_dir, paste0(base_name, ".", if(output_format == "html") "html" else "pdf"))

      cmd <- glue(
        "cd '{template_dir}' && ",
        "quarto render '{template_name}' ",
        "--to {output_format} ",
        "--execute-params '{abs_param_file}'"
      )

      result <- system(cmd, intern = FALSE)
      unlink(param_file)

      if (result != 0) {
        stop("Quarto render failed with exit code ", result)
      }

      # Move the rendered file to desired location
      if (file.exists(expected_output) && expected_output != output_file) {
        file.copy(expected_output, output_file, overwrite = TRUE)
        file.remove(expected_output)
      }
    }

    # Verify output file was created
    if (!file.exists(output_file)) {
      stop("Output file was not created: ", output_file)
    }

    TRUE

  }, error = function(e) {
    cli_warn("Failed to render {basename(output_file)}: {e$message}")
    FALSE
  })
}

# ============================================================================= #
# INTEGRATION HELPER FUNCTIONS ----
# ============================================================================= #

#' Check if reports need to be generated (for run_utils.R integration)
#'
#' Checks timestamp and file existence to determine if reports are up to date
#'
#' @param experiment_name Experiment identifier
#' @param report_data_file Path to report data RDS file
#'
#' @return List with status information
check_reports_status <- function(experiment_name, report_data_file = NULL) {

  report_data_file <- report_data_file %||%
    here("experiments", experiment_name, "outputs", "technical", "R",
         paste0(experiment_name, "_data.RDS"))

  # Check if report data exists
  if (!file.exists(report_data_file)) {
    return(list(status = "missing_data", needs_generation = TRUE))
  }

  # Check if overview.html exists
  overview_file <- here("experiments", experiment_name, "outputs", "overview.html")
  if (!file.exists(overview_file)) {
    return(list(status = "missing_reports", needs_generation = TRUE))
  }

  # Compare timestamps
  data_time <- file.info(report_data_file)$mtime
  overview_time <- file.info(overview_file)$mtime

  if (data_time > overview_time) {
    return(list(status = "outdated", needs_generation = TRUE))
  }

  list(status = "current", needs_generation = FALSE)
}

#' Generate reports wrapper for integration with analysis scripts
#'
#' Simplified interface for use in analysis_script.R templates
#'
#' @param experiment_name Experiment identifier
#' @param output_formats Vector of formats (default: c("html", "pdf"))
#' @param force_regenerate Force regeneration even if reports exist
run_report_generation <- function(experiment_name, output_formats = c("html", "pdf"), force_regenerate = FALSE) {

  cli_h1("Report Generation")

  # Check if regeneration needed
  if (!force_regenerate) {
    status <- check_reports_status(experiment_name)
    if (!status$needs_generation) {
      cli_inform("Reports are up to date (use --force to regenerate)")
      return(invisible(TRUE))
    }
    cli_inform("Reports need regeneration: {status$status}")
  }

  # Load report data
  report_data <- tryCatch({
    load_report_data(experiment_name)
  }, error = function(e) {
    cli_abort("Failed to load report data. Run analysis first: {e$message}")
  })

  # Generate all reports
  generated_files <- generate_all_reports(experiment_name, report_data, output_formats)

  if (length(generated_files) > 0) {
    cli_inform("Report generation completed successfully")
    cli_inform("Main overview: experiments/{experiment_name}/outputs/overview.html")
    invisible(TRUE)
  } else {
    cli_warn("No reports were generated")
    invisible(FALSE)
  }
}