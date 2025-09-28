#!/usr/bin/env Rscript

# RNAseq pipeline run script helpers
# Creation date: 2025-09-04

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(DESeq2)
  library(tximeta)
  library(glue)
  library(yaml)
  library(cli)
  library(here)
})

here::i_am("functions/run_utils.R")

#' List available experiments
#' @return character vector of experiment names
list_experiments <- function() {
  experiments_dir <- here("experiments")
  if (!dir.exists(experiments_dir)) {
    cli_abort("Experiments directory not found: {experiments_dir}")
  }

  experiments <- list.dirs(experiments_dir, full.names = FALSE, recursive = FALSE)
  experiments <- experiments[experiments != "TEMPLATE"]
  experiments <- experiments[experiments != ""]

  if (length(experiments) == 0) {
    cli_warn("No experiments found in {experiments_dir}")
  }

  return(experiments)
}

#' Validate experiment exists
#' @param experiment_name Name of the experiment
#' @return logical
validate_experiment <- function(experiment_name) {
  available_experiments <- list_experiments()

  if (!experiment_name %in% available_experiments) {
    cli_abort("Experiment '{experiment_name}' not found. Available experiments: {paste(available_experiments, collapse = ', ')}")
  }

  experiment_dir <- here("experiments", experiment_name)
  analysis_script <- file.path(experiment_dir, "analysis_script.R")

  if (!file.exists(analysis_script)) {
    cli_abort("Analysis script not found: {analysis_script}")
  }

  return(TRUE)
}

#' Run a single experiment with specified pipeline step
#' @param experiment_name Name of the experiment
#' @param pipeline_step One of: full, dds-only, analysis-only, reports-only
#' @param verbose Print detailed output
#' @return logical indicating success
run_single_experiment <- function(experiment_name, pipeline_step = "full", verbose = TRUE) {
  validate_experiment(experiment_name)

  valid_steps <- c("full", "dds-only", "analysis-only", "reports-only")
  if (!pipeline_step %in% valid_steps) {
    cli_abort("Invalid pipeline step '{pipeline_step}'. Must be one of: {paste(valid_steps, collapse = ', ')}")
  }

  experiment_dir <- here("experiments", experiment_name)
  analysis_script <- file.path(experiment_dir, "analysis_script.R")

  if (verbose) {
    cli_h1("Running {experiment_name} - Step: {pipeline_step}")
  }

  # Set environment variable for pipeline step
  Sys.setenv(PIPELINE_STEP = pipeline_step)

  # Change to experiment directory
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(experiment_dir)

  # Source the analysis script
  tryCatch({
    source(analysis_script, echo = verbose)
    if (verbose) {
      cli_alert_success("Completed {experiment_name} - {pipeline_step}")
    }
    return(TRUE)
  }, error = function(e) {
    cli_alert_danger("Failed {experiment_name} - {pipeline_step}: {e$message}")
    return(FALSE)
  })
}

#' Check if experiment outputs exist for a given pipeline step
#' @param experiment_name Name of the experiment
#' @param pipeline_step Pipeline step to check
#' @return logical indicating if outputs exist
check_experiment_outputs <- function(experiment_name, pipeline_step) {
  experiment_dir <- here("experiments", experiment_name, "outputs")

  if (!dir.exists(experiment_dir)) {
    return(FALSE)
  }

  result <- switch(as.character(pipeline_step),
    "dds-only" = {
      dds_file <- file.path(experiment_dir, paste0(experiment_name, "_dds.RDS"))
      file.exists(dds_file)
    },
    "analysis-only" = {
      results_file <- file.path(experiment_dir, "results", paste0(experiment_name, "_results.RDS"))
      file.exists(results_file)
    },
    "reports-only" = {
      reports_dir <- file.path(experiment_dir, "reports")
      dir.exists(reports_dir) && length(list.files(reports_dir, pattern = "\\.(pdf|html)$")) > 0
    },
    "full" = {
      check_experiment_outputs(experiment_name, "dds-only") &&
      check_experiment_outputs(experiment_name, "analysis-only") &&
      check_experiment_outputs(experiment_name, "reports-only")
    },
    # Default case
    {
      cli_warn("Unknown pipeline step: {pipeline_step}")
      FALSE
    }
  )

  return(result)
}

#' Get experiment status summary
#' @param experiment_name Name of the experiment (optional, if NULL checks all)
#' @return data.frame with experiment status
get_experiment_status <- function(experiment_name = NULL) {
  if (is.null(experiment_name)) {
    experiments <- list_experiments()
  } else {
    experiments <- experiment_name
  }

  status_df <- map_dfr(experiments, function(exp) {
    data.frame(
      experiment = exp,
      dds_exists = check_experiment_outputs(exp, "dds-only"),
      analysis_exists = check_experiment_outputs(exp, "analysis-only"),
      reports_exist = check_experiment_outputs(exp, "reports-only"),
      stringsAsFactors = FALSE
    )
  })

  return(status_df)
}

#' Print experiment status in a nice format
#' @param experiment_name Name of the experiment (optional, if NULL shows all)
print_experiment_status <- function(experiment_name = NULL) {
  status_df <- get_experiment_status(experiment_name)

  if (nrow(status_df) == 0) {
    cli_inform("No experiments found.")
    return(invisible(NULL))
  }

  cli_h2("Experiment Status")

  for (i in seq_len(nrow(status_df))) {
    exp <- status_df$experiment[i]
    dds_status <- if (status_df$dds_exists[i]) cli::col_green("✓") else cli::col_red("✗")
    analysis_status <- if (status_df$analysis_exists[i]) cli::col_green("✓") else cli::col_red("✗")
    reports_status <- if (status_df$reports_exist[i]) cli::col_green("✓") else cli::col_red("✗")

    cli_inform("{exp}: DDS {dds_status} | Analysis {analysis_status} | Reports {reports_status}")
  }
}

#' Parse command line arguments for run scripts
#' @param args Command line arguments vector
#' @return list with parsed arguments
parse_run_args <- function(args) {
  # Initialize defaults
  result <- list(
    experiment = NULL,
    pipeline_step = "full",
    verbose = TRUE,
    help = FALSE,
    status = FALSE,
    force = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[i]

    if (arg %in% c("-h", "--help")) {
      result$help <- TRUE
    } else if (arg %in% c("-s", "--status")) {
      result$status <- TRUE
    } else if (arg %in% c("-f", "--force")) {
      result$force <- TRUE
    } else if (arg %in% c("-q", "--quiet")) {
      result$verbose <- FALSE
    } else if (arg %in% c("--dds-only", "--dds")) {
      result$pipeline_step <- "dds-only"
    } else if (arg %in% c("--analysis-only", "--analysis")) {
      result$pipeline_step <- "analysis-only"
    } else if (arg %in% c("--reports-only", "--reports")) {
      result$pipeline_step <- "reports-only"
    } else if (arg %in% c("--full")) {
      result$pipeline_step <- "full"
    } else if (!startsWith(arg, "-")) {
      # Positional argument - experiment name
      if (is.null(result$experiment)) {
        result$experiment <- arg
      }
    }

    i <- i + 1
  }

  return(result)
}

#' Print usage information for run scripts
#' @param script_name Name of the script (for help text)
print_run_usage <- function(script_name = "run_experiment.R") {
  cli_h1("RNA-seq Pipeline Runner")
  cli_text("Usage: Rscript {script_name} [experiment_name] [options]")
  cli_text("")
  cli_h3("Arguments:")
  cli_text("  experiment_name    Name of experiment to run (required for single experiment)")
  cli_text("")
  cli_h3("Pipeline Steps:")
  cli_text("  --full            Run complete pipeline (default)")
  cli_text("  --dds-only        Generate DDS object only")
  cli_text("  --analysis-only   Run analysis and generate results")
  cli_text("  --reports-only    Generate reports only")
  cli_text("")
  cli_h3("Options:")
  cli_text("  -s, --status      Show experiment status")
  cli_text("  -f, --force       Force rerun even if outputs exist")
  cli_text("  -q, --quiet       Reduce output verbosity")
  cli_text("  -h, --help        Show this help message")
  cli_text("")
  cli_h3("Examples:")
  cli_text("  Rscript {script_name} rsRDM_all")
  cli_text("  Rscript {script_name} rsRDM_all --dds-only")
  cli_text("  Rscript {script_name} --status")
}
