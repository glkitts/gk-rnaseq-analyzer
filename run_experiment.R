#!/usr/bin/env Rscript
# Run Single Experiment Script
# Executes analysis pipeline for a single experiment with configurable steps

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(optparse)
  library(cli)
})

# Establish project root
here::i_am("run_experiment.R")

# Source run utilities
source(here("functions", "run_utils.R"))

# Define command line options
option_list <- list(
  make_option(c("-e", "--experiment"),
              type = "character",
              default = NULL,
              help = "Name of experiment to run (required)",
              metavar = "EXPERIMENT"),

  make_option(c("--dds-only"),
              action = "store_true",
              default = FALSE,
              help = "Generate DDS object only"),

  make_option(c("--analysis-only"),
              action = "store_true",
              default = FALSE,
              help = "Run analysis and generate results only"),

  make_option(c("--reports-only"),
              action = "store_true",
              default = FALSE,
              help = "Generate reports only"),

  make_option(c("-f", "--force"),
              action = "store_true",
              default = FALSE,
              help = "Force rerun even if outputs exist"),

  make_option(c("-q", "--quiet"),
              action = "store_true",
              default = FALSE,
              help = "Reduce output verbosity"),

  make_option(c("-s", "--status"),
              action = "store_true",
              default = FALSE,
              help = "Show experiment status")
)

# Parse arguments - allow positional argument for experiment name
opt_parser <- OptionParser(
  option_list = option_list,
  usage = "usage: %prog [experiment_name] [options]",
  description = paste(
    "\nRNA-seq Pipeline Single Experiment Runner",
    "\nExecutes analysis pipeline for a single experiment with configurable steps.",
    "\nThe experiment name can be provided as a positional argument or with --experiment.",
    sep = ""
  ),
  epilogue = paste(
    "\nExamples:",
    "  Rscript run_experiment.R rsRDM_all                   # Run full pipeline",
    "  Rscript run_experiment.R rsRDM_all --dds-only        # Generate DDS only",
    "  Rscript run_experiment.R -e rsRDM_all --reports-only # Generate reports only",
    "  Rscript run_experiment.R --status -e rsRDM_all       # Show experiment status",
    "",
    sep = "\n  "
  )
)

# Custom argument parsing to handle positional argument
args <- commandArgs(trailingOnly = TRUE)
positional_experiment <- NULL

# Extract positional argument (experiment name) if present and not starting with -
if (length(args) > 0 && !startsWith(args[1], "-")) {
  positional_experiment <- args[1]
  args <- args[-1]  # Remove positional arg for optparse
}

opt <- parse_args(opt_parser, args = args)

#' Determine pipeline step from options
#' @param opt Parsed options list
#' @return character pipeline step name
get_pipeline_step <- function(opt) {
  steps <- c(
    "dds-only" = opt$`dds-only`,
    "analysis-only" = opt$`analysis-only`,
    "reports-only" = opt$`reports-only`
  )

  active_steps <- names(steps)[steps]

  if (length(active_steps) > 1) {
    cli_abort("Multiple pipeline steps specified: {paste(active_steps, collapse = ', ')}. Please specify only one.")
  }

  if (length(active_steps) == 0) {
    return("full")
  }

  return(active_steps[1])
}

# Main execution
main <- function() {
  # Determine experiment name (positional takes precedence)
  experiment_name <- if (!is.null(positional_experiment)) {
    positional_experiment
  } else {
    opt$experiment
  }

  # Handle status request
  if (opt$status) {
    print_experiment_status(experiment_name)
    quit(status = 0)
  }

  # Validate experiment name is provided
  if (is.null(experiment_name)) {
    cli_abort("Error: Experiment name is required\n\nUse --help for usage information")
  }

  # Determine pipeline step
  pipeline_step <- get_pipeline_step(opt)
  verbose <- !opt$quiet
  force <- opt$force

  # Check if outputs already exist (unless forced)
  if (!force && check_experiment_outputs(experiment_name, pipeline_step)) {
    if (verbose) {
      cli_inform("Outputs already exist for {experiment_name} - {pipeline_step}")
      cli_inform("Use --force to rerun or specify a different pipeline step")
      cli_inform("Current status:")
      print_experiment_status(experiment_name)
    }
    quit(status = 0)
  }

  # Run the experiment
  success <- run_single_experiment(
    experiment_name = experiment_name,
    pipeline_step = pipeline_step,
    verbose = verbose
  )

  if (success) {
    if (verbose) {
      cli_h2("Final Status")
      print_experiment_status(experiment_name)
    }
    quit(status = 0)
  } else {
    quit(status = 1)
  }
}

# Run main function
main()
