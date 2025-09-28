#!/usr/bin/env Rscript
# Batch Processing Script for All Experiments
# Runs analysis pipeline across multiple experiments with selective regeneration

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(optparse)
  library(cli)
})

# Establish project root
here::i_am("run_all.R")

# Source run utilities
source(here("functions", "run_utils.R"))

# Define command line options
option_list <- list(
  make_option(c("-e", "--experiments"),
              type = "character",
              default = NULL,
              help = "Comma-separated list of experiments to run (default: all)",
              metavar = "EXP1,EXP2"),

  make_option(c("--dds-only"),
              action = "store_true",
              default = FALSE,
              help = "Generate DDS objects only"),

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
              help = "Show experiment status for all"),

  make_option(c("--stop-on-error"),
              action = "store_true",
              default = FALSE,
              help = "Stop if any experiment fails (default: continue)"),

  make_option(c("--dry-run"),
              action = "store_true",
              default = FALSE,
              help = "Show what would be run without executing")
)

# Parse arguments
opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nRNA-seq Pipeline Batch Runner\n\nRuns analysis pipeline across multiple experiments with selective regeneration options.",
  epilogue = paste(
    "\nExamples:",
    "  Rscript run_all.R                                    # Run all experiments (full pipeline)",
    "  Rscript run_all.R --reports-only                     # Regenerate all reports",
    "  Rscript run_all.R -e rsRDM_all,LG_pmB                # Run specific experiments",
    "  Rscript run_all.R -e rsRDM_all --dds-only            # DDS only for one experiment",
    "  Rscript run_all.R --status                           # Show status of all experiments",
    "  Rscript run_all.R --dry-run --analysis-only          # Preview what would run",
    "",
    sep = "\n  "
  )
)

opt <- parse_args(opt_parser)

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

#' Parse experiment list from option string
#' @param experiment_string Comma-separated experiment names
#' @return character vector of experiment names
parse_experiment_list <- function(experiment_string) {
  if (is.null(experiment_string)) {
    return(NULL)
  }

  experiments <- str_split(experiment_string, ",", simplify = TRUE) %>%
    str_trim() %>%
    .[. != ""]

  return(experiments)
}

#' Run experiments in batch (tidyverse style)
#' @param experiments Vector of experiment names (NULL for all)
#' @param pipeline_step Pipeline step to run
#' @param force Force rerun even if outputs exist
#' @param verbose Print detailed output
#' @param continue_on_error Continue if individual experiments fail
#' @param dry_run Show what would be run without executing
#' @return tibble with results for each experiment
run_batch_experiments <- function(experiments = NULL,
                                  pipeline_step = "full",
                                  force = FALSE,
                                  verbose = TRUE,
                                  continue_on_error = TRUE,
                                  dry_run = FALSE) {

  # Get list of experiments to run
  all_experiments <- list_experiments()

  if (is.null(experiments)) {
    target_experiments <- all_experiments
    if (verbose) {
      cli_inform("Running all {length(target_experiments)} experiments")
    }
  } else {
    # Validate specified experiments exist
    missing_experiments <- setdiff(experiments, all_experiments)
    if (length(missing_experiments) > 0) {
      cli_abort("Experiments not found: {paste(missing_experiments, collapse = ', ')}")
    }
    target_experiments <- experiments
    if (verbose) {
      cli_inform("Running {length(target_experiments)} specified experiments")
    }
  }

  if (length(target_experiments) == 0) {
    cli_warn("No experiments found to run")
    return(tibble())
  }

  # Create experiment plan
  step <- pipeline_step  # Store the parameter value to avoid name conflict
  experiment_plan <- tibble(
    experiment = target_experiments,
    pipeline_step = pipeline_step
  ) %>%
    mutate(
      outputs_exist = map_lgl(experiment, ~check_experiment_outputs(.x, step)),
      will_run = force | !outputs_exist,
      status = case_when(
        will_run ~ "planned",
        !will_run ~ "skip_existing",
        TRUE ~ "unknown"
      )
    )

  if (verbose) {
    cli_h1("Batch Processing Plan: {pipeline_step}")

    # Show summary
    planned_count <- sum(experiment_plan$status == "planned")
    skip_count <- sum(experiment_plan$status == "skip_existing")

    summary_parts <- c()
    if (planned_count > 0) summary_parts <- c(summary_parts, glue("{planned_count} will run"))
    if (skip_count > 0) summary_parts <- c(summary_parts, glue("{skip_count} will skip (outputs exist)"))

    cli_text("Plan: {paste(summary_parts, collapse = ', ')}")
    cli_text("")
  }

  if (dry_run) {
    cli_h2("Dry Run - Would Execute:")
    experiment_plan %>%
      filter(will_run) %>%
      pwalk(function(experiment, pipeline_step, ...) {
        cli_text("  {experiment} --{pipeline_step}")
      })
    return(experiment_plan)
  }

  # Execute experiments
  if (verbose && any(experiment_plan$will_run)) {
    cli_h2("Executing Experiments")
  }

  results <- experiment_plan %>%
    mutate(
      result = pmap_chr(list(experiment, will_run), function(exp, should_run) {
        if (!should_run) {
          if (verbose) cli_inform("⚬ Skipping {exp} - outputs exist")
          return("skipped")
        }

        if (verbose) {
          row_num <- which(experiment_plan$experiment == exp)
          total_to_run <- sum(experiment_plan$will_run)
          current_run <- sum(experiment_plan$will_run[1:row_num])
          cli_h3("Processing {exp} ({current_run}/{total_to_run})")
        }

        success <- tryCatch({
          run_single_experiment(
            experiment_name = exp,
            pipeline_step = step,
            verbose = verbose
          )
        }, error = function(e) {
          cli_alert_danger("Error in {exp}: {e$message}")
          FALSE
        })

        if (success) {
          if (verbose) cli_alert_success("✓ Completed {exp}")
          return("success")
        } else {
          if (verbose) cli_alert_danger("✗ Failed {exp}")
          if (!continue_on_error) {
            cli_abort("Stopping due to failure in {exp}")
          }
          return("failed")
        }
      })
    )

  # Print summary
  if (verbose) {
    cli_h2("Batch Processing Summary")

    successful_count <- sum(results$result == "success")
    failed_count <- sum(results$result == "failed")
    skipped_count <- sum(results$result == "skipped")

    if (successful_count > 0) cli_inform(cli::col_green(glue("✓ Successful: {successful_count}")))
    if (failed_count > 0) cli_inform(cli::col_red(glue("✗ Failed: {failed_count}")))
    if (skipped_count > 0) cli_inform(glue("⚬ Skipped: {skipped_count}"))

    # Show failed experiments if any
    failed_experiments <- results$experiment[results$result == "failed"]

    if (length(failed_experiments) > 0) {
      cli_text("Failed experiments: {paste(failed_experiments, collapse = ', ')}")
    }
  }

  return(results)
}

# Main execution
main <- function() {
  # Handle status request
  if (opt$status) {
    print_experiment_status()
    quit(status = 0)
  }

  # Determine pipeline step
  pipeline_step <- get_pipeline_step(opt)

  # Parse experiment list
  experiments <- parse_experiment_list(opt$experiments)

  # Run batch processing
  results <- run_batch_experiments(
    experiments = experiments,
    pipeline_step = pipeline_step,
    force = opt$force,
    verbose = !opt$quiet,
    continue_on_error = !opt$`stop-on-error`,
    dry_run = opt$`dry-run`
  )

  # Determine exit code
  if (opt$`dry-run`) {
    quit(status = 0)
  }

  failed_count <- sum(results$result == "failed")

  if (failed_count > 0) {
    quit(status = 1)
  } else {
    quit(status = 0)
  }
}

# Run main function
main()