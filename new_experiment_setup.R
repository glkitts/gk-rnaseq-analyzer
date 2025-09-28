#!/usr/bin/env Rscript
# New Experiment Setup Script
# Creates a new experiment folder with customized analysis_script.R

library(here)

# Get experiment name from command line argument
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("❌ Usage: Rscript new_experiment_setup.R <experiment_name>\n")
  cat("   Example: Rscript new_experiment_setup.R rsRDM_iron_study\n")
  quit(status = 1)
}

experiment_name <- args[1]

# Validate experiment name (basic checks)
if (grepl("[^a-zA-Z0-9_-]", experiment_name)) {
  cat("❌ Experiment name should only contain letters, numbers, underscores, and hyphens\n")
  quit(status = 1)
}

# Create experiment directory
experiment_dir <- here("experiments", experiment_name)

if (dir.exists(experiment_dir)) {
  cat(paste0("❌ Experiment directory already exists: ", experiment_dir, "\n"))
  cat("   Choose a different name or delete the existing directory.\n")
  quit(status = 1)
}

cat(paste0("CREATING: ", experiment_name, "\n"))

dir.create(experiment_dir, recursive = TRUE)

# Create output directory structure following collaborator-friendly design
output_dirs <- c(
  "outputs",
  "outputs/reports",
  "outputs/reports/html",
  "outputs/reports/pdf",
  "outputs/plots",
  "outputs/data_tables",
  "outputs/data_tables/individual_comparisons",
  "outputs/technical",
  "outputs/technical/R",
  "outputs/technical/metadata",
  "outputs/technical/logs"
)

for (dir in output_dirs) {
  dir.create(file.path(experiment_dir, dir), recursive = TRUE, showWarnings = FALSE)
}

# Read template analysis script
template_path <- here("experiments", "TEMPLATE", "analysis_script.R")

if (!file.exists(template_path)) {
  cat(paste0("❌ Template not found: ", template_path, "\n"))
  cat("   Make sure you have a TEMPLATE folder with analysis_script.R\n")
  quit(status = 1)
}

template_content <- readLines(template_path, warn = FALSE)

# Replace experiment name and date in comment lines
current_date <- format(Sys.Date(), "%Y-%m-%d")
template_content <- gsub(
  "# RNA-seq Analysis: experiment_template",
  paste0("# RNA-seq Analysis: ", experiment_name),
  template_content,
  fixed = TRUE
)
template_content <- gsub(
  "# Creation date: CREATION_DATE",
  paste0("# Creation date: ", current_date),
  template_content,
  fixed = TRUE
)

template_content <- gsub(
  "# Experiment name - CHANGE THIS!!",
  paste0("# Experiment name (must match experiment folder name!)"),
  template_content,
  fixed = TRUE
)


# Replace experiment name in the variable assignment
template_content <- gsub(
  'experiment_name <- "experiment_template"',
  paste0('experiment_name <- "', experiment_name, '"'),
  template_content,
  fixed = TRUE
)

# Write customized analysis script
analysis_script_path <- file.path(experiment_dir, "analysis_script.R")
writeLines(template_content, analysis_script_path)

cat("READY: Edit analysis_script.R and run analysis\n")
