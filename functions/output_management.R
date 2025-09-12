
# RNAseq analysis functions
# Creation date: 2025-09-04

library(tidyverse)
library(DESeq2)
library(tximeta)
library(glue)
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
    "outputs/reports", 
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
