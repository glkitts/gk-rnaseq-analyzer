#!/usr/bin/env Rscript

# RNA-seq Analysis: experiment_template
# Creation date: CREATION_DATE

# ============================================================================= #
# SETUP & CONFIG ----
# ============================================================================= #

library(tidyverse)
library(magrittr)
library(DESeq2)
library(tximeta)
library(yaml)
library(glue)
library(here)


# Experiment name (must match experiment folder name!)
experiment_name <- "experiment_template"

# Establish project root
here::i_am(paste0("experiments/", experiment_name, "/analysis_script.R"))

# Source enhanced functions

source(here("functions/dds_management.R"))
source(here("functions/differential_expression.R"))
source(here("functions/de_compilation.R"))
source(here("functions/report_generation.R"))
# source(here("functions/processing_helpers.R"))
source(here("functions/output_management.R"))

# ============================================================================= #
# EXPERIMENT CONFIG ----
# ============================================================================= #

output_label <- experiment_name
cat(glue("Running experiment: {experiment_name}\n"))

# Load global config using here()
config <- read_yaml(here("config", "global_settings.yaml"))

# Experiment-specific overrides (customize as needed)
# config$analysis$min_counts <- 1
# config$analysis$significance_threshold <- 0.01
# config$plotting$base_size <- 10

## Curated pathway subsets for DE compilation output
pathway_subsets.path <- here(config$paths$pathway_subsets)
pathway_subsets <- pathway_subsets.path %>% 
  readxl::excel_sheets() %>%
  purrr::set_names() %>%
  map(readxl::read_excel, path = pathway_subsets.path)

genesets <- read_rds(here(config$paths$genesets))

# Determine what processing steps to run
pipeline_step <- Sys.getenv("PIPELINE_STEP", "full")  # full, dds-only, analysis-only, reports-only
force_regenerate <- as.logical(Sys.getenv("FORCE_REGENERATE", "FALSE"))
cat(glue("Running {experiment_name} - Step: {pipeline_step}\n"))
if (force_regenerate) cat("Force regeneration enabled\n")

# ============================================================================= #
# COLDATA FILTERING LOGIC ----
# ============================================================================= #
# Customize this filtering logic for your experiment
# This section can be run interactively for testing

# Load master metadata
coldata.in <- readxl::read_excel(here(config$paths$metadata_file), sheet = "exp_metadata")

coldata <- coldata.in %>%
  filter(experiment_group %in% c(
    "iron"
    # "aki_rvvA",
    # "lb_rvvA"
  )) %>%
  filter(!strain %in% c(
    "dryhB"
    # "dfur",
    # "drvvA"
  )) %>%
  # filter(!sample %in% c("ID_dvxrB_3", "IR_dvxrB_2")) %>%
  mutate(
    strain = fct_relevel(strain, "WT"),
    condition = fct_relevel(condition, "iron_replete"),
    group = fct_inorder(sampleLabel),
    names = sample,
    files = file.path(
      here(config$paths$salmon_base),
      experiment_folder,
      sample_FolderName,
      "salmon_quants",
      "quant.sf"
    )
  ) %>%
  select(sample, strain, condition, group, files, temperature, names)

# coldata %>% view

# ============================================================================= #
# CONTRAST DEFINITIONS ----
# ============================================================================= #
# Define your comparisons - customize for your experiment

contrast_list <- list(
  c("group", "IR_drvvA", "IR_WT"),
  c("group", "ID_drvvA", "ID_WT"),
  c("group", "IR_dfur", "IR_WT"),
  c("group", "ID_dfur", "ID_WT"),
  c("group", "IR_dvxrB", "IR_WT"),
  c("group", "ID_dvxrB", "ID_WT")
  # Add more contrasts as needed
)

# ============================================================================= #
# LOAD ANNOTATIONS ----
# ============================================================================= #
# Loads default subsetted annotations
# Modify if different needed

# vc_anno <- read_rds(here(config$paths$annotation_db))
# 
# annotations <- vc_anno$all_unfiltered %>%
#   select(
#     Label,
#     Product,
#     Type,
#     Begin,
#     End,
#     Length,
#     any_of(c(
#       "Gene", "GeneNames", "PsortB", "SignalP_5.0"
#     )),
#     starts_with("GO terms"),
#     starts_with("yildiz_"),
#     any_of(c("Transcription Units")),
#     starts_with("TIGR"),
#     Locus_Tag_Old,
#     UniProtID,
#   ) %>% 
#   distinct(Label, .keep_all = T)
# 
# 
# annotations %>% write_rds(here(config$paths$annotation_subset))
annotations <- read_rds(here(config$paths$annotation_subset))

# ============================================================================= #
# DDS GENERATION ----
# ============================================================================= #

if (pipeline_step %in% c("full", "dds-only")) {
  # cat("Generating DDS object...\n")
  
  # Validate files exist
  missing_files <- coldata$files[!file.exists(coldata$files)]
  if (length(missing_files) > 0) {
    cat("ERROR: Missing salmon files:\n")
    cat(paste(missing_files, collapse = "\n"))
    stop("Missing salmon quantification files. Check paths.")
  }
  
  # cat(glue("Found {nrow(coldata)} samples for analysis\n"))
  
  # Create DDS object using enhanced functions
  dds <- create_dds_from_salmon(
    coldata = coldata, 
    design_formula = ~group,
    min_counts = config$analysis$min_counts,
    experiment_name = experiment_name
  )
  
  # Save DDS and metadata
  save_dds_with_metadata(dds, experiment_name, coldata, config)
  
}



# res <- generate_resLFC(dds, contrast_list[[1]])
# 
# res %>% create_result_subsets()
# # ============================================================================= #
# DIFFERENTIAL EXPRESSION ANALYSIS ----
# # ============================================================================= #

if (pipeline_step %in% c("full", "analysis-only")) {
  cat("Running differential expression analysis...\n")
  
  # Load DDS if not already in memory
  if (!exists("dds")) {
    dds <- load_dds(experiment_name)
  }
  
  # Create 'all' and 'DE' tibble per contrast
  res.l.all <- contrast_list %>%
    set_names(sapply(., format_comparison_name)) %>%
    lapply(function(contrast) {
      res.l <- generate_resLFC(
        dds = dds, 
        contrast = contrast, 
        shrink_method = config$analysis$lfcShrinkage_method
      ) %>%
        merge_annotations(annotations) %>%
        create_result_subsets(
          significance_threshold = config$analysis$significance_threshold,
          log2FC_threshold = config$analysis$log2FC_threshold
        )
    })
  
  csv.paths <- res.l.all %>%
    imap(function(res.l, comparison) {
      res.l %>%
        save_individual_csvs(
          experiment_name = experiment_name,
          comparison_name = comparison
          # types = c("all", "DE")
        )
    })
  
  # Master table data is included in consolidated report data (_data.RDS)
  d <- res.l.all %>%
    compile_master_table(dds) %>%
    merge_annotations(annotations)
  
  # Create DE compilation
  # Returns first sheet of compilation 
  # (usually sigHits_compiled)
  sigHits_compiled <- res.l.all %>%
    create_de_compilation(
      annotations = annotations,
      pathway_subsets = pathway_subsets,
      config = config,
      experiment_name = experiment_name,
      dds = dds,
      coldata = coldata
    ) %>% 
    save_de_compilation_excel(
      experiment_name = experiment_name
    )
  
  # =========================================================================== #
  # GSEA ANALYSIS ----
  # =========================================================================== #
  
  cat("Running GSEA analysis...\n")
  
  # Run GSEA for GO Biological Process
  if (!is.null(genesets$GObp_panther)) {
    gsea_gobp <- run_gsea_compilation(
      res.l.all = res.l.all,
      gene_sets = genesets$GObp_panther,
      gene_set_name = "GObp",
      subset_type = "all",
      experiment_name = experiment_name,
      save_excel = TRUE
    )
  }
  
  # Run GSEA for GO Molecular Function
  if (!is.null(genesets$GOmf_panther)) {
    gsea_gomf <- run_gsea_compilation(
      res.l.all = res.l.all,
      gene_sets = genesets$GOmf_panther,
      gene_set_name = "GOmf",
      subset_type = "all",
      experiment_name = experiment_name,
      save_excel = TRUE
    )
  }

  # Compile report data for publication-quality reports
  cat("Compiling report data...\n")

  # Collect GSEA results if available
  gsea_results <- list()
  if (exists("gsea_gobp")) gsea_results$GObp <- gsea_gobp
  if (exists("gsea_gomf")) gsea_results$GOmf <- gsea_gomf

  # Add experiment_name to config for report generation
  config$experiment_name <- experiment_name

  # Compile comprehensive report data
  report_data <- compile_report_data(
    dds = dds,
    res.l.all = res.l.all,
    contrast_list = contrast_list,
    coldata = coldata,
    config = config,
    gsea_results = if(length(gsea_results) > 0) gsea_results else NULL
  )

  # Save report data to technical directory
  save_report_data(experiment_name, report_data)
}

# ============================================================================= #
# REPORT GENERATION ----
# ============================================================================= #

if (pipeline_step %in% c("full", "reports-only")) {
  cat("Generating publication-quality reports...\n")

  # Generate complete report suite using new system
  success <- run_report_generation(
    experiment_name = experiment_name,
    output_formats = c("html", "pdf"),
    force_regenerate = force_regenerate
  )

  if (success) {
    cat(glue("\n✓ Report generation completed successfully!\n"))
    cat(glue("📊 Main overview: experiments/{experiment_name}/outputs/overview.html\n"))
    cat(glue("📁 All reports: experiments/{experiment_name}/outputs/reports/\n"))
  }
}
