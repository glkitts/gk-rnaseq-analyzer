#!/usr/bin/env Rscript

# RNAseq DE compilation generation functions
# Creation date: 2025-09-04

library(tidyverse)
library(magrittr)
library(DESeq2)
library(tximeta)
library(glue)
library(yaml)
library(cli)
library(openxlsx2)
library(here)

here::i_am("functions/de_compilation.R")


# ============================================================================= #
# CREATE COMPILATIONS ----
# ============================================================================= #

#' Create compiled master table using full_join approach
#'
#' @param res.l.all Named list of comparison results from workflow
#' @param dds DESeqDataSet object for normalized counts
#' @param config Configuration list with analysis parameters
#'
#' @return List with allHits_compiled and sigHits_compiled tibbles
compile_master_table <- function(res.l.all, dds) {
  # Rename comparison columns and full join for allHits_compiled
  allHits_compiled <- res.l.all %>%
    imap(function(comparison_res, comparison_name) {
      comparison_res$all %>%
        dplyr::rename(
          !!paste0(comparison_name, ".log2FC") := log2FoldChange,!!paste0(comparison_name, ".padj") := padj
        )
    }) %>%
    purrr::reduce(full_join)
  
  # Add normalized counts
  nc_data <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("original_names") %>%
    separate_wider_delim(
      original_names,
      delim = "|",
      names = c("Label", NA, NA),
      too_few = "align_start",
      too_many = "drop"
    ) %>%
    # select(-original_names) %>%
    distinct(Label, .keep_all = TRUE) %>%
    rename_with( ~ paste0("nc_", .x), -Label)
  
  allHits_compiled <- allHits_compiled %>%
    left_join(nc_data)
  
  return(allHits_compiled)
}

compile_significant_hits <- function(res.l.all) {
  sigHits_compiled <- res.l.all %>%
    imap(function(comparison_res, comparison_name) {
      comparison_res$sig %>%
        dplyr::rename(
          !!paste0(comparison_name, ".log2FC") := log2FoldChange,
          !!paste0(comparison_name, ".padj") := padj
          # !!paste0(comparison_name, ".baseMean") := baseMean
        )
    }) %>%
    purrr::reduce(full_join)
  
  return(sigHits_compiled)
  
}


create_pathway_subsets <- function(sigHits_compiled, genesets) {
  ## Just Label and log2FC to avoid join overlap
  d.trim <- sigHits_compiled %>%
    select(Label, contains("log2FC"))
  
  ## For each geneset, attach sigHits
  d.sig.genesets <- genesets %>%
    lapply(function(x) {
      x %>% left_join(d.trim)
    })
  
  return(d.sig.genesets)
}

#' Create README/Info sheet for Excel DEcompilation
#'
#' @param res.l.all Named list of comparison results
#' @param config Configuration list
#' @param experiment_name String identifier
#' @param dds DESeqDataSet object
#' @param coldata Sample metadata
#'
#' @return Data frame formatted for Excel README sheet
create_readme_sheet <- function(res.l.all,
                                config,
                                experiment_name,
                                dds = NULL,
                                coldata = NULL) {
  # Initialize info sections
  info_data <- tibble(Section = character(),
                      Parameter = character(),
                      Value = character())
  
  # ANALYSIS OVERVIEW
  info_data <- bind_rows(info_data,
                         tibble(
                           Section = "ANALYSIS OVERVIEW",
                           Parameter = c(
                             "Experiment Name",
                             "Analysis Date",
                             "Number of Comparisons",
                             "Total Genes Analyzed",
                             "R Version",
                             "DESeq2 Version"
                           ),
                           Value = c(
                             experiment_name,
                             format(Sys.Date(), "%Y-%m-%d"),
                             length(res.l.all),
                             if (!is.null(dds))
                               nrow(dds)
                             else
                               "N/A",
                             paste(R.version$major, R.version$minor, sep = "."),
                             as.character(packageVersion("DESeq2"))
                           )
                         ))
  
  # ANALYSIS THRESHOLDS
  info_data <- bind_rows(info_data,
                         tibble(
                           Section = "ANALYSIS THRESHOLDS",
                           Parameter = c(
                             "Significance Threshold (padj)",
                             "Log2 Fold Change Threshold (DE)",
                             "Minimum Counts Filter",
                             "LFC Shrinkage Method"
                           ),
                           Value = c(
                             config$analysis$significance_threshold,
                             config$analysis$log2FC_threshold,
                             config$analysis$min_counts,
                             config$analysis$lfcShrinkage_method
                           )
                         ))
  
  # SAMPLE INFORMATION
  if (!is.null(coldata)) {
    sample_summary <- coldata %>%
      group_by(across(any_of(c(
        "strain", "condition", "group"
      )))) %>%
      summarise(n_samples = n(), .groups = "drop") %>%
      unite("group_info", everything(), sep = " | ") %>%
      pull(group_info) %>%
      paste(collapse = "; ")
    
    info_data <- bind_rows(info_data,
                           tibble(
                             Section = "SAMPLE INFORMATION",
                             Parameter = c("Total Samples", "Sample Groups"),
                             Value = c(nrow(coldata), sample_summary)
                           ))
  }
  
  # COMPARISON RESULTS SUMMARY
  comparison_summary <- res.l.all %>%
    imap_dfr(function(comp_data, comp_name) {
      tibble(
        comparison = comp_name,
        total_genes = nrow(comp_data$all),
        sig_genes = nrow(comp_data$sig),
        de_genes = nrow(comp_data$DE)
      )
    })
  
  for (i in 1:nrow(comparison_summary)) {
    comp <- comparison_summary[i, ]
    info_data <- bind_rows(info_data,
                           tibble(
                             Section = "RESULTS SUMMARY",
                             Parameter = paste0(comp$comparison, " - Total/Sig/DE genes"),
                             Value = paste0(comp$total_genes, " / ", comp$sig_genes, " / ", comp$de_genes)
                           ))
  }
  
  # REFERENCE INFORMATION
  ref_info <- c(
    "Genome Assembly" = if (!is.null(config$reference$genome_assembly))
      config$reference$genome_assembly
    else
      "Not specified",
    "Genome Accession" = if (!is.null(config$reference$genome_accession))
      config$reference$genome_accession
    else
      "Not specified",
    "Organism" = if (!is.null(config$reference$organism))
      config$reference$organism
    else
      "Not specified",
    "Annotation Version" = if (!is.null(config$reference$annotation_version))
      config$reference$annotation_version
    else
      "Not specified",
    "Annotation Database" = basename(config$paths$annotation_db),
    "Salmon Base Path" = config$paths$salmon_base,
    "Metadata File" = basename(config$paths$metadata_file)
  )
  
  for (param_name in names(ref_info)) {
    info_data <- bind_rows(
      info_data,
      tibble(
        Section = "REFERENCE INFORMATION",
        Parameter = param_name,
        Value = ref_info[[param_name]]
      )
    )
  }
  
  # SHEET DESCRIPTIONS
  sheet_descriptions <- c(
    "sigHits_compiled: Genes significant (padj ≤ threshold) in ≥1 comparison",
    "pathway sheets: Genes from sigHits filtered by functional pathways",
    "comparison_all: All analyzed genes for individual comparisons",
    "comparison_DE: Differentially expressed genes (padj + |log2FC| thresholds)"
  )
  
  for (desc in sheet_descriptions) {
    info_data <- bind_rows(
      info_data,
      tibble(
        Section = "SHEET DESCRIPTIONS",
        Parameter = str_extract(desc, "^[^:]+"),
        Value = str_extract(desc, "(?<=: ).*")
      )
    )
  }
  
  # COLUMN DEFINITIONS
  column_definitions <- c(
    "Label: Gene identifier (VC####)",
    "Product: Gene product/function description",
    "log2FoldChange: Log2 fold change (treatment vs control)",
    "padj: Adjusted p-value (Benjamini-Hochberg)",
    "baseMean: Mean normalized expression across samples",
    "nc_: Normalized count columns for individual samples"
  )
  
  for (def in column_definitions) {
    info_data <- bind_rows(
      info_data,
      tibble(
        Section = "COLUMN DEFINITIONS",
        Parameter = str_extract(def, "^[^:]+"),
        Value = str_extract(def, "(?<=: ).*")
      )
    )
  }
  return(info_data)
}



create_de_compilation <- function(res.l.all,
                                  annotations,
                                  genesets = NULL,
                                  comparison_types = c("DE", "all"),
                                  config = NULL,
                                  experiment_name = NULL,
                                  dds = NULL,
                                  coldata = NULL) {
  # Start with sigHits as primary sheet
  excel.list <- list()
  
  sigHits_compiled <- res.l.all %>%
    compile_significant_hits() %>%
    merge_annotations(annotations) %>%
    relocate(contains(".padj"), .after = last_col())
  
  # Conditional naming and ordering based on number of comparisons
  if (length(res.l.all) > 1) {
    # Multiple comparisons: sigHits_compiled first
    excel.list$sigHits_compiled <- sigHits_compiled
  } else {
    # Single comparison: sigHits, then individual sheets, then pathways
    comparison_name <- names(res.l.all)[1]
    excel.list[[paste0(comparison_name, "_sigHits")]] <- sigHits_compiled
    
    # Add individual sheets right after sigHits
    single_comparison_sheets <- comparison_types %>%
      keep( ~ .x %in% names(res.l.all[[comparison_name]])) %>%
      set_names(paste0(comparison_name, "_", str_to_upper(.))) %>%
      map( ~ res.l.all[[comparison_name]][[.x]])
    
    excel.list <- c(excel.list, single_comparison_sheets)
  }
  
  # Add README sheet
  if (!is.null(config) && !is.null(experiment_name)) {
    readme_sheet <- create_readme_sheet(res.l.all, config, experiment_name, dds, coldata)
    excel.list$README <- readme_sheet
  }
  
  # Add pathway subsets if provided
  if (!is.null(genesets)) {
    d.sig.genesets <- sigHits_compiled %>%
      create_pathway_subsets(genesets)
    
    if (length(d.sig.genesets) > 0) {
      excel.list <- excel.list %>%
        append(d.sig.genesets)
    }
  }
  
  # Add individual comparison sheets (only for multiple comparisons)
  if (length(res.l.all) > 1) {
    individual_sheets <- res.l.all %>%
      imap(function(comparison_data, comparison_name) {
        comparison_types %>%
          keep( ~ .x %in% names(comparison_data)) %>%
          set_names(paste0(comparison_name, "_", str_to_upper(.))) %>%
          map( ~ comparison_data[[.x]])
      }) %>%
      flatten()
    
    excel.list <- c(excel.list, individual_sheets)
  }
  
  
  
  cli_inform("Excel structure: {length(excel.list)} sheets")
  return(excel.list)
}


