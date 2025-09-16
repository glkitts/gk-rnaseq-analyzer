#!/usr/bin/env Rscript

# RNAseq differential analysis functions
# Creation date: 2025-09-04

library(tidyverse)
library(DESeq2)
library(tximeta)
library(glue)
library(yaml)
library(cli)
library(here)

here::i_am("functions/differential_expression.R")

# ============================================================================= #
# BASIC DIFFERENTIAL EXPRESSION FUNCTIONS ----
# ============================================================================= #

#' Format comparison name from DESeq2 contrast
#'
#' Creates consistent comparison names for file naming and results organization
#'
#' @param contrast Character vector of length 3 (factor, treatment, control)
#' @param separator String, separator for name parts (default: "_v_")
#'
#' @return String with formatted comparison name
format_comparison_name <- function(contrast, separator = "_v_") {
  contrast_name <- paste0(contrast[2], "_v_", contrast[3])
  
  return(contrast_name)
}

#' Generate shrunken log2 fold change results from DESeq2
#'
#' Extracts results from DDS object, performs LFC shrinkage, and parses gene labels.
#' Handles rowname parsing to extract VC#### gene identifiers from salmon output.
#'
#' @param dds DESeqDataSet object (processed with DESeq())
#' @param contrast Character vector of length 3 (factor, treatment, control)
#' @param shrink_method String, LFC shrinkage method (default: "ashr")
#' @param parse_rownames Logical, whether to parse gene labels from rownames (default: TRUE)
#' @param ... Additional arguments passed to lfcShrink()
#'
#' @return Data frame with Label, log2FoldChange, padj, baseMean columns
generate_resLFC <- function(dds,
                            contrast,
                            shrink_method = "ashr",
                            parse_rownames = T,
                            ...) {
  
  res <- results(dds, contrast = contrast)
  
  res.LFC <- as.data.frame(lfcShrink(dds, contrast = contrast, type = shrink_method))
  
  if (parse_rownames == T) {
    res.LFC <- res.LFC %>%
      select(log2FoldChange, padj, baseMean) %>%
      rownames_to_column(var = "original_names") %>%
      separate_wider_delim(
        original_names,
        delim = "|",
        names = c('Label', NA, NA),
        too_few = "align_start",
        too_many = "drop",
        cols_remove = F
      ) %>%
      arrange(desc(baseMean)) %>%
      relocate(original_names, .after = last_col()) %>%
      distinct(Label, .keep_all = T) %>%
      arrange(log2FoldChange)
  }
  
  return(res.LFC)
}

#' Create filtered subsets from differential expression results
#'
#' Creates standard subsets: all genes, significant hits (padj), and 
#' differentially expressed genes (padj + fold change thresholds)
#'
#' @param resLFC Data frame from generate_resLFC()
#' @param significance_threshold Numeric, padj threshold (default: from config)
#' @param log2FC_threshold Numeric, absolute log2FC threshold (default: from config)
#'
#' @return Named list with 'all', 'sig', and 'DE' data frames
create_result_subsets <- function(resLFC,
                                  significance_threshold = config$analysis$significance_threshold,
                                  log2FC_threshold = config$analysis$log2FC_threshold) {
  res.l <- list()
  res.l$all <- resLFC
  res.l$sig <- resLFC %>%
    filter(padj <= significance_threshold)
  res.l$DE <- resLFC %>%
    filter(padj <= significance_threshold &
             abs(log2FoldChange) >= log2FC_threshold)
  
  return(res.l)
}

#' Merge differential expression results with gene annotations
#'
#' Performs left join with annotation data and optionally relocates columns
#' for consistent output formatting. Validates row count preservation.
#'
#' @param res_df Data frame with differential expression results
#' @param annotations Data frame with gene annotations (joined on Label)
#' @param relocate_cols Logical, whether to reorder columns (default: TRUE)
#' @param cols_to_relocate Character vector, specific columns to relocate (optional)
#'
#' @return Data frame with merged annotations
merge_annotations <- function(res_df,
                              annotations,
                              relocate_cols = TRUE,
                              cols_to_relocate = NULL) {
  # Get original row count for validation
  original_rows <- nrow(res_df)
  
  # Perform left join
  merged_res <- res_df %>%
    left_join(annotations, na_matches = "never")
  
  # Validate row count
  if (nrow(merged_res) != original_rows) {
    cli_abort("Row count changed during annotation merge: {original_rows} -> {nrow(merged_res)}")
  }
  
  # Handle column relocation
  if (relocate_cols) {
    if (!is.null(cols_to_relocate)) {
      merged_res <- merged_res %>%
        relocate(any_of(cols_to_relocate))
    } else {
      merged_res <- merged_res %>%
        relocate(
          any_of("Label"),
          any_of(c("Product")),
          contains(c("log2FoldChange", "log2FC")),
          contains("padj"),
          any_of(c(
            "PsortB", "SignalP_5.0", "Gene", "GeneNames", "Type"
          )),
          starts_with("GO terms"),
          starts_with("yildiz_"),
          any_of(c("Begin", "End", "Length")),
          any_of("Transcription Units"),
          starts_with("TIGR"),
          any_of(c("Locus_Tag_Old", "UniProtID")),
          starts_with("nc_")
        )
    }
  }
  
  return(merged_res)
}

# ============================================================================= #
# GENE SET ENRICHMENT ANALYSIS FUNCTIONS ----
# ============================================================================= #

#' Run gene set enrichment analysis (GSEA) on differential expression results
#'
#' Performs GSEA using fgsea package with ranked gene list based on log2FoldChange.
#' Uses fgsea defaults for pathway size filtering.
#'
#' @param res Data frame with differential expression results
#' @param genesets Named list of character vectors (gene sets)
#' @param geneset.type String, identifier for gene set type (e.g., "GObp", "GOmf")
#' @param gene_col Column name for gene identifiers (default: Label)
#' @param ranking_col Column name for ranking metric (default: log2FoldChange)
#' @param minSize Integer, minimum pathway size (default: 2)
#' @param maxSize Integer, maximum pathway size (default: uses fgsea default)
#'
#' @return Data frame with GSEA results from fgsea
enrich_genesets <- function(res,
                            genesets,
                            geneset.type = "GObp",
                            gene_col = Label,
                            ranking_col = log2FoldChange, 
                            minSize = config$analysis$gsea_min_pathway_size, 
                            maxSize = config$analysis$gsea_max_pathway_size) {
  # Input validation
  if (nrow(res) == 0) {
    cli::cli_abort("No data provided for GSEA")
    return(tibble())
  }
  
  # Create gene rankings using {{ }} in case altered ranking
  rankings_df <- res %>%
    filter(!is.na({{ ranking_col }}), !is.na({{ gene_col }})) %>%
    arrange(desc({{ ranking_col }})) %>%
    distinct({{ gene_col }}, .keep_all = TRUE) %>%
    select(gene = {{ gene_col }}, rank_value = {{ ranking_col }})
  
  ranked_genes <- rankings_df$rank_value %>% set_names(rankings_df$gene)
  
  # Use fgsea default maxSize if not specified
  if (is.null(maxSize)) {
    maxSize <- length(ranked_genes) - 1
  }
  
  gsea_results <- fgsea::fgsea(
    pathways = genesets, 
    stats = ranked_genes, 
    minSize = minSize,
    maxSize = maxSize) %>% 
    arrange(padj)
  
  return(gsea_results)
}

#' Unnest GSEA leading edge genes with annotations
#'
#' Expands GSEA results to show individual leading edge genes with their
#' differential expression statistics and annotations.
#'
#' @param gsea_results Data frame from enrich_genesets()
#' @param res Data frame with differential expression results for annotation
#' @param sigOnly Logical, filter to significant pathways only (default: FALSE)
#'
#' @return Data frame with leading edge genes and their annotations
unnest_genesets <- function(gsea_results, 
                            res, 
                            sigOnly = F) {
  
  if (sigOnly == T) {
    # Future implementation for filtering significant pathways
  }
  
  res_trimmed <- res %>% 
    select(any_of(c(
      "Label",
      "Product",
      "log2FoldChange")), 
      contains(c(
        "log2FC",
        "padj",
        "PsortB",
        "SignalP_5.0"
      )))
  
  gsea_unnested <- gsea_results %>% 
    select(pathway, NES, padj, size, leadingEdge) %>% 
    unnest(leadingEdge) %>% 
    dplyr::rename(Label = leadingEdge, gsea.padj = padj) %>% 
    arrange(gsea.padj, pathway, desc(abs(pick(contains(c("log2FC", "log2FoldChange")))))) %>% 
    left_join(res_trimmed)
  
  return(gsea_unnested)
}

# ============================================================================= #
# GSEA COMPILATION FUNCTIONS ----
# ============================================================================= #

#' Run GSEA compilation for all comparisons
#'
#' Creates Excel file with GSEA results and unnested leading edge genes
#' for all comparisons in the experiment
#'
#' @param res.l.all Named list of comparison results (from main workflow)
#' @param gene_sets Named list of gene sets for GSEA
#' @param gene_set_name String identifier (e.g., "GObp", "GOmf")
#' @param subset_type Which subset to use for ranking ("all", "sig", "DE")
#' @param experiment_name String, experiment identifier
#' @param save_excel Logical, whether to save Excel file (default: TRUE)
#'
#' @return List with gsea_results and unnested_results
run_gsea_compilation <- function(res.l.all,
                                 gene_sets,
                                 gene_set_name,
                                 subset_type = "all",
                                 experiment_name,
                                 save_excel = TRUE) {
  
  cli_inform("Running GSEA compilation: {gene_set_name} ({subset_type} subset)")
  
  # Run GSEA for each comparison
  gsea_results <- res.l.all %>%
    imap(function(comparison_data, comparison_name) {
      cli_inform("  Processing {comparison_name}")
      
      # Get the specified subset for ranking
      res_subset <- comparison_data[[subset_type]]
      
      if (is.null(res_subset) || nrow(res_subset) == 0) {
        cli_warn("No data in {subset_type} subset for {comparison_name}")
        return(NULL)
      }
      
      # Run GSEA
      gsea_result <- res_subset %>%
        enrich_genesets(
          genesets = gene_sets,
          geneset.type = gene_set_name,
          gene_col = Label,
          ranking_col = log2FoldChange
        )
      
      return(gsea_result)
    })
  
  # Create unnested results (leading edge genes with annotations)
  unnested_results <- gsea_results %>%
    imap(function(gsea_result, comparison_name) {
      if (nrow(gsea_result) == 0) return(NULL)
      
      # Get the original data for annotation
      original_data <- res.l.all[[comparison_name]][[subset_type]]
      
      unnested <- gsea_result %>%
        unnest_genesets(res = original_data, sigOnly = FALSE)
      
      return(unnested)
    })
  
  # Save Excel compilation if requested
  if (save_excel && length(gsea_results) > 0) {
    excel_file <- create_gsea_excel_compilation(
      gsea_results = gsea_results,
      unnested_results = unnested_results,
      gene_set_name = gene_set_name,
      subset_type = subset_type,
      experiment_name = experiment_name
    )
    
    cli_inform("Saved GSEA compilation: {basename(excel_file)}")
  }
  
  return(list(
    gsea_results = gsea_results,
    unnested_results = unnested_results
  ))
}

#' Create Excel compilation of GSEA results
#'
#' Uses existing Excel formatting infrastructure for consistent styling.
#' Creates multi-sheet workbook with GSEA results and unnested leading edge genes.
#'
#' @param gsea_results Named list of GSEA results
#' @param unnested_results Named list of unnested results
#' @param gene_set_name String identifier
#' @param subset_type String, subset used for analysis
#' @param experiment_name String, experiment identifier
#'
#' @return String path to created Excel file
create_gsea_excel_compilation <- function(gsea_results,
                                          unnested_results,
                                          gene_set_name,
                                          subset_type,
                                          experiment_name) {
  
  # Add all GSEA results sheets first using imap
  gsea_sheets <- gsea_results %>%
    imap(~ .x) %>%
    set_names(paste0(names(gsea_results), "_gsea"))
  
  # Add all unnested results sheets second using imap  
  unnested_sheets <- unnested_results %>%
    imap(~ .x) %>%
    set_names(paste0(names(unnested_results), "_unnested"))
  
  # Combine with GSEA sheets first, then unnested sheets
  excel_compilation <- c(gsea_sheets, unnested_sheets)
  
  # Create filename
  subset_suffix <- if (subset_type != "all") paste0("_", subset_type) else ""
  filename <- paste0(experiment_name, "_GSEA_", gene_set_name, subset_suffix, ".xlsx")
  
  # Use existing Excel saving infrastructure
  file_path <- save_gsea_excel_file(
    excel_compilation = excel_compilation,
    experiment_name = experiment_name,
    filename = filename
  )
  
  return(file_path)
}

#' Save GSEA Excel file using existing formatting infrastructure
#'
#' Wrapper around save_excel_enhanced for GSEA-specific files.
#' Reuses the Excel formatting system from DE compilation output.
#' Saves files in outputs/results/ directory.
#'
#' @param excel_compilation Named list of data frames
#' @param experiment_name String, experiment identifier
#' @param filename String, output filename
#'
#' @return String path to created file
save_gsea_excel_file <- function(excel_compilation, experiment_name, filename) {
  
  # Create output directory in results subfolder
  output_dir <- ensure_experiment_outputs(experiment_name)
  
  file_path <- file.path(output_dir, "results", filename)
  
  # Use existing Excel formatting function
  save_excel_enhanced(
    data_list = excel_compilation,
    file_path = file_path,
    experiment_name = experiment_name,
    freeze_panes = TRUE,
    add_filters = TRUE
  )
  
  return(file_path)
}