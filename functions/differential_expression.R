#!/usr/bin/env Rscript

# RNAseq differential analysis functions
# Creation date: 2025-09-04

library(tidyverse)
library(DESeq2)
library(tximeta)
library(glue)
library(yaml)
library(here)

here::i_am("functions/differential_expression.R")


format_comparison_name <- function(contrast, separator = "_v_") {
  contrast_name <- paste0(contrast[2], "_v_", contrast[3])
  
  return(contrast_name)
}


generate_resLFC <- function(dds,
                            contrast,
                            shrink_method = "ashr",
                            parse_rownames = T,
                            ...) {
  #' Extracts results from DDS and performs lfcShrink
  #' Parses VC#### from rownames and uses as Label column
  #'
  #' @return tibble with gene labels as rownames
  
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
            "PsortB", 
            "SignalP_5.0",
            "Gene", 
            "GeneNames", 
            "Type"
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



