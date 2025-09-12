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

generate_resLFC <- function(dds, contrast, shrink_method = "ashr", parse_rownames = T, ...) {
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
  res.l$DE <- resLFC %>%
    filter(padj <= significance_threshold &
             abs(log2FoldChange) >= log2FC_threshold)
  
  return(res.l)
}

merge_annotations <- function(res_df,
                              annotations,
                              join_col = "Label",
                              relocate_cols = T,
                              cols_to_relocate = NULL) {
  merged_res <- res_df %>%
    left_join(annotations, by = "Label")
  
  if (relocate_cols == T) {
    if (!is.null(cols_to_relocate)) {
      merged_res <- merged_res %>%
        relocate(cols_to_relocate)
      
      else {
        merged_res <- merged_res %>%
          relocate(
            Label,
            any_of(c("Product", "Type")),
            contains(c("log2FoldChange", "log2FC")),
            contains(c("padj")),
            any_of(c(
              "Gene", "GeneNames", "PsortB", "SignalP_5.0"
            )),
            starts_with("GO terms"),
            starts_with("yildiz_"),
            any_of(c("Begin", "End", "Length", ))
            any_of(c("Transcription Units")),
            starts_with("TIGR"),
            Locus_Tag_Old,
            UniProtID,
          )
      }
    }
  }
  return(merged_res)
}

save_individual_csvs <- function(res.l,
                                 experiment_name,
                                 comparison_name,
                                 output_dir = NULL,
                                 types = c("all", "DE")) {
  ensure_experiment_outputs(experiment_name)
  
  if (is.null(output_dir)) {
    output_dir <- here("experiments",
                       experiment_name,
                       "outputs/results/individual_comparisons/")
  }
  
  filtered_res.l <- res.l[names(res.l) %in% types]
  
  # Return file paths instead of data
  file_paths <- filtered_res.l %>%
    imap_chr(function(res, name) {
      output_path <- paste0(output_dir, comparison_name, "_", name, ".csv")
      
      res %>%
        readr::write_excel_csv(file = output_path)
      
      return(output_path)
    })
  
  return(file_paths)
}
