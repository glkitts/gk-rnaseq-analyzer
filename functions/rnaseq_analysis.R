
# RNAseq analysis functions
# Creation date: 2025-09-04

library(tidyverse)
library(DESeq2)
library(tximeta)
library(glue)
library(yaml)

# =============================================================================
# DDS CREATION AND MANAGEMENT ----
# =============================================================================

#' Create DDS from salmon quantification with enhanced error handling
#' @param coldata Data frame with sample metadata and file paths
#' @param design_formula Formula for DESeq2 design (e.g., ~group)
#' @param min_counts Minimum count threshold for filtering
#' @param experiment_name Name for logging and validation
create_dds_from_salmon <- function(coldata, design_formula = ~group, 
                                  min_counts = 2, experiment_name = NULL) {
  
  # Validate inputs
  required_cols <- c("files", "names")
  missing_cols <- required_cols[!required_cols %in% colnames(coldata)]
  if (length(missing_cols) > 0) {
    stop(glue("Missing required columns in coldata: {paste(missing_cols, collapse=', ')}"))
  }
  
  cat(glue("🔬 Creating DDS for {nrow(coldata)} samples...\n"))
  
  # Create tximeta object
  gse <- tximeta(coldata, skipMeta = TRUE)
  
  # Create DESeq2 object
  dds <- DESeqDataSet(gse, design = design_formula)
  
  # Run DESeq2
  cat("🧬 Running DESeq2...\n")
  dds <- DESeq(dds)
  
  # Filter low count transcripts
  dds <- filter_low_counts(dds, min_counts = min_counts, 
                          experiment_name = experiment_name)
  
  return(dds)
}

#' Enhanced low count filtering with reporting
#' @param dds DESeqDataSet object
#' @param min_counts Minimum count threshold
#' @param experiment_name Experiment name for logging
filter_low_counts <- function(dds, min_counts = 2, experiment_name = NULL) {
  nrow_before <- nrow(dds)
  keep <- rowSums(counts(dds)) > min_counts
  dds_filtered <- dds[keep, ]
  nrow_after <- nrow(dds_filtered)
  nrow_removed <- nrow_before - nrow_after
  
  cat(glue("🔍 Filtering: {nrow_removed} of {nrow_before} transcripts removed (≤{min_counts} counts)\n"))
  cat(glue("✅ {nrow_after} transcripts retained for analysis\n"))
  
  # Log filtering stats if experiment name provided
  if (!is.null(experiment_name)) {
    log_filtering_stats(nrow_before, nrow_after, nrow_removed, 
                       min_counts, experiment_name)
  }
  
  return(dds_filtered)
}


# =============================================================================
# SAVE AND LOAD DDS ----
# =============================================================================

save_dds_with_metadata <- function(dds, experiment_name, coldata = NULL, config = NULL) {
  # Save DESeq2 object with comprehensive metadata for reproducibility
  #
  # This function saves the processed DESeq2 object and extracts key metadata
  # including analysis parameters, quality metrics, and dataset characteristics.
  # All files are saved in the experiment's outputs directory.
  #
  # Parameters:
  #   dds: DESeqDataSet object, must be processed with DESeq() function
  #   experiment_name: String, experiment identifier used for file naming
  #   coldata: Data frame, original sample metadata for additional context (optional)
  #   config: List, analysis configuration parameters from global_settings.yaml (optional)
  #
  # Returns:
  #   List containing:
  #     - dds_path: Full path to saved DDS RDS file
  #     - metadata_path: Full path to saved metadata YAML file
  #
  # Files Created:
  #   - {experiment_name}_dds.RDS: The DESeq2 object
  #   - {experiment_name}_metadata.yaml: Analysis metadata and quality metrics
  #
  # Example:
  #   result <- save_dds_with_metadata(dds, "my_experiment", coldata, config)
  #   # Creates: outputs/my_experiment_dds.RDS and outputs/my_experiment_metadata.yaml
  
  # Validate inputs
  if (!is(dds, "DESeqDataSet")) {
    stop("dds must be a DESeqDataSet object")
  }
  if (is.null(experiment_name) || experiment_name == "") {
    stop("experiment_name is required")
  }
  
  # Ensure all output directories exist
  outputs_dir <- ensure_experiment_outputs(experiment_name)
  
  # Construct file paths
  dds_path <- file.path(outputs_dir, paste0(experiment_name, "_dds.RDS"))
  metadata_path <- file.path(outputs_dir, paste0(experiment_name, "_metadata.yaml"))
  
  # Extract metadata from DDS object
  metadata <- list(
    experiment = list(
      name = experiment_name,
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      deseq2_version = as.character(packageVersion("DESeq2"))
    ),
    
    dataset = list(
      n_samples = ncol(dds),
      n_transcripts = nrow(dds),
      design_formula = as.character(design(dds))[2],
      sample_names = colnames(dds)
    ),
    
    analysis_params = list(
      min_counts = if (!is.null(config)) config$analysis$min_counts else 2,
      significance_threshold = if (!is.null(config)) config$analysis$significance_threshold else 0.05,
      fold_change_threshold = if (!is.null(config)) config$analysis$fold_change_threshold else 1.0
    )
  )
  
  # Add quality metrics if available
  if (!is.null(sizeFactors(dds))) {
    metadata$quality_metrics <- list(
      size_factors_range = round(range(sizeFactors(dds)), 3),
      library_sizes_range = range(colSums(counts(dds)))
    )
  }
  
  if (!is.null(dispersions(dds))) {
    if (is.null(metadata$quality_metrics)) {
      metadata$quality_metrics <- list()
    }
    metadata$quality_metrics$dispersion_median <- round(median(dispersions(dds), na.rm = TRUE), 4)
  }
  
  # Save DDS object
  cat(sprintf("Saving DDS: %s\n", basename(dds_path)))
  saveRDS(dds, dds_path)
  
  # Save metadata
  cat(sprintf("Saving metadata: %s\n", basename(metadata_path)))
  yaml::write_yaml(metadata, metadata_path)
  
  # Report success with key information
  file_size_mb <- round(file.size(dds_path) / 1024^2, 1)
  cat(sprintf("DDS saved successfully (%s MB) | Samples: %d | Transcripts: %d\n", 
              file_size_mb, ncol(dds), nrow(dds)))
  
  return(list(
    dds_path = dds_path, 
    metadata_path = metadata_path
  ))
}


#' 
#' 
#' #' Load DDS with validation
#' #' @param experiment_name Experiment name
#' load_dds_with_validation <- function(experiment_name) {
#'   dds_path <- file.path("outputs", "dds", paste0(experiment_name, "_dds.RDS"))
#'   metadata_path <- file.path("outputs", "dds", paste0(experiment_name, "_metadata.yaml"))
#'   
#'   if (!file.exists(dds_path)) {
#'     stop(glue("❌ DDS file not found: {dds_path}\nRun DDS generation first with --dds-only"))
#'   }
#'   
#'   dds <- readRDS(dds_path)
#'   
#'   # Load and display metadata if available
#'   if (file.exists(metadata_path)) {
#'     metadata <- read_yaml(metadata_path)
#'     cat(glue("📊 Loaded DDS: {metadata$n_samples} samples, {metadata$n_transcripts} transcripts\n"))
#'     cat(glue("🕐 Generated: {metadata$timestamp}\n"))
#'   } else {
#'     cat(glue("📊 Loaded DDS: {ncol(dds)} samples, {nrow(dds)} transcripts\n"))
#'   }
#'   
#'   return(dds)
#' }
#' 
#' # =============================================================================
#' # ANNOTATION PROCESSING
#' # =============================================================================
#' 
#' #' Prepare annotations in standardized format
#' #' @param annotations Raw annotation object (from RDS file)
#' prepare_annotations <- function(annotations) {
#'   # Your current annotation processing logic, standardized
#'   if (is.list(annotations) && "pathways_functions" %in% names(annotations)) {
#'     sub_anots <- annotations$pathways_functions
#'   } else {
#'     sub_anots <- annotations
#'   }
#'   
#'   # Standardize columns and add helper columns
#'   sub_anots <- sub_anots %>%
#'     mutate(
#'       label_rows = if_else(
#'         is.na(yildiz_name),
#'         Label,
#'         paste0(Label, "  *(", yildiz_name, ")*")
#'       )
#'     ) %>%
#'     distinct(Label, .keep_all = TRUE)
#'   
#'   return(sub_anots)
#' }
#' 
#' # =============================================================================
#' # ENHANCED DIFFERENTIAL EXPRESSION ANALYSIS
#' # =============================================================================
#' 
#' #' Enhanced version of your deseq_makeOutput function
#' #' Drop-in replacement with better organization and error handling
#' enhanced_deseq_makeOutput <- function(contrast_list, dds, sub_anots, 
#'                                      merge_anno = TRUE, experiment_name = NULL, 
#'                                      config = NULL) {
#'   
#'   cat(glue("🧬 Running differential expression for {length(contrast_list)} comparisons...\n"))
#'   
#'   # Generate output names from contrast list
#'   contrast_names <- contrast_list %>%
#'     map_chr(function(x) paste0(x[2], "_vs_", x[3]))
#'   
#'   contrast_list <- set_names(contrast_list, contrast_names)
#'   
#'   # Generate results for each contrast
#'   cat("📊 Generating results...\n")
#'   res_list <- contrast_list %>%
#'     imap(function(contrast, name) {
#'       cat(glue("   • {name}\n"))
#'       deseq_result_generator(dds, contrast)
#'     })
#'   
#'   # Process results with annotations
#'   cat("🏷️  Adding annotations...\n")
#'   res_out <- res_list %>%
#'     map(function(x) {
#'       deseq_sigFC_subsetter(x, anno = sub_anots, merge_anno = merge_anno)
#'     }) %>%
#'     unlist(recursive = FALSE)
#'   
#'   # Create Excel-formatted results
#'   cat("📋 Compiling Excel results...\n")
#'   res_excel <- res_list %>%
#'     imap(function(x, name) {
#'       x %>%
#'         deseq_sigFC_subsetter(anno = sub_anots, merge_anno = merge_anno) %>%
#'         .$sig %>%
#'         deseq_sigCompiler(name)
#'     }) %>%
#'     reduce(full_join)
#'   
#'   # Add annotations to Excel results
#'   if (merge_anno && !is.null(sub_anots)) {
#'     res_excel <- res_excel %>%
#'       left_join(sub_anots, by = "Label")
#'   }
#'   
#'   # Organize Excel results
#'   res_excel <- res_excel %>%
#'     relocate(
#'       Label,
#'       any_of(c("Product")),
#'       contains("log2FC"),
#'       contains("log2FoldChange"),
#'       any_of(c("Gene", "GeneNames", "PsortB", "SignalP_5.0")),
#'       starts_with("GO terms"),
#'       starts_with("yildiz_"), 
#'       any_of(c("Transcription Units")),
#'       starts_with("TIGR"),
#'       contains("padj")
#'     ) %>%
#'     list() %>%
#'     set_names("sigHits_compiled") %>%
#'     append(res_out)
#'   
#'   # Filter to relevant results
#'   res_excel <- res_excel %>%
#'     keep(str_detect(names(.), "sigHits_compiled|all|sig_fc"))
#'   
#'   # Add gene subset analysis if genesets available
#'   if (exists("genesets", envir = .GlobalEnv)) {
#'     gene_subsets <- res_excel$sigHits_compiled %>%
#'       rseq_geneSubsettor(genesets)
#'     res_excel <- append(res_excel, gene_subsets, after = 1)
#'   }
#'   
#'   # Final output structure
#'   output <- list(
#'     out = res_out,
#'     xlsx = res_excel,
#'     dds = dds,
#'     list = res_list,
#'     metadata = list(
#'       experiment_name = experiment_name,
#'       timestamp = Sys.time(),
#'       n_comparisons = length(contrast_list),
#'       contrast_names = contrast_names
#'     )
#'   )
#'   
#'   cat(glue("✅ Differential expression analysis complete!\n"))
#'   
#'   return(output)
#' }
#' 
#' # =============================================================================
#' # RESULT MANAGEMENT  
#' # =============================================================================
#' 
#' #' Save experiment results in organized structure
#' #' @param results Results list from enhanced_deseq_makeOutput
#' #' @param experiment_name Experiment name
#' #' @param config Configuration list
#' save_experiment_results <- function(results, experiment_name, config = NULL) {
#'   
#'   # Create results directories
#'   results_dir <- file.path("outputs", "results")
#'   individual_dir <- file.path(results_dir, "individual_comparisons")
#'   dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)
#'   
#'   # Save combined Excel file
#'   excel_path <- file.path(results_dir, paste0(experiment_name, "_results.xlsx"))
#'   writexl::write_xlsx(results$xlsx, excel_path)
#'   
#'   # Save individual CSV files
#'   iwalk(results$xlsx, function(data, name) {
#'     if (str_detect(name, "all$|sig$|sig_fc$")) {
#'       csv_path <- file.path(individual_dir, paste0(name, ".csv"))
#'       write_csv(data, csv_path)
#'     }
#'   })
#'   
#'   # Save R object for programmatic access
#'   rds_path <- file.path(results_dir, paste0(experiment_name, "_results.RDS"))
#'   saveRDS(results, rds_path)
#'   
#'   cat(glue("💾 Results saved:\n"))
#'   cat(glue("   📊 Excel: {excel_path}\n"))
#'   cat(glue("   📁 CSVs: {individual_dir}/\n"))
#'   cat(glue("   🔧 RDS: {rds_path}\n"))
#'   
#'   return(list(excel_path = excel_path, rds_path = rds_path))
#' }
#' 
#' #' Load saved experiment results
#' #' @param experiment_name Experiment name
#' rnaseq_load_results <- function(experiment_name) {
#'   rds_path <- file.path("outputs", "results", paste0(experiment_name, "_results.RDS"))
#'   
#'   if (!file.exists(rds_path)) {
#'     stop(glue("❌ Results file not found: {rds_path}\nRun analysis first with --analysis-only"))
#'   }
#'   
#'   results <- readRDS(rds_path)
#'   cat(glue("📊 Loaded results: {results$metadata$n_comparisons} comparisons\n"))
#'   
#'   return(results)
#' }
#' 
# =============================================================================
# LOGGING AND UTILITIES
# =============================================================================

#' Log filtering statistics
log_filtering_stats <- function(nrow_before, nrow_after, nrow_removed,
                               min_counts, experiment_name) {

  log_dir <- file.path("outputs", "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  log_entry <- list(
    timestamp = Sys.time(),
    experiment = experiment_name,
    filtering = list(
      transcripts_before = nrow_before,
      transcripts_after = nrow_after,
      transcripts_removed = nrow_removed,
      min_counts_threshold = min_counts,
      percent_retained = round(nrow_after / nrow_before * 100, 1)
    )
  )

  log_path <- file.path(log_dir, paste0(experiment_name, "_filtering.yaml"))
  write_yaml(log_entry, log_path)
}

#' # =============================================================================
#' # WRAPPER FUNCTIONS FOR YOUR EXISTING CODE
#' # =============================================================================
#' 
#' # Keep your existing function names for compatibility
#' deseq_result_generator <- function(dds, contrast) {
#'   comparison_name <- paste0(contrast[2], "_vs_", contrast[3])
#'   res <- results(dds, contrast = contrast)
#'   res.LFC <- as.data.frame(lfcShrink(dds, contrast = contrast, type = "ashr"))
#'   return(res.LFC)
#' }
#' 
#' deseq_sigFC_subsetter <- function(res.LFC, anno = NULL, merge_anno = FALSE) {
#'   res.l <- list()
#'   
#'   res.l$all <- res.LFC %>%
#'     select(log2FoldChange, padj, baseMean) %>%
#'     rownames_to_column() %>%
#'     separate(rowname, sep = '\\|', into = c('Label', NA, NA)) %>%
#'     arrange(log2FoldChange) %>% 
#'     distinct(Label, .keep_all = TRUE)
#'   
#'   res.l$sig <- res.l$all %>% filter(padj <= 0.05)
#'   res.l$sig_fc <- res.l$sig %>% filter(abs(log2FoldChange) >= 1)
#'   
#'   if (merge_anno && !is.null(anno)) {
#'     res.l <- res.l %>%
#'       map(function(x) left_join(x, anno, by = "Label"))
#'   }
#'   
#'   return(res.l)
#' }
#' 
#' deseq_sigCompiler <- function(d, name) { 
#'   new_names <- c("log2FC", "padj", "baseMean")
#'   new_names <- c("Label", paste0(name, "_", new_names))
#'   
#'   d <- d %>% 
#'     filter(padj <= 0.05) %>%
#'     arrange(padj) %>% 
#'     distinct(Label, .keep_all = TRUE) %>% 
#'     select(Label, log2FoldChange, padj, baseMean) %>% 
#'     mutate(log2FoldChange = signif(log2FoldChange, digits = 3))
#'   
#'   colnames(d) <- new_names
#'   return(d)
#' }
#' 
#' rseq_geneSubsettor <- function(res, genesets) {
#'   a <- res %>% select(Label, contains("log2FC"))
#' 
#'   gset_res <- genesets %>%
#'     map(function(x) left_join(x, a, by = "Label"))
#' 
#'   return(gset_res)
#' }