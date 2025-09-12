# Enhanced RNA-seq Pipeline - Detailed Implementation Guide

## Implementation Order & Development Strategy

This guide provides detailed specifications for building enhanced functions and automation on top of your existing workflow, ordered by development priority and testing approach.

---

## Phase 1: Enhanced Functions & Basic Templates (Week 1)

### 1. `setup.R` - Project Setup & Structure Creation
**Purpose:** One-time project initialization and directory creation
**When to run:** First time only

#### **Functionality:**
- Check R version and required packages
- Create directory structure for hybrid system  
- Generate global settings template
- Create experiment template files
- Test basic functionality

#### **Key Functions to Implement:**
```r
check_r_version(required = "4.1.0")
install_required_packages(cran_packages, bioc_packages)
create_hybrid_directory_structure()
create_global_settings_template()
create_experiment_template()
test_enhanced_functions()
```

#### **Outputs:**
- Complete directory structure
- `config/global_settings.yaml` template
- `experiments/TEMPLATE/analysis_script.R` template
- Setup validation report

---

### 2. `config/global_settings.yaml` - Global Resources Only
**Purpose:** Paths, styling, and analysis defaults (no experiment-specific configs)
**When to create:** During setup, modify for your environment

#### **Structure:**
```yaml
paths:
  metadata_file: "shared_resources/metadata/rnaseq_yildiz_metadata.csv"
  salmon_base: "../../../vc-refs/rnaseq_counts" 
  annotation_db: "shared_resources/annotations/vc_genome_anno.RDS"

analysis:
  significance_threshold: 0.05
  fold_change_threshold: 1.0
  min_counts: 2

styling:
  base_size: 8
  figure_width: 7
  figure_height: 5
  pathway_colors: "tol"

output:
  generate_pdf: true
  generate_html: false  # Optional interactive reports
  save_plots_separately: true
```

---

### 3. `functions/enhanced_analysis.R` - Drop-in Function Replacements
**Purpose:** Enhanced versions of your current functions with same interfaces
**Dependencies:** DESeq2, tidyverse

#### **Core Functions to Implement:**

```r
enhanced_deseq_makeOutput <- function(contrast_list, dds, sub_anots, 
                                     merge_anno = TRUE, experiment_name = NULL, 
                                     config = NULL) {
  # Drop-in replacement for your current deseq_makeOutput function
  # Same interface, enhanced functionality:
  # - Better error handling and progress reporting
  # - Professional plot generation and auto-saving
  # - Organized output file management  
  # - Enhanced Excel export with better formatting
  
  config <- config %||% load_global_config()
  
  # Your existing logic enhanced:
  contrast_names <- generate_contrast_names(contrast_list)
  res_list <- generate_results_list(contrast_list, dds)
  res_out <- process_results_with_annotations(res_list, sub_anots, merge_anno)
  res_excel <- compile_excel_results(res_list, sub_anots, merge_anno)
  
  # Enhanced outputs:
  if (!is.null(experiment_name)) {
    save_organized_results(res_excel, experiment_name)
    generate_summary_plots(res_list, experiment_name, config)
  }
  
  return(res_excel)
}

filter_low_counts <- function(dds, min_counts = 2, report = TRUE, experiment_name = NULL) {
  # Enhanced version of your filtering with automatic reporting
  nrow_before <- nrow(dds)
  keep <- rowSums(counts(dds)) > min_counts
  dds_filtered <- dds[keep, ]
  nrow_after <- nrow(dds_filtered)
  nrow_removed <- nrow_before - nrow_after
  
  if (report) {
    cat(glue("Filtering: {nrow_removed} of {nrow_before} transcripts removed (< {min_counts} counts)\n"))
    if (!is.null(experiment_name)) {
      write_filtering_report(nrow_before, nrow_after, nrow_removed, experiment_name)
    }
  }
  
  return(dds_filtered)
}

save_dds <- function(dds, experiment_name) {
  # Save DDS with metadata and validation
  dds_path <- file.path("outputs", "dds", paste0(experiment_name, "_dds.RDS"))
  dir.create(dirname(dds_path), recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(dds, dds_path)
  
  # Save metadata
  metadata <- list(
    experiment = experiment_name,
    timestamp = Sys.time(),
    n_samples = ncol(dds),
    n_transcripts = nrow(dds),
    design_formula = design(dds)
  )
  
  metadata_path <- file.path("outputs", "dds", paste0(experiment_name, "_metadata.yaml"))
  yaml::write_yaml(metadata, metadata_path)
  
  cat(glue("DDS saved: {dds_path}\n"))
  return(dds_path)
}

load_dds <- function(experiment_name) {
  # Load DDS with validation
  dds_path <- file.path("outputs", "dds", paste0(experiment_name, "_dds.RDS"))
  
  if (!file.exists(dds_path)) {
    stop(glue("DDS file not found: {dds_path}. Run DDS generation first."))
  }
  
  dds <- readRDS(dds_path)
  cat(glue("DDS loaded: {dds_path}\n"))
  return(dds)
}
```

---

### 4. `functions/enhanced_plotting.R` - Professional Visualization Functions
**Purpose:** Enhanced plotting with publication-ready styling
**Dependencies:** ggplot2, ggrepel, ggtext, viridis, patchwork

#### **Core Functions to Implement:**

```r
setup_plot_theme <- function(config = NULL) {
  # Set up professional plot theme based on config
  config <- config %||% load_global_config()
  
  base_size <- config$styling$base_size
  
  theme_set(
    ggthemes::theme_few(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold"),
      axis.title = element_markdown(),
      legend.title = element_markdown(),
      strip.text = element_text(face = "bold")
    )
  )
}

enhanced_volcano_plot <- function(res_df, config = NULL, pathway_col = "yildiz_pathway",
                                 label_top_genes = TRUE, n_labels = 25) {
  # Enhanced version of your volcano plot
  config <- config %||% load_global_config()
  
  # Your current volcano plot logic enhanced with:
  # - Better pathway coloring
  # - Professional styling  
  # - Automatic labeling optimization
  # - Consistent sizing and formatting
  
  setup_plot_theme(config)
  
  # Process data for plotting (your current logic)
  plot_data <- prepare_volcano_data(res_df, pathway_col, n_labels)
  
  # Create plot with enhanced styling
  p <- plot_data %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj), 
               color = path_cols, fill = path_cols, label = Label)) +
    create_volcano_layers() +
    apply_volcano_styling(config) +
    add_volcano_labels(label_top_genes, n_labels)
  
  return(p)
}

enhanced_pathway_plot <- function(pathway_data, config = NULL, max_pathways = 20) {
  # Enhanced GO biological process plot
  config <- config %||% load_global_config()
  
  # Your current pathway plotting logic enhanced with:
  # - Dynamic height calculation
  # - Better color schemes
  # - Professional formatting
  
  plot_height <- calculate_pathway_plot_height(pathway_data, max_pathways)
  
  p <- pathway_data %>%
    slice_head(n = max_pathways) %>%
    ggplot(aes(x = log2FC, y = pathway, color = -log10(padj), size = n)) +
    geom_point(aes(fill = after_scale(alpha(color, 0.8))), shape = 21) +
    apply_pathway_styling(config)
  
  return(list(plot = p, height = plot_height))
}

save_plot_professional <- function(plot, filename, width = NULL, height = NULL, 
                                  config = NULL, experiment_name = NULL) {
  # Professional plot saving with consistent formatting
  config <- config %||% load_global_config()
  width <- width %||% config$styling$figure_width
  height <- height %||% config$styling$figure_height
  
  if (!is.null(experiment_name)) {
    plot_dir <- file.path("outputs", "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(plot_dir, filename)
  }
  
  # Save both PDF and PNG
  ggsave(paste0(filename, ".pdf"), plot, width = width, height = height, device = "pdf")
  ggsave(paste0(filename, ".png"), plot, width = width, height = height, device = "png", dpi = 300)
  
  cat(glue("Plots saved: {filename}.pdf, {filename}.png\n"))
}
```

---

### 5. `functions/output_management.R` - File Organization System
**Purpose:** Automated, consistent output organization
**Dependencies:** fs, glue

#### **Core Functions to Implement:**

```r
setup_experiment_output <- function(experiment_name) {
  # Create organized directory structure for experiment
  base_dir <- file.path("outputs")
  
  directories <- c(
    "dds",
    "results",
    "results/individual_comparisons", 
    "plots",
    "plots/volcano_plots",
    "plots/pathway_plots",
    "plots/qc_plots",
    "reports",
    "logs"
  )
  
  walk(directories, ~{
    dir.create(file.path(base_dir, .x), recursive = TRUE, showWarnings = FALSE)
  })
  
  cat(glue("Output directories created for {experiment_name}\n"))
  return(base_dir)
}

save_organized_results <- function(results_list, experiment_name) {
  # Save results in organized structure
  
  # Combined Excel file
  excel_path <- file.path("outputs", "results", paste0(experiment_name, "_results.xlsx"))
  writexl::write_xlsx(results_list, excel_path)
  
  # Individual CSV files for each comparison
  comparison_dir <- file.path("outputs", "results", "individual_comparisons")
  
  iwalk(results_list, ~{
    if (str_detect(.y, "all$|sig$|sig_fc$")) {
      csv_path <- file.path(comparison_dir, paste0(.y, ".csv"))
      write_csv(.x, csv_path)
    }
  })
  
  cat(glue("Results saved: {excel_path} and individual CSV files\n"))
}

generate_experiment_readme <- function(experiment_name, results_summary = NULL) {
  # Auto-generate experiment summary README
  readme_content <- glue("
# {experiment_name} Analysis Summary

**Generated:** {Sys.time()}

## Experiment Overview
- **DDS Object:** outputs/dds/{experiment_name}_dds.RDS
- **Results:** outputs/results/{experiment_name}_results.xlsx
- **Reports:** outputs/reports/

## Comparisons Analyzed
{generate_comparison_list(results_summary)}

## Output Files
- `dds/`: DESeq2 dataset objects
- `results/`: Excel and CSV result files
- `plots/`: All plots organized by type
- `reports/`: PDF and HTML reports
- `logs/`: Processing logs and metadata

## Regeneration Commands
```bash
# Regenerate reports only
Rscript ../../run_experiment.R {experiment_name} --reports-only

# Regenerate analysis + reports  
Rscript ../../run_experiment.R {experiment_name} --analysis-only

# Full pipeline
Rscript ../../run_experiment.R {experiment_name}
```
")
  
  readme_path <- file.path("outputs", "README.md")
  writeLines(readme_content, readme_path)
  
  return(readme_path)
}
```

---

### 6. `templates/enhanced_comparison.qmd` - Professional Report Template
**Purpose:** Enhanced version of your current template with professional styling
**Dependencies:** Quarto, knitr, all analysis functions

#### **Key Enhancements:**
- Professional PDF styling with better typography
- Enhanced volcano plots with pathway coloring
- Improved GO pathway analysis visualization
- Better table formatting and layout
- Optional interactive HTML output
- Consistent figure sizing and quality

#### **Template Structure:**
```yaml
---
title: "RNA-seq Analysis: `r params$contrast[2]` vs `r params$contrast[3]`"
subtitle: "`r params$experiment_name`"
author: "Enhanced RNA-seq Pipeline"
date: "`r Sys.Date()`"
format:
  pdf:
    toc: true
    number-sections: true
    colorlinks: true
    geometry:
      - top=0.7in
      - right=0.7in  
      - bottom=0.7in
      - left=0.7in
    fig-cap-location: bottom
  html:
    toc: true
    code-fold: true
    interactive: true
params:
  contrast: !r c("group", "treatment", "control")
  dds_path: "path/to/dds.RDS"
  results_data: null
  experiment_name: "experiment"
  config: null
---
```

---

### 7. `experiments/TEMPLATE/analysis_script.R` - Experiment Template
**Purpose:** Template for creating new experiment analysis scripts
**When to use:** Copy and modify for each new experiment

#### **Template Structure:**
```r
#!/usr/bin/env Rscript
library(tidyverse)
library(DESeq2)
library(tximeta)

# Load enhanced functions
source("../../functions/enhanced_analysis.R")
source("../../functions/enhanced_plotting.R") 
source("../../functions/output_management.R")

# Experiment configuration
experiment_name <- "EXPERIMENT_NAME"  # Change this
config <- yaml::read_yaml("../../config/global_settings.yaml")

# Set up output directory
setup_experiment_output(experiment_name)

# Determine processing step
pipeline_step <- Sys.getenv("PIPELINE_STEP", "full")
cat(glue("Running {experiment_name} - Step: {pipeline_step}\n"))

# STEP 1: DDS GENERATION (your current logic)
if (pipeline_step %in% c("full", "dds-only")) {
  cat("🔬 Generating DDS...\n")
  
  # YOUR CURRENT FILTERING LOGIC HERE (customize per experiment)
  coldata.in <- read.csv(config$paths$metadata_file)
  
  coldata <- coldata.in %>%
    filter(experiment_group %in% c("iron")) %>%  # Customize filters
    filter(!strain %in% c("dryhB")) %>%          # Customize filters
    # Add your custom filtering logic here
    mutate(
      strain = fct_relevel(strain, "WT"),
      condition = fct_relevel(condition, "iron_replete"),
      group = fct_inorder(sampleLabel),
      names = sample,
      files = file.path(
        config$paths$salmon_base,
        experiment_folder,
        sample_FolderName,
        "salmon_quants",
        "quant.sf"
      )
    )
  
  # DDS creation (your current approach)
  gse <- tximeta(coldata, skipMeta = TRUE)
  dds <- DESeqDataSet(gse, design = ~group)
  dds <- DESeq(dds)
  dds <- filter_low_counts(dds, min_counts = 2, experiment_name = experiment_name)
  save_dds(dds, experiment_name)
}

# STEP 2: ANALYSIS (enhanced functions)
if (pipeline_step %in% c("full", "analysis-only")) {
  if (!exists("dds")) dds <- load_dds(experiment_name)
  
  # YOUR CONTRAST LIST (customize per experiment)
  contrast_list <- list(
    c("group", "IR_drvvA", "IR_WT"),
    c("group", "ID_drvvA", "ID_WT")
    # Add your contrasts here
  )
  
  # Enhanced analysis (replaces your deseq_makeOutput)
  results <- enhanced_deseq_makeOutput(
    contrast_list = contrast_list,
    dds = dds,
    sub_anots = readRDS(config$paths$annotation_db),
    experiment_name = experiment_name,
    config = config
  )
}

# STEP 3: REPORTS (automated generation)
if (pipeline_step %in% c("full", "analysis-only", "reports-only")) {
  if (!exists("results")) {
    results <- load_experiment_results(experiment_name)
    dds <- load_dds(experiment_name)
  }
  
  generate_all_reports(results, dds, experiment_name, config)
  generate_experiment_readme(experiment_name, results)
}

cat(glue("✅ {experiment_name} analysis complete!\n"))
```

---

## Phase 2: Batch Processing & Automation (Week 2)

### 8. `run_experiment.R` - Single Experiment Runner
**Purpose:** Run individual experiments with step control
**Dependencies:** Command line argument parsing

#### **Key Features:**
- Step-selective processing (DDS, analysis, reports)
- Error handling and progress reporting  
- Working directory management
- Environment variable communication

### 9. `run_all.R` - Batch Processing Controller  
**Purpose:** Process multiple experiments with selective regeneration
**Dependencies:** All experiment functions

#### **Command Options:**
```bash
Rscript run_all.R --full                    # Complete rebuild all experiments
Rscript run_all.R --dds-only               # Regenerate DDS for all
Rscript run_all.R --analysis-only          # Analysis + reports for all  
Rscript run_all.R --reports-only           # Reports only for all
Rscript run_all.R --experiment NAME        # Single experiment processing
Rscript run_all.R --list                   # List available experiments
```

### 10. `functions/report_automation.R` - Automated Report Generation
**Purpose:** Automated rendering of reports across experiments
**Dependencies:** Quarto, all analysis functions

---

## Phase 3: Professional Styling & Advanced Features (Week 3)

### 11. Enhanced Templates with Professional Styling
- Publication-ready PDF formatting
- Interactive HTML versions
- Enhanced visualizations
- Better table formatting

### 12. Advanced Visualization Features  
- Pathway enrichment enhancements
- Multi-experiment comparison plots
- Quality control dashboards  
- Custom plot themes

### 13. Experiment Summary Generation
- Cross-experiment analysis
- Shared gene identification
- Multi-experiment reports

---

## Implementation Workflow

### **Week 1 Testing:**
```bash
# 1. Set up system
Rscript setup.R

# 2. Test enhanced functions with existing experiment
cp -r experiments/TEMPLATE experiments/test_rsRDM_all
# Edit analysis_script.R with your rsRDM_all filtering logic

# 3. Run single experiment
Rscript run_experiment.R test_rsRDM_all

# 4. Compare outputs with your current workflow
```

### **Week 2 Scaling:**
```bash
# 1. Convert all existing experiments to new system
# 2. Test batch processing
Rscript run_all.R --reports-only

# 3. Validate consistency across experiments
```

### **Week 3 Enhancement:**
```bash  
# 1. Polish templates and styling
# 2. Add advanced features
# 3. Performance optimization
```

This implementation guide provides a systematic approach to building the enhanced system while preserving your current workflow flexibility.