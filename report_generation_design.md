# RNA-seq Report Generation Design Specification

*Design Date: December 12, 2025*

## Overview

This document outlines the design for automated report generation in the RNA-seq Enhanced Analysis Pipeline. The system will generate both experimental summary reports and individual differential expression comparison reports using Quarto, with outputs in HTML and PDF formats.

## Core Requirements

### Report Types
1. **Experimental Summary Report** - Overview of entire experiment (samples, QC, PCA, all comparisons)
2. **Individual DE Comparison Reports** - Detailed analysis per comparison (volcano plots, GSEA, top genes)
3. **Overview Landing Page** - Main entry point with links to all reports

### Output Formats
- **HTML** - Primary format with interactive elements and linking
- **PDF** - Publication-ready static reports
- **Dynamic sizing** - Plot dimensions adjust based on data (e.g., number of significant pathways)

### Integration
- Called from `analysis_script.R` as pipeline step
- Linked navigation between reports
- Portable experiment directories for collaboration

## Architecture Decision: File-Based Data Flow

**Selected Approach:** Use saved files as primary data source rather than passing complex objects as Quarto parameters.

### Benefits
- Avoids Quarto parameter serialization issues with S4 objects (DDS)
- Enables selective data loading
- Maintains consistency with existing file-based pipeline
- Makes reports debuggable and reproducible independently
- Allows efficient batch processing

### Data Flow
```
analysis_script.R → Save report-ready data → Quarto reports load files → Generate reports
```

## Data Storage Decision: Single RDS File

**Selected Approach:** Single RDS file containing all report data rather than multiple separate files.

### File Structure
```r
# Single file: {experiment_name}_report_data.RDS
report_data <- list(
  dds = dds,                          # QC data, counts, sample info
  res.l.all = res.l.all,              # All DE results with subsets & annotations
  contrast_list = contrast_list,       # Comparison definitions
  coldata = coldata,                  # Sample metadata
  config_snapshot = config,           # Analysis parameters used
  gsea_results = list(                # GSEA results if generated
    GObp = gsea_gobp_results,
    GOmf = gsea_gomf_results
  ),
  metadata = list(                    # Additional metadata
    experiment_name = experiment_name,
    analysis_date = Sys.Date(),
    n_comparisons = length(contrast_list)
  )
)
```

### Advantages
- **Smaller size** - RDS compression more efficient across related objects
- **Atomic loading** - One command loads everything needed
- **Simpler management** - Only one file to track/transfer
- **Consistency guarantee** - All objects from same analysis run
- **Portable** - Single file contains complete analysis state

## Directory Organization: Option 1 (Collaborator-Friendly)

**Selected Structure:**
```
experiments/rsRDM_all/
├── analysis_script.R                  # Keep at root - easy to find/edit
├── README.md                         # Brief experiment description
└── outputs/
    ├── overview.html                  # 👈 NEW: Main landing page
    ├── rsRDM_all_DEcompilation.xlsx   # Main results - front and center
    ├── reports/                       # Visual summaries & detailed reports
    │   ├── html/
    │   │   ├── rsRDM_all_summary.html        # Experimental summary
    │   │   ├── rsRDM_all_IR_drvvA_v_IR_WT_DE.html  # Individual comparisons
    │   │   └── rsRDM_all_ID_drvvA_v_ID_WT_DE.html
    │   └── pdf/                       # PDF versions
    ├── plots/                         # Individual plot files (if enabled)
    ├── data_tables/                   # Individual CSVs (renamed from results/)
    │   └── individual_comparisons/
    └── technical/                     # Technical files tucked away
        ├── R/                         # RDS files, report data
        │   └── rsRDM_all_report_data.RDS
        ├── metadata/                  # YAML files, coldata
        └── logs/                      # Processing logs
```

### Key Features
- **analysis_script.R stays at root** - no code changes needed, here() works fine
- **Important results front and center** - Excel file and overview.html prominent
- **Technical details hidden** - RDS files, logs in technical/ subfolder
- **Clear navigation** - overview.html serves as entry point

## Report Structure & Linking

### Overview Landing Page (`outputs/overview.html`)
**Purpose:** Main entry point for collaborators
**Content:**
- Experiment description and key parameters
- Quick summary statistics (samples, comparisons, significant genes)
- Links to:
  - Experimental summary report
  - Individual comparison reports
  - Main Excel file
- Generated from simple template with experiment metadata

### Experimental Summary Report (`reports/html/rsRDM_all_summary.html`)
**Content:**
- Sample metadata table
- QC plots (PCA, sample clustering, library size distributions)
- Heatmap of top variable genes
- Summary table of all comparisons (genes tested, significant hits)
- Links to individual comparison reports

### Individual DE Comparison Reports (`reports/html/rsRDM_all_{comparison}_DE.html`)
**Content:**
- Volcano plot with labeled top genes
- GSEA enrichment plots (if available)
- Top differentially expressed genes table
- Gene set enrichment tables
- Comparison metadata (treatment vs control, samples)

### Linking Strategy
**Consistent naming convention:**
- Overview: `overview.html`
- Experimental summary: `{experiment_name}_summary.html`
- DE comparisons: `{experiment_name}_{comparison_name}_DE.html`

**Navigation flow:**
```
overview.html → 
├── experimental_summary.html → individual_comparison_DE.html
├── individual_comparison_DE.html (direct links)
└── DEcompilation.xlsx (direct link)
```

## Dynamic Sizing Implementation

### Plot Dimension Calculation
Create helper functions for optimal sizing:
```r
# In report templates or helper functions
calculate_pathway_plot_height <- function(n_pathways) {
  base_height <- 3
  height_per_pathway <- 0.3
  max_height <- 12
  min(max_height, base_height + (n_pathways * height_per_pathway))
}

calculate_volcano_dimensions <- function(n_genes, n_labeled) {
  # Logic for optimal volcano plot sizing based on data density
}
```

### Conditional Content
Use conditional chunks in Quarto templates:
```r
#| eval: !expr nrow(gsea_results) > 0
# Only generate GSEA section if results exist

#| fig-height: !expr calculate_pathway_plot_height(nrow(significant_pathways))
# Dynamic plot sizing based on data
```

## Integration with analysis_script.R

### Pipeline Step Integration
```r
# Report generation as pipeline step
if (pipeline_step %in% c("full", "reports-only")) {
  cat("Generating reports...\n")
  
  # Load data if not in memory (for reports-only mode)
  if (!exists("report_data")) {
    report_data_path <- file.path(outputs_dir, "technical", "R", 
                                 paste0(experiment_name, "_report_data.RDS"))
    if (file.exists(report_data_path)) {
      report_data <- read_rds(report_data_path)
    } else {
      # Compile report data if not saved yet
      report_data <- compile_report_data(dds, res.l.all, contrast_list, 
                                       coldata, config, gsea_results)
      write_rds(report_data, report_data_path)
    }
  }
  
  # Generate all reports
  generate_all_reports(experiment_name, report_data, 
                      output_formats = config$reports$formats)
}
```

### Data Preparation
```r
# Save report data after analysis completion
if (pipeline_step %in% c("full", "analysis-only")) {
  # ... existing analysis code ...
  
  # Compile everything needed for reports
  cat("Preparing report data...\n")
  report_data <- list(
    dds = dds,
    res.l.all = res.l.all,
    contrast_list = contrast_list,
    coldata = coldata,
    config = config,
    gsea_results = list()
  )
  
  # Add GSEA results if they exist
  if (exists("gsea_gobp")) report_data$gsea_results$GObp <- gsea_gobp
  if (exists("gsea_gomf")) report_data$gsea_results$GOmf <- gsea_gomf
  
  # Save single report data file
  report_data_path <- file.path(outputs_dir, "technical", "R", 
                               paste0(experiment_name, "_report_data.RDS"))
  write_rds(report_data, report_data_path)
  cat(glue("Report data saved: {basename(report_data_path)}\n"))
}
```

## Template Organization

### Template Structure
```
templates/
├── reports/
│   ├── overview.qmd                    # Landing page template
│   ├── experimental_summary.qmd        # Main experiment overview
│   ├── comparison_report.qmd            # Individual DE comparison
│   ├── _common.qmd                     # Shared functions and setup
│   └── assets/
│       ├── styles.css                  # Custom styling
│       └── report_functions.R          # Report-specific helper functions
```

### Template Parameters
Keep parameters simple and focused:

**Overview Parameters:**
- `experiment_name`
- `data_dir` (path to technical/R/)

**Experimental Summary Parameters:**
- `experiment_name`
- `comparison_names` (vector)
- `data_dir`

**Comparison Report Parameters:**
- `experiment_name`
- `comparison_name`
- `data_dir`

### Data Loading in Templates
```r
# At top of each template
report_data <- read_rds(file.path(params$data_dir, 
                                 paste0(params$experiment_name, "_report_data.RDS")))

# Extract needed components
dds <- report_data$dds
res.l.all <- report_data$res.l.all
# etc.
```

## Implementation Functions

### Core Report Generation Functions
```r
# In functions/report_generation.R

generate_overview_page(experiment_name, report_data, output_formats = c("html"))
# Creates main landing page with links and summary stats

generate_experiment_summary_report(experiment_name, report_data, output_formats = c("html", "pdf"))
# Creates comprehensive experiment overview with QC and summary

generate_comparison_report(experiment_name, comparison_name, report_data, output_formats = c("html", "pdf"))
# Creates detailed individual comparison report

generate_all_reports(experiment_name, report_data, output_formats = c("html", "pdf"))
# Orchestrates generation of all report types

compile_report_data(dds, res.l.all, contrast_list, coldata, config, gsea_results = NULL)
# Compiles all data needed for reports into single list structure
```

## Implementation Timeline

### Phase 1: Core Infrastructure
1. Create report data compilation function
2. Build overview page template and generation
3. Create basic experimental summary template
4. Implement data loading patterns

### Phase 2: Detailed Reports
1. Build individual comparison report template
2. Implement dynamic sizing functions
3. Add GSEA visualization components
4. Create linking between reports

### Phase 3: Polish & Integration
1. Integrate with analysis_script.R pipeline
2. Add professional styling and formatting
3. Test with multiple experiment types
4. Optimize for collaborator usability

## Success Criteria

- **Reports generate automatically** from pipeline with no manual intervention
- **Dynamic sizing works** - plots adjust appropriately to data
- **Navigation works seamlessly** - all links functional between reports
- **Collaborator-friendly** - biologists can easily find and use results
- **Portable directories** - experiments can be compressed and shared
- **Multiple formats** - HTML and PDF both functional

## Benefits of This Design

1. **Minimal storage overhead** - single RDS file approach
2. **Robust data handling** - no complex object serialization issues
3. **Maintainable code** - clear separation of concerns
4. **Collaborator-friendly** - important results prominent, technical details hidden
5. **Flexible and extensible** - easy to add new report types
6. **Consistent with existing pipeline** - leverages current file-based architecture