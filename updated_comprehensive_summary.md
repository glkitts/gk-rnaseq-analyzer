# RNA-seq Enhanced Analysis Pipeline - Design Summary

## Project Overview

**Goal:** Enhance existing flexible RNA-seq workflow with professional automation, batch processing, and publication-ready outputs while preserving current filtering flexibility and Excel-based metadata management.

**Philosophy:** "Best of Both Worlds" - Keep scientific flexibility, add professional automation layer.

**Current Challenge:** Functional Excel-based workflow with flexible R filtering but lacks automated report generation, organized outputs, and batch processing capabilities for multiple experiments.

## Architecture Decision: Hybrid Enhancement System

### **Layer 1: Resource Management**
- **YAML for global settings only** - Paths, styling, annotation databases
- **Excel/CSV metadata system** - Existing master metadata file (preserved)
- **No rigid experiment configs** - Maintain current filtering flexibility

### **Layer 2: Enhanced Analysis Functions**  
- **Drop-in function replacements** - Enhanced versions of existing functions
- **Preserve current workflow** - Same interfaces, better outputs
- **Professional visualizations** - Publication-ready plots with consistent styling
- **Organized file management** - Automated output structure

### **Layer 3: Automated Processing**
- **Self-contained experiments** - Each experiment gets organized folder structure
- **Flexible batch processing** - Selective regeneration options
- **Enhanced report templates** - Professional PDF/HTML outputs  
- **Multi-level processing control** - DDS generation, analysis, reports

## Key Design Decisions

### **Workflow Preservation**
✅ **RECOMMENDED:** Keep current Excel metadata + R filtering approach
- Maintain existing flexible filtering logic per experiment
- No learning curve for new configuration formats
- Preserve ability to filter on any combination of columns
- Support complex experimental designs

### **Enhancement Strategy** 
✅ **RECOMMENDED:** Drop-in function replacements + automation
- Enhanced versions of existing functions (same interfaces)
- Professional report templates with better styling
- Automated output organization and batch processing
- Gradual adoption path - replace components one at a time

### **Processing Architecture**
✅ **RECOMMENDED:** Three-step pipeline with selective execution
- **Step 1:** DDS Generation (expensive, only when needed)
- **Step 2:** Analysis & Results (moderate cost, rerun for annotation changes)
- **Step 3:** Report Generation (cheap, rerun for template changes)

### **Batch Processing Options**
✅ **RECOMMENDED:** Multi-level selective regeneration
- Full pipeline rebuild (new reference transcriptome)
- DDS regeneration (sample filtering changes)
- Analysis only (annotation updates)
- Reports only (template/styling changes)
- Single experiment testing (validate changes)

## Directory Structure

```
rnaseq-enhanced-pipeline/
├── config/
│   └── global_settings.yaml           # Global resources & styling only
├── functions/
│   ├── enhanced_analysis.R            # Drop-in replacements for current functions
│   ├── enhanced_plotting.R            # Professional publication plots
│   ├── output_management.R            # Automated file organization
│   └── report_automation.R            # Automated report generation
├── templates/
│   ├── enhanced_comparison.qmd        # Enhanced version of current template
│   └── experiment_summary.qmd         # Multi-comparison overview
├── experiments/                       # Self-contained experimental analyses
│   ├── rsRDM_all/                    # Each experiment gets own folder
│   │   ├── analysis_script.R         # Custom filtering + contrast logic
│   │   ├── outputs/
│   │   │   ├── dds/                  # DDS objects
│   │   │   ├── results/              # Excel files, tables
│   │   │   ├── plots/                # All plots organized by type
│   │   │   ├── reports/              # PDF/HTML reports
│   │   │   └── README.md             # Auto-generated experiment summary
│   │   └── logs/                     # Processing logs and metadata
│   ├── rsRDM_rvvFur/
│   │   ├── analysis_script.R         # Different filtering logic
│   │   └── outputs/...
│   └── rvvA_all/
│       ├── analysis_script.R
│       └── outputs/...
├── shared_resources/
│   ├── metadata/
│   │   └── rnaseq_yildiz_metadata.csv # Master metadata file
│   └── annotations/
│       └── vc_genome_anno.RDS         # Genome annotations
├── run_all.R                         # Batch processing controller
├── run_experiment.R                  # Single experiment runner
└── setup.R                          # One-time setup script
```

## Core Workflow

### **Setup (One-time):**
```bash
# Install dependencies and create directory structure
Rscript setup.R

# Update config/global_settings.yaml with your paths and preferences
```

### **Create New Experiment:**
```bash
# Copy template analysis script and modify filtering logic (5 minutes)
cp experiments/TEMPLATE/analysis_script.R experiments/my_experiment/
# Edit filtering logic and contrast list in analysis_script.R
```

### **Single Experiment Processing:**
```bash
# Full pipeline for one experiment
Rscript run_experiment.R my_experiment

# Or step-by-step:
Rscript run_experiment.R my_experiment --dds-only      # Generate DDS only
Rscript run_experiment.R my_experiment --analysis-only  # Analysis + reports
Rscript run_experiment.R my_experiment --reports-only   # Reports only
```

### **Batch Processing (Key Feature):**
```bash
# Template updates → regenerate all reports
Rscript run_all.R --reports-only

# Annotation updates → rerun analysis for all experiments
Rscript run_all.R --analysis-only

# New reference transcriptome → full rebuild
Rscript run_all.R --full

# Test changes on one experiment first
Rscript run_all.R --experiment rsRDM_all --reports-only
```

## Key Enhancements Over Current System

### **Workflow Management**
- ✅ **Flexible batch processing** with selective regeneration options
- ✅ **Self-contained experiments** with organized output structure
- ✅ **Professional automation** while preserving current filtering flexibility
- ✅ **Multi-step pipeline** with intelligent dependency management
- ✅ **Easy template/function updates** applied across all experiments

### **Analysis Quality**
- ✅ **Enhanced visualization functions** with publication-ready styling
- ✅ **Drop-in function replacements** with better error handling
- ✅ **Consistent output organization** across all experiments
- ✅ **Professional report templates** with enhanced formatting
- ✅ **Automated quality control** and validation reporting

### **Output Generation**
- ✅ **Publication-ready PDF reports** with professional styling
- ✅ **Interactive HTML reports** with plotly integration (optional)
- ✅ **Enhanced visualizations** - pathway-colored volcano plots, professional GO analysis
- ✅ **Organized file structure** - plots, tables, reports automatically organized
- ✅ **Multi-format outputs** - PDF, PNG plots plus Excel/CSV data tables

### **Maintenance & Updates**
- ✅ **Template updates** automatically applied to all experiments
- ✅ **Function enhancements** immediately available across projects  
- ✅ **Annotation database updates** easily propagated
- ✅ **Consistent styling** and formatting across all outputs
- ✅ **Version control friendly** - track changes to templates and functions

## Processing Examples

### **Daily Development Workflow:**
```bash
# 1. Modify template or function
# 2. Test on one experiment
Rscript run_all.R --experiment rsRDM_all --reports-only
# 3. When satisfied, apply to all
Rscript run_all.R --reports-only
```

### **New Annotation Database:**
```bash  
# 1. Update shared_resources/annotations/vc_genome_anno.RDS
# 2. Regenerate analysis for all experiments (skip expensive DDS generation)
Rscript run_all.R --analysis-only
# ✅ All experiments get updated annotations, pathway analysis, gene names
```

### **Template Styling Update:**
```bash
# 1. Modify templates/enhanced_comparison.qmd
# 2. Regenerate all reports
Rscript run_all.R --reports-only
# ✅ All experiments get updated formatting in minutes
```

## Configuration Examples

### **Global Settings (`config/global_settings.yaml`):**
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
  generate_html: true
  save_plots_separately: true
```

### **Experiment Analysis Script (Example - Keep Current Logic):**
```r
#!/usr/bin/env Rscript

# Load enhanced functions (drop-in replacements)
source("../../functions/enhanced_analysis.R")

# YOUR CURRENT FILTERING LOGIC (keep exactly as is)
coldata <- coldata.in %>%
  filter(experiment_group %in% c("iron")) %>%
  filter(!strain %in% c("dryhB")) %>%
  filter(temperature == 37) %>%  # Your custom filtering
  # ... your current logic

# YOUR CURRENT CONTRAST LIST (keep exactly as is)
contrast_list <- list(
  c("group", "IR_drvvA", "IR_WT"),
  c("group", "ID_drvvA", "ID_WT")
)

# ENHANCED PROCESSING (replaces your deseq_makeOutput)
results <- enhanced_deseq_makeOutput(contrast_list, dds, annotations, experiment_name)

# AUTOMATED REPORT GENERATION (new capability)
generate_all_reports(results, dds, experiment_name)
```

## Migration Strategy

### **From Current Workflow:**
1. **Set up enhanced system** - Run setup.R, configure paths
2. **Test with one experiment** - Copy current analysis, use enhanced functions
3. **Validate outputs** - Compare results with current workflow
4. **Gradual adoption** - Replace functions one at a time as validated
5. **Scale to batch processing** - Apply enhancements across all experiments

### **Enhancement Timeline:**
- **Week 1:** Enhanced functions + basic templates
- **Week 2:** Batch processing + output organization  
- **Week 3:** Professional styling + advanced features

## Success Metrics

The enhanced system succeeds when it provides:
1. **Preserved flexibility** - All current filtering approaches still work
2. **Professional outputs** - Publication-ready reports automatically generated
3. **Time savings** - Batch processing and automated organization
4. **Easy maintenance** - Template/function updates applied system-wide
5. **Better organization** - Consistent, professional file structure
6. **Enhanced quality** - Better visualizations and reporting

This design maintains complete scientific flexibility while adding professional automation and batch processing capabilities.