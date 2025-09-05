# RNA-seq Enhanced Analysis Pipeline - Complete Overview

## Project Philosophy

**Goal:** Enhance existing flexible RNA-seq workflow with professional automation, batch processing, and publication-ready outputs while preserving current filtering flexibility and Excel-based metadata management.

**Philosophy:** "Best of Both Worlds" - Keep scientific flexibility, add professional automation layer.

**Design Principle:** KISS (Keep It Simple, Stupid) - Build exactly what you need, no more, no less.

---

## Architecture Overview

### **Hybrid Enhancement System**

**What We Keep (Current Strengths):**
- ✅ **Excel/CSV metadata management** - Single master file with all samples  
- ✅ **Flexible R-based filtering** - Custom logic per experiment using any combination of columns  
- ✅ **Dynamic contrast generation** - Each experiment defines its own comparisons  
- ✅ **Full control over experimental design** - No rigid structure constraints  

**What We Add (Professional Automation):**
- ✅ **YAML for global resources only** - Paths, styling, annotation databases (no experiment configs)  
- ✅ **Enhanced analysis functions** - Drop-in replacements for current functions  
- ✅ **Professional report templates** - Publication-ready outputs with enhanced visualizations  
- ✅ **Organized output management** - Self-contained experiment folders with consistent structure  
- ✅ **Batch processing system** - Selective regeneration across all experiments  

---

## Directory Structure

```
gk-rnaseq-analyzer.Rproj
├── config/
│   └── global_settings.yaml           # Global resources & styling only
├── functions/
│   ├── rnaseq_analysis.R              # Core analysis functions (CRITICAL)
│   ├── rnaseq_plotting.R              # Visualization functions  
│   └── rnaseq_output_management.R     # Report generation and file organization
├── experiments/                       # Self-contained experimental analyses
│   ├── TEMPLATE/
│   │   └── analysis_script.R          # Template for new experiments
│   ├── rsRDM_all/                     # Example experiment
│   │   ├── analysis_script.R          # Custom filtering + contrast logic
│   │   └── outputs/
│   │       ├── rsRDM_all_dds.RDS      # DESeq2 dataset object
│   │       ├── rsRDM_all_metadata.yaml # Analysis metadata
│   │       ├── results/               # Excel files and CSV tables  
│   │       │   ├── rsRDM_all_results.xlsx
│   │       │   └── individual_comparisons/
│   │       ├── plots/                 # All plots (no subfolders)
│   │       ├── reports/               # PDF and HTML reports
│   │       ├── logs/                  # Processing logs
│   │       └── README.md              # Auto-generated experiment summary
│   └── rsRDM_rvvFur/                  # Another experiment
│       ├── analysis_script.R          # Different filtering logic
│       └── outputs/...
├── shared_resources/
│   ├── metadata/
│   │   └── rnaseq_yildiz_metadata.csv # Master metadata file
│   └── annotations/
│       └── vc_genome_anno.RDS         # Genome annotations
├── new_experiment_setup.R             # Experiment creation script
└── README.md                          # This documentation
```

---

## Setup Instructions

### **One-Time Setup**

1. **Initialize project structure:**
   ```bash
   Rscript setup.R
   ```

2. **Configure global settings:**
   Edit `config/global_settings.yaml` with your paths:
   ```yaml
   paths:
     metadata_file: "shared_resources/metadata/rnaseq_yildiz_metadata.csv"
     salmon_base: "../../../vc-refs/rnaseq_counts"
     annotation_db: "shared_resources/annotations/vc_genome_anno.RDS"
   ```

3. **Place your resources:**
   - Master metadata CSV in `shared_resources/metadata/`
   - Annotation RDS file in `shared_resources/annotations/`

4. **Create function files:**
   Create the three function files based on the Function Reference Guide

---

## Core Workflow

### **Creating New Experiments**

```bash
# Create new experiment with organized structure
Rscript new_experiment_setup.R my_experiment_name

# Results in:
# - experiments/my_experiment_name/analysis_script.R (customized)
# - experiments/my_experiment_name/outputs/ (full directory structure)
```

### **Customizing Analysis Script**

Each experiment has its own `analysis_script.R` with three customization sections:

#### **1. Coldata Filtering (Lines ~45-65)**
```r
# Customize this filtering logic for your experiment
coldata <- coldata.in %>%
  filter(experiment_group %in% c("iron")) %>%
  filter(!strain %in% c("dryhB")) %>%
  # Add your specific filtering logic here
  mutate(
    strain = fct_relevel(strain, "WT"),
    condition = fct_relevel(condition, "iron_replete"),
    group = fct_inorder(sampleLabel)
  )
```

#### **2. Contrast Definitions (Lines ~70-75)**
```r
# Define your comparisons
contrast_list <- list(
  c("group", "IR_drvvA", "IR_WT"),
  c("group", "ID_drvvA", "ID_WT")
  # Add more contrasts as needed
)
```

#### **3. Design Formula (Line ~85)**
```r
design_formula = ~group  # Customize for your experimental design
```

### **Running Analysis**

#### **Full Pipeline:**
```bash
cd experiments/my_experiment_name
Rscript analysis_script.R
```

#### **Selective Processing:**
```bash
# Generate DDS only (expensive, only when needed)
PIPELINE_STEP=dds-only Rscript analysis_script.R

# Analysis + reports (after DDS exists)
PIPELINE_STEP=analysis-only Rscript analysis_script.R

# Reports only (fast, for template changes)
PIPELINE_STEP=reports-only Rscript analysis_script.R
```

---

## Configuration Management

### **Global Settings (config/global_settings.yaml)**

```yaml
project:
  name: "RNA-seq Analysis Pipeline"
  author: "Your Name"
  
paths:
  metadata_file: "shared_resources/metadata/rnaseq_yildiz_metadata.csv"
  salmon_base: "../../../vc-refs/rnaseq_counts"
  annotation_db: "shared_resources/annotations/vc_genome_anno.RDS"

analysis:
  significance_threshold: 0.05
  fold_change_threshold: 1.0
  min_counts: 2

plotting:
  base_size: 8
  figure_width: 7
  figure_height: 5
  dpi: 300

reports:
  generate_pdf: true
  generate_html: false
  save_plots_separately: true
```

### **Experiment-Specific Overrides**

In each `analysis_script.R`, override global settings as needed:

```r
# Load global config
config <- read_yaml(here("config", "global_settings.yaml"))

# Override for this experiment
config$analysis$min_counts <- 1        # Lower threshold
config$plotting$base_size <- 10        # Larger text
```

---

## Function Overview

### **Core Analysis (rnaseq_analysis.R)**
- `create_dds_from_salmon()` - **CRITICAL** - Replace current DDS workflow
- `enhanced_deseq_makeOutput()` - **CRITICAL** - Drop-in replacement for current function
- `prepare_annotations()` - **CRITICAL** - Standardize annotation data
- `save_dds_with_metadata()` / `load_dds_with_validation()` - File management
- `save_experiment_results()` / `load_experiment_results()` - Result organization

### **Visualization (rnaseq_plotting.R)**
- `setup_plot_theme()` - Consistent theming
- `enhanced_volcano_plot()` - Publication-ready volcano plots
- `enhanced_pathway_plot()` - Professional pathway enrichment plots
- `save_plot_professional()` - Multi-format plot saving

### **Output Management (rnaseq_output_management.R)**
- `generate_comparison_reports()` - Individual comparison reports
- `generate_experiment_summary_report()` - Multi-comparison overview
- `generate_experiment_readme()` - Auto-generated documentation

---

## Generated Outputs

### **Per Experiment Structure:**
```
experiments/my_experiment/outputs/
├── my_experiment_dds.RDS              # DESeq2 object
├── my_experiment_metadata.yaml        # Analysis metadata
├── results/
│   ├── my_experiment_results.xlsx     # Combined Excel workbook
│   ├── my_experiment_results.RDS      # Programmatic access
│   └── individual_comparisons/        # Individual CSV files
├── plots/                             # All plots (no subfolders)
│   ├── treatment_vs_control_volcano.pdf
│   ├── treatment_vs_control_volcano.png
│   └── pathway_enrichment.pdf
├── reports/                           # Generated reports
│   ├── DGE_treatment_vs_control.pdf
│   └── experiment_summary.pdf
├── logs/                              # Processing logs
└── README.md                          # Auto-generated summary
```

---

## Development Strategy

### **Phase 1: Core Functions (Priority)**
1. **Create `rnaseq_analysis.R`** with essential functions
2. **Test with real dataset** using existing experiment
3. **Validate drop-in compatibility** with current workflow
4. **Create minimal placeholder** plotting and output functions

### **Phase 2: Enhanced Features**
1. **Add professional plotting** functions
2. **Implement report generation**
3. **Add comprehensive error handling**
4. **Test batch processing capabilities**

### **Phase 3: Polish & Automation**
1. **Optimize performance** for large datasets
2. **Add advanced visualization** features
3. **Create comprehensive documentation**
4. **Implement quality control** dashboards

---

## Migration Strategy

### **From Current Workflow:**

#### **Week 1: Foundation**
1. **Set up enhanced system** - Run `setup.R`, configure paths
2. **Create test experiment** - Use `new_experiment_setup.R` with existing data
3. **Build core functions** - Start with minimal versions of critical functions
4. **Validate compatibility** - Ensure same results as current workflow

#### **Week 2: Enhancement**
1. **Add enhanced features** - Better plotting, error handling, organization
2. **Test selective processing** - Validate DDS-only, analysis-only, reports-only modes
3. **Create additional experiments** - Test with different datasets
4. **Refine templates** - Based on real usage patterns

#### **Week 3: Automation**
1. **Implement report generation** - Professional PDF/HTML outputs
2. **Add batch processing** - Cross-experiment capabilities
3. **Create comprehensive documentation** - User guides and troubleshooting
4. **Performance optimization** - For larger datasets

---

## Interactive Development Workflow

### **Testing Coldata Filtering:**
```r
# Open analysis_script.R in RStudio
# Run setup sections (lines 1-35)
# Run filtering section interactively (lines 45-65)

# Test your filtering:
View(coldata)
table(coldata$strain, coldata$condition)

# Once satisfied, run full pipeline
```

### **Function Development:**
```r
# Start with existing functions as drop-in replacements
# Test each function individually with real data
# Add enhanced features incrementally
# Maintain backward compatibility
```

---

## Key Features

### **Safety & Reliability**
- **Auto-detection** of experiment names prevents accidental overwrites
- **Comprehensive validation** with helpful error messages  
- **Metadata tracking** for full reproducibility
- **Selective processing** to avoid expensive re-computation

### **Professional Outputs**
- **Publication-ready reports** with consistent formatting
- **Enhanced visualizations** with pathway coloring and professional styling
- **Organized file structure** for easy navigation and sharing
- **Multiple output formats** (PDF, HTML, Excel, CSV, RDS)

### **Flexible Workflow**
- **Preserve current filtering approach** - No rigid configuration requirements
- **Custom contrast definitions** per experiment
- **Experiment-specific overrides** for analysis parameters
- **Incremental adoption** - Replace components as validated

### **Batch Processing**
- **Template updates** automatically applied to all experiments
- **Function enhancements** immediately available across projects  
- **Annotation updates** easily propagated
- **Selective regeneration** (reports-only, analysis-only, etc.)

---

## Usage Examples

### **Daily Development:**
```bash
# 1. Modify template or function
# 2. Test on one experiment
cd experiments/test_experiment
PIPELINE_STEP=reports-only Rscript analysis_script.R

# 3. When satisfied, apply to all experiments
# (batch processing commands to be implemented)
```

### **New Annotation Database:**
```bash
# 1. Update shared_resources/annotations/vc_genome_anno.RDS
# 2. Regenerate analysis for all experiments (skip expensive DDS generation)
cd experiments/experiment1 && PIPELINE_STEP=analysis-only Rscript analysis_script.R
cd experiments/experiment2 && PIPELINE_STEP=analysis-only Rscript analysis_script.R
# etc.
```

### **Template Styling Update:**
```bash
# 1. Modify templates or plotting functions
# 2. Regenerate all reports
cd experiments/experiment1 && PIPELINE_STEP=reports-only Rscript analysis_script.R
cd experiments/experiment2 && PIPELINE_STEP=reports-only Rscript analysis_script.R
# etc.
```

---

## Troubleshooting

### **Common Issues**

#### **Missing Function Errors**
```
Error: could not find function "create_dds_from_salmon"
```
**Solution:** Create the function files (`rnaseq_analysis.R`, etc.) or start with minimal placeholder functions.

#### **File Path Issues**
```
Error: cannot open file 'config/global_settings.yaml'
```
**Solution:** Ensure you're working from the R project root and `here::i_am()` is correctly set.

#### **Metadata Loading Issues**
```
Error: object 'experiment_folder' not found
```
**Solution:** Check that your metadata CSV has the expected column names and your filtering logic matches your actual data structure.

### **Development Tips**

1. **Start simple** - Get basic workflow working before adding enhancements
2. **Test incrementally** - Validate each function with real data
3. **Use interactive development** - Run sections of analysis_script.R interactively
4. **Maintain compatibility** - Keep existing function interfaces when possible
5. **Document changes** - Update templates based on real usage patterns

---

## Success Metrics

The enhanced system succeeds when it provides:

1. **Preserved flexibility** - All current filtering approaches still work
2. **Professional outputs** - Publication-ready reports automatically generated
3. **Time savings** - Batch processing and automated organization
4. **Easy maintenance** - Template/function updates applied system-wide
5. **Better organization** - Consistent, professional file structure
6. **Enhanced quality** - Better visualizations and reporting

This hybrid approach maintains complete scientific flexibility while adding professional automation and organization capabilities.