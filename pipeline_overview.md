# Enhanced RNA-seq Analysis Pipeline - Complete Overview

## Project Philosophy

**Goal:** Enhance existing flexible RNA-seq workflow with professional automation, batch processing, and publication-ready outputs while preserving current filtering flexibility and Excel-based metadata management.

**Design Principle:** "Best of Both Worlds" - Keep scientific flexibility, add professional automation layer.

## Architecture: Hybrid Development + Automation System

The pipeline uses a **two-phase architecture** that supports both interactive development and automated batch processing:

### **Phase 1: Interactive Development**
- Work within experiment directories for custom filtering logic development
- Test and iterate on experimental design interactively
- Edit configuration scripts with experiment-specific parameters

### **Phase 2: Automated Execution** 
- Run from project root using orchestrator scripts
- Batch processing across multiple experiments
- Selective pipeline step execution
- Consistent path resolution and output management

## Core Components

### **1. Execution Orchestrators (Project Root)**

#### `run_experiment.R` - Single Experiment Controller
**Purpose:** Execute individual experiments with step control and working directory management

**Usage:**
```bash
# From project root only
Rscript run_experiment.R <experiment_name> [--step]
Rscript run_experiment.R rsRDM_all --dds-only
Rscript run_experiment.R rsRDM_all --analysis-only  
Rscript run_experiment.R rsRDM_all --reports-only
Rscript run_experiment.R rsRDM_all  # full pipeline
```

**Responsibilities:**
- Load enhanced functions from `functions/` directory
- Set environment variables for pipeline step control
- Change to experiment directory and source `analysis_script.R`
- Handle errors and restore working directory
- Coordinate with global configuration

#### `run_all.R` - Batch Processing Controller
**Purpose:** Process multiple experiments with selective regeneration

**Usage:**
```bash
# Batch processing options
Rscript run_all.R --full                    # Complete rebuild all experiments
Rscript run_all.R --dds-only               # Regenerate DDS for all experiments
Rscript run_all.R --analysis-only          # Analysis + reports for all
Rscript run_all.R --reports-only           # Reports only for all experiments
Rscript run_all.R --experiment rsRDM_all   # Single experiment processing
Rscript run_all.R --list                   # List available experiments
```

**Responsibilities:**
- Discover available experiments
- Coordinate pipeline steps across multiple experiments
- Report batch processing results
- Handle experiment dependencies and failures

### **2. Configuration Scripts (Experiment Directories)**

#### `experiments/{name}/analysis_script.R` - Experiment-Specific Logic
**Purpose:** Define experiment-specific filtering, contrasts, and parameters

**Key Sections:**
```r
# COLDATA FILTERING (customize per experiment)
coldata <- coldata.in %>%
  filter(experiment_group %in% c("iron")) %>%
  filter(!strain %in% c("dryhB")) %>%
  # Add your custom filtering logic here

# CONTRAST DEFINITIONS (customize per experiment)  
contrast_list <- list(
  c("group", "IR_drvvA", "IR_WT"),
  c("group", "ID_drvvA", "ID_WT")
)

# PIPELINE EXECUTION (uses loaded functions)
# Step 1: DDS Generation
# Step 2: Differential Expression
# Step 3: Report Generation
```

**Workflow:**
- Gets sourced by `run_experiment.R` from within experiment directory
- Receives pipeline step and experiment name via environment variables
- Uses enhanced functions loaded by orchestrator
- Focuses on experiment-specific logic only

### **3. Enhanced Function Library**

#### `functions/dds_management.R` - DDS Creation & Storage
**Key Functions:**
- `create_dds_from_salmon()` - Enhanced DDS creation with validation
- `save_dds_with_metadata()` - Save DDS with comprehensive metadata
- `load_dds_with_validation()` - Load DDS with helpful error messages
- `filter_low_counts()` - Enhanced filtering with reporting

#### `functions/differential_expression.R` - DESeq2 Analysis
This will contain functions similar to some of those in gk_deseq2.R. Specifically this script contains differential analysis helper functions. It will have functions helping with:

- generating res LFC for given contrasts
- creating significant (padj < 0.05) and differential expression (padj < 0.05 AND abs(log2FC) >= 1) subsets from differential expression lists
- joining resLFC data with pathway gene-sets (similar to rseq_geneSubsettor in gk_deseq2.R)
- merging gene annotation data onto resLFC and subsetted resLFC files
- creating the 'sigHits_compiled' tibble
  - this is a tibble where the log2FC values (and possibly other columns) from each individual RNAseq comparison for a given experiment are present IF the gene was significant (padj < 0.05) for that given comparison. If a gene was not significant for a given comparison, the log2FC value is left blank. This essentially creates a compilation where all significant hits can be viewed for all comparisons at once. Annotation data is also included. 
  - might also upgrade this to create a DEhits_compiled version
  - similar to deseq_sigCompiler from gk_deseq2.R
- Creating the master DEcompilation (a list of tibbles that will eventually be exported to excel)
  - each tibble in list is a separate excel workbook
  - included sheets (in approximate order)
    - sigHits_compiled is usually first sheet
    - individual pathway subsets (like sigHits_compiled but for chosen sets of genes i.e. pathways curated in a different file)
    - individual comparisons from an experiment - unfiltered (essentially just resLFC) and differential expressed (padj < 0.05 AND abs(log2FC) >= 1)
      - These individual unfiltered comparisons will likely also be saved on their own as .csv files elsewhere
- gene-set-enrichment-analysis function
  - likely very similar or copy of genesetEnrich from gk_deseq2.R

**Key Functions:** TBD / in-progress

- `generate_resLFC()` 
  - updated version of `deseq_result_generator()` from gk_deseq2.R
  - receives contrast and dds
  - Generates differential expresion results
  - performs lfcShrink using ashr as default
  - returns res.lfc
- In-progress

#### `functions/plotting_functions.R` - Publication Visualizations
**Key Functions:** TBD

#### `functions/output_management.R` - File Organization
**Key Functions:**

- `ensure_experiment_outputs()` - Create organized directory structure
- `save_experiment_results()` - Save results in multiple formats
- `load_experiment_results()` - Load saved results with validation

#### `functions/report_generation.R` - Automated Reports
**Key Functions:**
- `generate_comparison_reports()` - Individual comparison reports
- `generate_experiment_summary_report()` - Multi-comparison overview
- `generate_experiment_readme()` - Auto-generated documentation

### **4. Global Configuration**

#### `config/global_settings.yaml` - Shared Resources Only
**Purpose:** Paths, styling, and analysis defaults (no experiment-specific configs)

```yaml
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

## Directory Structure

```
enhanced-rnaseq-pipeline/
├── config/
│   └── global_settings.yaml              # Global resources & styling only
├── functions/                            # Modular function library
│   ├── dds_management.R                  # DDS creation, saving, loading
│   ├── differential_expression.R         # DESeq2 analysis functions
│   ├── annotation_processing.R           # Annotation handling
│   ├── plotting_functions.R              # Publication-ready visualizations
│   ├── output_management.R               # File organization system
│   └── report_generation.R               # Automated report generation
├── experiments/                          # Self-contained experiment analyses
│   ├── TEMPLATE/
│   │   └── analysis_script.R             # Template for new experiments
│   ├── rsRDM_all/                        # Example experiment
│   │   ├── analysis_script.R             # Custom filtering + contrast logic
│   │   └── outputs/
│   │       ├── rsRDM_all.dds.RDS         # DESeq2 dataset object
│   │       ├── rsRDM_all.metadata.yaml   # Analysis metadata
│   │       ├── rsRDM_all.coldata.csv     # Sample metadata
│   │       ├── results/                  # Excel files and CSV tables
│   │       │   ├── rsRDM_all_results.xlsx
│   │       │   ├── rsRDM_all_results.RDS
│   │       │   └── individual_comparisons/
│   │       ├── plots/                    # All plots (organized)
│   │       ├── reports/                  # PDF and HTML reports
│   │       ├── logs/                     # Processing logs
│   │       └── README.md                 # Auto-generated experiment summary
│   ├── rsRDM_rvvFur/                     # Another experiment
│   └── rvvA_all/                         # Additional experiments...
├── shared_resources/
│   ├── metadata/
│   │   └── rnaseq_yildiz_metadata.csv    # Master metadata file
│   └── annotations/
│       └── vc_genome_anno.RDS            # Genome annotations
├── templates/                            # Report templates
│   ├── enhanced_comparison.qmd           # Individual comparison template
│   └── experiment_summary.qmd            # Multi-comparison overview
├── run_experiment.R                      # Single experiment orchestrator
├── run_all.R                            # Batch processing orchestrator
├── new_experiment_setup.R                # Experiment creation utility
└── setup.R                             # One-time project setup
```

## Processing Pipeline

The pipeline uses a **three-step architecture** with selective execution:

### **Step 1: DDS Generation** (Expensive - Only When Needed)
- Load and validate sample metadata
- Create tximeta object from salmon quantification
- Run DESeq2 analysis with filtering
- Save DDS with comprehensive metadata

**When to Run:**
- New experiments
- Sample filtering changes
- Reference transcriptome updates

### **Step 2: Analysis & Results** (Moderate Cost - For Annotation Changes)
- Load saved DDS object
- Run differential expression analysis
- Process results with annotations
- Generate organized Excel/CSV outputs
- Create summary visualizations

**When to Run:**
- Annotation database updates
- New contrast definitions
- Analysis parameter changes

### **Step 3: Report Generation** (Fast - For Template Changes)
- Load saved results
- Generate publication-ready reports
- Create individual comparison PDFs
- Generate experiment summary documentation

**When to Run:**
- Template or styling updates
- Plot customization
- Report format changes

## Core Workflows

### **Setup (One-Time)**
```bash
# Initialize project structure and dependencies
Rscript setup.R

# Configure global settings
# Edit config/global_settings.yaml with your paths
```

### **Create New Experiment**
```bash
# Create experiment template
Rscript new_experiment_setup.R my_experiment_name

# Customize experiment configuration
# Edit experiments/my_experiment_name/analysis_script.R
# - Update coldata filtering logic
# - Define contrast list
# - Set any experiment-specific parameters
```

### **Interactive Development Workflow**
```bash
# Work within experiment directory for development
cd experiments/my_experiment_name

# Edit analysis_script.R interactively in RStudio
# Test filtering logic:
# View(coldata)
# table(coldata$strain, coldata$condition)

# When ready, run from project root
cd ../..
Rscript run_experiment.R my_experiment_name
```

### **Production Processing**

#### Single Experiment
```bash
# Full pipeline
Rscript run_experiment.R rsRDM_all

# Selective processing
Rscript run_experiment.R rsRDM_all --dds-only      # DDS generation only
Rscript run_experiment.R rsRDM_all --analysis-only  # Analysis + reports
Rscript run_experiment.R rsRDM_all --reports-only   # Reports only
```

#### Batch Processing
```bash
# Process all experiments
Rscript run_all.R --reports-only           # Regenerate all reports
Rscript run_all.R --analysis-only          # Rerun analysis for all
Rscript run_all.R --full                   # Complete rebuild

# Process specific experiment via batch system
Rscript run_all.R --experiment rsRDM_all --reports-only

# List available experiments
Rscript run_all.R --list
```

### **Common Update Scenarios**

#### Template/Styling Updates
```bash
# 1. Modify templates/enhanced_comparison.qmd or plotting functions
# 2. Regenerate all reports (fast)
Rscript run_all.R --reports-only
```

#### Annotation Database Updates  
```bash
# 1. Update shared_resources/annotations/vc_genome_anno.RDS
# 2. Regenerate analysis for all experiments (skips expensive DDS generation)
Rscript run_all.R --analysis-only
```

#### Function Enhancement Testing
```bash
# 1. Modify functions in functions/ directory
# 2. Test on one experiment first
Rscript run_all.R --experiment rsRDM_all --reports-only
# 3. When satisfied, apply to all
Rscript run_all.R --reports-only
```

## Key Features

### **Flexibility Preservation**
- **Excel metadata management** - Existing master metadata file workflow
- **Custom filtering logic** - Any R filtering approach per experiment
- **Dynamic contrast generation** - Each experiment defines its own comparisons
- **Experimental design freedom** - No rigid configuration constraints

### **Professional Automation**
- **Enhanced analysis functions** - Drop-in replacements with better functionality
- **Publication-ready reports** - Professional PDF/HTML outputs
- **Organized file structure** - Consistent output organization
- **Batch processing capabilities** - Selective regeneration across experiments

### **Development Support**
- **Interactive development** - Work within experiment directories
- **Modular function library** - Focused, testable components
- **Template system** - Consistent report generation
- **Error handling** - Helpful validation and debugging

### **Maintenance Benefits**
- **Template updates** - Automatically applied to all experiments
- **Function enhancements** - Immediately available across projects
- **Annotation propagation** - Easy database updates
- **Version control friendly** - Clear separation of code and configuration

## Migration Strategy

### **Phase 1: Setup & Core Functions (Week 1)**
1. Run `setup.R` to create directory structure
2. Configure `global_settings.yaml` with your paths
3. Implement core functions in `functions/` directory
4. Test with one experiment conversion

### **Phase 2: Template & Automation (Week 2)**
1. Create enhanced report templates
2. Implement batch processing scripts
3. Convert additional experiments
4. Validate consistency across experiments

### **Phase 3: Enhancement & Optimization (Week 3)**
1. Add advanced visualization features
2. Optimize batch processing performance
3. Create comprehensive documentation
4. Fine-tune professional styling

## Success Metrics

The enhanced system succeeds when it provides:

1. **Preserved Flexibility** - All current filtering approaches still work
2. **Professional Outputs** - Publication-ready reports automatically generated
3. **Time Savings** - Batch processing and automated organization
4. **Easy Maintenance** - Template/function updates applied system-wide
5. **Better Organization** - Consistent, professional file structure
6. **Enhanced Quality** - Better visualizations and reporting

This architecture maintains complete scientific flexibility while adding professional automation and sophisticated batch processing capabilities.