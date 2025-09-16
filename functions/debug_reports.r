#!/usr/bin/env Rscript

# Debug script for report generation
# Run this to test components before full report generation

library(here)
library(cli)

# 1. Check if templates directory exists
templates_dir <- here("templates")
if (!dir.exists(templates_dir)) {
  cli_alert_danger("Templates directory missing: {templates_dir}")
  cli_alert_info("Create it with: dir.create('{templates_dir}')")
} else {
  cli_alert_success("Templates directory exists")
}

# 2. Check for required template files
required_templates <- c("experimental_summary.qmd", "comparison_report.qmd")
for (template in required_templates) {
  template_path <- file.path(templates_dir, template)
  if (!file.exists(template_path)) {
    cli_alert_danger("Template missing: {template}")
    cli_alert_info("Create {template} in templates/ directory")
  } else {
    cli_alert_success("Template found: {template}")
  }
}

# 3. Check if Quarto is available
cli_alert_info("Checking Quarto installation...")
quarto_check <- system("quarto --version", intern = TRUE, ignore.stderr = TRUE)
if (length(quarto_check) > 0) {
  cli_alert_success("Quarto found: {quarto_check[1]}")
} else {
  cli_alert_danger("Quarto CLI not found")
  cli_alert_info("Install from: https://quarto.org")
}

# 4. Check required packages
required_packages <- c("tidyverse", "DESeq2", "ComplexHeatmap", "DT", "knitr", "scales", "ggrepel", "viridis")
optional_packages <- c("plotly", "kableExtra")

cli_alert_info("Checking required packages...")
missing_required <- character(0)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_required <- c(missing_required, pkg)
  }
}

if (length(missing_required) > 0) {
  cli_alert_danger("Missing required packages: {paste(missing_required, collapse = ', ')}")
  cli_alert_info("Install with: install.packages(c({paste(paste0('\"', missing_required, '\"'), collapse = ', ')}))")
} else {
  cli_alert_success("All required packages available")
}

cli_alert_info("Checking optional packages...")
missing_optional <- character(0)
for (pkg in optional_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_optional <- c(missing_optional, pkg)
  }
}

if (length(missing_optional) > 0) {
  cli_alert_warning("Missing optional packages: {paste(missing_optional, collapse = ', ')}")
  cli_alert_info("These packages enable enhanced features but are not required")
} else {
  cli_alert_success("All optional packages available")
}

# 5. Test basic template rendering (if everything is available)
if (dir.exists(templates_dir) && 
    file.exists(file.path(templates_dir, "experimental_summary.qmd")) &&
    length(quarto_check) > 0 &&
    length(missing_required) == 0) {
  
  cli_alert_info("Testing basic template rendering...")
  
  # Create minimal test data
  test_data <- list(
    experiment = list(
      metadata = list(
        n_samples = 6,
        n_transcripts = 1000,
        n_comparisons = 2,
        comparison_names = c("test1", "test2"),
        design_formula = "~group"
      ),
      coldata = data.frame(
        sample_id = paste0("sample_", 1:6),
        group = rep(c("A", "B"), each = 3),
        strain = "WT"
      ),
      analysis_date = Sys.Date()
    ),
    plotting = list(
      comparison_summary = data.frame(
        comparison = c("test1", "test2"),
        total_genes = c(1000, 1000),
        sig_genes = c(50, 75),
        de_genes = c(25, 40),
        max_log2fc = c(2.5, 3.1),
        min_padj = c(0.001, 0.0005)
      )
    )
  )
  
  # Try to save test data and render
  tryCatch({
    test_dir <- tempdir()
    test_data_file <- file.path(test_dir, "test_data.RDS")
    saveRDS(test_data, test_data_file)
    
    output_file <- file.path(test_dir, "test_summary.html")
    
    # Simple quarto command
    cmd <- glue::glue(
      "quarto render '{file.path(templates_dir, 'experimental_summary.qmd')}' ",
      "--output '{output_file}' ",
      "--to html ",
      "-P experiment_name:'TEST' ",
      "-P temp_data_file:'{test_data_file}' ",
      "-P output_format:'html'"
    )
    
    result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    
    if (file.exists(output_file)) {
      cli_alert_success("Test rendering successful!")
      cli_alert_info("Test file created: {output_file}")
    } else {
      cli_alert_warning("Test rendering failed - output file not created")
    }
    
    # Clean up
    unlink(c(test_data_file, output_file))
    
  }, error = function(e) {
    cli_alert_warning("Test rendering failed: {e$message}")
  })
} else {
  cli_alert_info("Skipping test rendering - missing requirements")
}

cli_rule("Debug Summary")
cli_alert_info("If all checks pass, try running reports again")
cli_alert_info("If issues persist, check Quarto documentation: https://quarto.org/docs/")
