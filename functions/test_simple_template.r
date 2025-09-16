#!/usr/bin/env Rscript

# Test simple template rendering
# Save this as functions/test_simple_template.R and run it

library(here)
library(cli)

test_simple_template <- function() {
  
  cli_alert_info("Testing simple Quarto template...")
  
  # Check if simple template exists
  simple_template <- here("templates", "simple_experimental_summary.qmd")
  
  if (!file.exists(simple_template)) {
    cli_alert_danger("Simple template not found: {simple_template}")
    cli_alert_info("Create templates/simple_experimental_summary.qmd first")
    return(FALSE)
  }
  
  # Create test output directory
  test_dir <- here("test_output")
  if (!dir.exists(test_dir)) {
    dir.create(test_dir)
  }
  
  # Try basic rendering without parameters
  cli_alert_info("Testing basic rendering...")
  
  output_file <- file.path(test_dir, "test_basic.html")
  
  cmd <- paste(
    "quarto render",
    shQuote(simple_template),
    "--output", shQuote(output_file),
    "--to html"
  )
  
  tryCatch({
    result <- system(cmd, intern = TRUE)
    
    if (file.exists(output_file)) {
      cli_alert_success("Basic rendering works!")
      cli_alert_info("Test file: {output_file}")
      
      # Try with parameters
      cli_alert_info("Testing with parameters...")
      
      output_file2 <- file.path(test_dir, "test_params.html")
      
      cmd2 <- paste(
        "quarto render",
        shQuote(simple_template),
        "--output", shQuote(output_file2),
        "--to html",
        "-P experiment_name:TEST_EXPERIMENT"
      )
      
      result2 <- system(cmd2, intern = TRUE)
      
      if (file.exists(output_file2)) {
        cli_alert_success("Parameter passing works!")
        cli_alert_info("Test file: {output_file2}")
        return(TRUE)
      } else {
        cli_alert_warning("Parameter test failed")
        return(FALSE)
      }
      
    } else {
      cli_alert_danger("Basic rendering failed")
      cli_alert_info("Command: {cmd}")
      return(FALSE)
    }
    
  }, error = function(e) {
    cli_alert_danger("Error during testing: {e$message}")
    return(FALSE)
  })
}

# Run the test
if (!interactive()) {
  result <- test_simple_template()
  if (result) {
    cli_alert_success("Template testing passed!")
  } else {
    cli_alert_danger("Template testing failed!")
  }
}