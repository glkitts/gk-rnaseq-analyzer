# Report-specific helper functions
# Based on gk_functions.R reference standards

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(plotly)
})

# ============================================================================= #
# CONSTANTS AND THEME SETUP ----
# ============================================================================= #

# Font and sizing constants (from gk_functions.R)
base_family <- "Helvetica"
base_size <- 8
base_line_size <- base_size / 22
base_rect_size <- base_size / 22
half_line <- base_size / 2
base.geomsize <- base_size / 8

# Text and point sizing
txt.mm_to_pts <- 25.4 / 72.27
geom.text_size <- base_size * 0.8 * txt.mm_to_pts
geom.point_size <- 1.2

# Plot styling constants
plt.linetype <- "dashed"
plt.color <- "black"
plt.alpha <- 0.2

# ============================================================================= #
# THEME FUNCTIONS ----
# ============================================================================= #

#' Publication-quality theme based on gk_functions.R
theme_report <- function(base_size = 8,
                        base_family = "Helvetica",
                        grid_lines = "none",
                        grid_color = "gray92") {

  # Calculate base parameters
  half_line <- base_size / 2
  base_line_size <- base_size / 22

  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # Plot elements
      plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0),
      plot.subtitle = element_text(size = rel(1), color = "gray30"),
      plot.tag = element_text(face = "bold", size = rel(1.25)),
      plot.margin = margin(half_line, half_line, half_line, half_line),

      # Panel and background
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = base_line_size),

      # Grid lines
      panel.grid.major.x = if (grid_lines %in% c("major", "both")) {
        element_line(color = grid_color, linewidth = base_line_size * 0.5)
      } else {
        element_blank()
      },
      panel.grid.major.y = if (grid_lines %in% c("major", "both", "y_major")) {
        element_line(color = grid_color, linewidth = base_line_size * 0.5)
      } else {
        element_blank()
      },
      panel.grid.minor = if (grid_lines %in% c("minor", "both")) {
        element_line(color = grid_color, linewidth = base_line_size * 0.25)
      } else {
        element_blank()
      },

      # Axes
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = base_line_size),
      axis.text = element_text(color = "black", size = rel(0.9)),
      axis.title = element_text(color = "black", size = rel(1)),

      # Legend
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.background = element_rect(color = NA, fill = "white"),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = rel(0.9)),
      legend.title = element_text(size = rel(1)),

      # Strip for facets
      strip.background = element_rect(fill = "gray95", color = "gray80"),
      strip.text = element_text(size = rel(0.9), face = "plain")
    )
}

# ============================================================================= #
# DATA PROCESSING FUNCTIONS ----
# ============================================================================= #

#' Handle padj values of 0 for volcano plots (from gk_deseq2.R)
#' @param res Results data frame with padj column
#' @return Data frame with padj=0 values replaced
replace_padj_0 <- function(res) {
  vp.minpv <- res %>%
    filter(padj != 0, !is.na(padj)) %>%
    pull(padj) %>%
    min()

  pval_0 <- res %>%
    filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
    filter(padj == 0) %>%
    select(any_of(c("gene_id", "Label")))

  if (nrow(pval_0) > 0) {
    pval_0$padj2 <- 1 * 10^-rnorm(nrow(pval_0), mean = -log10(vp.minpv), sd = 1)

    res.out <- res %>%
      filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
      left_join(pval_0, by = names(pval_0)[1]) %>%
      mutate(padj = if_else(padj == 0 & !is.na(padj2), padj2, padj)) %>%
      select(-any_of("padj2"))
  } else {
    res.out <- res %>%
      filter(!is.na(padj) & !is.na(log2FoldChange))
  }

  return(res.out)
}

#' Calculate symmetric axis limits for volcano plots
#' @param x Vector of log2FoldChange values
#' @param buffer_factor Factor to add buffer space (default 1.1)
#' @return Symmetric limits c(-max_abs, max_abs)
calculate_symmetric_limits <- function(x, buffer_factor = 1.1) {
  max_abs <- max(abs(x), na.rm = TRUE) * buffer_factor
  return(c(-max_abs, max_abs))
}

#' Determine genes to label in volcano plot
#' @param res_sig Significant results data frame
#' @param n_up Number of upregulated genes to label (default 5)
#' @param n_down Number of downregulated genes to label (default 5)
#' @return Vector of gene names/labels to highlight
get_volcano_labels <- function(res_sig, n_up = 5, n_down = 5) {
  # Try different label columns in order of preference
  label_cols <- c("yildiz_lbl", "label_rows", "Label", "gene_name", "gene_id")
  label_col <- NULL

  for (col in label_cols) {
    if (col %in% colnames(res_sig)) {
      label_col <- col
      break
    }
  }

  if (is.null(label_col)) {
    warning("No suitable label column found")
    return(character(0))
  }

  top_up <- res_sig %>%
    filter(log2FoldChange > 0) %>%
    arrange(padj) %>%
    head(n_up) %>%
    pull(!!sym(label_col))

  top_down <- res_sig %>%
    filter(log2FoldChange < 0) %>%
    arrange(padj) %>%
    head(n_down) %>%
    pull(!!sym(label_col))

  return(c(top_up, top_down))
}

# ============================================================================= #
# PLOTTING FUNCTIONS ----
# ============================================================================= #

#' Create publication-quality volcano plot
#' @param res_all All results data frame
#' @param res_sig Significant results data frame
#' @param comparison_name Name of comparison for title
#' @param sig_threshold Significance threshold (default 0.05)
#' @param fc_threshold Fold change threshold for labeling (default 1)
#' @param interactive Create interactive plotly version (default FALSE)
#' @return ggplot object or plotly object
create_volcano_plot <- function(res_all, res_sig, comparison_name,
                               sig_threshold = 0.05, fc_threshold = 1,
                               interactive = FALSE) {

  # Process data
  volcano_data <- res_all %>%
    replace_padj_0() %>%
    mutate(
      log_padj = -log10(padj),
      significant = ifelse(padj < sig_threshold, "Significant", "Not Significant"),
      direction = case_when(
        padj < sig_threshold & log2FoldChange > fc_threshold ~ "Up",
        padj < sig_threshold & log2FoldChange < -fc_threshold ~ "Down",
        padj < sig_threshold ~ "Significant",
        TRUE ~ "Not Significant"
      )
    )

  # Get genes to label
  genes_to_label <- get_volcano_labels(res_sig)

  # Determine label column
  label_cols <- c("yildiz_lbl", "label_rows", "Label", "gene_name", "gene_id")
  label_col <- NULL
  for (col in label_cols) {
    if (col %in% colnames(volcano_data)) {
      label_col <- col
      break
    }
  }

  if (!is.null(label_col)) {
    volcano_data <- volcano_data %>%
      mutate(
        plot_label = ifelse(!!sym(label_col) %in% genes_to_label, !!sym(label_col), ""),
        hover_label = !!sym(label_col)
      )
  } else {
    volcano_data <- volcano_data %>%
      mutate(
        plot_label = "",
        hover_label = row_number()
      )
  }

  # Calculate symmetric limits
  x_limits <- calculate_symmetric_limits(volcano_data$log2FoldChange)

  # Create color scheme
  colors <- c("Up" = "#E31A1C", "Down" = "#1F78B4",
              "Significant" = "#33A02C", "Not Significant" = "gray70")

  # Create base plot
  p <- ggplot(volcano_data, aes(x = log2FoldChange, y = log_padj)) +
    geom_point(aes(color = direction, text = hover_label),
               alpha = 0.6, size = 1) +
    scale_color_manual(values = colors, name = "Direction") +
    scale_x_continuous(limits = x_limits) +
    geom_hline(yintercept = -log10(sig_threshold),
               linetype = "dashed", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 0,
               linetype = "solid", color = "black", alpha = 0.3) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold),
               linetype = "dotted", color = "gray50", alpha = 0.5) +
    labs(
      title = paste("Volcano Plot:", comparison_name),
      x = "log₂ Fold Change",
      y = "-log₁₀(Adjusted P-value)",
      subtitle = paste(nrow(res_sig), "significant genes out of",
                      format(nrow(volcano_data), big.mark = ","), "tested")
    ) +
    theme_report()

  # Add labels for top genes
  if (length(genes_to_label) > 0 && !is.null(label_col)) {
    p <- p +
      geom_text_repel(
        aes(label = plot_label),
        size = geom.text_size,
        max.overlaps = 15,
        box.padding = 0.3,
        point.padding = 0.3,
        min.segment.length = 0,
        color = "black"
      )
  }

  if (interactive) {
    p <- ggplotly(p, tooltip = c("x", "y", "text", "colour")) %>%
      layout(
        title = list(text = paste("Volcano Plot:", comparison_name),
                    font = list(size = 14)),
        showlegend = TRUE
      )
  }

  return(p)
}

#' Calculate dynamic height for pathway plots
#' @param n_pathways Number of pathways to display
#' @param base_height Base height in inches (default 3)
#' @param height_per_pathway Height per pathway in inches (default 0.3)
#' @param max_height Maximum height in inches (default 12)
#' @return Calculated height in inches
calculate_pathway_plot_height <- function(n_pathways, base_height = 3,
                                        height_per_pathway = 0.3, max_height = 12) {
  min(max_height, base_height + (n_pathways * height_per_pathway))
}

#' Create GSEA enrichment plot
#' @param gsea_results GSEA results data frame
#' @param comparison_name Name of comparison
#' @param max_pathways Maximum pathways to show (default 20)
#' @param interactive Create interactive version (default FALSE)
#' @return ggplot object or NULL if no significant pathways
create_gsea_plot <- function(gsea_results, comparison_name, max_pathways = 20,
                           interactive = FALSE) {

  if (is.null(gsea_results) || nrow(gsea_results) == 0) {
    return(NULL)
  }

  # Filter and prepare data
  plot_data <- gsea_results %>%
    filter(padj <= 0.05) %>%
    arrange(desc(abs(NES))) %>%
    head(max_pathways) %>%
    mutate(
      pathway_clean = str_replace(pathway, "\\(GO", " \\(GO"),
      pathway_clean = str_trunc(pathway_clean, 60),
      pathway_clean = fct_reorder(pathway_clean, NES),
      significant = padj <= 0.05
    )

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  # Create plot
  p <- ggplot(plot_data, aes(x = NES, y = pathway_clean)) +
    geom_point(aes(color = -log10(padj), size = lengths(leadingEdge)),
               alpha = 0.8) +
    scale_color_gradient(low = "blue", high = "red",
                        name = "-log₁₀(padj)") +
    scale_size_continuous(range = c(2, 8), name = "Gene Count") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.3) +
    labs(
      title = paste("Gene Set Enrichment:", comparison_name),
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway",
      subtitle = paste(nrow(plot_data), "significant pathways shown")
    ) +
    theme_report() +
    theme(
      axis.text.y = element_text(size = rel(0.8)),
      legend.position = "right"
    )

  if (interactive) {
    p <- ggplotly(p, tooltip = c("x", "y", "colour", "size")) %>%
      layout(
        title = list(text = paste("Gene Set Enrichment:", comparison_name)),
        showlegend = TRUE
      )
  }

  return(p)
}

# ============================================================================= #
# UTILITY FUNCTIONS ----
# ============================================================================= #

#' Format numbers for display in reports
#' @param x Numeric vector
#' @param digits Number of significant digits (default 3)
#' @return Formatted character vector
format_numbers <- function(x, digits = 3) {
  ifelse(is.na(x), "—", signif(x, digits))
}

#' Create summary statistics table
#' @param res_all All results
#' @param res_sig Significant results
#' @param comparison_name Comparison name
#' @return Data frame with summary statistics
create_summary_stats <- function(res_all, res_sig, comparison_name) {
  tibble(
    Comparison = comparison_name,
    `Genes Tested` = nrow(res_all),
    `Significant Genes` = nrow(res_sig),
    `Upregulated` = sum(res_sig$log2FoldChange > 0, na.rm = TRUE),
    `Downregulated` = sum(res_sig$log2FoldChange < 0, na.rm = TRUE),
    `% Significant` = round(100 * nrow(res_sig) / nrow(res_all), 1)
  )
}