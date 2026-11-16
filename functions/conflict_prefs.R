# Function-name conflict preferences
#
# Bioconductor packages (DESeq2, tximeta -> S4Vectors, IRanges, GenomicRanges,
# matrixStats, Biobase, ...) export functions with the same names as tidyverse
# verbs. Whichever package is attached LAST wins, so e.g. `count()` silently
# became matrixStats::count and errored on tibbles. Load order cannot fix this
# reliably (re-attaching a package does not move it up the search path).
#
# The `conflicted` package resolves each name explicitly, independent of load
# order. Preferences can be declared before a package is attached.
# Any NEW conflict not listed here becomes an informative error
# ("[conflicted] x found in 2 packages") instead of a silent mask - add a line
# below when that happens.
#
# Source this AFTER library() calls, e.g.
#   source(here::here("functions/conflict_prefs.R"))

suppressPackageStartupMessages(library(conflicted))

# matrixStats row/col functions: prefer the MatrixGenerics S4 generics that
# DESeq2/SummarizedExperiment use (they dispatch to matrixStats for matrices).
# Applied first so the explicit preferences below (e.g. rowRanges) override it.
conflicted::conflict_prefer_all("MatrixGenerics", c("matrixStats", "Biobase"), quiet = TRUE)

conflicted::conflicts_prefer(
  # tidyverse verbs over Bioconductor / plotly / stats versions
  dplyr::filter, dplyr::select, dplyr::count, dplyr::rename, dplyr::slice,
  dplyr::first, dplyr::desc, dplyr::combine, dplyr::collapse, dplyr::lag,
  purrr::reduce, purrr::set_names, purrr::discard,
  tidyr::expand, tidyr::extract,
  readr::col_factor,
  lubridate::second, lubridate::`second<-`, lubridate::`%within%`,
  # base versions of generics re-exported by Bioconductor/generics
  base::union, base::intersect, base::setdiff, base::setequal,
  base::as.factor, base::as.difftime, base::Position,
  # other non-tidyverse clashes: keep the package that won before conflicted
  magrittr::subtract, glue::trim, Biobase::anyMissing,
  SummarizedExperiment::rowRanges,
  plotly::layout, ComplexHeatmap::add_heatmap, heatmaply::normalize,
  .quiet = TRUE
)
