#!/usr/bin/env Rscript

# RNAseq DE compilation generation functions
# Creation date: 2025-09-04

library(tidyverse)
library(magrittr)
library(DESeq2)
library(tximeta)
library(glue)
library(yaml)
library(here)

here::i_am("functions/processing_helpers.R")

