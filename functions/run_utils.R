#!/usr/bin/env Rscript

# RNAseq pipeline run script helpers
# Creation date: 2025-09-04

library(tidyverse)
library(magrittr)
library(DESeq2)
library(tximeta)
library(glue)
library(yaml)
library(cli)
library(here)

here::i_am("functions/run_utils.R")

