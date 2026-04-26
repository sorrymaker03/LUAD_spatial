# LUAD_spatial
Benchmarking four cell-type annotation methods (Seurat label transfer, CARD, RCTD, and SPOTlight) on 45 LUAD spatial transcriptomics. 

---

## Overview

This repository contains the complete R code for benchmarking four cell-type annotation methods on 45 lung adenocarcinoma (LUAD) spatial transcriptomics (10x Visium) samples. The four methods evaluated are Seurat-based label transfer (RPCA), CARD, RCTD (spacexr), and SPOTlight.

The analysis includes cell-type composition comparison, correlation and RMSE assessment, spatial autocorrelation (Moran's I and LISA), neighborhood enrichment analysis, module score evaluation, and gene expression gradient analysis.

---

## Data sources

| Data type | Source | Samples |
|-----------|--------|---------|
| LUAD spatial data | Human NSCLC lesions (10x Visium) | 20 |
| LUAD spatial data | GSE307534 | 25 |
| scRNA-seq reference | Human Cell Atlas Lung Network | 114,549 cells, 10 cell types |

---

## Requirements

**R version:** 4.3

**Required packages:**

```r
library(circlize)
library(CARD)
library(ComplexHeatmap)
library(dplyr)
library(FNN)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(patchwork)
library(reshape2)
library(rstatix)
library(Seurat)
library(SingleCellExperiment)
library(spacexr)
library(SpatialExperiment)
library(spdep)
library(SPOTlight)
library(SummarizedExperiment)
library(tidyr)
library(viridis)
