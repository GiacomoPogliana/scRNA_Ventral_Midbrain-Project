# scRNA_Ventral_Midbrain-Project

## Author: 
Giacomo Pogliana

## Overview
This R-based pipeline is designed to process single-cell RNA sequencing (scRNA-seq) data from the ventral midbrain of postnatal day 21 (P21) mice. It automates quality control, cell filtering, doublet removal, data normalization, dimensionality reduction, cell-type clustering, canonical marker annotation, and differential expression analysis. The pipeline specifically targets midbrain neuronal lineages (including dopaminergic and GABAergic neurons) and surrounding glial populations (astrocytes, microglia, and oligodendrocytes) at the weaning developmental stage.

## Prerequisites & Dependencies
To run this pipeline, ensure R (>= 4.2) and the following bioinformatics packages are installed and accessible in your R environment:
- **Seurat** (v4 / v5 - Data structure, normalization, dimensional reduction, and clustering)

- **tidyverse** (Data manipulation and visualization via `dplyr` and `ggplot2`)

- **SingleCellExperiment** (Data container infrastructure for single-cell assays)

- **DoubletFinder** (Computational doublet identification and removal)

- **patchwork** (Multi-panel plot layout engine)

- **ComplexHeatmap** (Advanced heatmap generation for gene expression panels)

- **rmarkdown & knitr** (Automated HTML report compilation)


## Directory Structure and Input Files
The script expects a specific directory layout to function correctly. You must execute the pipeline from the main repository root containing raw 10x Genomics matrices or HDF5 files.

**Required Directory Layout:**
```text
.
├── data/
│   ├── raw/                  # Raw 10x Genomics matrices (barcodes, features, matrix)
│   └── processed/            # Serialized Seurat objects (.rds)
├── scripts/                  # R pipeline scripts
│   ├── 01_qc_filtering.R
│   ├── 02_normalization_clustering.R
│   ├── 03_cell_type_annotation.R
│   └── 04_differential_expression.R
├── docs/                     # Compiled HTML analysis reports
└── figures/                  # High-resolution exported plots
```

## Usage
The analysis pipeline can be executed sequentially via `Rscript` command-line arguments or sourced directly within RStudio. Quality control thresholds and clustering resolution parameters are configured via command-line flags.

**Available Flag Options:**
- `--min_genes` : Minimum detected genes per cell (default: `500`)

- `--max_genes` : Maximum detected genes per cell to exclude doublets (default: `6000`)

- `--max_mito`  : Maximum allowed percentage of mitochondrial reads (default: `5.0`)

- `--resolution`: Louvain/Leiden graph-based clustering resolution (default: `0.5`)

## Pipeline Workflow
For each sample dataset, the pipeline executes the following processing steps:

1. **Data Loading & Preprocessing:** Reads raw 10x Genomics count matrices and constructs a unified `Seurat` object with associated experiment metadata.

2. **Quality Control & Filtering:** Calculates mitochondrial gene percentages (`percent.mt`) and ribosomal ratios; filters out damaged cells, empty droplets, and high-mitochondrial outliers.

3. **Doublet Detection:** Executes `DoubletFinder` to predict and remove artificial transcriptomic doublets.

4. **Normalization & Variance Stabilization:** Applies `SCTransform` (or log-normalization) while regressing out confounding technical variation (`percent.mt`, `nCount_RNA`) and identifies highly variable genes (HVGs).

5. **Dimensionality Reduction & Clustering:** Performs Principal Component Analysis (PCA), calculates UMAP coordinates from significant PCs, and executes graph-based clustering.

6. **Cell-Type Annotation:** Assigns cell identities using canonical marker genes for ventral midbrain populations:
   * **Dopaminergic Neurons (DA):** *Th*, *Ddc*, *Slc6a3*, *Slc18a2*
   * **GABAergic / Glutamatergic Neurons:** *Gad1*, *Gad2*, *Slc17a6*
   * **Astrocytes:** *Aqp4*, *Gfap*, *S100b*
   * **Oligodendrocytes & OPCs:** *Mbp*, *Olig1*, *Olig2*, *Pdgfra*
   * **Microglia:** *Cx3cr1*, *P2ry12*, *Aif1*

7. **Differential Expression Analysis:** Conducts Wilcoxon rank-sum testing (`FindAllMarkers`) to extract cluster-specific gene signatures and differential expression profiles across subpopulations.

## Outputs
Within the `data/processed/`, `figures/`, and `docs/` directories, the pipeline generates:
- Quality control violin plots before and after cell filtering (`qc_post_filtering.png`).

- Fully annotated and serialized Seurat RDS objects (`midbrain_p21_annotated.rds`).

- High-resolution UMAP embeddings grouped by cluster and cell-type annotations (`umap_cell_types.png`).

- Marker gene expression dot plots and heatmaps (`marker_heatmap.png`).

- CSV tables listing cluster-specific differentially expressed genes with log2 fold-changes and adjusted p-values (`deg_cluster_markers.csv`).

- A compiled, fully reproducible HTML analysis report (`scRNA_Ventral_Midbrain_Report.html`).
