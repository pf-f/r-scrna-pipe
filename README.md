# r-scrna

A comprehensive automated single-cell RNA-seq analysis pipeline in R, providing modular functionality for quality control, normalization, clustering, annotation, differential expression, enrichment analysis, and optional advanced analyses including pseudotime, cell communication, SCENIC, RNA velocity, and CNV detection.

## Features

- **Core Pipeline**: Complete analysis workflow from raw count matrix to annotated clusters
- **Modular Design**: Optional add-on modules for advanced analyses
- **Flexible Configuration**: YAML-based parameter management
- **Multi-format Support**: Accepts CellRanger, raw count matrices, and Seurat RDS files
- **Species Support**: Human and mouse (with extensibility for other species)
- **Batch Integration**: Harmony and CCA-based integration methods
- **Auto-annotation**: Multiple annotation methods (marker database, SingleR, manual)
- **Comprehensive Reports**: RMarkdown-generated HTML/PDF reports
- **Reproducible**: Version-controlled configuration and intermediate file saving

## Pipeline Overview

```
00.check_input()      → Validate input files
01.qc_filter(obj)      → QC filtering (returns qc.rds)
02.normalize_hvg(qc.rds)  → Normalization & HVG (returns norm.rds)
03.dim_cluster(norm.rds)  → Dimensionality reduction & clustering (returns cluster.rds)
04.auto_annotate(cluster.rds) → Cell annotation (returns label.rds)
05.diffexp(label.rds)     → Differential expression (DEG.csv)
06.enrichment(DEG.csv)   → GO/KEGG enrichment
07.pseudotime(label.rds)  → Optional: Pseudotime analysis (Monocle3)
08.cellchat(label.rds)    → Optional: Cell communication
09.scenic(label.rds)      → Optional: SCENIC analysis
10.rna_velocity(label.rds) → Optional: RNA velocity
11.cnv(label.rds)        → Optional: CNV analysis
09.report()              → Generate final report
```

## Installation

### Option 1: Using Conda (Recommended)

```bash
# Create environment from YAML file
conda env create -f environment.yml

# Activate environment
conda activate r-scrna

# Install additional R packages if needed
R
> install.packages(c(
    "Seurat", "SingleR", "CellChat", 
    "monocle3", "SCENIC", "clusterProfiler",
    "DESeq2", "patchwork", "ComplexHeatmap"
  ))
> quit()
```

### Option 2: Using R Only

```r
# Install required packages
install.packages(c("remotes"))
remotes::install_github("satijalab/seurat-object")
remotes::install_github("satijalab/seurat")
remotes::install_github("sqjin/CellChat")

# Install from CRAN/Bioconductor
install.packages(c("Seurat", "dplyr", "ggplot2", "yaml"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "SingleR", "clusterProfiler"))
```

## Quick Start

### 1. Prepare Configuration

Edit `config/config.yaml` to set your analysis parameters:

```yaml
# Set species
species:
  name: "human"
  genome: "hg38"
  mt_pattern: "^MT-"

# Set QC thresholds
qc:
  min_features_per_cell: 200
  max_mt_pct: 20.0

# Enable optional modules
modules:
  pseudotime:
    enabled: false
  cellchat:
    enabled: false
```

### 2. Run the Pipeline

```r
# From R or Rscript
Rscript scripts/main_scRNA.R --config config/config.yaml

# Or run interactively in R
source("scripts/main_scRNA.R")
```

### 3. Using as R Package

```r
# Install package
devtools::install("..")

# Use in R
library(r.scrna)
run_pipeline(config_path = "config/config.yaml")
```

## Project Structure

```
r-scrna/
├── config/
│   ├── config.yaml               # Main configuration file
│   └── marker_db/               # Cell marker databases
├── R/                            # Modular R functions
│   ├── 00_input_checks.R
│   ├── 01_qc_filter.R
│   ├── 02_normalize_hvg.R
│   ├── 03_dim_cluster.R
│   ├── 04_auto_annotate.R
│   ├── 05_diffexp.R
│   ├── 06_enrichment.R
│   ├── 07_pseudotime.R          # Optional
│   ├── 08_cellchat.R            # Optional
│   ├── 09_scenic.R              # Optional
│   ├── 10_rna_velocity.R         # Optional
│   ├── 11_cnv.R               # Optional
│   └── utils.R
├── scripts/
│   └── main_scRNA.R            # Main entry point
├── Rmd/
│   ├── 01_qc_report.Rmd
│   ├── 02_clustering_report.Rmd
│   ├── 03_annotation_report.Rmd
│   ├── 04_de_report.Rmd
│   ├── 05_enrichment_report.Rmd
│   └── 06_full_report.Rmd
├── data/
├── outputs/
│   ├── intermediate/             # Intermediate RDS files
│   └── reports/                 # HTML/PDF reports
├── examples/
│   └── example_config.yaml
├── environment.yml               # Conda environment
├── requirements.txt              # R package list
├── README.md
└── LICENSE
```

## Configuration

Key parameters in `config/config.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `species.name` | "human" | Species type (human/mouse) |
| `species.genome` | "hg38" | Reference genome build |
| `qc.min_features_per_cell` | 200 | Minimum features per cell |
| `qc.max_mt_pct` | 20.0 | Max mitochondrial percentage |
| `normalization.method` | "LogNormalize" | Normalization method |
| `normalization.n_top_genes` | 3000 | Number of HVGs |
| `clustering.resolutions` | [0.2,0.4,0.5,0.6,0.8,1.0] | Clustering resolutions |
| `clustering.npcs` | 30 | Number of PCs for clustering |
| `de.test_use` | "wilcox" | Differential expression test |
| `annotation.method` | "auto" | Annotation method |

## Modules

### Core Modules (Always Run)

1. **Input Checks** - Validate input data format and structure
2. **QC & Filtering** - Quality control and cell/gene filtering
3. **Normalization & HVG** - Data normalization and highly variable gene detection
4. **Dimensionality Reduction** - PCA, UMAP, t-SNE
5. **Clustering** - Graph-based clustering at multiple resolutions
6. **Annotation** - Automatic or manual cell type annotation
7. **Differential Expression** - Find marker genes and DE genes
8. **Enrichment Analysis** - GO, KEGG pathway enrichment

### Optional Modules (Enable in config)

9. **Pseudotime Analysis** - Monocle3, Slingshot trajectory inference
10. **CellChat** - Cell-cell communication analysis
11. **SCENIC** - Gene regulatory network inference
12. **RNA Velocity** - Splicing dynamics and RNA velocity
13. **CNV Analysis** - Copy number variation detection

## Output Files

All outputs are saved to `outputs/` directory by default:

- `intermediate/` - RDS objects from each analysis step
- `reports/` - HTML/PDF reports with visualizations
- `DEG.csv` - Differential expression results
- `markers.csv` - Cluster marker genes
- `enrichment/` - GO/KEGG enrichment results

## License

This project is licensed under the MIT License - see the LICENSE file for details.
