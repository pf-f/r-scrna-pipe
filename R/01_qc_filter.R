#' @title Module 01: Quality Control and Filtering
#' @description Perform QC filtering on Seurat object
#' @author r-scrna development team

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

#' Main QC filtering function
#' @param obj Seurat object
#' @param config Configuration object
#' @return Filtered Seurat object
run_qc <- function(obj, config) {
  log_message("Running QC filtering...")
  
  # QC parameters
  min_cells <- config$qc$min_cells_per_gene
  min_features <- config$qc$min_features_per_cell
  max_features <- config$qc$max_features_per_cell
  max_mt_pct <- config$qc$max_mt_pct
  max_ribo_pct <- config$qc$max_ribo_pct
  
  # Calculate QC metrics
  species_settings <- get_species_settings(config)
  mt_pattern <- species_settings$mt_pattern
  
  obj[["percent.mt"] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  if (!is.null(config$species$ribosomal_pattern)) {
    ribo_pattern <- config$species$ribosomal_pattern
    obj[["percent.ribo"] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
  }
  
  # Doublet detection (optional)
  if (config$qc$doublet_detection) {
    log_message("Doublet detection not implemented in V1")
  }
  
  # Apply filters
  log_message("Applying QC filters...")
  
  # Filter cells
  obj_filtered <- subset(obj, subset = nFeature_RNA > min_features & 
                                        nFeature_RNA < max_features & 
                                        percent.mt < max_mt_pct)
  
  if (!is.null(config$qc$max_ribo_pct)) {
    obj_filtered <- subset(obj_filtered, subset = percent.ribo < max_ribo_pct)
  }
  
  # Filter genes
  obj_filtered <- subset(obj_filtered, features = rowSums(x = GetAssayData(obj, slot = "counts")) > min_cells)
  
  # Log filtering statistics
  cells_before <- ncol(obj)
  cells_after <- ncol(obj_filtered)
  genes_before <- nrow(obj)
  genes_after <- nrow(obj_filtered)
  
  log_message("Cells before: ", cells_before, ", after: ", cells_after, 
             " (removed: ", cells_before - cells_after, ")")
  log_message("Genes before: ", genes_before, ", after: ", genes_after,
             " (removed: ", genes_before - genes_after, ")")
  
  return(obj_filtered)
}

#' Calculate detailed QC metrics
#' @param obj Seurat object
#' @param config Configuration object
#' @return Data frame with QC metrics
calculate_qc_metrics <- function(obj, config) {
  log_message("Calculating QC metrics...")
  
  qc_df <- obj@meta.data
  
  # Calculate additional metrics
  qc_df$counts_per_cell <- Matrix::colSums(GetAssayData(obj, slot = "counts"))
  qc_df$features_per_cell <- Matrix::colSums(GetAssayData(obj, slot = "counts") > 0)
  
  # Create violin plots
  p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.1)
  
  # Create scatter plots
  p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    ggtitle("MT % vs UMI count")
  p3 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    ggtitle("Genes vs UMI count")
  
  combined_plot <- (p1 | p2 / p3)
  
  return(list(qc_df = qc_df, plot = combined_plot))
}

#' Filter cells based on QC thresholds
#' @param obj Seurat object
#' @param config Configuration object
#' @return Filtered Seurat object
filter_cells <- function(obj, config) {
  min_features <- config$qc$min_features_per_cell
  max_features <- config$qc$max_features_per_cell
  max_mt_pct <- config$qc$max_mt_pct
  max_ribo_pct <- config$qc$max_ribo_pct
  
  subset_string <- "nFeature_RNA > min_features"
  
  if (!is.null(max_features)) {
    subset_string <- paste0(subset_string, " & nFeature_RNA < max_features")
  }
  
  subset_string <- paste0(subset_string, " & percent.mt < ", max_mt_pct)
  
  if (!is.null(max_ribo_pct) && "percent.ribo" %in% colnames(obj@meta.data)) {
    subset_string <- paste0(subset_string, " & percent.ribo < ", max_ribo_pct)
  }
  
  obj_filtered <- subset(obj, subset = eval(parse(text = subset_string)))
  
  log_message("Applied subset: ", subset_string)
  
  return(obj_filtered)
}

#' Generate QC report plots
#' @param obj Seurat object
#' @param config Configuration object
#' @return List of ggplot objects
generate_qc_plots <- function(obj, config) {
  log_message("Generating QC plots...")
  
  # Violin plots
  p_violin <- VlnPlot(obj, 
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     ncol = 3,
                     pt.size = 0.1) +
    plot_annotation(title = "QC Metrics")
  
  # Scatter plots
  p_scatter1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    geom_smooth(method = "lm") +
    ggtitle("MT % vs UMI count")
  
  p_scatter2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm") +
    ggtitle("Genes vs UMI count")
  
  # Combine plots
  qc_plots <- list(
    violin = p_violin,
    scatter_mt = p_scatter1,
    scatter_features = p_scatter2
  )
  
  return(qc_plots)
}
