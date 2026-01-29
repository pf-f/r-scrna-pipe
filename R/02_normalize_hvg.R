#' @title Module 02: Normalization and Highly Variable Genes
#' @description Normalize data and detect HVGs
#' @author r-scrna development team

library(Seurat)
library(dplyr)
library(ggplot2)

#' Normalize data using specified method
#' @param obj Seurat object
#' @param config Configuration object
#' @return Normalized Seurat object
normalize_data <- function(obj, config) {
  method <- config$normalization$method
  log_factor <- config$normalization$scale_factor
  
  log_message("Normalizing data using method: ", method)
  
  if (method == "LogNormalize") {
    obj <- NormalizeData(obj, 
                          normalization.method = "LogNormalize", 
                          scale.factor = log_factor)
  } else if (method == "SCTransform") {
    if (!requireNamespace("SCTransform", quietly = TRUE)) {
      stop("SCTransform package not installed. Install with: install.packages('SCTransform')")
    }
    obj <- SCTransform(obj, vars.to.regress = config$normalization$vars_to_regress)
    log_message("Applied SCTransform normalization")
  } else {
    stop("Unsupported normalization method: ", method)
  }
  
  log_message("Normalization complete")
  return(obj)
}

#' Find highly variable genes
#' @param obj Normalized Seurat object
#' @param config Configuration object
#' @return Seurat object with HVGs
find_hvg <- function(obj, config) {
  n_top_genes <- config$normalization$n_top_genes
  
  log_message("Finding highly variable genes (top ", n_top_genes, ")...")
  
  if (config$normalization$method == "SCTransform") {
    # SCTransform already identifies HVGs
    hvg <- VariableFeatures(obj)
    log_message("Found ", length(hvg), " HVGs via SCTransform")
  } else {
    obj <- FindVariableFeatures(obj, 
                                 selection.method = "vst", 
                                 nfeatures = n_top_genes)
    hvg <- VariableFeatures(obj)
    log_message("Found ", length(hvg), " HVGs via vst method")
  }
  
  # Plot HVGs
  p1 <- VariableFeaturePlot(obj, ncol = 2)
  p2 <- LabelPoints(plot = p1, points = head(hvg, 10), repel = TRUE)
  hvg_plot <- p1 + p2 + plot_annotation(title = paste0("Top ", n_top_genes, " Variable Genes"))
  
  return(list(obj = obj, hvg = hvg, plot = hvg_plot))
}

#' Scale data
#' @param obj Seurat object with HVGs
#' @param config Configuration object
#' @return Scaled Seurat object
scale_data <- function(obj, config) {
  vars_to_regress <- config$normalization$vars_to_regress
  
  log_message("Scaling data...")
  
  if (length(vars_to_regress) > 0) {
    obj <- ScaleData(obj, 
                     features = rownames(obj), 
                     vars.to.regress = vars_to_regress)
    log_message("Regressing out: ", paste(vars_to_regress, collapse = ", "))
  } else {
    obj <- ScaleData(obj, features = rownames(obj))
  }
  
  log_message("Scaling complete")
  return(obj)
}

#' Main function for normalization and HVG detection
#' @param obj Seurat object
#' @param config Configuration object
#' @return List with normalized object and plots
run_normalization <- function(obj, config) {
  log_message("Module 02: Normalization and HVG Detection")
  
  # Step 1: Normalize
  norm_obj <- normalize_data(obj, config)
  
  # Step 2: Find HVGs
  hvg_result <- find_hvg(norm_obj, config)
  obj <- hvg_result$obj
  
  # Step 3: Scale
  scaled_obj <- scale_data(obj, config)
  
  # Create plots
  plots <- list()
  plots[["hvg"]] <- hvg_result$plot
  
  log_message("Module 02 complete")
  return(list(obj = scaled_obj, hvg = hvg_result$hvg, plots = plots))
}

#' Calculate variance statistics
#' @param obj Seurat object
#' @return Data frame with variance stats
calculate_variance_stats <- function(obj) {
  vst_results <- hvg_result <- VariableFeatures(obj)
  
  stats <- data.frame(
    mean_expr = rowMeans(as.matrix(GetAssayData(obj, slot = "data"))),
    var_expr = apply(as.matrix(GetAssayData(obj, slot = " "data")), 1, var),
    std_expr = apply(as.matrix(GetAssayData(obj, slot = "data")), 1, sd)
  )
  stats$gene <- rownames(stats)
  
  return(stats)
}
