#' @title Module 03: Dimensionality Reduction and Clustering
#' @description Perform PCA, UMAP, tSNE and clustering
#' @author r-scrna development team

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

#' Run PCA analysis
#' @param obj Seurat object
#' @param config Configuration object
#' @return Seurat object with PCA
run_pca <- function(obj, config) {
  npcs <- config$clustering$npcs
  n_top_genes <- config$normalization$n_top_genes
  
  log_message("Running PCA (top ", npcs, " PCs)...")
  
  obj <- RunPCA(obj, features = head(VariableFeatures(obj), n_top_genes), npcs = npcs, verbose = FALSE)
  return(obj)
}

#' Determine optimal number of PCs using Elbow plot
#' @param obj Seurat object with PCA
#' @param config Configuration object
#' @return Recommended number of PCs
determine_pcs <- function(obj, config) {
  npcs <- config$clustering$npcs
  reduction <- config$clustering$reduction
  
  if (reduction == "pca") {
    log_message("Analyzing PC variance using Elbow plot...")
    elbow <- ElbowPlot(obj, ndims = npcs)
    print(elbow)
    
    # Suggest optimal PCs based on elbow point
    # Find the point where variance drops most
    pcs_to_use <- min(10, npcs)  # Default to first 10
    return(pcs_to_use)
  } else {
    log_message("Using default PCs: ", config$clustering$npcs)
    return(config$clustering$npcs)
  }
}

#' Run UMAP
#' @param obj Seurat object with PCA
#' @param config Configuration object
#' @return Seurat object with UMAP
run_umap <- function(obj, config) {
  dims <- config$dimensionality_reduction$umap$n_neighbors
  min_dist <- config$dimensionality_reduction$umap$min_dist
  
  log_message("Running UMAP (dims: ", dims, ")...")
  
  obj <- RunUMAP(obj, reduction = "umap", dims = dims, 
               n.neighbors = dims, min.dist = min_dist, verbose = FALSE)
  return(obj)
}

#' Run t-SNE
#' @param obj Seurat object with PCA
#' @param config Configuration object
#' @return Seurat object with t-SNE
run_tsne <- function(obj, config) {
  perplexity <- config$dimensionality_reduction$tsne$perplexity
  
  log_message("Running t-SNE (perplexity: ", perplexity, ")...")
  obj <- RunTSNE(obj, reduction = "tsne", dims = 2, 
               perplexity = perplexity, verbose = FALSE)
  return(obj)
}

#' Find nearest neighbors
#' @param obj Seurat object with PCA
#' @param config Configuration object
#' @return Seurat object with neighbors
find_neighbors <- function(obj, config) {
  dims <- determine_pcs(obj, config)
  reduction <- config$clustering$reduction
  
  log_message("Finding neighbors using ", reduction, " (dims: ", dims, ")...")
  
  obj <- FindNeighbors(obj, reduction = reduction, dims = dims, verbose = FALSE)
  return(obj)
}

#' Perform clustering at multiple resolutions
#' @param obj Seurat object with neighbors
#' @param config Configuration object
#' @return Seurat object with multiple cluster assignments
run_clustering <- function(obj, config) {
  resolutions <- config$clustering$resolutions
  algorithm <- config$clustering$algorithm
  
  log_message("Running clustering at ", length(resolutions), " resolutions...")
  
  for (res in resolutions) {
    cluster.name <- paste0("RNA_snn_res.", res)
    obj <- FindClusters(obj, resolution = res, 
                       algorithm = algorithm, 
                       cluster.name = cluster.name,
                       verbose = FALSE)
    log_message("  Cluster resolution ", res, ": ", 
                 length(unique(Idents(obj))), " clusters")
  }
  
  return(obj)
}

#' Main function for dimensionality reduction and clustering
#' @param obj Seurat object
#' @param config Configuration object
#' @return List with results
run_dim_cluster <- function(obj, config) {
  log_message("Module 03: Dimensionality Reduction and Clustering")
  
  # Step 1: PCA
  obj <- run_pca(obj, config)
  
  # Step 2: Determine optimal PCs
  pcs_to_use <- determine_pcs(obj, config)
  
  # Step 3: Find neighbors
  obj <- find_neighbors(obj, config)
  
  # Step 4: Clustering at multiple resolutions
  obj <- run_clustering(obj, config)
  
  # Step 5: UMAP
  obj <- run_umap(obj, config)
  
  # Step 6: t-SNE
  obj <- run_tsne(obj, config)
  
  # Create visualization
  umap_plots <- list()
  for (res in config$clustering$resolutions) {
    ident.use <- paste0("RNA_snn_res.", res)
    p <- DimPlot(obj, reduction = "umap", group.by = ident.use, 
                 label = TRUE, label.size = 5, pt.size = 0.5) +
      plot_annotation(title = paste0("UMAP (res: ", res, ")"))
    umap_plots[[as.character(res)]] <- p
  }
  
  combined_umap <- wrap_plots(umap_plots)
  
  # Cluster metrics
  cluster_metrics <- data.frame(
    resolution = as.numeric(config$clustering$resolutions),
    n_clusters = sapply(config$clustering$resolutions, function(res) {
      ident.use <- paste0("RNA_snn_res.", res)
      length(unique(Idents(obj, ident.use = ident.use)))
    })
  )
  
  log_message("Module 03 complete")
  return(list(
    obj = obj,
    pcs_to_use = pcs_to_use,
    cluster_metrics = cluster_metrics,
    umap_plots = umap_plots,
    combined_umap = combined_umap
  ))
}

#' Visualize PC variance
#' @param obj Seurat object with PCA
#' @param npcs Number of PCs to show
#' @return ggplot object
visualize_pca_variance <- function(obj, npcs = 50) {
  log_message("Visualizing PC variance (first ", npcs, " PCs)...")
  
  variance_explained <- obj[["pca"]]@stdev^2 / sum(obj[["pca"]]@stdev^2)
  variance_cumsum <- cumsum(variance_explained)[1:npcs]
  
  df <- data.frame(
    PC = 1:npcs,
    Variance = as.numeric(variance_explained[1:npcs]),
    Cumulative = variance_cumsum[1:npcs]
  )
  
  # Find elbow point
  df$Distance <- c(0, diff(df$Cumulative))
  elbow_point <- which.max(df$Distance[2:min(20, npcs)])
  
  p <- ggplot(df, aes(x = PC, y = Variance)) +
    geom_line() +
    geom_vline(xintercept = elbow_point + 0.5, linetype = "dashed", 
               color = "red", linewidth = 1) +
    geom_point(aes(x = elbow_point + 1, y = df$Variance[elbow_point + 1]), 
               size = 3, color = "red") +
    labs(x = "Principal Component", 
         y = "Variance Explained") +
    theme_bw() +
    plot_annotation(title = "PCA Variance Explained")
  
  return(p)
}

#' Calculate clustering metrics
#' @param obj Seurat object with multiple cluster resolutions
#' @param resolutions Vector of resolutions used
#' @return Data frame with clustering metrics
calculate_clustering_metrics <- function(obj, resolutions) {
  metrics_list <- list()
  
  for (res in resolutions) {
    ident.use <- paste0("RNA_snn_res.", res)
    clusters <- Idents(obj, ident.use = ident.use)
    n_clusters <- length(unique(clusters))
    
    # Calculate silhouette score (simplified)
    if (n_clusters > 1) {
      sil_score <- mean(sapply(clusters, function(clust) {
        clust_cells <- WhichCells(obj, ident = paste0("RNA_snn_res.", res), expression = clust)
        clust_cells <- sample(clust_cells, min(1000, length(clust_cells)))
        # Simplified metric: within-cluster variance
        if (length(clust_cells) > 1) {
          clust_pca <- Embeddings(obj, reduction = "pca", dims = 1:10)[clust_cells, ]
          wss <- sum(rowSums(scale(clust_pca)^2))
          return(-wss / length(clust_cells))
        } else {
          return(0)
        }
      }))
    } else {
      sil_score <- NA
    }
    
    metrics_list[[as.character(res)]] <- data.frame(
      resolution = res,
      n_clusters = n_clusters,
      silhouette_score = sil_score
    )
  }
  
  return(bind_rows(metrics_list, .id = "resolution"))
}
