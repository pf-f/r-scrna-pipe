#' @title Module 05: Differential Expression Analysis
#' @description Perform differential expression analysis between groups
#' @author r-scrna development team

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

#' Find differential expressed genes between two groups
#' @param obj Annotated Seurat object
#' @param ident.1 First identity
#' @param ident.2 Second identity (optional, compares to all others if NULL)
#' @param config Configuration object
#' @return Data frame with DEG results
run_de_test <- function(obj, ident.1, ident.2 = NULL, config) {
  min_pct <- config$de$min_pct
  logfc_threshold <- config$de$logfc_threshold
  test_use <- config$de$test_use
  only_pos <- config$de$only_pos
  
  log_message("Running DE analysis: ident.1 = ", ident.1, 
             "ident.2 = ", ifelse(is.null(ident.2), "all other clusters", ident.2))
  
  # Find markers
  de_result <- FindMarkers(
    obj, 
    ident.1 = ident.1, 
    ident.2 = ident.2,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    only.pos = only_pos,
    test.use = test_use
  )
  
  log_message("Found ", nrow(de_result), " differentially expressed genes")
  return(de_result)
}

#' Compare multiple groups
#' @param obj Annotated Seurat object
#' @param groups_to_compare List of groups to compare
#' @param config Configuration object
#' @return List of DEG data frames
compare_groups <- function(obj, groups_to_compare, config) {
  de_list <- list()
  
  for (i in seq_along(groups_to_compare)) {
    for (j in seq_along(groups_to_compare)) {
      if (i >= j) next  # Skip duplicates and self-comparisons
      
      group1 <- groups_to_compare[[i]]
      group2 <- groups_to_compare[[j]]
      
      de_result <- run_de_test(obj, group1, group2, config)
      de_result$group1 <- group1
      de_result$group2 <- group2
      de_result$comparison <- paste0(group1, "_vs_", group2)
      
      de_list[[length(de_list) + 1]] <- de_result
    }
  }
  
  log_message("Completed ", length(de_list), " group comparisons")
  return(de_list)
}

#' Compare all groups vs rest
#' @param obj Annotated Seurat object
#' @param config Configuration object
#' @return List of DEG data frames
compare_all_vs_rest <- function(obj, config) {
  de_result <- FindAllMarkers(
    obj, 
    only.pos = config$de$only_pos,
    min.pct = config$de$min_pct,
    logfc.threshold = config$de$logfc_threshold,
    test.use = config$$de$test_use
  )
  
  de_result$comparison <- "all_vs_rest"
  log_message("Found ", nrow(de_result), " DEGs (all clusters vs rest)")
  return(de_result)
}

#' Run DE analysis for specified groups
#' @param obj Annotated Seurat object
#' @param groups_to_compare List of groups to compare
#' @param config Configuration object
#' @return List of DEG data frames
run_de_analysis <- function(obj, groups_to_compare = NULL, config) {
  log_message("Module 05: Differential Expression Analysis")
  
  if (is.null(groups_to_compare) || length(groups_to_compare) == 0) {
    log_message("No groups specified, running all vs rest")
    de_results <- compare_all_vs_rest(obj, config)
    de_results <- list(de_results)
  } else {
    log_message("Comparing specific groups: ", paste(groups_to_compare, collapse = ", "))
    de_results <- compare_groups(obj, groups_to_compare, config)
  }
  
  # Add statistics
  de_results <- lapply(de_results, function(de_df) {
    de_df <- de_df %>%
      dplyr::mutate(
        neg_log10_padj = if("p_val_adj" %in% colnames(de_df), -log10(de_df$p_val_adj)),
        regulated = dplyr::case_when(
          avg_log2FC > 0 ~ "Up",
          avg_log2FC < 0 ~ "Down",
          TRUE ~ "NS"
        )
      )
    return(de_df)
  })
  
  log_message("DE analysis complete")
  return(de_results)
}

#' Generate volcano plot
#' @param de_df Data frame with DE results
#' @param config Configuration object
#' @return ggplot object
generate_volcano <- function(de_df, config) {
  log_message("Generating volcano plot...")
  
  # Define thresholds
  logfc_threshold <- config$de$logfc_threshold
  pvalue_threshold <- config$de$pvalue_cutoff
  
  de_df <- de_df %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        avg_log2FC > logfc_threshold & p_val_adj < pvalue_threshold ~ "Significant",
        TRUE ~ "NS"
      )
    )
  
  # Create volcano plot
  p <- ggplot(de_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significance, size = -log10(p_val_adj)), alpha = 0.6) +
    scale_color_manual(values = c("Up" = "#E41A36", 
                                        "Down" = "#3366CC",
                                        "NS" = "#7BAFD0")) +
    geom_vline(xintercept = -logfc_threshold, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = "grey50") +
    geom_text_repel(aes(label = ifelse(avg_log2FC > logfc_threshold & abs(avg_log2FC) > 2, 
                                                   as.character(de_df$gene)[1:10]), 
                                       vjust = -1,
                                       color = ifelse(avg_log2FC > 0, "red", "blue")),
                       fontface = "italic",
                       size = 3) +
    labs(x = "Log2 Fold Change", 
         y = "-Log10 Adjusted P-value",
         color = "Regulation") +
    theme_bw() +
    plot_annotation(title = "Volcano Plot")
  
  return(p)
}

#' Generate heatmap of top DE genes
#' @param obj Seurat object
#' @param de_df Data frame with DE results
#' @param n_top Number of top genes to show
#' @param config Configuration object
#' @return ggplot object
generate_de_heatmap <- function(obj, de_df, n_top = 20, config) {
  log_message("Generating DE heatmap (top ", n_top, " genes per cluster)...")
  
  # Get top genes per cluster
  top_genes <- de_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc = avg_log2FC) %>%
    dplyr::slice_head(n = n_top)
  
  genes_to_plot <- unique(top_genes$gene)
  
  # Get expression data
  expr_data <- GetAssayData(obj, slot = "scale.data")
  expr_subset <- expr_data[genes_to_plot, ]
  
  # Order by DE score and cluster
  de_df$rank <- as.numeric(rownames(de_df))
  de_df <- de_df %>% arrange(rank)
  
  # Subset to top genes
  de_top <- de_df %>% slice_head(n = n_top * length(unique(de_df$cluster)))
  
  # Create annotation for heatmap
  anno_col <- data.frame(
    Cluster = de_top$cluster,
    Log2FC = de_top$avg_log2FC,
    P_adj = de_top$p_val_adj,
    Regulation = de_top$regulated
  )
  
  # Scale data for heatmap
  if (length(genes_to_plot) > 0) {
    expr_to_plot <- expr_subset[rev(genes_to_plot), ]
    expr_scaled <- t(scale(t(expr_to_plot)))
  } else {
    expr_scaled <- matrix(nrow = 0, ncol = 0)
    anno_col <- data.frame()
  }
  
  # Create heatmap
  p <- pheatmap(
    as.matrix(expr_scaled),
    annotation_col = anno_col,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    color = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#F44E3B")),
    main = paste0("Top ", n_top, " DE Genes"),
    fontsize_row = 6,
    fontsize_col = 8,
    border_color = "white",
    cellwidth = 20,
    cellheight = 0.2,
    scale = "row"
  )
  
  return(p)
}

#' Generate MA plot
#' @param de_df Data frame with DE results
#' @param config Configuration object
#' @return ggplot object
generate_ma_plot <- function(de_df, config) {
  log_message("Generating MA plot...")
  
  # Create MA plot
  p <- ggplot(de_df, aes(x = (avg_log2FC + 1)/2, y = avg_log2FC)) +
    geom_point(aes(color = ifelse(avg_log2FC > 0, "red", "blue"), alpha = 0.6)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(xintercept = log2(config$de$logfc_threshold) / 2, 
               linetype = "dashed", color = "grey50") +
    labs(x = "M (Mean expression)", 
         y = "A (Log2 Fold Change)") +
    theme_bw() +
    plot_annotation(title = "MA Plot")
  
  return(p)
}

#' Main function for DE analysis
#' @param obj Annotated Seurat object
#' @param config Configuration object
#' @param groups_to_compare List of groups to compare (optional)
#' @return List with DE results and plots
run_diffexp <- function(obj, groups_to_compare = NULL, config) {
  log_message("Module 05: Differential Expression Analysis")
  
  # Get current identity
  current_idents <- levels(Idents(obj))
  log_message("Available groups: ", paste(current_idents, collapse = ", "))
  
  # Run DE analysis
  de_results <- run_de_analysis(obj, groups_to_compare, config)
  
  # Generate visualizations
  if (length(de_results) > 0) {
    de_df <- de_results[[1]]  # Use first comparison for visualization
    
    volcano_p <- generate_volcano(de_df, config)
    de_heatmap <- generate_de_heatmap(obj, de_df, config = n_top = 20)
    ma_p <- generate_ma_plot(de_df, config)
    
    plots <- list(
      volcano = volcano_p,
      heatmap = de_heatmap,
      ma = ma_p
    )
  } else {
    plots <- list()
  }
  
  # Summary statistics
  if (length(de_results) > 0) {
    summary <- data.frame(
      metric = c("Total DEGs (up+down)", "Upregulated", "Downregulated"),
      value = c(
        nrow(de_df[avg_log2FC > 0]) + nrow(de_df[avg_log2FC < 0]),
        nrow(de_df[avg_log2FC > 0]),
        nrow(de_df[avg_log2FC < 0])
      )
    )
  } else {
    summary <- data.frame(
      metric = c("Total DEGs (up+down)", "Upregulated", "Downregulated"),
      value = c(0, 0, 0)
    )
  }
  
  log_message("Module 05 complete")
  return(list(
    de_results = de_results,
    plots = plots,
    summary = summary
  ))
}

#' Export DEG results to CSV
#' @param de_results List of DE data frames
#' @param config Configuration object
#' @return Invisibly NULL
export_deg_csv <- function(de_results, config) {
  if (config$output$format %in% "csv") {
    output_dir <- file.path(config$output$dir, "enrichment")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    for (i in seq_along(de_results)) {
      de_df <- de_results[[i]]
      filename <- file.path(output_dir, paste0("DEG_", de_df$comparison[1], ".csv"))
      write.csv(de_df, file = filename, row.names = FALSE)
      log_message("Exported: ", filename)
    }
  }
  }
  
  return(invisible(NULL))
}
