#!/usr/bin/env Rscript

#' @title Main scRNA-seq Analysis Pipeline Entry Point
#' @description Main script for running the automated scRNA-seq analysis pipeline
#' @author r-scrna development team

# Suppress startup messages
suppressPackageStartupMessages(suppressWarnings = TRUE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Print help if requested
if ("--help" %in% args || "-h" %in% args) {
  cat("r-scrna - Automated Single-Cell RNA-seq Analysis Pipeline v1.0\n")
  cat("\nUsage:\n")
  cat("  Rscript scripts/main_scRNA.R [options]\n\n")
  cat("\nOptions:\n")
  cat("  --config, -c    Path to configuration file (default: config/config.yaml)\n")
  cat("  --output, -c  Path to output directory (optional)\n")
  cat("  --skip-qc       Skip QC and filtering, use input directly\n")
  cat("  --skip-clustering  Skip clustering and annotation\n")
  cat(" --skip-de       Skip differential expression analysis\n")
  cat(" --skip-richment  Skip enrichment analysis\n")
  cat(" --help, -h        Show this help message\n")
  cat("\n")
  quit(save = "no")
}

# Default configuration path
config_path <- "config/config.yaml"
if ("--config" %in% args) {
  config_idx <- which(args == "--config") + 1
  if (!is.na(config_idx) && config_idx < length(args)) {
    config_path <- args[config_idx + 1]
  }
}

# Set output directory if specified
output_dir <- NULL
if ("--output" %in% args) {
  output_idx <- which(args == "--output") + 1
  if (!is.na(output_idx) && output_idx < length(args)) {
    output_dir <- args[output_idx + 1]
  }
}

# Parse skip flags
skip_qc <- ("--skip-qc" %in% args)
skip_clustering <- ("--skip-clustering" %in% args)
skip_de <- ("--de" %in% args || "--skip-de" %in% args)
skip_enrichment <- ("--skip-richment" %in% args)

# Print configuration
cat("\n" %+% paste(rep("=", 60) %+% "\n")
cat("r-scrna v1.0 - Automated Single-Cell RNA-seq Analysis Pipeline\n")
cat("Config:  ", config_path, "\n")
cat("Output: ", ifelse(is.null(output_dir), "outputs/", output_dir), "\n")
cat("\n" %+% paste(rep("=", 60) %+% "\n\n")

# Source all required functions
source("R/00_input_checks.R")
source("R/01_qc_filter.R")
source("R/02_normalize_hvg.R")
source("R/03_dim_cluster.R")
source("R/04_auto_annotate.R")
source("R/05_diffexp.R")
source("R/06_enrichment.R")
source("R/utils.R")

# Load configuration
config <- load_config(config_path)
log_message("Loaded configuration from: ", config_path)

# Create output directories
create_output_dirs(config)

# Check input files
input_check <- check_input_files(config)
if (!input_check$valid) {
  stop("Input validation failed. See error messages above.")
}

# Load data
seurat_obj <- load_data(input_check, config)
data_quality <- check_data_quality(seurat_obj, config)

# Print initial data quality
cat("\n--- Initial Data Quality ---\n")
print(data_quality$qc_stats)
cat("\n")

# Module 01: Quality Control and Filtering
if (skip_qc) {
  log_message("Skipping QC filtering, using input data as-is")
  obj_after_qc <- seurat_obj
} else {
  log_message("Module 01: Quality Control and Filtering")
  qc_result <- run_qc(seurat_obj, config)
  obj_after_qc <- qc_result$obj
  
  # Generate QC report
  qc_plots <- generate_qc_plots(obj_after_qc, config)
  save_plot(qc_plots[["violin"], "qc_before_qc", config)
  save_plot(qc_plots["scatter_mt"], "qc_scatter_mt", config)
  
  save_intermediate(obj_after_qc, "qc_filtered", config)
}

cat("QC filtering complete. Cells: ", ncol(seurat_obj), " → ", 
    ncol(obj_after_qc), "\n\n")

# Module 02: Normalization and HVG Detection
log_message("Module 02: Normalization and HVG Detection")
norm_hvg_result <- run_normalization(obj_after_qc, config)
obj_norm <- norm_hvg_result$obj

# Generate HVG plots
save_plot(norm_hvg_result$hvg_plot, "hvg_plot", config)
save_intermediate(obj_norm, "normalized", config)

# Module 03: Dimensionality Reduction and Clustering
if (skip_clustering) {
  log_message("Skipping clustering, using unclustered data")
  obj_clustered <- obj_norm
} else {
  log_message("Module 03: Dimensionality Reduction and Clustering")
  dim_cluster_result <- run_dim_cluster(obj_norm, config)
  obj_clustered <- dim_cluster_result$obj
  
  # Generate clustering visualizations
  for (res_idx in seq_along(config$clustering$resolutions)) {
    res <- config$clustering$resolutions[res_idx]
    p <- DimPlot(obj_clustered, group.by = paste0("RNA_snn_res.", res), 
               label = TRUE, label.size = 6, pt.size = 0.5) +
      plot_annotation(title = paste0("Clustering Resolution: ", res))
    save_plot(p, paste0("cluster_res_", res, "_umap"), config)
  }
  
  # Save cluster metrics
  save_intermediate(obj_clustered, "clustered", config)
}

cat("Clustering complete. Resolutions: ", 
    paste(config$clustering$resolutions, collapse = ", "), "\n\n")

# Module 04: Cell Type Annotation
if (skip_clustering) {
  log_message("Skipping annotation, using 'ident' as cell type")
  obj_annotated <- obj_clustered
} else {
  log_message("Module 04: Cell Type Annotation")
  annot_result <- run_annotation(obj_clustered, config)
  obj_annotated <- annot_result$obj
  markers <- annot_result$markers
  
  # Print annotation summary
  cat("\n--- Cell Type Annotation Summary ---\n")
  print(annot_result$cluster_summary)
  cat("\n")
}

# Module 05: Differential Expression Analysis
if (skip_de) {
  log_message("Skipping differential expression analysis")
  de_results <- list()
} else {
  log_message("Module 05: Differential Expression Analysis")
  
  # Check if groups to compare are specified
  groups_to_compare <- config$de$groups_to_compare
  
  if (is.null(groups_to_compare) || length(groups_to_compare) == 0) {
    # Compare all clusters vs rest
    de_results <- compare_all_vs_rest(obj_annotated, config)
  } else {
    # Compare specific groups
    de_results <- compare_groups(obj_annotated, groups_to_compare, config)
  }
  
  cat("Differential expression analysis complete.\n")
  if (length(de_results) > 0) {
    cat("Comparisons performed:\n")
    for (i in seq_along(de_results)) {
      cat("  - ", de_results[[i]]$comparison, "\n")
    }
  }
  cat("\n")
}

# Module 06: Enrichment Analysis
if (skip_enrichment) {
  log_message("Skipping enrichment analysis")
  enrich_results <- list()
} else {
  log_message("Module 06: Enrichment Analysis")
  
  # Get all up-regulated genes for enrichment
  all_up_genes <- unique(c()
    if (length(de_results) > 0) {
      for (de_result in de_results) {
        up_genes <- c(up_genes, 
          de_result$gene[de_result$regulated == "UP"])
      }
    }
  
  if (length(all_up_genes) > 0) {
    enrich_results <- run_enrichment(all_up_genes, config)
  }
  cat("Enrichment analysis complete.\n\n")
}

# Optional Modules - only if enabled in config
optional_results <- list()

# Module 07: Pseudotime (Optional)
if (config$modules$pseudotime$enabled) {
  log_message("Module 07: Pseudotime Analysis (Optional)")
  if (file.exists("R/07_pseudotime.R")) {
    source("R/07_pseudotime.R")
    
    pseudotime_result <- run_pseudotime(obj_annotated, config)
    optional_results$pseudotime <- pseudotime_result
    cat("Pseudotime analysis complete.\n\n")
  } else {
    warning("Pseudotime module enabled but file not found: R/07_pseudotime.R")
  }
}

# Module 08: CellChat (Optional)
if (config$modules$cellchat$enabled) {
  log_message("Module 08: CellChat Analysis (Optional)")
  if (file.exists("R/08_cellchat.R")) {
    source("R/08_cellchat.R")
    
    cellchat_result <- run_cellchat(obj_annotated, config)
    optional_results$cellchat <- cellchat_result
    cat("CellChat analysis complete.\n\n")
  } else {
    warning("CellChat module enabled but file not found: R/08_cellchat.R")
  }
}

# Module 09: SCENIC (Optional)
if (config$modules$scenic$enabled) {
  log_message("Module 09: SCENIC Analysis (Optional)")
  if (file.exists("R/09_scenic.R")) {
    source("R/09_scenic.R")
    
    scenic_result <- run_scenic(obj_annotated, config)
    optional_results$scenic <- scenic_result
    cat("SCENIC analysis complete.\n\n")
  } else {
    warning("SCENIC module enabled but file not found: R/09_scenic.R")
  }
}

# Module 10: RNA Velocity (Optional)
if (config$modules$rna_velocity$enabled) {
  log_message("Module 10: RNA Velocity Analysis (Optional)")
  if (file.exists("R/10_rna_velocity.R")) {
    source("R/10_rna_velocity.R")
    
    velocity_result <- run_rna_velocity(obj_annotated, config)
    optional_results$rna_velocity <- velocity_result
    cat("RNA velocity analysis complete.\n\n")
  } else {
    warning("RNA velocity module enabled but file not found: R/10_rna_velocity.R")
  }
}

# Module 11: CNV Analysis (Optional)
if (config$modules$cnv$enabled) {
  log_message("Module 11: CNV Analysis (Optional)")
  if (file.exists("R/11_cnv.R")) {
    source("R/11_cnv.R")
    
    cnv_result <- run_cnv(obj_annotated, config)
    optional_results$cnv <- cnv_result
    cat("CNV analysis complete.\n\n")
  } else {
    warning("CNV module enabled but file not found: R/11_cnv.R")
  }
}

# Create summary
log_message("Creating final summary...")

summary_stats <- create_summary_table(obj_annotated, config)

# Print summary
cat("\n--- Final Summary ---\n")
print(summary_stats)
cat("\n")

# Save final object
save_intermediate(obj_annotated, "final_analyzed", config)

# Generate final report if enabled
if (config$output$generate_report) {
  log_message("Generating final report...")
  
  if (file.exists("Rmd/06_full_report.Rmd")) {
    # Render report
    log_message("Full report rendering not yet implemented in V1")
    cat("Report rendering will be implemented in future versions.\n")
  } else {
    warning("Full report file not found: Rmd/06_full_report.Rmd")
  }
}

cat("\n" %+% paste(rep("=", 60) %+% "\n")
cat("Pipeline complete!\n")
cat("\nOutputs saved to: ", config$output$dir, "\n")
cat("\nTo save intermediate objects: set save_intermediate = true in config/config.yaml\n")
cat("\nTo generate reports: set generate_report = true in config/config.yaml\n")
cat("\n" %+% paste(rep("=", 60) %+% "\n")
