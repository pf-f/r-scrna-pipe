#' @title Module 00: Input Data Validation and Loading
#' @description Check input files and load data into Seurat object
#' @author r-scrna development team

library(Seurat)
library(data.table)
library(R.utils)

#' Validate input files exist
#' @param config Configuration object
#' @return List with validation status and data path
check_input_files <- function(config) {
  input_dir <- config$input$data_dir
  data_type <- tolower(config$input$data_type)
  
  log_message("Checking input files...")
  
  if (!dir.exists(input_dir)) {
    stop("Input directory not found: ", input_dir)
  }
  
  result <- list(
    valid = FALSE,
    data_path = NULL,
    data_type = NULL,
    metadata_path = NULL
  )
  
  # Check different input types
  if (data_type == "cellranger") {
    # Check for filtered_feature_bc_matrix directory
    matrix_dir <- list.files(input_dir, pattern = "filtered_feature_bc_matrix", 
                           full.names = TRUE, recursive = TRUE)
    if (length(matrix_dir) > 0) {
      result$data_path <- file.path(dirname(matrix_dir[1]), basename(matrix_dir[1]))
      result$data_type <- "cellranger"
      result$valid <- TRUE
    } else {
      stop("CellRanger output directory (filtered_feature_bc_matrix) not found")
    }
    
  } else if (data_type == "counts") {
    # Check for count matrix file
    count_files <- list.files(input_dir, pattern = "\\.(txt|tsv|csv|gz|rds)$", 
                         full.names = TRUE)
    if (length(count_files) > 0) {
      result$data_path <- file.path(input_dir, count_files[1])
      result$data_type <- "counts"
      result$valid <- TRUE
    } else {
      stop("No count matrix file found in: ", input_dir)
    }
    
  } else if (data_type == "seurat_rds") {
    # Check for RDS file
    rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(rds_files) > 0) {
      result$data_path <- file.path(input_dir, rds_files[1])
      result$data_type <- "seurat_rds"
      result$valid <- TRUE
    } else {
      stop("No Seurat RDS file found in: ", input_dir)
    }
    
  } else {
    stop("Unsupported data type: ", data_type)
  }
  
  # Check for metadata file
  metadata_path <- config$input$sample_info
  if (!is.null(metadata_path) && !is.na(metadata_path)) {
    full_path <- file.path(input_dir, metadata_path)
    if (file.exists(full_path)) {
      result$metadata_path <- full_path
      log_message("Found metadata file: ", full_path)
    }
  }
  
  log_message("Input validation: ", result$valid)
  return(result)
}

#' Load data based on type
#' @param check_result Result from check_input_files
#' @param config Configuration object
#' @return Seurat object
load_data <- function(check_result, config) {
  data_type <- check_result$data_type
  data_path <- check_result$data_path
  
  log_message("Loading data from: ", data_path)
  
  seurat_obj <- NULL
  
  if (data_type == "cellranger") {
    # Load CellRanger output
    counts <- Read10X(data.dir = data_path)
    project_name <- config$project$name
    seurat_obj <- CreateSeuratObject(
      counts = counts,
      project = project_name,
      min.cells = config$qc$min_cells_per_gene,
      min.features = config$qc$min_features_per_cell
    )
    
  } else if (data_type == "counts") {
    # Load count matrix
    if (grepl("\\.rds$", data_path)) {
      counts <- readRDS(data_path)
    } else {
      counts <- fread(data_path, sep = "\t", header = TRUE)
      rownames(counts) <- counts[[1]]
      counts <- as.matrix(counts[, -1, with = FALSE])
    }
    
    project_name <- config$project$name
    seurat_obj <- CreateSeuratObject(
      counts = counts,
      project = project_name,
      min.cells = config$qc$min_cells_per_gene,
      min.features = config$qc$min_features_per_cell
    )
    
  } else if (data_type == "seurat_rds") {
    # Load Seurat RDS file
    seurat_obj <- readRDS(data_path)
    log_message("Loaded Seurat object from RDS")
  }
  
  # Load metadata if available
  if (!is.null(check_result$metadata_path)) {
    metadata <- read.csv(check_result$metadata_path, row.names = 1)
    
    # Match metadata cells with Seurat object
    common_cells <- intersect(rownames(metadata), colnames(seurat_obj))
    if (length(common_cells) > 0) {
      seurat_obj <- seurat_obj[, common_cells]
      seurat_obj <- AddMetaData(seurat_obj, metadata[common_cells, , drop = FALSE])
      log_message("Added metadata for ", length(common_cells), " cells")
    }
  }
  
  log_message("Loaded Seurat object: ", ncol(seurat_obj), " cells, ", nrow(seurat_obj), " genes")
  
  return(seurat_obj)
}

#' Check data quality before filtering
#' @param obj Seurat object
#' @param config Configuration object
#' @return List with QC statistics
check_data_quality <- function(obj, config) {
  log_message("Checking initial data quality...")
  
  species_settings <- get_species_settings(config)
  mt_pattern <- species_settings$mt_pattern
  
  # Calculate mitochondrial percentage
  obj[["percent.mt"] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  
  # Calculate ribosomal percentage if pattern specified
  if (!is.null(config$species$ribosomal_pattern)) {
    ribo_pattern <- config$species$ribosomal_pattern
    obj[["percent.ribo"] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
  }
  
  # Get initial stats
  qc_stats <- data.frame(
    metric = c("total_cells", "total_genes", "mean_features", "median_features",
                "mean_counts", "median_counts", "mean_mt_pct", "median_mt_pct"),
    value = c(
      ncol(obj),
      nrow(obj),
      round(mean(obj@meta.data$nFeature_RNA), 2),
      round(median(obj@meta.data$nFeature_RNA), 2),
      round(mean(obj@meta.data$nCount_RNA), 2),
      round(median(obj@meta.data$nCount_RNA), 2),
      round(mean(obj@meta.data$percent.mt), 2),
      round(median(obj@meta.data$percent.mt), 2)
    )
  )
  
  log_message("Initial data quality:")
  print(qc_stats)
  
  return(list(obj = obj, qc_stats = qc_stats))
}
