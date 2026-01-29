#' @title Shared Utility Functions
#' @description Collection of utility functions used across multiple modules
#' @author r-scrna development team
#' @date 2024-01-23

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(yaml)
library(RColorBrewer)

#' Load configuration file
#' @param config_path Path to YAML config file
#' @return List of configuration parameters
load_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  config <- yaml::read_yaml(config_path)
  
  # Validate required config sections
  required_sections <- c("project", "input", "species", "qc", 
                           "normalization", "clustering", "output")
  missing_sections <- setdiff(required_sections, names(config))
  if (length(missing_sections) > 0) {
    warning("Missing config sections: ", paste(missing_sections, collapse = ", "))
  }
  
  return(config)
}

#' Log messages with timestamp
#' @param message Message to log
#' @param level Log level (INFO, WARNING, ERROR)
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s - %s\n", timestamp, level, message))
}

#' Create output directory structure
#' @param config Configuration object
create_output_dirs <- function(config) {
  output_dir <- config$output$dir
  dirs_to_create <- c(
    file.path(output_dir, config$output$intermediate_dir),
    file.path(output_dir, config$output$report_dir),
    file.path(output_dir, "enrichment"),
    file.path(output_dir, "figures")
  )
  
  for (dir_path in dirs_to_create) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      log_message("Created directory: ", dir_path)
    }
  }
}

#' Save intermediate object
#' @param object R object to save
#' @param filename Name for the saved file (without extension)
#' @param config Configuration object
save_intermediate <- function(object, filename, config) {
  if (!config$output$save_intermediate) {
    return(invisible(NULL))
  }
  
  save_dir <- file.path(config$output$dir, config$output$intermediate_dir)
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  file_path <- file.path(save_dir, paste0(filename, ".rds"))
  saveRDS(object, file = file_path)
  log_message("Saved intermediate: ", file_path)
}

#' Load intermediate object
#' @param filename Name of the saved file (without extension)
#' @param config Configuration object
#' @return Loaded R object
load_intermediate <- function(filename, config) {
  save_dir <- file.path(config$output$dir, config$output$intermediate_dir)
  file_path <- file.path(save_dir, paste0(filename, ".rds"))
  
  if (!file.exists(file_path)) {
    stop("Intermediate file not found: ", file_path)
  }
  
  log_message("Loading intermediate: ", file_path)
  return(readRDS(file_path))
}

#' Check if a Seurat object is valid
#' @param obj Seurat object to validate
#' @return TRUE if valid, stop with error message if invalid
validate_seurat <- function(obj) {
  if (!inherits(obj, "Seurat")) {
    stop("Object is not a Seurat object")
  }
  
  if (ncol(obj) == 0) {
    stop("Seurat object has no cells")
  }
  
  if (nrow(obj) == 0) {
    stop("Seurat object has no genes")
  }
  
  return(TRUE)
}

#' Get species-specific settings
#' @param config Configuration object
#' @return List of species settings
get_species_settings <- function(config) {
  species <- tolower(config$species$name)
  
  settings <- switch(species,
    "human" = list(
      mt_pattern = config$species$mt_pattern %||% "^MT-",
      org_db = "org.Hs.eg.db",
      cellchat_db = "CellChatDB.human",
      scenic_org = "hgnc",
      genome = config$species$genome
    ),
    "mouse" = list(
      mt_pattern = config$species$mt_pattern %||% "^mt-",
      org_db = "org.Mm.eg.db",
      cellchat_db = "CellChatDB.mouse",
      scenic_org = "mgi",
      genome = config$species$genome
    ),
    stop("Unsupported species: ", species)
  )
  
  return(settings)
}

#' Generate consistent color palette
#' @param n Number of colors needed
#' @return Vector of hex colors
get_colors <- function(n) {
  if (n <= 8) {
    return(brewer.pal(n, "Set2"))
  } else if (n <= 12) {
    return(brewer.pal(n, "Set3"))
  } else if (n <= 20) {
    return(colorRampPalette(brewer.pal(20, "Spectral"))(n))
  } else {
    return(colorRampPalette(brewer.pal(min(20, n), "Paired"))(n))
  }
}

#' Save plot to file
#' @param plot ggplot object
#' @param filename Output filename
#' @param config Configuration object
#' @param width Plot width in inches
#' @param height Plot height in inches
save_plot <- function(plot, filename, config, width = NULL, height = NULL) {
  width <- width %||% config$visualization$figure_width
  height <- height %||% config$visualization$figure_height
  dpi <- config$visualization$figure_dpi
  
  for (fmt in config$output$plot_formats) {
    filepath <- file.path(config$output$dir, "figures", 
                          paste0(filename, ".", fmt))
    ggsave(filepath, plot = plot, width = width, height = height, dpi = dpi)
    log_message("Saved plot: ", filepath)
  }
}

#' Calculate cluster proportions
#' @param obj Seurat object with cluster annotations
#' @param meta_var Metadata variable containing cluster information
#' @return Data frame with cluster proportions
calc_cluster_props <- function(obj, meta_var = "ident") {
  cluster_col <- obj@meta.data[[meta_var]]
  
  if (is.null(cluster_col)) {
    cluster_col <- Idents(obj)
  }
  
  prop_table <- table(cluster_col) / length(cluster_col)
  prop_df <- data.frame(
    cluster = names(prop_table),
    proportion = as.numeric(prop_table),
    n_cells = as.numeric(table(cluster_col))
  )
  
  return(prop_df)
}

#' Merge multiple Seurat objects
#' @param obj_list List of Seurat objects
#' @param add.cell.ids Prefix for cell barcodes (optional)
#' @return Merged Seurat object
merge_seurat_objects <- function(obj_list, add.cell.ids = NULL, project = "Merged") {
  if (length(obj_list) < 2) {
    warning("Less than 2 objects, no merging performed")
    return(obj_list[[1]])
  }
  
  merged <- Reduce(function(x, y) {
    merge(x, y, add.cell.ids = add.cell.ids, project = project)
  }, obj_list)
  
  log_message("Merged ", length(obj_list), " Seurat objects")
  return(merged)
}

#' Check if required packages are installed
#' @param packages Vector of package names
#' @param stop_if_missing Whether to stop if packages are missing
check_packages <- function(packages, stop_if_missing = TRUE) {
  missing <- character()
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  
  if (length(missing) > 0) {
    msg <- paste("Missing packages:", paste(missing, collapse = ", "))
    if (stop_if_missing) {
      stop(msg)
    } else {
      warning(msg)
    }
  }
  
  return(invisible(NULL))
}

#' Parallel processing setup
#' @param config Configuration object
#' @return NULL (sets up global options)
setup_parallel <- function(config) {
  if (!config$parallel$enabled) {
    return(invisible(NULL))
  }
  
  library(future)
  n_cores <- config$parallel$n_cores
  max_size <- config$parallel$future_global_max_size
  
  if (.Platform$OS.type == "unix") {
    plan("multicore", workers = n_cores)
  } else {
    plan("multisession", workers = n_cores)
  }
  
  options(future.globals.maxSize = max_size * 1024 * 1024^2)
  log_message("Parallel processing enabled: ", n_cores, " cores")
  
  return(invisible(NULL))
}

#' Create a summary table of analysis results
#' @param obj Final Seurat object
#' @param config Configuration object
#' @return Data frame with summary statistics
create_summary_table <- function(obj, config) {
  summary_df <- data.frame(
    metric = c(
      "Total Cells",
      "Total Genes",
      "Clusters",
      "Mean Features per Cell",
      "Median MT Percentage"
    ),
    value = c(
      ncol(obj),
      nrow(obj),
      length(unique(Idents(obj))),
      round(mean(obj@meta.data$nFeature_RNA), 2),
      round(median(obj@meta.data$percent.mt), 2)
    )
  )
  
  return(summary_df)
}
