#' @title Module 06: Enrichment Analysis
#' @description Perform GO and KEGG pathway enrichment analysis
#' @author r-scrna development team

library(Seurat)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#' Load annotation database
#' @param config Configuration object
#' @return Annotation database object
load_annotation_db <- function(config) {
  species_settings <- get_species_settings(config)
  org_db <- species_settings$org_db
  
  log_message("Loading annotation database: ", org_db)
  
  tryCatch({
    if (org_db == "org.Hs.eg.db") {
      library(org.Hs.eg.db)
      db <- org.Hs.eg.db::org.Hs.eg.db
    } else if (org_db == "org.Mm.eg.db") {
      library(org.Mm.eg.db)
      db <- org.Mm.eg.db::org.Mm.eg.db::org.Mm.eg.db
    } else {
      stop("Unsupported annotation database: ", org_db)
    }
  }, error = function(e) {
    stop("Failed to load annotation database: ", e$message)
  })
  
  return(db)
}

#' Convert gene symbols to Entrez IDs
#' @param genes Vector of gene symbols
#' @param config Configuration object
#' @return Data frame with Entrez IDs
convert_to_entrez <- function(genes, config) {
  species_settings <- get_species_settings(config)
  org_db <- species_settings$org_db
  
  log_message("Converting ", length(genes), " genes to Entrez IDs...")
  
  tryCatch({
    if (org_db == "org.Hs.eg.db") {
      gene_df <- bitr(genes, fromType = "SYMBOL", 
                     toType = c("ENTREZID", "ENSEMBL"),
                     OrgDb = org.Hs.eg.db)
    } else if (org_db == "org.Mm.eg.db") {
      gene_df <- bitr(genes, fromType = "SYMBOL", 
                     toType = c("ENTREZID", "ENSEMBL"),
                     OrgDb = org.Mm.eg.db)
    } else {
      stop("Unsupported annotation database: ", org_db)
    }
  }, error = function(e) {
    warning("Failed to convert some genes: ", e$message)
    return(NULL)
  })
  
  return(gene_df)
}

#' Run GO enrichment analysis
#' @param genes Vector of gene symbols
#' @param config Configuration object
#' @return enrichplot result object
run_go_enrich <- function(genes, config) {
  db <- load_annotation_db(config)
  pvalue_cutoff <- config$enrichment$pvalue_cutoff
  p_adjust_method <- config$enrichment$p_adjust_method
  qvalue_cutoff <- config$enrichment$qvalue_cutoff
  min_gssize <- config$enrichment$min_gssize
  max_gssize <- config$enrichment$max_gssize
  
  log_message("Running GO enrichment for ", length(genes), " genes...")
  
  tryCatch({
    ego <- enrichGO(
      gene         = genes,
      OrgDb        = db,
      ont           = "BP",
      pvalueCutoff  = pvalue_cutoff,
      pAdjustMethod = p_adjust_method,
      qvalueCutoff  = qvalue_cutoff,
      minGSSize      = min_gssize,
      maxGSSize      = max_gssize,
      readable      = TRUE,
      keyType       = "SYMBOL"
    )
    
    log_message("GO enrichment complete: ", nrow(ego@result), " GO terms enriched")
    return(ego)
  }, error = function(e) {
    stop("GO enrichment failed: ", e$message)
  })
}

#' Run KEGG enrichment analysis
#' @param genes Vector of gene symbols
#' @param config Configuration object
#' @return enrichplot result object
run_kegg_enrich <- function(genes, config) {
  db <- load_annotation_db(config)
  pvalue_cutoff <- config$enrichment$pvalue_cutoff
  p_adjust_method <- config$enrichment$p_adjust_method
  qvalue_cutoff <- config$enrichment$qvalue_cutoff
  min_gssize <- config$enrichment$min_gssize
  max_gssize <- config$enrichment$max_gssize
  
  log_message("Running KEGG enrichment for ", length(genes), " genes...")
  tryCatch({
    kegg <- enrichKEGG(
      gene         = genes,
      organism     = db,
      pvalueCutoff = pvalue_cutoff,
      pAdjustMethod = p_adjust_method,
      qvalueCutoff  = qvalue_cutoff,
      minGSSize      = min_gssize,
      maxGSSize      = max_gssize
    )
    
    log_message("KEGG enrichment complete: ", nrow(kegg@result), " KEGG pathways enriched")
    return(kegg)
  }, error = function(e) {
    stop("KEGG enrichment failed: ", e$message)
  })
}

#' Generate enrichment summary
#' @param enrich_result enrichplot result object
#' @param config Configuration object
#' @return Data frame with summary statistics
summarize_enrichment <- function(enrich_result, config) {
  log_message("Generating enrichment summary...")
  
  result_df <- as.data.frame(enrich_result@result)
  
  # Count by category
  summary_stats <- result_df %>%
    dplyr::group_by(ONTOLOGY) %>%
    dplyr::summarise(
      n_terms = n(),
      n_significant = sum(adjusted.p.value < config$enrichment$qvalue_cutoff),
      min_qvalue = min(p.adjust, na.rm = TRUE)
    )
  )
  
  log_message("Enrichment summary:")
  print(summary_stats)
  return(summary_stats)
}

#' Create GO bar plot
#' @param enrich_result enrichplot result object
#' @param n_terms Number of top terms to show
#' @param config Configuration object
#' @return ggplot object
create_go_barplot <- function(enrich_result, n_terms = 20, config) {
  log_message("Creating GO bar plot...")
  top_terms <- enrich_result@result %>%
    dplyr::arrange(p.adjust) %>%
    head(n_terms)
  
  # Create dot plot
  p <- dotplot(enrich_result, showCategory = TRUE, 
                 term = top_terms$escription, 
                 size = p.adjust, 
                 color = "p.adjust",
                 title = paste0("GO Enrichment (Top ", n_terms, " Terms)"),
                 label_format = "{p.adjust:.2e} ~ {p.adjust:.2e}",
                 geneRatio = "geneRatio/%s" %in% n) + ")")
  
  return(p)
}

#' Create KEGG bar plot
#' @param enrich_result enrichplot result object
#' @param n_terms Number of top terms to show
#' @param config Configuration object
#' @return ggplot object
create_kegg_barplot <- function(enrich_result, n_terms = 20, config) {
  log_message("Creating KEGG bar plot...")
  top_terms <- enrich_result@result %>%
    dplyr::arrange(p.adjust) %>%
    head(n_terms)
  
  # Create bar plot
  p <- barplot(enrich_result, showCategory = "KEGG", 
               drop = TRUE, 
               showCategory = TRUE,
               x = top_terms$geneRatio, 
               color = "p.adjust", 
               title = paste0("KEGG Enrichment (Top ", n_terms, " Pathways)"),
               label_format = "{p.adjust:.2e}",
               font.size = 10,
               order = "count")
  )
  
  return(p)
}

#' Create GO dot plot
#' @param enrich_result enrichplot result object
#' @param n_terms Number of top terms to show
#' @param config Configuration object
#' @return ggplot object
create_go_dotplot <- function(enrich_result, n_terms = 15, config) {
  log_message("Creating GO dot plot...")
  top_terms <- enrich_result@result %>%
    dplyr::arrange(p.adjust) %>%
    head(n_terms)
  
  # Create dot plot
  p <- dotplot(enrich_result, showCategory = TRUE, 
                 x = "GeneRatio", 
                 category = "ONTOLOGY",
                 size = p.adjust, 
                 color = "p.adjust", 
                 title = paste0("GO Enrichment (Top ", n_terms, " Terms)"),
                 label_format = "{p.adjust:.2e} ~ {p.adjust:.2e}",
                 label_size = 4)
  )
  return(p)
}

#' Create KEGG dot plot
#' @param enrich_result enrichplot result object
#' @param n_terms Number of top terms to show
#' @param config Configuration object
#' @return ggplot object
create_kegg_dotplot <- function(enrich_result, n_terms = 15, config) {
  log_message("Creating KEGG dot plot...")
  top_terms <- enrich_result@result %>%
    dplyr::arrange(p.adjust) %>%
    head(n_terms)
  
  # Create dot plot
  p <- dotplot(enrich_result, showCategory = "KEGG", 
                 x = "GeneRatio", 
                 category = "Pathway", 
                 size = p.adjust, 
                 color = "p.adjust", 
                 title = paste0("KEGG Enrichment (Top ", n_terms, " Pathways)"),
                 label_format = "{p.adjust:.2e} ~ {p.adjust:.2e}",
                 label_size = 4
  )
  return(p)
}

#' Create enrichment comparison plot
#' @param go_result GO enrichment result
#' @param kegg_result KEGG enrichment result
#' @param config Configuration object
#' @return Combined ggplot object
create_enrichment_comparison <- function(go_result, kegg_result, config) {
  log_message("Creating enrichment comparison plot...")
  
  n_terms <- 10
  
  # Get top terms
  top_go <- go_result@result %>%
    dplyr::arrange(p.adjust) %>%
    head(n_terms)
  
  top_kegg <- kegg_result@result %>%
    dplyr::arrange(p.adjust) %>%
    head(n_terms)
  
  # Create comparison data frame
  go_terms <- top_go$description
  kegg_terms <- top_kegg$description
  
  combined_terms <- unique(c(go_terms, kegg_terms))
  
  n_combined <- length(combined_terms)
  
  # Create comparison table
  comparison_df <- data.frame(
    term = combined_terms,
    go_pvalue = NA,
    go_qvalue = NA,
    kegg_pvalue = NA,
    go_reg = NA,
    kegg_reg = NA
  )
  
  for (i in seq_len(n_combined)) {
    term <- comparison_df$term[i]
    if (term %in% go_terms) {
      go_row <- top_go[top_go$description == term, ]
      comparison_df$go_pvalue[i] <- go_row$p.value
      comparison_df$go_reg[i] <- ifelse(go_row$regulated == "UP", "Enriched", "Enriched")
    }
    if (term %in% kegg_terms) {
      kegg_row <- top_kegg[top_kegg$description == term, ]
      comparison_df$kegg_pvalue[i] <- kegg_row$p.value
      comparison_df$kegg_reg[i] <- ifelse(kegg_row$regulated == "UP", "Enriched", "Enriched")
    }
  }
  
  # Fill missing values
  comparison_df$go_pvalue[is.na(comparison_df$go_pvalue)] <- NA
  comparison_df$kegg_pvalue[is.na(comparison_df$kegg_pvalue)] <- NA
  
  # Create comparison plot
  p <- ggplot(comparison_df, aes(x = term, y = -log10(p.value))) +
    geom_point(aes(color = "Type", shape = "Type"), size = 2, alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    facet_wrap(~Type, scales = "free_y") +
    labs(x = "", y = "-Log10(p.value)", 
           title = "GO vs KEGG Enrichment Comparison") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Main function for enrichment analysis
#' @param genes Vector of genes or list of gene vectors
#' @param config Configuration object
#' @return List with enrichment results
run_enrichment <- function(genes, config) {
  if (!config$enrichment$enabled) {
    log_message("Enrichment disabled, skipping...")
    return(list())
  }
  
  log_message("Module 06: Enrichment Analysis")
  
  # Convert genes to Entrez IDs
  if (is.list(genes)) {
    gene_df <- sapply(genes, function(g) {
      convert_to_entrez(g, config)
    })
    all_genes <- unique(unlist(gene_df$ENTREZID[!is.na(gene_df$ENTREZID)])
  } else {
    gene_df <- convert_to_entrez(genes, config)
    all_genes <- gene_df$ENTREZID[!is.na(gene_df$ENTREZID)]
  }
  
  # Run GO and KEGG enrichment
  go_result <- run_go_enrich(all_genes, config)
  kegg_result <- run_kegg_enrich(all_genes, config)
  
  # Summarize
  go_summary <- summarize_enrichment(go_result, config)
  kegg_summary <- summarize_enrichment(kegg_result, config)
  
  # Create visualizations
  n_terms <- 10
  
  plots <- list()
  plots[["go_bar"]] <- create_go_barplot(go_result, n_terms, config)
  plots[["go_dot"]] <- create_go_dotplot(go_result, n_terms, config)
  plots[["kegg_bar"]] <- create_kegg_barplot(kegg_result, n_terms, config)
  plots[["kegg_dot"]] <- create_kegg_dotplot(kegg_result, n_terms, config)
  plots[["comparison"]] <- create_enrichment_comparison(go_result, kegg_result, config)
  
  # Save results
  output_dir <- file.path(config$output$dir, "enrichment")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save GO results
  write.csv(as.data.frame(go_result@result), 
          file = file.path(output_dir, "go_results.csv"), 
          row.names = FALSE)
  )
  
  # Save KEGG results
  write.csv(as.data.frame(kegg_result@result), 
          file = file.path(output_dir, "kegg_results.csv"), 
          row.names = FALSE
  )
  
  # Save summaries
  write.csv(go_summary, 
          file = file.path(output_dir, "go_summary.csv"), 
          row.names = FALSE
  )
  write.csv(kegg_summary, 
          file = file.path(output_dir, "kegg_summary.csv"), 
          row.names = FALSE
  )
  
  # Save plots
  save_plot(plots[["go_bar"]], "go_barplot", config)
  save_plot(plots[["go_dot"]], "go_dotplot", config)
  save_plot(plots[["kegg_bar"]], "kegg_barplot", config)
  save_plot(plots[["kegg_dot"]], "kegg_dotplot", config)
  save_plot(plots[["comparison"]], "enrichment_comparison", config)
  
  log_message("Module 06 complete")
  return(list(
    go = go_result,
    kegg = kegg_result,
    go_summary = go_summary,
    kegg_summary = kegg_summary,
    plots = plots
  ))
}
