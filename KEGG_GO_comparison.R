# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringi)
library(GOplot)
library(stringr)
library(tidyverse)
library(pheatmap)

# Set global option
R.utils::setOption("clusterProfiler.download.method", 'auto')

# Function: Perform GO and KEGG analysis
perform_enrichment_analysis <- function(input_file, prefix) {
  # Read input data
  input_diff = read.table(input_file, sep = "\t", header = TRUE, check.names = FALSE)
  input_gene <- input_diff[,1]
  input_gene = unique(as.vector(input_gene))
  
  # Convert gene symbols to Entrez IDs
  entrezIDs = mget(input_gene, org.Hs.egSYMBOL2EG, ifnotfound = NA)
  entrezIDs = as.character(entrezIDs)
  gene = entrezIDs[entrezIDs != "NA"]
  gene = gsub("c\\(\"(\\d+)\".*", "\\1", gene)
  
  # Parameter settings
  pvalueFilter = 0.05
  qvalueFilter = 0.05
  colorSel = ifelse(qvalueFilter > 0.05, "pvalue", "qvalue")
  
  # GO analysis
  kk_go = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
                   qvalueCutoff = 1, ont = "all", readable = TRUE)
  GO = as.data.frame(kk_go)
  GO = GO[(GO$pvalue < pvalueFilter & GO$qvalue < qvalueFilter),]
  
  # KEGG analysis
  kk_kegg = enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
  KEGG = as.data.frame(kk_kegg)
  KEGG$geneID = as.character(sapply(KEGG$geneID, function(x)
    paste(input_gene[match(strsplit(x,"/")[[1]], as.character(entrezIDs))],
          collapse = "/")))
  KEGG = KEGG[(KEGG$pvalue < pvalueFilter & KEGG$qvalue < qvalueFilter),]
  
  # Save results
  write.table(GO, file = paste0(prefix, "_GO.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(KEGG, file = paste0(prefix, "_KEGG.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(list(GO = GO, KEGG = KEGG))
}

# Function: Prepare enrichment analysis data
prepare_enrichment_data <- function(go_result, kegg_result) {
  # Process GO results
  go_df <- go_result %>%
    as.data.frame() %>%
    mutate(
      ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF")),
      Description = gsub("\n", " ", Description),
      Description = trimws(Description)
    ) %>%
    dplyr::select(Description, ONTOLOGY, GeneRatio, pvalue, Count)
  
  # Process KEGG results
  kegg_df <- kegg_result %>%
    as.data.frame() %>%
    mutate(
      ONTOLOGY = "KEGG",
      Description = gsub("\n", " ", Description),
      Description = trimws(Description)
    ) %>%
    dplyr::select(Description, ONTOLOGY, GeneRatio, pvalue, Count)
  
  # Combine data
  combined_df <- bind_rows(go_df, kegg_df) %>%
    mutate(
      GeneRatio = sapply(GeneRatio, function(x) {
        nums <- as.numeric(strsplit(x, "/")[[1]])
        nums[1] / nums[2]
      })
    ) %>%
    group_by(ONTOLOGY) %>%
    slice_head(n = 5) %>%  # Get top 5 per category
    ungroup()
  
  return(combined_df)
}

# Function: Create enrichment bubble plot
create_enrichment_plot <- function(df, title) {
  # Define color palette
  pal <- c(
    'BP' = '#7faced',
    'CC' = '#98d8a0',
    'KEGG' = '#baa4de',
    'MF' = '#f7a799'
  )
  
  # Create plot
  p <- ggplot(df, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
    geom_point(aes(size = Count, color = -log10(pvalue))) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(2, 8)) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    labs(
      title = title,
      x = "Gene Ratio",
      y = "Description",
      size = "Gene Count",
      color = "-log10(p-value)"
    )
  
  return(p)
}

# Function: Create comparison bubble plot
create_comparison_plot <- function(survival_data, response_data) {
  # Find shared pathways
  common_paths <- intersect(survival_data$Description, response_data$Description)
  
  # Prepare comparison data
  comparison_data <- bind_rows(
    survival_data %>% 
      filter(Description %in% common_paths) %>%
      mutate(Group = "Survival"),
    response_data %>%
      filter(Description %in% common_paths) %>%
      mutate(Group = "Response")
  )
  
  # Create comparison bubble plot
  p <- ggplot(comparison_data, 
              aes(x = Group, y = Description, size = Count, color = -log10(pvalue))) +
    geom_point() +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 14)
    ) +
    labs(title = "Comparison of Shared Pathways",
         size = "Gene Count",
         color = "-log10(p-value)")
  
  return(p)
}

# Function: Create similarity heatmap
create_similarity_heatmap <- function(survival_data, response_data) {
  # Get all unique pathways
  all_paths <- unique(c(survival_data$Description, response_data$Description))
  
  # Create matrix
  similarity_matrix <- matrix(0, nrow = length(all_paths), ncol = 2)
  rownames(similarity_matrix) <- all_paths
  colnames(similarity_matrix) <- c("Survival", "Response")
  
  # Fill matrix
  for(path in all_paths) {
    surv_idx <- which(survival_data$Description == path)
    if(length(surv_idx) > 0) {
      similarity_matrix[path, "Survival"] <- -log10(survival_data$pvalue[surv_idx])
    }
    resp_idx <- which(response_data$Description == path)
    if(length(resp_idx) > 0) {
      similarity_matrix[path, "Response"] <- -log10(response_data$pvalue[resp_idx])
    }
  }
  
  # Generate heatmap
  pheatmap(similarity_matrix,
           main = "Pathway Enrichment Similarity",
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           fontsize_row = 12,
           fontsize_col = 12)
}

# Main program
# Run survival analysis
survival_results <- perform_enrichment_analysis("survival_uniCox.txt", "survival")
survival_enrichment <- prepare_enrichment_data(survival_results$GO, survival_results$KEGG)
survival_plot <- create_enrichment_plot(survival_enrichment, "Survival Prediction Enrichment")
ggsave("survival_enrichment_plot.pdf", survival_plot, width = 12, height = max(10, nrow(survival_enrichment) * 0.3))

# Run response prediction analysis
response_results <- perform_enrichment_analysis("response_uniCox.txt", "response")
response_enrichment <- prepare_enrichment_data(response_results$GO, response_results$KEGG)
response_plot <- create_enrichment_plot(response_enrichment, "Response Prediction Enrichment")
ggsave("response_enrichment_plot.pdf", response_plot, width = 12, height = max(10, nrow(response_enrichment) * 0.3))

# Create comparison plot
comparison_plot <- create_comparison_plot(survival_enrichment, response_enrichment)
ggsave("pathway_comparison_plot.pdf", comparison_plot, width = 12,
       height = max(10, length(intersect(survival_enrichment$Description,
                                         response_enrichment$Description)) * 0.3))

# Create similarity heatmap
pdf("pathway_similarity_heatmap.pdf", width = 10, height = 12)
create_similarity_heatmap(survival_enrichment, response_enrichment)
dev.off()

# Output shared and specific pathways
shared_paths <- intersect(survival_enrichment$Description, response_enrichment$Description)
survival_specific <- setdiff(survival_enrichment$Description, response_enrichment$Description)
response_specific <- setdiff(response_enrichment$Description, survival_enrichment$Description)

write.table(data.frame(Pathway = shared_paths), 
            "shared_pathways.txt", row.names = FALSE, quote = FALSE)
write.table(data.frame(Pathway = survival_specific), 
            "survival_specific_pathways.txt", row.names = FALSE, quote = FALSE)
write.table(data.frame(Pathway = response_specific), 
            "response_specific_pathways.txt", row.names = FALSE, quote = FALSE)