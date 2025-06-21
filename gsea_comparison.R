# Load required packages
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggridges)
library(ggsci)
library(pheatmap)
library(tidyverse)
library(ggrepel)

# Function: Perform GSEA analysis
perform_gsea <- function(input_file, prefix) {
  # Read and process input file
  rt <- read.table(input_file, header=TRUE, sep="\t", check.names=FALSE)
  rt <- rt[order(rt[,"logFC"], decreasing=TRUE),]
  logFC <- as.vector(rt[,"logFC"])
  names(logFC) <- as.vector(rt[,"ID"])
  
  # Read Hallmark gene set file
  gmt <- read.gmt("h.all.v2024.1.Hs.symbols.gmt")
  
  # Perform GSEA analysis
  kk <- GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
  kkTab <- as.data.frame(kk)
  kkTab <- kkTab[kkTab$pvalue < 0.05,]
  
  # Save GSEA results
  write.table(kkTab, file=paste0(prefix, "_GSEA.result.Hallmark.txt"),
              sep="\t", quote=FALSE, row.names = FALSE)
  
  # Compute absolute NES and sort
  kkTab$absNES <- abs(kkTab$NES)
  kkTab <- kkTab[order(kkTab$absNES, decreasing=TRUE),]
  rownames(kkTab) <- kkTab$Description
  
  return(list(kk=kk, kkTab=kkTab))
}

# Function: Create improved ridgeplot
create_ridgeplot <- function(gsea_obj, title, color_scheme) {
  ridgeplot(gsea_obj, showCategory=10) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 16),
      plot.title = element_text(size = 16),
      axis.title = element_text(size = 16),
      panel.grid.major = element_line(color = "black"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Enrichment Score",
      title = title,
      y = ""
    ) +
    scale_fill_gradient(
      low = color_scheme[1],
      high = color_scheme[2]
    )
}

# Function: Create comparison heatmap (all pathways)
create_comparison_heatmap <- function(survival_gsea, response_gsea) {
  # Get all significant pathways
  all_pathways <- unique(c(
    rownames(survival_gsea$kkTab),
    rownames(response_gsea$kkTab)
  ))
  
  # Create comparison matrix
  comparison_matrix <- matrix(0, nrow=length(all_pathways), ncol=4)
  rownames(comparison_matrix) <- all_pathways
  colnames(comparison_matrix) <- c("Survival_NES", "Survival_pval",
                                   "Response_NES", "Response_pval")
  
  # Fill in data
  for(pathway in all_pathways) {
    if(pathway %in% rownames(survival_gsea$kkTab)) {
      comparison_matrix[pathway, "Survival_NES"] <-
        survival_gsea$kkTab[pathway, "NES"]
      comparison_matrix[pathway, "Survival_pval"] <-
        -log10(survival_gsea$kkTab[pathway, "pvalue"])
    }
    if(pathway %in% rownames(response_gsea$kkTab)) {
      comparison_matrix[pathway, "Response_NES"] <-
        response_gsea$kkTab[pathway, "NES"]
      comparison_matrix[pathway, "Response_pval"] <-
        -log10(response_gsea$kkTab[pathway, "pvalue"])
    }
  }
  
  # Plot heatmap
  pheatmap(comparison_matrix,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           scale = "row",
           show_rownames = TRUE,
           fontsize_row = 14,
           fontsize_col = 14,
           angle_col = 45,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}

# Function: Create bidirectional barplot comparison
create_comparison_barplot <- function(survival_gsea, response_gsea, topN = 10) {
  # Prepare data
  survival_top <- head(survival_gsea$kkTab, topN)
  response_top <- head(response_gsea$kkTab, topN)
  
  # Combine data
  plot_data <- rbind(
    data.frame(
      Pathway = survival_top$Description,
      NES = survival_top$NES,
      Type = "Survival",
      Significance = -log10(survival_top$pvalue)
    ),
    data.frame(
      Pathway = response_top$Description,
      NES = response_top$NES,
      Type = "Response",
      Significance = -log10(response_top$pvalue)
    )
  )
  
  # Create plot
  ggplot(plot_data, aes(x=reorder(Pathway, NES), y=NES, fill=Type)) +
    geom_bar(stat="identity", position="dodge") +
    coord_flip() +
    scale_fill_manual(values=c("Survival"="#E65100", "Response"="#1565C0")) +
    theme_minimal() +
    labs(x="", y="Normalized Enrichment Score", title="") +
    theme(
      axis.text.x  = element_text(size=16,  color="black"),
      axis.text.y  = element_text(size=18,  color="black"),
      axis.title.x = element_text(size=18,  color="black"),
      axis.title.y = element_text(size=18,  color="black"),
      plot.title   = element_text(size=20, color="black"),
      legend.title = element_text(size=16,  color="black"),
      legend.text  = element_text(size=16,  color="black"),
      legend.position="bottom"
    )
}

# Get pathways used in barplot
get_barplot_paths <- function(survival_gsea, response_gsea, topN = 10) {
  survival_top <- head(survival_gsea$kkTab, topN)
  response_top <- head(response_gsea$kkTab, topN)
  unique(c(as.character(survival_top$Description), as.character(response_top$Description)))
}

# Function: Only plot heatmap for barplot pathways
create_comparison_heatmap_barplot <- function(survival_gsea, response_gsea, barplot_paths) {
  # Create comparison matrix
  comparison_matrix <- matrix(0, nrow=length(barplot_paths), ncol=4)
  rownames(comparison_matrix) <- barplot_paths
  colnames(comparison_matrix) <- c("Survival_NES", "Survival_pval",
                                   "Response_NES", "Response_pval")
  
  # Fill in data
  for(pathway in barplot_paths) {
    if(pathway %in% survival_gsea$kkTab$Description) {
      idx <- which(survival_gsea$kkTab$Description == pathway)
      comparison_matrix[pathway, "Survival_NES"] <- survival_gsea$kkTab$NES[idx]
      comparison_matrix[pathway, "Survival_pval"] <- -log10(survival_gsea$kkTab$pvalue[idx])
    }
    if(pathway %in% response_gsea$kkTab$Description) {
      idx <- which(response_gsea$kkTab$Description == pathway)
      comparison_matrix[pathway, "Response_NES"] <- response_gsea$kkTab$NES[idx]
      comparison_matrix[pathway, "Response_pval"] <- -log10(response_gsea$kkTab$pvalue[idx])
    }
  }
  
  # Plot heatmap
  pheatmap(comparison_matrix,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           scale = "row",
           show_rownames = TRUE,
           fontsize_row = 18,
           fontsize_col = 16,
           angle_col = 45,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}

# Common pathway analysis function
# Scatter plot
create_common_pathway_scatter <- function(survival_gsea, response_gsea, common_paths) {
  # Prepare data
  common_data <- data.frame(
    Pathway = common_paths,
    Survival_NES = survival_gsea$kkTab[common_paths, "NES"],
    Response_NES = response_gsea$kkTab[common_paths, "NES"],
    Survival_pvalue = survival_gsea$kkTab[common_paths, "pvalue"],
    Response_pvalue = response_gsea$kkTab[common_paths, "pvalue"]
  )
  
  ggplot(common_data, aes(x = Survival_NES, y = Response_NES)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    geom_point(aes(size = -log10(Survival_pvalue),
                   color = -log10(Response_pvalue)),
               alpha = 0.7) +
    geom_text_repel(aes(label = Pathway),
                    size = 6,
                    max.overlaps = 20) +
    theme_minimal() +
    labs(x = "Survival NES",
         y = "Response NES",
         title = "Common Pathways Comparison",
         size = "-log10(Survival p-value)",
         color = "-log10(Response p-value)") +
    scale_color_gradient(low = "blue", high = "red") +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          legend.position = "right")
}

# Bidirectional NES barplot for common pathways
create_common_pathway_comparison <- function(survival_gsea, response_gsea, common_paths) {
  # Prepare data
  comparison_data <- data.frame(
    Pathway = rep(common_paths, 2),
    NES = c(survival_gsea$kkTab[common_paths, "NES"],
            response_gsea$kkTab[common_paths, "NES"]),
    Type = rep(c("Survival", "Response"), each = length(common_paths)),
    Pvalue = -log10(c(survival_gsea$kkTab[common_paths, "pvalue"],
                      response_gsea$kkTab[common_paths, "pvalue"]))
  )
  
  # Order by average absolute NES
  avg_nes <- tapply(abs(comparison_data$NES), comparison_data$Pathway, mean)
  pathway_order <- names(sort(avg_nes, decreasing = TRUE))
  comparison_data$Pathway <- factor(comparison_data$Pathway, levels = pathway_order)
  
  ggplot(comparison_data, aes(x = Pathway, y = NES, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    scale_fill_manual(values = c("Survival" = "#E65100", "Response" = "#1565C0")) +
    theme_minimal() +
    labs(x = "",
         y = "Normalized Enrichment Score",
         title = "Common Pathways NES Comparison") +
    theme(axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom")
}

# Heatmap for common pathways
create_common_pathway_heatmap <- function(survival_gsea, response_gsea, common_paths) {
  # Prepare data
  heatmap_data <- matrix(0, nrow = length(common_paths), ncol = 4)
  rownames(heatmap_data) <- common_paths
  colnames(heatmap_data) <- c("Survival_NES", "Survival_Sig",
                              "Response_NES", "Response_Sig")
  
  # Fill in data
  heatmap_data[, 1] <- survival_gsea$kkTab[common_paths, "NES"]
  heatmap_data[, 2] <- -log10(survival_gsea$kkTab[common_paths, "pvalue"])
  heatmap_data[, 3] <- response_gsea$kkTab[common_paths, "NES"]
  heatmap_data[, 4] <- -log10(response_gsea$kkTab[common_paths, "pvalue"])
  
  # Column annotation
  annotation_col <- data.frame(
    Type = c("NES", "Significance", "NES", "Significance"),
    Group = c("Survival", "Survival", "Response", "Response")
  )
  rownames(annotation_col) <- colnames(heatmap_data)
  
  # Define annotation colors
  annotation_colors <- list(
    Type = c(NES = "red", Significance = "blue"),
    Group = c(Survival = "#E65100", Response = "#1565C0")
  )
  
  # Plot heatmap
  pheatmap(heatmap_data,
           main = "Common Pathways Analysis",
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           scale = "none",
           show_rownames = TRUE,
           fontsize_row = 12,
           fontsize_col = 12,
           angle_col = 45,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}

# ------------------------------ Main Program ------------------------------

# Run survival analysis GSEA
survival_gsea <- perform_gsea("survival_TCGA.diffall.edgeR.txt", "survival")
pdf("survival_GSEA.ridgeplot.continuous7.pdf", width=12, height=8)
print(create_ridgeplot(survival_gsea$kk,
                       "Survival GSEA Hallmark Pathways",
                       c("#FFF3E0", "#E65100")))
dev.off()

# Run response prediction GSEA
response_gsea <- perform_gsea("response_TCGA.diffall.edgeR.txt", "response")
pdf("response_GSEA.ridgeplot.continuous7.pdf", width=12, height=8)
print(create_ridgeplot(response_gsea$kk,
                       "Response GSEA Hallmark Pathways",
                       c("#E3F2FD", "#1565C0")))
dev.off()

# Create comparison heatmap (all pathways)
pdf("GSEA_comparison_heatmap.pdf", width=12, height=10)
create_comparison_heatmap(survival_gsea, response_gsea)
dev.off()

# Create comparison barplot (adjust barplot_topN as needed)
barplot_topN <- 10
pdf("GSEA_comparison_barplot.pdf", width=12, height=10)
print(create_comparison_barplot(survival_gsea, response_gsea, topN = barplot_topN))
dev.off()

# Use barplot pathways for heatmap
barplot_paths <- get_barplot_paths(survival_gsea, response_gsea, topN = barplot_topN)
pdf("GSEA_comparison_heatmap_barplotpaths.pdf", width=12, height=8)
create_comparison_heatmap_barplot(survival_gsea, response_gsea, barplot_paths)
dev.off()

# Output common pathway analysis
survival_paths <- rownames(survival_gsea$kkTab)
response_paths <- rownames(response_gsea$kkTab)
common_paths <- intersect(survival_paths, response_paths)

# Save common pathway information
write.table(
  data.frame(
    Pathway = common_paths,
    Survival_NES = survival_gsea$kkTab[common_paths, "NES"],
    Survival_pvalue = survival_gsea$kkTab[common_paths, "pvalue"],
    Response_NES = response_gsea$kkTab[common_paths, "NES"],
    Response_pvalue = response_gsea$kkTab[common_paths, "pvalue"]
  ),
  "common_pathways_comparison.txt",
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

# Create scatter plot comparison for common pathways
pdf("common_pathways_scatter.pdf", width = 12, height = 10)
print(create_common_pathway_scatter(survival_gsea, response_gsea, common_paths))
dev.off()

# Create bidirectional NES comparison barplot for common pathways
pdf("common_pathways_nes_comparison.pdf", width = 12, height = max(8, length(common_paths) * 0.3))
print(create_common_pathway_comparison(survival_gsea, response_gsea, common_paths))
dev.off()

# Create heatmap for common pathways
pdf("common_pathways_heatmap.pdf", width = 12, height = max(8, length(common_paths) * 0.3))
create_common_pathway_heatmap(survival_gsea, response_gsea, common_paths)
dev.off()