# Load necessary libraries
library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(dplyr)

# # Directory
# setwd("Your/Working/Directory") # Set your working directory as needed

# Read expression matrix file
data <- read.table("TCGA_LIHC_TPM_mRNA.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
# Convert to matrix
dimnames <- list(rownames(data), colnames(data))
mat <- matrix(as.numeric(as.matrix(data)), nrow = nrow(data), dimnames = dimnames)

# Read gene set file
geneSet <- getGmt("immune.gmt", geneIdType = SymbolIdentifier())

# ssgsea analysis (GSVA >= 1.52 syntax)
gsvapar <- gsvaParam(mat, geneSet, kcdf = 'Gaussian', absRanking = TRUE)
ssgseaScore <- gsva(gsvapar)

# Define ssGSEA score normalization function
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
# Normalize ssGSEA score
rt <- normalize(ssgseaScore)
ssgseaOut <- rbind(id = colnames(rt), rt)
write.table(ssgseaOut, file = "immScore.txt", sep = "\t", quote = FALSE, col.names = FALSE)
rt <- t(rt)
# Keep only first 12 characters in sample names
rownames(rt) <- substr(rownames(rt), 1, 12)

# Read survival data
cli <- read.table("clinical_survival_bio_CT.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
group <- cli$group

# Get common samples
sameSample <- intersect(row.names(rt), row.names(cli))
# Extract common samples
rt <- rt[sameSample, , drop = FALSE]
cli <- cli[sameSample, "group", drop = FALSE]
rt1 <- cbind(rt, cli)

# Prepare data for plotting
rt_long <- melt(rt1, id.vars = c("group"))
colnames(rt_long) <- c("Risk", "Type", "Score")
rt_long$Risk <- factor(rt_long$Risk, levels = c("1", "0"), labels = c("Low risk", "High risk"))

# Boxplot
p <- ggboxplot(rt_long, x = "Type", y = "Score", color = "Risk",
               ylab = "Score", add = "none", xlab = "", palette = c("Orange2", "DarkOrchid"))
p <- p + rotate_x_text(50)
p <- p + theme(
  panel.background = element_rect(fill = "transparent", color = NA),
  plot.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent", color = NA),
  legend.box.background = element_rect(fill = "transparent", color = NA)
)
p <- p + stat_compare_means(aes(group = Risk), symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                                  symbols = c("***", "**", "*", "")),
                            label = "p.signif", size = 8, vjust = 0.9, label.y = max(rt_long$Score) + 0.1)
pdf(file = "imm.pdf", width = 8, height = 5)
print(p)
dev.off()

# Identify significantly different cell types (Wilcoxon test)
significant_cells <- rt_long %>%
  group_by(Type) %>%
  summarise(
    p_value = wilcox.test(Score ~ Risk)$p.value
  ) %>%
  filter(p_value < 0.05) %>%
  pull(Type)

print("Significantly different cell types:")
print(significant_cells)

# Improved violin plot function
create_improved_violin_plot <- function(data, cell_type) {
  y_range <- range(data$Score[data$Type == cell_type])
  y_max <- y_range[2]
  y_position <- y_max + (y_range[2] - y_range[1]) * 0.1
  
  colors <- c(
    "Non-response" = '#F3B1A0',
    "Response" = "#57C3F3"
  )
  
  p <- ggplot(data %>% filter(Type == cell_type), 
              aes(x = Risk, y = Score, fill = Risk)) +
    geom_violin(
      alpha = 0.7,
      trim = FALSE,
      scale = "width"
    ) +
    geom_boxplot(
      width = 0.2,
      alpha = 0.5,
      position = position_dodge(0.9),
      fill = "white",
      outlier.shape = NA
    ) +
    geom_jitter(
      size = 0.8,
      alpha = 0.4,
      width = 0.1,
      color = "gray30"
    ) +
    scale_fill_manual(
      values = colors,
      name = "Group",
      labels = c("Non-response", "Response")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 12,
        face = "bold",
        margin = margin(b = 15)
      ),
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        size = 10,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 10,
        color = "black"
      ),
      axis.title.y = element_text(
        size = 10,
        margin = margin(r = 10)
      ),
      legend.position = "bottom",
      panel.grid.major = element_line(
        color = "gray90",
        size = 0.2
      ),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    ) +
    labs(
      title = cell_type,
      y = "Score"
    )
  
  # Add statistical significance annotation
  p <- p + stat_compare_means(
    aes(group = Risk),
    method = "wilcox.test",
    label = "p.signif",
    size = 12,
    vjust = -0.5,
    label.y = y_position,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    )
  )
  
  return(p)
}

# Create and arrange violin plots
plots_list <- lapply(significant_cells, function(cell) {
  create_improved_violin_plot(rt_long, cell)
})

combined_plot <- ggarrange(
  plotlist = plots_list,
  ncol = 3,
  nrow = ceiling(length(significant_cells)/3),
  common.legend = TRUE,
  legend = "bottom"
)

final_plot <- annotate_figure(
  combined_plot,
  top = text_grob(
    "Immune Cell Type Distributions by Risk Group",
    face = "bold",
    size = 14
  )
)

ggsave(
  "improved_immune_cell_violin_plots.pdf",
  final_plot,
  width = 15,
  height = 5 * ceiling(length(significant_cells)/3) + 1,
  units = "in",
  dpi = 300
)