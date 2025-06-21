# 设置工作环境
rm(list=ls())
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)
library(ggplot2)



# 读取数据
cli <- read.table("clinical_survival_bio_CT.txt", header=T, sep="\t", check.names=F, row.names=1)
data <- read.table("TCGA_LIHC_count_mRNA.txt", header=T, sep="\t", check.names=F, row.names=1)

# 数据预处理
cat("数据预处理开始...\n")
# 转化为matrix
dimnames <- list(rownames(data), colnames(data))
data <- matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

# 去除低表达的基因
data <- data[rowMeans(data)>1,]
data <- t(data)

# 样本名处理
rownames(data) <- substr(rownames(data),1,12)
rownames(data) <- gsub('[.]', '-', rownames(data))

# 获取共同样本
sameSample <- intersect(row.names(data), row.names(cli))
data <- data[sameSample,]
cli <- cli[sameSample,]
group <- cli$group
data <- t(data)

# 样本分组
conNum <- length(group[group==1])  # responder
treatNum <- length(group[group==0]) # non-responder
Type <- c(rep(1,conNum), rep(2,treatNum))

# 数据排序
data1 <- data[,group == 1]
data2 <- data[,group == 0]
data <- cbind(data1,data2)

# 设计矩阵
Type <- factor(Type)
design <- model.matrix(~0+Type)
rownames(design) <- colnames(data)
colnames(design) <- levels(Type)

# 差异分析
cat("开始差异表达分析...\n")
DGElist <- DGEList(counts = data, group = Type)
keep_gene <- rowSums(cpm(DGElist) > 1) >= 2
cat(sprintf("保留基因数量: %d\n", sum(keep_gene)))

DGElist <- DGElist[keep_gene, , keep.lib.sizes = FALSE]
DGElist <- calcNormFactors(DGElist)
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)

# 拟合模型
fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)

# 保存所有差异结果
write.table(data.frame(ID=rownames(nrDEG_edgeR),nrDEG_edgeR),
            file="TCGA.diffall.edgeR.txt",
            sep="\t", row.names=F, quote=F)

# 差异基因筛选
padj <- 0.05
logFC <- 1
nrDEG_edgeR$regulated <- ifelse(nrDEG_edgeR$logFC>logFC & nrDEG_edgeR$PValue<padj,
                                "up",
                                ifelse(nrDEG_edgeR$logFC<(-logFC) & nrDEG_edgeR$PValue<padj,
                                       "down", "normal"))
# 统计上调和下调基因数量
up_in_responder <- sum(nrDEG_edgeR$regulated == "up")
down_in_responder <- sum(nrDEG_edgeR$regulated == "down")

cat("Responder组上调基因数:", up_in_responder, "\n")
cat("Responder组下调基因数:", down_in_responder, "\n")

# Non-responder组上调和下调基因数与Responder组相反
up_in_non_responder <- down_in_responder
down_in_non_responder <- up_in_responder

cat("Non-responder组上调基因数:", up_in_non_responder, "\n")
cat("Non-responder组下调基因数:", down_in_non_responder, "\n")
# 火山图
cat("绘制火山图...\n")
p1.1 <- ggplot(data=nrDEG_edgeR, 
               aes(x=logFC, y=-log10(PValue), color=regulated)) + 
  geom_point(alpha=0.5, size=1.8) + 
  theme_bw(base_size=20) + 
  xlab("log2FC") + 
  ylab("-log10(Pvalue)") +
  scale_colour_manual(values = c('blue','black','red'))

ggsave("volcano_edgr.pdf", p1.1, width=10, height=8)

# 筛选显著差异基因
outDiff <- nrDEG_edgeR[(nrDEG_edgeR$FDR < padj &
                          (nrDEG_edgeR$logFC>logFC | nrDEG_edgeR$logFC<(-logFC))),]
write.table(data.frame(ID=rownames(outDiff),outDiff),
            file="TCGA.diff.edgeR.txt",
            sep="\t", row.names=F, quote=F)

# 对数转换
hmExp <- log2(hmExp + 0.01)

# 准备注释信息
Type <- data.frame(Type=c(rep("Normal",conNum),rep("Tumor",treatNum)))
rownames(Type) <- colnames(hmExp)





# 加载必要的包
library(ggplot2)
library(ggrepel)  # 用于添加文本标签
library(ggrastr)  # 用于处理大量点的渲染
library(viridis)  # 现代配色方案
library(gridExtra) # 用于组合多个图
library(ggnewscale) # 用于添加多个颜色图层

# 选择top10最显著的基因（基于P值和fold change）
top_genes <- nrDEG_edgeR[order(-abs(nrDEG_edgeR$logFC), nrDEG_edgeR$PValue)[1:10], ]

# 创建优化后的火山图
p <- ggplot(data=nrDEG_edgeR, aes(x=logFC, y=-log10(PValue))) +
  # 添加参考线
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgray", alpha=0.5) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="darkgray", alpha=0.5) +
  
  # 主要散点
  geom_point(aes(color=regulated, size=-log10(PValue)), alpha=0.7) +
  scale_color_manual(values=c("royalblue3", "gray80", "indianred3")) +
  
  # 添加基因标签
  geom_text_repel(
    data = top_genes,
    aes(label = rownames(top_genes)),
    size = 10,
    box.padding = 0.5,
    point.padding = 0.3,
    force = 10,
    segment.color = "grey50",
    segment.alpha = 0.5,
    max.overlaps = Inf
  ) +
  
  # 添加统计信息
  annotate("text", 
           x = max(nrDEG_edgeR$logFC), 
           y = max(-log10(nrDEG_edgeR$PValue)),
           hjust = 1, 
           vjust = 1,
           size = 10,
           label = paste0("Total: ", nrow(nrDEG_edgeR), "\n",
                          "Up: ", sum(nrDEG_edgeR$regulated=="up"), "\n",
                          "Down: ", sum(nrDEG_edgeR$regulated=="down"))) +
  
  # 主题设置
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray95"),  # 淡化主网格线
    panel.grid.minor = element_blank(),                 # 删除次网格线
    # plot.title = element_text(size = 26),
    plot.subtitle = element_text(size = 22)
  ) +
  
  # 标签设置
  labs(
    x = "log2(Fold Change)",
    y = "-log10(P-value)",
    title = "",
    subtitle = "",
    color = "Regulation",
    size = "Significance"
  )+
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), # 完全去除网格线
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "gray95", size = 1),  # 保留XY轴实线
    axis.text.x = element_text(size = 26),  # 横轴刻度字体
    axis.text.y = element_text(size = 26),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.position = "right",
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 26)   # 纵轴刻度字体
  )

# 保存图片
ggsave("volcano_enhanced_labeled.pdf", p, width = 12, height = 10, dpi = 300)

