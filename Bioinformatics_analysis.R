# -------------------------load packages-------------------------------------------------------------------
library(openxlsx)
library(stringr)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(DiffBind)
library(TxDb.Mmusculus.UCSC.mm9.knownGene) 
library(purrr)
library(patchwork)


library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggh4x)
library(sva)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

library(readxl)
library(writexl)
library(ChIPseeker)
library(org.Mm.eg.db)
library(dplyr)

library(openxlsx)
library(tidyverse)
library(ropls)
library(pheatmap)
library(ggrepel)
library(limma)

# -------------------------1.RNAseq---------------------------------------------------------------- -----------------------------------------------------------------
# 1.LNAFv Bulk RNAseq：Preprocessing  --------------------------------------------------------------------
sample_LNAFv<-data.frame(sampleid=c("4408", "4409", "4407", "4406",
                                     "4237", "4402", "4403", "4239",
                                     "4401", "4404", "717", "719",
                                     "143", "144", "1385", "723"),
                          groups=c(rep("Control",4),
                                   rep("StageA",4),
                                   rep("StageB",4),
                                   rep("StageC",4)))
rownames(sample_LNAFv)<-sample_LNAFv$sampleid
count_LNAFv <- read.table('RNAseq/RNAseq.LNAFv.count.txt', header = T, sep = "\t" ,row.names = NULL,check.names = FALSE)
gene_LNAFv<-count_LNAFv[,c("gene_id","gene_name")]
gene_LNAFv$gene_name<-str_to_title(gene_LNAFv$gene_name)
rownames(count_LNAFv)<-count_LNAFv$gene_id
count_LNAFv<-count_LNAFv[,-c(1:2)]
dds_LNAFv <- DESeqDataSetFromMatrix(countData = count_LNAFv,colData = sample_LNAFv,design = ~ groups)
min_count<-30
sample_LNAFv$groups<-factor(sample_LNAFv$groups)
genelist <- DGEList(counts = count_LNAFv, group = relevel(sample_LNAFv$groups, ref = "Control"))
keep <- filterByExpr(genelist, min.count = min_count)
total_gene<-gene_LNAFv[keep,]
dds_LNAFv_PCA <- dds_LNAFv[keep,]
count_LNAFv_PCA<-count_LNAFv[keep,]
vsd <- vst(dds_LNAFv_PCA, blind = TRUE)
vst_matrix <- assay(vsd)
gene_mean <- rowMeans(vst_matrix)
gene_var <- apply(vst_matrix, 1, var)
log2_mean <- log2(gene_mean + 1e-6) 
log2_var <- log2(gene_var + 1e-6)
trend_model <- loess(log2_var ~ log2_mean, span = 0.75)
predicted_log2_var <- predict(trend_model, newdata = data.frame(log2_mean = log2_mean)) 
residuals <- log2_var - predicted_log2_var
hv_genes <- names(residuals[residuals > 1])
hv_matrix <- vst_matrix[hv_genes,]
pca <- prcomp(t(hv_matrix), scale. = FALSE)
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)
pca.dat <- cbind(sample_LNAFv, pca.dat)
pca.dat$groups <- factor(pca.dat$groups)
PCA<-ggplot(pca.dat, aes(PC1, PC2, color = groups)) +geom_point(size=6)+scale_color_manual(values = c("Control" = "black","StageA" = "#EBC678","StageB" = "#EB7815","StageC" = "red"))+labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),y = paste0('PC2: ', pca.var.percent[2], ' %'))+geom_text_repel(aes(label = sampleid),hjust = 0.5,size=5) + theme_minimal()+theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text = element_text(size = 10, color = "black"),axis.ticks = element_line(color="black", linewidth=0.8),axis.title = element_text(size=12))
PCA
samples.to.be.excluded <- c("4239","4404","4406","723")
count_LNAFv<-count_LNAFv[,!colnames(count_LNAFv)%in%samples.to.be.excluded]
sample_LNAFv<-sample_LNAFv[!sample_LNAFv$sampleid %in% samples.to.be.excluded,]
intersect(sample_LNAFv$sampleid,colnames(count_LNAFv))
count_LNAFv<-count_LNAFv[,intersect(sample_LNAFv$sampleid,colnames(count_LNAFv))]
sample_LNAFv<-sample_LNAFv[intersect(sample_LNAFv$sampleid,colnames(count_LNAFv)),]

sample_LNAFv$groups<-factor(sample_LNAFv$groups)
dds_LNAFv <- DESeqDataSetFromMatrix(countData = count_LNAFv,colData = sample_LNAFv,design = ~ groups)
min_count<-30
sample_LNAFv$groups<-factor(sample_LNAFv$groups)
genelist <- DGEList(counts = count_LNAFv, group = relevel(sample_LNAFv$groups, ref = "Control"))
keep <- filterByExpr(genelist, min.count = min_count)#
total_gene<-gene_LNAFv[keep,]
dds_LNAFv <- dds_LNAFv[keep,]
count_LNAFv<-count_LNAFv[keep,]

vsd_LNAFv <- vst(dds_LNAFv, blind = TRUE) 
vst_matrix <- assay(vsd_LNAFv)
gene_mean <- rowMeans(vst_matrix) 
gene_var <- apply(vst_matrix, 1, var)
log2_mean <- log2(gene_mean + 1e-6)
log2_var <- log2(gene_var + 1e-6)
trend_model <- loess(log2_var ~ log2_mean, span = 0.75)
predicted_log2_var <- predict(trend_model, newdata = data.frame(log2_mean = log2_mean))
residuals <- log2_var - predicted_log2_var 
hv_genes <- names(residuals[residuals > 1])
hv_matrix <- vst_matrix[hv_genes, ] 
pca <- prcomp(t(hv_matrix), scale. = FALSE)
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)
pca.dat <- cbind(sample_LNAFv, pca.dat)
pca.dat$groups <- factor(pca.dat$groups)
PCA<-ggplot(pca.dat, aes(PC1, PC2, color = groups)) +geom_point(size=6) +scale_color_manual(values = c("Control" = "black", "StageA" = "#EBC678", "StageB" = "#EB7815", "StageC" = "red"))+labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),y = paste0('PC2: ', pca.var.percent[2], ' %'))+theme_minimal()+theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text = element_text(size = 10, color = "black"), axis.ticks = element_line(color="black", linewidth=0.8),axis.title = element_text(size=12))

# differentially expressed genes -----------------------------------------
dds_LNAFv_diff <- suppressMessages(DESeq(dds_LNAFv, parallel = FALSE))
dds_LNAFv_StageA_Control <- results(dds_LNAFv_diff, alpha = 0.05, contrast = c("groups","StageA","Control")) 
dds_LNAFv_StageB_Control <- results(dds_LNAFv_diff, alpha = 0.05, contrast = c("groups","StageB","Control")) 
dds_LNAFv_StageC_Control <- results(dds_LNAFv_diff, alpha = 0.05, contrast = c("groups","StageC","Control")) 
dds_LNAFv_StageA_StageB <- results(dds_LNAFv_diff, alpha = 0.05, contrast = c("groups","StageA","StageB")) 
dds_LNAFv_StageA_StageC <- results(dds_LNAFv_diff, alpha = 0.05, contrast = c("groups","StageA","StageC")) 
dds_LNAFv_StageB_StageC <- results(dds_LNAFv_diff, alpha = 0.05, contrast = c("groups","StageB","StageC")) 


dds_LNAFv_StageA_Control <- as.data.frame(dds_LNAFv_StageA_Control)
dds_LNAFv_StageA_Control <- na.omit(dds_LNAFv_StageA_Control)
dds_LNAFv_StageB_Control <- as.data.frame(dds_LNAFv_StageB_Control)
dds_LNAFv_StageB_Control <- na.omit(dds_LNAFv_StageB_Control)
dds_LNAFv_StageC_Control <- as.data.frame(dds_LNAFv_StageC_Control)
dds_LNAFv_StageC_Control <- na.omit(dds_LNAFv_StageC_Control)
dds_LNAFv_StageA_StageB <- as.data.frame(dds_LNAFv_StageA_StageB)
dds_LNAFv_StageA_StageB <- na.omit(dds_LNAFv_StageA_StageB)
dds_LNAFv_StageA_StageC <- as.data.frame(dds_LNAFv_StageA_StageC)
dds_LNAFv_StageA_StageC <- na.omit(dds_LNAFv_StageA_StageC)
dds_LNAFv_StageB_StageC <- as.data.frame(dds_LNAFv_StageB_StageC)
dds_LNAFv_StageB_StageC <- na.omit(dds_LNAFv_StageB_StageC)

dds_LNAFv_StageA_Control$gene_id<-rownames(dds_LNAFv_StageA_Control)
dds_LNAFv_StageA_Control<-merge(gene_LNAFv,dds_LNAFv_StageA_Control,by="gene_id")
rownames(dds_LNAFv_StageA_Control)<-dds_LNAFv_StageA_Control$gene_id
dds_LNAFv_StageB_Control$gene_id<-rownames(dds_LNAFv_StageB_Control)
dds_LNAFv_StageB_Control<-merge(gene_LNAFv,dds_LNAFv_StageB_Control,by="gene_id")
rownames(dds_LNAFv_StageB_Control)<-dds_LNAFv_StageB_Control$gene_id
dds_LNAFv_StageC_Control$gene_id<-rownames(dds_LNAFv_StageC_Control)
dds_LNAFv_StageC_Control<-merge(gene_LNAFv,dds_LNAFv_StageC_Control,by="gene_id")
rownames(dds_LNAFv_StageC_Control)<-dds_LNAFv_StageC_Control$gene_id
dds_LNAFv_StageA_StageB$gene_id<-rownames(dds_LNAFv_StageA_StageB)
dds_LNAFv_StageA_StageB<-merge(gene_LNAFv,dds_LNAFv_StageA_StageB,by="gene_id")
rownames(dds_LNAFv_StageA_StageB)<-dds_LNAFv_StageA_StageB$gene_id
dds_LNAFv_StageA_StageC$gene_id<-rownames(dds_LNAFv_StageA_StageC)
dds_LNAFv_StageA_StageC<-merge(gene_LNAFv,dds_LNAFv_StageA_StageC,by="gene_id")
rownames(dds_LNAFv_StageA_StageC)<-dds_LNAFv_StageA_StageC$gene_id
dds_LNAFv_StageB_StageC$gene_id<-rownames(dds_LNAFv_StageB_StageC)
dds_LNAFv_StageB_StageC<-merge(gene_LNAFv,dds_LNAFv_StageB_StageC,by="gene_id")
rownames(dds_LNAFv_StageB_StageC)<-dds_LNAFv_StageB_StageC$gene_id

dds_LNAFv_StageA_Control[["DifferentialExpression"]] <- ifelse(dds_LNAFv_StageA_Control$padj < 0.01 & abs(dds_LNAFv_StageA_Control$log2FoldChange) > 0.5,
                                                          ifelse(dds_LNAFv_StageA_Control$log2FoldChange > 0, "Up-regulated", "Down-regulated"), "Non-significant")
dds_LNAFv_StageB_Control[["DifferentialExpression"]] <- ifelse(dds_LNAFv_StageB_Control$padj < 0.01 & abs(dds_LNAFv_StageB_Control$log2FoldChange) > 0.5,
                                                               ifelse(dds_LNAFv_StageB_Control$log2FoldChange > 0, "Up-regulated", "Down-regulated"), "Non-significant")
dds_LNAFv_StageC_Control[["DifferentialExpression"]] <- ifelse(dds_LNAFv_StageC_Control$padj < 0.01 & abs(dds_LNAFv_StageC_Control$log2FoldChange) > 0.5,
                                                               ifelse(dds_LNAFv_StageC_Control$log2FoldChange > 0, "Up-regulated", "Down-regulated"), "Non-significant")
dds_LNAFv_StageA_StageB[["DifferentialExpression"]] <- ifelse(dds_LNAFv_StageA_StageB$padj < 0.01 & abs(dds_LNAFv_StageA_StageB$log2FoldChange) > 0.5,
                                                               ifelse(dds_LNAFv_StageA_StageB$log2FoldChange > 0, "Up-regulated", "Down-regulated"), "Non-significant")
dds_LNAFv_StageA_StageC[["DifferentialExpression"]] <- ifelse(dds_LNAFv_StageA_StageC$padj < 0.01 & abs(dds_LNAFv_StageA_StageC$log2FoldChange) > 0.5,
                                                               ifelse(dds_LNAFv_StageA_StageC$log2FoldChange > 0, "Up-regulated", "Down-regulated"), "Non-significant")
dds_LNAFv_StageB_StageC[["DifferentialExpression"]] <- ifelse(dds_LNAFv_StageB_StageC$padj < 0.01 & abs(dds_LNAFv_StageB_StageC$log2FoldChange) > 0.5,
                                                               ifelse(dds_LNAFv_StageB_StageC$log2FoldChange > 0, "Up-regulated", "Down-regulated"), "Non-significant")

write.xlsx(dds_LNAFv_StageA_Control,"RNAseq/result/LNAFv/dds_LNAFv_StageA_Control.xlsx")
write.xlsx(dds_LNAFv_StageB_Control,"RNAseq/result/LNAFv/dds_LNAFv_StageB_Control.xlsx")
write.xlsx(dds_LNAFv_StageC_Control,"RNAseq/result/LNAFv/dds_LNAFv_StageC_Control.xlsx")
write.xlsx(dds_LNAFv_StageA_StageB,"RNAseq/result/LNAFv/dds_LNAFv_StageA_StageB.xlsx")
write.xlsx(dds_LNAFv_StageA_StageC,"RNAseq/result/LNAFv/dds_LNAFv_StageA_StageC.xlsx")
write.xlsx(dds_LNAFv_StageB_StageC,"RNAseq/result/LNAFv/dds_LNAFv_StageB_StageC.xlsx")

# Genes differentially expressed between any two groups were group --------
dds_LNAFv_StageA_Control<-read.xlsx("RNAseq/result/LNAFv/dds_LNAFv_StageA_Control.xlsx")
dds_LNAFv_StageB_Control<-read.xlsx("RNAseq/result/LNAFv/dds_LNAFv_StageB_Control.xlsx")
dds_LNAFv_StageC_Control<-read.xlsx("RNAseq/result/LNAFv/dds_LNAFv_StageC_Control.xlsx")
dds_LNAFv_StageA_StageB<-read.xlsx("RNAseq/result/LNAFv/dds_LNAFv_StageA_StageB.xlsx")
dds_LNAFv_StageA_StageC<-read.xlsx("RNAseq/result/LNAFv/dds_LNAFv_StageA_StageC.xlsx")
dds_LNAFv_StageB_StageC<-read.xlsx("RNAseq/result/LNAFv/dds_LNAFv_StageB_StageC.xlsx")
sample_meta <- as.data.frame(colData(vsd_LNAFv))
get_group_mean <- function(gene_id){
  v <- vst_matrix[gene_id,]
  tibble(
    gene_id = gene_id,
    Ctrl = mean(v[sample_meta$groups=="Control"]),
    A    = mean(v[sample_meta$groups=="StageA"]),
    B    = mean(v[sample_meta$groups=="StageB"]),
    C    = mean(v[sample_meta$groups=="StageC"])
  )
}
gene_list <- rownames(vst_matrix)
group_mean <- bind_rows(lapply(gene_list, get_group_mean))


sig_any <- unique(c(
  dds_LNAFv_StageA_Control$gene_id[!is.na(dds_LNAFv_StageA_Control$padj) & dds_LNAFv_StageA_Control$padj<0.01 & abs(dds_LNAFv_StageA_Control$log2FoldChange)>0.5],
  dds_LNAFv_StageB_Control$gene_id[!is.na(dds_LNAFv_StageB_Control$padj) & dds_LNAFv_StageB_Control$padj<0.01 & abs(dds_LNAFv_StageB_Control$log2FoldChange)>0.5],
  dds_LNAFv_StageC_Control$gene_id[!is.na(dds_LNAFv_StageC_Control$padj) & dds_LNAFv_StageC_Control$padj<0.01 & abs(dds_LNAFv_StageC_Control$log2FoldChange)>0.5],
  dds_LNAFv_StageA_StageB$gene_id[!is.na(dds_LNAFv_StageA_StageB$padj) & dds_LNAFv_StageA_StageB$padj<0.01 & abs(dds_LNAFv_StageA_StageB$log2FoldChange)>0.5],
  dds_LNAFv_StageA_StageC$gene_id[!is.na(dds_LNAFv_StageA_StageC$padj) & dds_LNAFv_StageA_StageC$padj<0.01 & abs(dds_LNAFv_StageA_StageC$log2FoldChange)>0.5],
  dds_LNAFv_StageB_StageC$gene_id[!is.na(dds_LNAFv_StageB_StageC$padj) & dds_LNAFv_StageB_StageC$padj<0.01 & abs(dds_LNAFv_StageB_StageC$log2FoldChange)>0.5]
  
))
length(sig_any)
group_mean1 <- group_mean %>% filter(gene_id %in% sig_any)
group_mean_df <- group_mean1 %>%
  mutate(
    gene_mean = (Ctrl + A + B + C)/4,
    gene_sd   = sqrt( ((Ctrl-gene_mean)^2 + (A-gene_mean)^2 + (B-gene_mean)^2 + (C-gene_mean)^2 ) / 3 ),
    z_ctrl = (Ctrl - gene_mean)/gene_sd,
    z_A    = (A - gene_mean)/gene_sd,
    z_B    = (B - gene_mean)/gene_sd,
    z_C    = (C - gene_mean)/gene_sd
  ) %>%
  mutate(
    Cluster = case_when(
      (z_ctrl > z_A - 0.2) & (z_A > z_B - 0.2) & (z_B > z_C - 0.2) ~ "C1",
      (z_ctrl > z_A - 0.2) & (z_A > z_B - 0.2) & (z_C > z_B - 0.2) ~ "C2",
      (z_A > z_ctrl - 0.2) & (z_A > z_B - 0.2) & (z_B > z_C - 0.2) ~ "C3",
      (z_ctrl < z_A + 0.2) & (z_A < z_B + 0.2) & (z_C < z_B + 0.2) ~ "C4",
      (z_C > z_ctrl - 0.2) & (z_C > z_A - 0.2) & (z_C > z_B - 0.2) ~ "C5"
    )
  )

group_mean_df <- merge(gene_LNAFv, group_mean_df, by = "gene_id")
group_mean_df <- group_mean_df %>% filter(Cluster %in% c("C1","C2","C3","C4","C5")) %>%
  mutate(Cluster = factor(Cluster,levels = c("C1","C2","C3","C4","C5"))) %>%arrange(Cluster)


write.xlsx(group_mean_df, "RNAseq/result/LNAFv/LNAFv_gene_Cluster_class.xlsx")

# S4_B--------------------------------------------------------
group_mean_df<-read.xlsx("RNAseq/result/LNAFv/LNAFv_gene_Cluster_class.xlsx")
gene_order <- pull(group_mean_df, gene_id)
vst_sub <- vst_matrix[gene_order, ]
meta <- as.data.frame(colData(vsd_LNAFv))
meta$groups <- factor(meta$groups, levels = c("Control","StageA","StageB","StageC"))
meta <- meta %>% arrange(groups)
vst_sub <- vst_sub[, rownames(meta)]
n1 <- nrow(filter(group_mean_df, Cluster=="C1"))
n2 <- nrow(filter(group_mean_df, Cluster%in%c("C1","C2")))
n3 <- nrow(filter(group_mean_df, Cluster%in%c("C1","C2","C3")))
n4 <- nrow(filter(group_mean_df, Cluster%in%c("C1","C2","C3","C4")))
gaps_row <- c(n1, n2, n3, n4)
cnt_ctrl  <- sum(meta$groups=="Control")
cnt_A     <- sum(meta$groups=="StageA")
cnt_B     <- sum(meta$groups=="StageB")
gaps_col  <- c(cnt_ctrl, cnt_ctrl+cnt_A, cnt_ctrl+cnt_A+cnt_B)
row_anno_df <- group_mean_df %>% dplyr::select(Cluster)
rownames(row_anno_df) <- gene_order
col_anno_df <- meta %>% dplyr::select(groups)
ann_row_col <- list(
  Cluster = c(
    C1 = "#992288",
    C2 = "#ff77aa",
    C3 = "#ffddcc",
    C4 = "#223377",
    C5 = "#0099aa"
  )
)
ann_col_col <- list(
  groups = c(
    Control = "#b8d8f0",
    StageA = "#2b67b2",
    StageB = "#b0e288",
    StageC = "#379937"
  )
)
ann_col_col <- list(
  groups = c(
    Control = "black",
    StageA = "#EBC678",
    StageB = "#EB7815",
    StageC = "red"
  )
)
mat_scaled <- t(scale(t(vst_sub)))
#key genes
target_ids <- c("ENSMUSG00000052135", "ENSMUSG00000041616", "ENSMUSG00000029019","ENSMUSG00000020988","ENSMUSG00000073609")
target_map <- c("ENSMUSG00000052135"="Foxo6",
                "ENSMUSG00000041616"="Nppa",
                "ENSMUSG00000029019"="Nppb",
                "ENSMUSG00000020988"="L2hgdh",
                "ENSMUSG00000073609"="D2hgdh")


target_pos  <- which(rownames(vst_sub) %in% target_ids)
target_labs <- target_map[rownames(vst_sub)[target_pos]]

cluster_vec <- group_mean_df$Cluster
names(cluster_vec) <- gene_order
cluster_vec <- cluster_vec[rownames(vst_sub)]
ha_row <- rowAnnotation(
  Cluster = cluster_vec,
  col = list(Cluster = c(C1="#992288", C2="#ff77aa", C3="#ffddcc",
                         C4="#223377", C5="#0099aa")),
  show_annotation_name = FALSE
)
ha_col <- HeatmapAnnotation(
  group = meta$groups,
  col = list(group = c(Control="black", StageA="#EBC678",
                       StageB="#EB7815", StageC="red")),
  annotation_name_gp = gpar(fontsize = 11)
)
ha_mark <- rowAnnotation(
  genes = anno_mark(
    at = target_pos,
    labels = target_labs,
    labels_gp = gpar(fontsize = 11),
    link_gp = gpar(col = "black", lwd = 1)
  )
)
S4_B <- Heatmap(
  mat_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, -1, 0, 1, 2),
                   c("blue", "#135DA8", "#F4F6F5", "#CC5533", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = ha_row,
  top_annotation = ha_col,
  right_annotation = ha_mark,
  split = factor(cluster_vec, levels = c("C1","C2","C3","C4","C5")),
  column_split = meta$groups,
  gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  rect_gp = gpar(col = NA),
  heatmap_legend_param = list(direction = "vertical")
)
S4_B




# S4_D--------------------------------------------------------

C1 <- group_mean_df %>% filter(Cluster=="C1") %>% pull(gene_name)
C2 <- group_mean_df %>% filter(Cluster=="C2") %>% pull(gene_name)
C3 <- group_mean_df %>% filter(Cluster=="C3") %>% pull(gene_name)
C4 <- group_mean_df %>% filter(Cluster=="C4") %>% pull(gene_name)
C5 <- group_mean_df %>% filter(Cluster=="C5") %>% pull(gene_name)


LNAFv_C1_GO <- enrichGO(gene = C1, OrgDb = org.Mm.eg.db,keyType = 'SYMBOL', ont = 'BP',pAdjustMethod = 'fdr', pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
LNAFv_C2_GO <- enrichGO(gene = C2, OrgDb = org.Mm.eg.db,keyType = 'SYMBOL', ont = 'BP',pAdjustMethod = 'fdr', pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
LNAFv_C3_GO <- enrichGO(gene = C3, OrgDb = org.Mm.eg.db,keyType = 'SYMBOL', ont = 'BP',pAdjustMethod = 'fdr', pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
LNAFv_C4_GO <- enrichGO(gene = C4, OrgDb = org.Mm.eg.db,keyType = 'SYMBOL', ont = 'BP',pAdjustMethod = 'fdr', pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
LNAFv_C5_GO <- enrichGO(gene = C5, OrgDb = org.Mm.eg.db,keyType = 'SYMBOL', ont = 'BP',pAdjustMethod = 'fdr', pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
LNAFv_C1_GO <- as.data.frame(LNAFv_C1_GO)
LNAFv_C2_GO <- as.data.frame(LNAFv_C2_GO)
LNAFv_C3_GO <- as.data.frame(LNAFv_C3_GO)
LNAFv_C4_GO <- as.data.frame(LNAFv_C4_GO)
LNAFv_C5_GO <- as.data.frame(LNAFv_C5_GO)


write.xlsx(LNAFv_C1_GO, "RNAseq/result/LNAFv/LNAFv_C1_GO.xlsx")
write.xlsx(LNAFv_C2_GO, "RNAseq/result/LNAFv/LNAFv_C2_GO.xlsx")
write.xlsx(LNAFv_C3_GO, "RNAseq/result/LNAFv/LNAFv_C3_GO.xlsx")
write.xlsx(LNAFv_C4_GO, "RNAseq/result/LNAFv/LNAFv_C4_GO.xlsx")
write.xlsx(LNAFv_C5_GO, "RNAseq/result/LNAFv/LNAFv_C5_GO.xlsx")



LNAFv_C1_GO<-read.xlsx( "RNAseq/result/LNAFv/LNAFv_C1_GO.xlsx")
LNAFv_C2_GO<-read.xlsx( "RNAseq/result/LNAFv/LNAFv_C2_GO.xlsx")
LNAFv_C3_GO<-read.xlsx( "RNAseq/result/LNAFv/LNAFv_C3_GO.xlsx")
LNAFv_C4_GO<-read.xlsx( "RNAseq/result/LNAFv/LNAFv_C4_GO.xlsx")
LNAFv_C5_GO<-read.xlsx( "RNAseq/result/LNAFv/LNAFv_C5_GO.xlsx")

energy<-data.frame(
  Cluster = c(rep("C1", 3), rep("C3", 1), rep("C4", 8)),
  Description = c(
    "fatty acid beta-oxidation",
    "cellular respiration",
    "branched-chain amino acid catabolic process",
    "response to insulin",
    "glycolytic process",
    "fatty acid oxidation",
    "lipid storage",
    "regulation of lipid biosynthetic process",
    "glucose transmembrane transport",
    "regulation of gluconeogenesis",
    "regulation of cellular ketone metabolic process",
    "response to insulin"),
  category=c(rep("ENERGY METABOLISM", 12)))

immune<-data.frame(
  Cluster = c(rep("C2", 3), rep("C3", 2),rep("C5", 4)),
  Description = c(
    "positive regulation of lymphocyte activation",
    "tumor necrosis factor production",
    "B cell activation",
    "regulation of myeloid cell differentiation",
    "antigen processing and presentation of peptide antigen via MHC class II",
    "macrophage migration",
    "immune response-activating signaling pathway",
    "regulation of T cell activation",
    "B cell activation"),
  category=c(rep("IMMUNE AND INFLAMMATORY RESPONSE", 9)))

cellular<-data.frame(
  Cluster = c( rep("C3", 6), rep("C4", 4), rep("C5", 3)),
  Description = c(
    
    "response to glucocorticoid",
    "response to oxidative stress",
    "response to hypoxia",
    "response to unfolded protein",
    "response to cAMP",
    "stress-activated MAPK cascade",
    "regulation of autophagy",
    "response to oxidative stress",
    "response to hypoxia",
    "response to unfolded protein",
    "response to hypoxia",
    "response to oxidative stress",
    "response to unfolded protein"),
  category=c(rep("CELLULAR STRESS AND SIGNALING", 13)))

hypertrophy <-data.frame(
  Cluster = c(rep("C3", 3), rep("C4", 4), rep("C5", 5)),
  Description = c(
    "positive regulation of supramolecular fiber organization",
    "regulation of collagen biosynthetic process",
    "cardiac muscle hypertrophy",
    "cellular response to transforming growth factor beta stimulus",
    "fibroblast proliferation",
    "cardiac muscle hypertrophy",
    "positive regulation of supramolecular fiber organization",
    "extracellular matrix organization",
    "cellular response to transforming growth factor beta stimulus",
    "positive regulation of supramolecular fiber organization",
    "cardiac muscle hypertrophy",
    "regulation of collagen biosynthetic process"),
  category=c(rep("HYPERTROPHY AND EXTRACELLULAR REMODELING", 12)))
target_go_df<-bind_rows(energy,
                        immune,
                        cellular,
                        hypertrophy)
write.xlsx(target_go_df, "RNAseq/result/LNAFv/LNAFv_target_go_df.xlsx")

target_go_df<-read.xlsx( "RNAseq/result/LNAFv/LNAFv_target_go_df.xlsx")

LNAFv_C1_GO$Cluster<-"C1"
LNAFv_C1_GO$Cluster<-"C2"
LNAFv_C1_GO$Cluster<-"C3"
LNAFv_C1_GO$Cluster<-"C4"
LNAFv_C1_GO$Cluster<-"C5"

all_go_results <- bind_rows(LNAFv_C1_GO %>% mutate(Cluster="C1"),
                            LNAFv_C2_GO %>% mutate(Cluster="C2"),
                            LNAFv_C3_GO %>% mutate(Cluster="C3"),
                            LNAFv_C4_GO %>% mutate(Cluster="C4"),
                            LNAFv_C5_GO %>% mutate(Cluster="C5"))
all_go <- target_go_df %>%
  left_join(all_go_results, by = c("Cluster", "Description")) %>%
  mutate(
    Description = str_replace(Description, "^(.)", toupper),
    EnrichGeneCount = as.integer(gsub("/.*", "", GeneRatio)),
    TotalGeneCount  = as.integer(gsub(".*/", "", GeneRatio)),
    GeneRatio_num   = EnrichGeneCount / TotalGeneCount,
    Cluster = factor(Cluster, levels = c("C1","C2","C3","C4","C5"), ordered = TRUE),
    log10_P = -log10(p.adjust)
  )%>%
  mutate(
    category = case_when(
      category == "ENERGY METABOLISM" ~ "ENERGY\nMETABOLISM",
      category == "IMMUNE AND INFLAMMATORY RESPONSE" ~ "IMMUNE AND\nINFLAMMATORY\nRESPONSE",
      category == "CELLULAR STRESS AND SIGNALING" ~ "CELLULAR STRESS\nAND SIGNALING",
      category == "HYPERTROPHY AND EXTRACELLULAR REMODELING" ~ "HYPERTROPHY AND\nEXTRACELLULAR\nREMODELING",
      TRUE ~ as.character(category)
    ),
    # 重新设定因子顺序，保持与原顺序一致
    category = factor(category, levels = c(
      "ENERGY\nMETABOLISM",
      "IMMUNE AND\nINFLAMMATORY\nRESPONSE",
      "CELLULAR STRESS\nAND SIGNALING",
      "HYPERTROPHY AND\nEXTRACELLULAR\nREMODELING"
    ))
  )

all_go <- all_go %>%
  arrange(category, Cluster) %>%
  mutate(
    Term_Ordered = factor(Description, levels = rev(unique(Description)), ordered = TRUE)
  )
write.xlsx(all_go, "RNAseq/result/LNAFv/LNAFv_gene_Cluster_class_visualize.xlsx")

cat_colors <- c("#2E8B57", "#CD5C5C", "#9370DB", "#CD853F")
S4_D <- ggplot(all_go, aes(x = Cluster, y = Term_Ordered)) +
  geom_point(aes(color = log10_P, size = GeneRatio_num), alpha = 0.8) +
  scale_color_gradientn(
    colors = c("#E8D789", "#FF9A44", "#D83A25", "#4A0404"),
    values = c(0, 0.1, 0.5, 1),
    limits = c(min(all_go$log10_P), max(all_go$log10_P)),
    name = "-log10(P)"
  ) +
  
  scale_size_continuous(
    range = c(3, 8),
    labels = scales::percent,
    name = "GeneRatio"
  ) +
  guides(
    color = guide_colorbar(
      barwidth = unit(0.4, "cm"), 
      barheight = unit(2, "cm"),  
      frame.colour = "black",  
      frame.linewidth = 1, 
      ticks.colour = "black",
      ticks.linewidth = 0.5)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = " ", y = " ", title = " ") +
  facet_grid2(category ~ ., scales = "free_y", space = "free_y", switch = "y",
              strip = strip_themed(
                background_y = elem_list_rect(fill = cat_colors, color = NA),
                text_y = elem_list_text(colour = "white", face = "bold",
                                        size = 14, hjust = 0.5,
                                        lineheight = 1.2)   # 👈 增加行距，让换行更清晰
              )) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 16, hjust = 1, colour = "black"),
    axis.text.x = element_text(face = "bold", size = 14, colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 10)),
    panel.grid.major = element_line(linetype = "dashed", color = "gray80", size = 0.3),
    panel.grid.minor = element_line(linetype = "dashed", color = "gray90", size = 0.2),
    legend.key.width = unit(1, "cm"),
    legend.spacing.x = unit(1.5, "cm"),
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")
  )
S4_D









# S4_C -------------------------------------------------------------------
line_chart<-read.xlsx("RNAseq/result/LNAFv/LNAFv_gene_Cluster_class_visualize_line_chart.xlsx")
group_mean_df<-read.xlsx("RNAseq/result/LNAFv/LNAFv_gene_Cluster_class.xlsx")
library(tidyverse)
library(ggplot2)
line_chart_split <- line_chart %>%
  separate_rows(geneID, sep = "/") %>%
  rename(gene_name = geneID)
table(line_chart_split$Cluster)
plot_data <- line_chart_split %>%
  inner_join(group_mean_df, by = "gene_name") %>%
  select(Cluster = Cluster.x, gene_name, Ctrl, A, B, C)
plot_data <- plot_data %>%
  distinct(gene_name, Cluster, .keep_all = TRUE)

plot_data_norm <- plot_data %>%
  mutate(
    min_expr   = pmin(Ctrl, A, B, C, na.rm = TRUE),
    max_expr   = pmax(Ctrl, A, B, C, na.rm = TRUE),
    range_expr = max_expr - min_expr,
    n_Ctrl = if_else(range_expr > 0, (Ctrl - min_expr) / range_expr, 0.5),
    n_A    = if_else(range_expr > 0, (A    - min_expr) / range_expr, 0.5),
    n_B    = if_else(range_expr > 0, (B    - min_expr) / range_expr, 0.5),
    n_C    = if_else(range_expr > 0, (C    - min_expr) / range_expr, 0.5)
  )
cluster_stats <- plot_data_norm %>%
  pivot_longer(
    cols = c(n_Ctrl, n_A, n_B, n_C),
    names_to = "Condition",
    values_to = "Expression"
  ) %>%
  mutate(
    Condition = factor(
      case_when(
        Condition == "n_Ctrl" ~ "Control",
        Condition == "n_A"    ~ "StageA",
        Condition == "n_B"    ~ "StageB",
        Condition == "n_C"    ~ "StageC"
      ),
      levels = c("Control", "StageA", "StageB", "StageC")
    ),
    x_num = as.numeric(Condition) 
  ) %>%
  group_by(Cluster, Condition, x_num) %>%
  summarise(
    mean_expr = mean(Expression, na.rm = TRUE),
    q25 = quantile(Expression, 0.25, na.rm = TRUE),
    q75 = quantile(Expression, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cluster_stats$Cluster <- factor(
  cluster_stats$Cluster,
  levels = c("C1", "C2&C5", "C3", "C4", "C3&C4&C5", "C3&C4&C5_1", "C5")
)

S4_C <- ggplot(cluster_stats, aes(x = x_num, y = mean_expr, group = Cluster, color = Cluster)) +
  
  geom_ribbon(
    aes(ymin = q25, ymax = q75, fill = Cluster), 
    alpha = 0.2, 
    color = NA  
  ) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  facet_grid(Cluster ~ ., scales = "free_y") +
  scale_x_continuous(
    breaks = 1:4, 
    labels = c("Control", "StageA", "StageB", "StageC")
  ) +
  scale_color_manual(values = c(
    "C1" = "black",
    "C2&C5" = "#EC75A1",
    "C3" = "#F39C12",
    "C4" = "#0E6EB8",
    "C3&C4&C5" = "#902280",
    "C3&C4&C5_1" = "#E60013",
    "C5" = "gray"
  )) +
  scale_fill_manual(values = c(
    "C1" = "black",
    "C2&C5" = "#EC75A1",
    "C3" = "#F39C12",
    "C4" = "#0E6EB8",
    "C3&C4&C5" = "#902280",
    "C3&C4&C5_1" = "#E60013",
    "C5" = "gray"
  )) +
  
  # 7. 主题设置（和你的原代码完全一致）
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.border = element_rect(color = "grey70", fill = NA,linetype = "dashed", linewidth = 0.5),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold", size = 12, hjust = 1),
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(), 
    #axis.text.y = element_text(color = "black", size = 8), # 可选：设置Y轴数字颜色和大小
    axis.text.x = element_text(color = "black", size = 10),
    legend.position = "none",
    panel.spacing = unit(0, "lines")
  ) +
  labs(x = NULL, y = "Normalized expression (0-1)")
S4_C



# 2.LNAFg Bulk RNAseq：Preprocessing -----------------------------------------------------------
sample_LNAFg<-data.frame(sampleid=c(c("1171", "1174", "630", "631", "80", "81", 
                                      "1155", "171", "168", "174", "167", 
                                      "374", "352", "1608", "1616", "1655")),
                                       groups=c(rep("Control",6),
                                                rep("StageC",5),
                                                rep("Foxo6kd",5)))
sample_LNAFg$sampleid<-as.character(sample_LNAFg$sampleid)
rownames(sample_LNAFg)<-sample_LNAFg$sampleid
count_LNAFg <- read.table('RNAseq/RNAseq.LNAFg.count.txt', header = T, sep = "\t")
gene_LNAFg<-count_LNAFg[,c("gene_id","gene_name")]
gene_LNAFg$gene_name<-str_to_title(gene_LNAFg$gene_name)
names(count_LNAFg) <- gsub("^X", "", names(count_LNAFg))
rownames(count_LNAFg)<-count_LNAFg$gene_id
count_LNAFg<-count_LNAFg[,-c(1:2)]
colnames(count_LNAFg)<-as.character(colnames(count_LNAFg))
count_LNAFg<-count_LNAFg[,intersect(sample_LNAFg$sampleid,colnames(count_LNAFg))]
sample_LNAFg<-sample_LNAFg[intersect(sample_LNAFg$sampleid,colnames(count_LNAFg)),]
dds_LNAFg <- DESeqDataSetFromMatrix(
  countData = count_LNAFg,
  colData = sample_LNAFg,
  design = ~ groups
)
min_count<-30
sample_LNAFg$groups<-factor(sample_LNAFg$groups)
genelist <- DGEList(counts = count_LNAFg, group = relevel(sample_LNAFg$groups, ref = "Control"))
keep <- filterByExpr(genelist, min.count = min_count)
total_gene<-gene_LNAFg[keep,]
dds_LNAFg <- dds_LNAFg[keep,]
count_LNAFg<-count_LNAFg[keep,]

sample_LNAFg_PCA<-sample_LNAFg%>%filter(sample_LNAFg$groups %in% c("Control","StageC"))
count_LNAFg_PCA<-count_LNAFg[,sample_LNAFg_PCA$sampleid]
sample_LNAFg_PCA$groups <- droplevels(sample_LNAFg_PCA$groups)
dds_LNAFg_PCA <- DESeqDataSetFromMatrix(
  countData = count_LNAFg_PCA,
  colData = sample_LNAFg_PCA,
  design = ~ groups
)
min_count<-30
sample_LNAFg_PCA$groups<-factor(sample_LNAFg_PCA$groups)
genelist <- DGEList(counts = count_LNAFg_PCA, group = relevel(sample_LNAFg_PCA$groups, ref = "Control"))
keep <- filterByExpr(genelist, min.count = min_count)#
dds_LNAFg_PCA <- dds_LNAFg_PCA[keep,]
count_LNAFg_PCA<-count_LNAFg_PCA[keep,]
# differentially expressed genes -----------------------------------------
dds_LNAFg <- suppressMessages(DESeq(dds_LNAFg, parallel = FALSE))
dds_LNAFg_Foxo6kd_Control <- results(dds_LNAFg, alpha = 0.05, contrast = c("groups","Foxo6kd","Control")) 
dds_LNAFg_Foxo6kd_Control <- as.data.frame(dds_LNAFg_Foxo6kd_Control)
dds_LNAFg_Foxo6kd_Control <- na.omit(dds_LNAFg_Foxo6kd_Control)
dds_LNAFg_Foxo6kd_Control$gene_id<-rownames(dds_LNAFg_Foxo6kd_Control)
dds_LNAFg_Foxo6kd_Control<-merge(gene_LNAFg,dds_LNAFg_Foxo6kd_Control,by="gene_id")
rownames(dds_LNAFg_Foxo6kd_Control)<-dds_LNAFg_Foxo6kd_Control$gene_id

dds_LNAFg_StageC_Control <- results(dds_LNAFg, alpha = 0.05, contrast = c("groups","StageC","Control")) 
dds_LNAFg_StageC_Control <- as.data.frame(dds_LNAFg_StageC_Control)
dds_LNAFg_StageC_Control <- na.omit(dds_LNAFg_StageC_Control)
dds_LNAFg_StageC_Control$gene_id<-rownames(dds_LNAFg_StageC_Control)
dds_LNAFg_StageC_Control<-merge(gene_LNAFg,dds_LNAFg_StageC_Control,by="gene_id")
rownames(dds_LNAFg_StageC_Control)<-dds_LNAFg_StageC_Control$gene_id



dds_LNAFg_Foxo6kd_Control[["DifferentialExpression"]] <- ifelse(dds_LNAFg_Foxo6kd_Control$padj < 0.01 & abs(dds_LNAFg_Foxo6kd_Control$log2FoldChange) > 0.5,
                                                          ifelse(dds_LNAFg_Foxo6kd_Control$log2FoldChange > 0, "Up-regulated", "Down-regulated"), "Non-significant")
dds_LNAFg_StageC_Control[["DifferentialExpression"]] <- ifelse(dds_LNAFg_StageC_Control$padj < 0.01 & abs(dds_LNAFg_StageC_Control$log2FoldChange) > 0.5,
                                                         ifelse(dds_LNAFg_StageC_Control$log2FoldChange > 0, "Up-regulated", "Down-regulated"), "Non-significant")

write.xlsx(dds_LNAFg_Foxo6kd_Control, "RNAseq/result/LNAFg/dds_LNAFg_Foxo6kd_Control.xlsx")
write.xlsx(dds_LNAFg_StageC_Control, "RNAseq/result/LNAFg/dds_LNAFg_StageC_Control.xlsx")

# S4_A   ---------------------------------------------------------
vsd_LNAFv <- vst(dds_LNAFv, blind = FALSE)
vst_mat_LNAFv <- assay(vsd_LNAFv)
vsd_LNAFg <- vst(dds_LNAFg_PCA, blind = FALSE)
vst_mat_LNAFg <- assay(vsd_LNAFg)
count_LNAFv <- counts(dds_LNAFv)
count_LNAFg <- counts(dds_LNAFg_PCA)
common_genes <- intersect(rownames(count_LNAFv), rownames(count_LNAFg))
count_comb <- cbind(count_LNAFv[common_genes, ], count_LNAFg[common_genes, ])
meta1 <- as.data.frame(colData(dds_LNAFv)) %>%mutate(batch = "Set1")
meta2 <- as.data.frame(colData(dds_LNAFg_PCA)) %>%mutate(batch = "Set2")
meta_comb <- bind_rows(meta1, meta2)
meta_comb <- meta_comb[colnames(count_comb), ]
dds_comb <- DESeqDataSetFromMatrix(
  countData = count_comb,
  colData   = meta_comb,
  design    = ~ batch + groups
)
vst_comb <- vst(dds_comb, blind = FALSE)
expr_vst <- assay(vst_comb)
mod <- model.matrix(~ groups, data = meta_comb)
expr_combat <- ComBat(
  dat = expr_vst,
  batch = meta_comb$batch,
  mod = mod,
  par.prior = TRUE
)
pca_combat <- prcomp(t(expr_combat), scale. = FALSE)
pca_combat_df <- as.data.frame(pca_combat$x[, c("PC1", "PC2")]) %>%
  mutate(
    groups = meta_comb$groups,
    batch  = meta_comb$batch
  )
percent_var_combat <- round(100 * pca_combat$sdev^2 / sum(pca_combat$sdev^2))
S4_A<-ggplot(pca_combat_df, aes(PC1, PC2, color = groups, shape = batch))+
  geom_point(size=6) +
  scale_color_manual(values = c("Control" = "black", 
                                "StageA" = "#EBC678", 
                                "StageB" = "#EB7815", 
                                "StageC" = "#e7141b"))+
  xlab(paste0("PC1: ", percent_var_combat[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var_combat[2], "% variance")) +
  scale_y_continuous(limits = c(-40, 40)) +
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10, color = "black"), # 显示坐标轴数字
        axis.ticks = element_line(color="black", linewidth=0.8), # 显示刻度短线
        axis.title = element_text(size=12))
S4_A







# S14 -------------------------------------------------------------------------
vsd_LNAFg <- vst(dds_LNAFg, blind = TRUE)
vst_matrixrix <- assay(vsd_LNAFg)
sample_Control_Foxo6kd<-sample_LNAFg%>%filter(sample_LNAFg$groups %in% c("Control","Foxo6kd"))
vst_matrixrix<-vst_matrixrix[,sample_Control_Foxo6kd$sampleid]
diff_dds_LNAFg_Foxo6kd_Control<-dds_LNAFg_Foxo6kd_Control%>%filter(!grepl("Non",dds_LNAFg_Foxo6kd_Control$DifferentialExpression))%>%arrange(desc(log2FoldChange))
vst_sub <- vst_matrixrix[diff_dds_LNAFg_Foxo6kd_Control$gene_id,]
rownames(vst_sub)<-diff_dds_LNAFg_Foxo6kd_Control$gene_name
immune_ifn <- c(
  "Ifitm1", "Irf1", "Irf7", "Bst2", "Ly6e",
  "Nmi", "Rtp4", "Tlr2", "Cd14", "Cd40",
  "Nfkbia", "Cxcl9", "Tnfsf10", "Il10ra", "Atpif1"
)
ecm_adhesion <- c(
  "Itgb3", "Icam1", "Edn1",
  "Gpc3", "Emp3"
)
fatty_acid_oxidation <- c(
  "Acsl1",  "Acadm", "Acads", "Acadsb",
  "Hadha", "Hadh", "Acot1", "Acot2"
)
glycolysis <- c(
  "Hk1", "Hk2", "Pdha1"
)
mitochondrial <- c(
  "Tfam", "Hspd1", "Afg3l2", "Pmpca", "Mrps2",
  "Ogdh", "Dlat", "Dld", "Ndufs1", "Sdha"
)
mark_gene <- c(immune_ifn, ecm_adhesion, fatty_acid_oxidation, glycolysis, mitochondrial)
gene_go_color <- data.frame(
  gene = mark_gene,
  text_color = c(
    rep("#E63946",15),
    rep("#457B9D",5),
    rep("#F1A208",8),
    rep("#805B45",3),
    rep("#9B5DE5",10)
  ))

info <- sample_LNAFg[match(colnames(vst_sub), sample_LNAFg$sampleid), ]
stopifnot(identical(info$sampleid, colnames(vst_sub)))
info <- info[order(info$groups), ]
vst_sub <- vst_sub[, info$sampleid]
col_anno_df <- data.frame(
  group = factor(info$groups, levels = c("Control", "Foxo6kd")),
  row.names = colnames(vst_sub)
)
ann_col_col <- list(group = c(Control = "#000000", Foxo6kd = "#1B287C"))
up_genes   <- dds_LNAFg_Foxo6kd_Control$gene_name[dds_LNAFg_Foxo6kd_Control$DifferentialExpression == "Up-regulated"]
down_genes <- dds_LNAFg_Foxo6kd_Control$gene_name[dds_LNAFg_Foxo6kd_Control$DifferentialExpression == "Down-regulated"]
dir_vec <- c(
  rep("up",   length(up_genes)),
  rep("down", length(down_genes))
)
names(dir_vec) <- c(up_genes, down_genes)

up_genes   <- up_genes[up_genes %in% rownames(vst_sub)]
down_genes <- down_genes[down_genes %in% rownames(vst_sub)]

panels_up <- list(
  "Immune / IFN" = immune_ifn,
  "ECM"          = ecm_adhesion
)
panels_down <- list(
  "Fatty acid oxidation" = fatty_acid_oxidation,
  "Glycolysis"           = glycolysis,
  "Mitochondrial"        = mitochondrial
)
all_panel_genes <- unique(c(unlist(panels_up), unlist(panels_down)))

large_gap_up    <- 180
large_gap_down  <- 100
head_gap_up     <- 100
head_gap_down   <- 80

panel_small_gaps <- list(
  "Immune / IFN"         = c(10, 30, 50, 40,20,30,20,50,70,20,30,60,20,30),
  "ECM"                  = c(10, 30, 50, 25),
  "Fatty acid oxidation" = c(50, 30, 40, 30,18,25,35),
  "Glycolysis"           = c(20, 50),
  "Mitochondrial"        = c(20, 30, 50, 30,20,30,40,40,30)
)

split_n <- function(x, n) {
  if (length(x) == 0) return(list())
  n <- min(n, length(x))
  idx <- cut(seq_along(x), breaks = n, labels = FALSE)
  unname(split(x, idx))
}

build_section <- function(panels, direction, other_pool,
                          small_gap_list, large_gap, head_gap) {
  panel_segs <- lapply(names(panels), function(pname) {
    g <- intersect(panels[[pname]], rownames(vst_sub))
    g <- g[dir_vec[g] == direction]
    n_seg <- length(small_gap_list[[pname]]) + 1
    if (length(g) < 2 || n_seg <= 1) return(list(g))
    split_n(g, n_seg)
  })
  names(panel_segs) <- names(panels)
  panel_segs <- panel_segs[sapply(panel_segs, function(s) sum(lengths(s))) > 0]
  total_small <- sum(unlist(small_gap_list[names(panel_segs)]))
  n_large <- max(0, length(panel_segs) - 1)
  need    <- head_gap + total_small + n_large * large_gap
  
  if (need > length(other_pool) && need > 0) {
    ratio <- length(other_pool) / need
    head_gap <- max(1, floor(head_gap * ratio))
    small_gap_list <- lapply(small_gap_list, function(gaps) {
      if (length(gaps) == 0) return(gaps)
      pmax(1, floor(gaps * ratio))
    })
    large_gap <- max(1, floor(large_gap * ratio))
  }
  gene_seq  <- character(0)
  block_seq <- character(0)
  title_seq <- character(0)
  other_idx <- 1
  
  add_block <- function(genes, block_name, title) {
    if (length(genes) == 0) return(invisible(NULL))
    gene_seq  <<- c(gene_seq, genes)
    block_seq <<- c(block_seq, rep(block_name, length(genes)))
    title_seq <<- c(title_seq, rep(title, length(genes)))
  }
  
  if (head_gap > 0 && other_idx <= length(other_pool)) {
    end <- min(other_idx + head_gap - 1, length(other_pool))
    add_block(other_pool[other_idx:end], paste0(direction, "_head"), "")
    other_idx <- end + 1
  }
  
  for (pi in seq_along(panel_segs)) {
    pname <- names(panel_segs)[pi]
    segs  <- panel_segs[[pi]]
    small_gaps <- small_gap_list[[pname]]
    
    for (si in seq_along(segs)) {
      seg <- segs[[si]]
      add_block(seg, paste0(direction, "_", pname, "_seg", si), pname)
      
      if (si <= length(small_gaps) && other_idx <= length(other_pool)) {
        gap_size <- small_gaps[si]
        end <- min(other_idx + gap_size - 1, length(other_pool))
        add_block(other_pool[other_idx:end],
                  paste0(direction, "_small_", pi, "_", si), "")
        other_idx <- end + 1
      }
    }
    
    if (pi < length(panel_segs) && other_idx <= length(other_pool)) {
      end <- min(other_idx + large_gap - 1, length(other_pool))
      add_block(other_pool[other_idx:end],
                paste0(direction, "_large_", pi), "")
      other_idx <- end + 1
    }
  }
  
  if (other_idx <= length(other_pool)) {
    add_block(other_pool[other_idx:length(other_pool)],
              paste0(direction, "_tail"), "")
  }
  
  list(gene_seq  = gene_seq,
       block_seq = block_seq,
       title_seq = title_seq)
}

up_pool   <- up_genes
down_pool <- down_genes
up_other   <- setdiff(up_pool,   all_panel_genes)
down_other <- setdiff(down_pool, all_panel_genes)

up_part   <- build_section(panels_up,   "up",   up_other,
                           panel_small_gaps, large_gap_up,   head_gap_up)
down_part <- build_section(panels_down, "down", down_other,
                           panel_small_gaps, large_gap_down, head_gap_down)


gene_seq   <- c(up_part$gene_seq,   down_part$gene_seq)
block_seq  <- c(up_part$block_seq,  down_part$block_seq)
title_seq  <- c(up_part$title_seq,  down_part$title_seq)

n_up <- length(up_part$gene_seq)
n_down <- length(down_part$gene_seq)
up_down_split <- factor(
  c(rep("Up-regulated", n_up), rep("Down-regulated", n_down)),
  levels = c("Up-regulated", "Down-regulated")
)

mat_show <- t(scale(t(vst_sub[gene_seq, ])))

gene_pos_final    <- which(gene_seq %in% mark_gene)
gene_labels_final <- gene_seq[gene_pos_final]
row_anno_right <- rowAnnotation(
  mark_gene = anno_mark(
    at = gene_pos_final,
    labels = gene_labels_final,
    labels_gp = gpar(col = "black", fontsize = 9),
    extend = unit(3, "mm")
  ),
  show_legend = FALSE
)

S14 <- Heatmap(
  matrix = mat_show,
  name = "Z-score",
  col = colorRamp2(c(-2, -1, 0, 1, 2),
                   c("blue", "#135DA8", "#F4F6F5", "#CC5533", "red")),
  
  cluster_rows = FALSE,
  row_order = seq_len(nrow(mat_show)),
  row_split = up_down_split,       # << 只有上调、下调2个大切片
  row_gap = unit(1, "mm"),         # << 仅此一处gap，上下调中间，自行调大小
  row_title = NULL,
  cluster_row_slices = FALSE,
  
  cluster_columns = FALSE,
  column_split = factor(col_anno_df$group, levels = c("Control", "Foxo6kd")),
  column_gap = unit(1, "mm"),
  cluster_column_slices = FALSE,
  column_title = NULL,
  
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = HeatmapAnnotation(df = col_anno_df, col = ann_col_col,
                                     show_annotation_name = FALSE),
  right_annotation = row_anno_right,
  use_raster = FALSE,
  border = FALSE
)

S14








# -------------------------2.Chipseq:H3k4me3 & H3K27me3---------------------------------------------------------------- -----------------------------------------------------------------
H3K4me3_peak <- read.table("Chipseq/H3K4me3_peak_count.txt", header = FALSE)#shell文件得到的
colnames(H3K4me3_peak) <- c("chr","start","end","name","score","strand",
                            "signalValue","pValue","qValue","rep1_raw","rep2_raw")
width_kb <- (H3K4me3_peak$end - H3K4me3_peak$start)/1000
rep1_total <- as.integer(readLines("Chipseq/H3K4me3_Rep1_total.txt"))
rep2_total <- as.integer(readLines("Chipseq/H3K4me3_Rep2_total.txt"))

H3K4me3_peak$rep1_RPKM <- (H3K4me3_peak$rep1_raw / width_kb) * (1e6 / rep1_total)
H3K4me3_peak$rep2_RPKM <- (H3K4me3_peak$rep2_raw / width_kb) * (1e6 / rep2_total)
H3K4me3_peak$mean_RPKM <- (H3K4me3_peak$rep1_RPKM + H3K4me3_peak$rep2_RPKM)/2
print(paste0("Pearson correlation between replicates RPKM: ", round(cor(H3K4me3_peak$rep1_RPKM, H3K4me3_peak$rep2_RPKM),3)))
threshold<-quantile(H3K4me3_peak$mean_RPKM, probs = c(0.85), na.rm=TRUE)
H3K4me3_peak$strength_group <- cut(H3K4me3_peak$mean_RPKM,
                                   breaks = c(-Inf, threshold, Inf),
                                   labels = c("Low", "High"))
table(H3K4me3_peak$strength_group)
H3K4me3_low<-H3K4me3_peak%>%filter(strength_group %in% c("Low"))
write.table(H3K4me3_low,"Chipseq/result/H3K4me3_low.bed", sep="\t", quote=F, row.names=F, col.names=F)

print(paste0("Pearson correlation between RPKM and signalValue: ", round(cor(H3K4me3_peak$mean_RPKM, H3K4me3_peak$signalValue),3)))
#K4  annotatePeak
H3K4me3 <- readPeakFile("Chipseq/RenLab-H3K4me3-heart-peaks.bed")##下载的bam文件
H3K4me3_peakAnno <- annotatePeak(H3K4me3,
                                 tssRegion = c(-4000, 4000),
                                 TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                                 annoDb = "org.Mm.eg.db")
H3K4me3_peakAnno_gene<-H3K4me3_peakAnno@anno
H3K4me3_peakAnno_gene<-as.data.frame(H3K4me3_peakAnno_gene)
colnames(H3K4me3_peakAnno_gene)[6:11]<-c("name","score","strand1","signalValue","pValue","qValue")
H3K4me3_peakAnno_gene$start <- H3K4me3_peakAnno_gene$start - 1L
merge_df <- H3K4me3_peak %>%
  dplyr::select(chr,start,end,rep1_raw,rep2_raw,width_kb,rep1_RPKM,rep2_RPKM,mean_RPKM,strength_group)
H3K4me3_peakAnno_gene <- merge(H3K4me3_peakAnno_gene,
                               merge_df,
                               by.x = c("seqnames","start","end"),
                               by.y = c("chr","start","end"))
front_cols <- c("seqnames", "start", "end", "strength_group", "annotation","distanceToTSS", "SYMBOL", "geneChr", "geneStart", "geneEnd")
H3K4me3_peakAnno_gene <- H3K4me3_peakAnno_gene %>%
  dplyr::select(all_of(front_cols), everything())
write.xlsx(H3K4me3_peakAnno_gene,"Chipseq/result/H3K4me3_peakAnno_gene.xlsx")

#K27  annotatePeak
H3K27me3 <- readPeakFile('F:/西湖大学/林佳一/文章专用/2026.8西湖大学/方法及代码/Chipseq/data/RenLab-H3K27me3-heart.peaks')
H3K27me3_peakAnno <- annotatePeak(H3K27me3, tssRegion = c(-4000, 4000), TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                                  annoDb = "org.Mm.eg.db")
H3K27me3_peakAnno_gene<-H3K27me3_peakAnno@anno
H3K27me3_peakAnno_gene<-as.data.frame(H3K27me3_peakAnno_gene)
colnames(H3K27me3_peakAnno_gene)[6:11]<-c("name","score","strand1","signalvalue","pvalue","qvalue")
write.xlsx(H3K27me3_peakAnno_gene,"Chipseq/result/H3K27me3_peakAnno_gene.xlsx")

# K4 low peak and K27 peak:bivalent peak ------------------------------------------------------------
bivalent_low <- readPeakFile('Chipseq/overlap_abs200.merged_K4_low_K27.bed')
bivalent_low_peakAnno <- annotatePeak(bivalent_low, tssRegion = c(-4000, 4000), TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                                      annoDb = "org.Mm.eg.db")
bivalent_low_peakAnno_gene<-bivalent_low_peakAnno@anno
bivalent_low_peakAnno_gene<-as.data.frame(bivalent_low_peakAnno_gene)
write.xlsx(bivalent_low_peakAnno_gene,"Chipseq/result/bivalent_low_peakAnno_gene.xlsx")

# F3 -------------------------------------------------------------
dds_LNAFg_StageC_Control<-read.xlsx("RNAseq/result/LNAFg/dds_LNAFg_StageC_Control.xlsx")
Up_StageC_Control<-dds_LNAFg_StageC_Control%>%filter(grepl("Up",DifferentialExpression))%>%pull(gene_name)
Down_StageC_Control<-dds_LNAFg_StageC_Control%>%filter(grepl("Down",DifferentialExpression))%>%pull(gene_name)
total_gene<-dds_LNAFg_StageC_Control$gene_name

H3K27me3_peakAnno_gene<-read.xlsx("Chipseq/result/H3K27me3_peakAnno_gene.xlsx")
H3K27me3_gene<-H3K27me3_peakAnno_gene%>%filter(grepl("Promoter",H3K27me3_peakAnno_gene$annotation))%>%pull(SYMBOL)

H3K4me3_peakAnno_gene<-read.xlsx("Chipseq/result/H3K4me3_peakAnno_gene.xlsx")
H3K4me3_gene<-H3K4me3_peakAnno_gene%>%filter(grepl("Promoter",H3K4me3_peakAnno_gene$annotation))%>%pull(SYMBOL)

H3K4me3_gene_High<-H3K4me3_peakAnno_gene%>%filter(grepl("High",strength_group))%>%filter(grepl("Promoter",annotation))%>%pull(SYMBOL)
H3K4me3_gene_Low<-H3K4me3_peakAnno_gene%>%filter(grepl("Low",strength_group))%>%filter(grepl("Promoter",annotation))%>%pull(SYMBOL)


HK<-read.csv("Chipseq/Housekeeping_GenesMouse.csv")
HK_list<-list()
for(i in 1:nrow(HK)){
  list<-strsplit(HK[i,1],";")
  HK_list[[i]]<-unlist(list)
}
HK_list<-do.call(rbind,HK_list)
HK_gene<-unique(HK_list[,2])


bivalent_low_peakAnno_gene_list<-str_to_title(bivalent_low_peakAnno_gene_list)
bivalent_low_peakAnno_gene_list<-unique(bivalent_low_peakAnno_gene_list)

Up_StageC_Control<-str_to_title(Up_StageC_Control)
Down_StageC_Control<-str_to_title(Down_StageC_Control)
total_gene<-str_to_title(total_gene)
H3K27me3_gene<-str_to_title(H3K27me3_gene)
H3K4me3_gene<-str_to_title(H3K4me3_gene)
H3K4me3_gene_High<-str_to_title(H3K4me3_gene_High)
H3K4me3_gene_Low<-str_to_title(H3K4me3_gene_Low)
HK_gene<-str_to_title(HK_gene)
bivalent_low_peakAnno_gene_list<-str_to_title(bivalent_low_peakAnno_gene_list)

Up_StageC_Control<-unique(Up_StageC_Control)
Down_StageC_Control<-unique(Down_StageC_Control)
total_gene<-unique(total_gene)
H3K27me3_gene<-unique(H3K27me3_gene)
H3K4me3_gene<-unique(H3K4me3_gene)
H3K4me3_gene_High<-unique(H3K4me3_gene_High)
H3K4me3_gene_Low<-unique(H3K4me3_gene_Low)
HK_gene<-unique(HK_gene)

run_permutation <- function(target_genes, diff_genes, gene_set_name, direction,
                            total_genes = total_gene, hk_genes = HK_gene,
                            n_sim = 100000, seed = 123) {
  set.seed(seed)
  observed_overlap <- length(intersect(diff_genes, target_genes))
  sample_size <- length(diff_genes)
  overlap_random <- replicate(n_sim, {
    sum(sample(total_genes, sample_size) %in% target_genes)
  })
  mean_overlap <- mean(overlap_random)
  enrichment_factor <- observed_overlap/ mean_overlap
  hits <- sum(abs(overlap_random - mean_overlap) >= abs(observed_overlap - mean_overlap))
  p_value <- (hits + 1) / (n_sim + 1)
  overlap_hk <- replicate(n_sim, {
    sum(sample(hk_genes, sample_size) %in% target_genes)
  })
  mean_overlap_hk <- mean(overlap_hk)
  format_p <- function(p, thresh = 0.001) {
    if (p < thresh) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3)
  }
  p_label <- format_p(p_value)
  label_txt <- paste(
    paste("me3:", observed_overlap),
    paste("random:", round(mean_overlap, 2)),
    paste("HK:", round(mean_overlap_hk, 2)),
    paste("FC:", round(enrichment_factor, 2)),
    paste("p:", p_label),
    paste("gene:", sample_size),
    sep = "\n"
  )
  plot_df <- data.frame(Overlap_Count = overlap_random)
  p <- ggplot(plot_df, aes(x = Overlap_Count)) +
    geom_histogram(binwidth = 1, fill = "#F2F0F0", color = "#C1C2C7", boundary = 0) +
    geom_density(aes(y = after_stat(count)), fill = "#C8B1A1", color = "#082038",
                 adjust = 2.5, size = 1.1) +
    geom_vline(xintercept = mean_overlap_hk, color = "green", size = 1.2, lineend = "round") +
    geom_vline(xintercept = observed_overlap, color = "#690F26", size = 1.2, lineend = "round") +
    geom_vline(xintercept = mean_overlap, color = "pink", size = 1.2, linetype = "dashed", lineend = "round") +
    labs(
      title = paste(paste(length(target_genes), gene_set_name, "genes"),
                    paste(sample_size, direction, "regulated genes"), sep = " | "),
      x = "Number of overlapping genes", y = "Frequency"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      plot.title = element_text(size = 8)
    ) +
    annotate("text", x = Inf, y = Inf, label = paste("FC:", round(enrichment_factor, 2),
                                                     "\np:", p_label, sep = ""),
             hjust = 1.1, vjust = 1.1, size = 5, fontface = "bold")
  
  
  list(
    observed = observed_overlap,
    mean_random = mean_overlap,
    mean_HK = mean_overlap_hk,
    FC = enrichment_factor,
    p_value = p_value,
    plot = p,
    label = label_txt
  )
}

gene_sets <- list(
  H3K4me3       = H3K4me3_gene,
  H3K4me3_High  = H3K4me3_gene_High,
  H3K4me3_Low   = H3K4me3_gene_Low,
  H3K27me3      = H3K27me3_gene)
diff_sets <- list(
  Up   = Up_StageC_Control,
  Down = Down_StageC_Control
)
results <- list()
for (gs_name in names(gene_sets)) {
  for (dir in names(diff_sets)) {
    cat("Running:", gs_name, dir, "\n")
    res <- run_permutation(
      target_genes = gene_sets[[gs_name]],
      diff_genes   = diff_sets[[dir]],
      gene_set_name = gs_name,
      direction    = dir
    )
    results[[paste(gs_name, dir, sep = "_")]] <- res
  }
}
plot_order <- c("H3K4me3_Up", "H3K4me3_Down",
                "H3K4me3_High_Up", "H3K4me3_High_Down",
                "H3K4me3_Low_Up", "H3K4me3_Low_Down",
                "H3K27me3_Up", "H3K27me3_Down")
F3 <- wrap_plots(
  plots = lapply(results[plot_order], `[[`, "plot"),
  ncol = 2
)
F3
# S8-------------------------------------------------------------
bivalent_low_peakAnno_gene<-read.xlsx("Chipseq/result/bivalent_low_peakAnno_gene.xlsx")
bivalent_low_peakAnno_gene_list<-bivalent_low_peakAnno_gene%>%pull(SYMBOL)
bivalent_low_peakAnno_gene_list<-str_to_title(bivalent_low_peakAnno_gene_list)
bivalent_low_peakAnno_gene_list<-unique(bivalent_low_peakAnno_gene_list)
gene_sets1 <- list(bivalent_low  = bivalent_low_peakAnno_gene_list)
results1 <- list()
for (gs_name in names(gene_sets1)) {
  for (dir in names(diff_sets)) {
    cat("Running:", gs_name, dir, "\n")
    res <- run_permutation(
      target_genes = gene_sets1[[gs_name]],
      diff_genes   = diff_sets[[dir]],
      gene_set_name = gs_name,
      direction    = dir
    )
    results1[[paste(gs_name, dir, sep = "_")]] <- res
  }
}
plot_order <- c("bivalent_low_Up", "bivalent_low_Down")
combined_K4_K27_bivalent_plot <- wrap_plots(
  plots = lapply(results1[plot_order], `[[`, "plot"),
  ncol = 2
)
combined_K4_K27_bivalent_plot
# S8---------------------------------------------------------------------
vst_mat <- vst(dds_LNAFg, blind = TRUE)
vst_df <- assay(vst_mat)
Control_samples <- sample_LNAFg %>% filter(groups == "Control") %>% pull(sampleid)
Control_vst <- vst_df[, Control_samples, drop = FALSE]
gene_Control_mean <- rowMeans(Control_vst) # 每个基因Control平均VST表达
gene_anno <- dds_LNAFg_StageC_Control %>%
  dplyr::select(gene_id, gene_name, log2FoldChange, padj, DifferentialExpression) %>%
  mutate(
    Control_mean = gene_Control_mean[gene_id],
    enriched_K4 = ifelse(gene_name %in% H3K4me3_gene, "K4_pos", "K4_neg"),
    enriched_K27 = ifelse(gene_name %in% H3K27me3_gene, "K27_pos", "K27_neg"),
  ) %>%
  drop_na(Control_mean)
gene_anno_sort <- gene_anno %>% arrange(Control_mean)
n_bin <- 100
gene_anno_sort$bin <- cut(1:nrow(gene_anno_sort), breaks = n_bin, labels = 1:n_bin)
bin_K4 <- gene_anno_sort %>%group_by(bin) %>%summarise(mean_expr = mean(Control_mean),K4_frac = mean(enriched_K4 == "K4_pos"),n_gene = n(),.groups = "drop")
plot_k4<-ggplot(bin_K4, aes(x = mean_expr, y = K4_frac)) +geom_point(size=2.5) +geom_smooth(method = "lm", color="#E63946", se = TRUE) +labs(x = "Expression in Control",y = "Fraction of genes with H3K4me3",title = "H3K4me3 occupancy increases with Control expression") +theme_classic()
bin_K27 <- gene_anno_sort %>%group_by(bin) %>%summarise(mean_expr = mean(Control_mean),K27_frac = mean(enriched_K27 == "K27_pos"),n_gene = n(),.groups = "drop")
plot_k27<-ggplot(bin_K27, aes(x = mean_expr, y = K27_frac)) +geom_point(size=2.5) +geom_smooth(method = "lm", color="#E63946", se = TRUE) +labs(x = "Expression in Control",y = "Fraction of genes with H3K27me3",title = "H3K27me3 occupancy decreases with Control expression") +theme_classic()
S8 <- plot_k4 + plot_k27
S8

# -------------------------3.ATACseq ---------------------------------------------------------------- ----------------------------------------------------------------
dbObj1 <- dba(sampleSheet = "ATACseq/samplesheet.csv")
dbObj <- dba.count(dbObj1, bUseSummarizeOverlaps=TRUE,summits=FALSE)
save(dbObj,file="ATACseq/dbObj.RData")

load("E:/上传代码/ATACseq/dbObj.RData")
dbObj$samples$Treatment <- factor(dbObj$samples$Treatment, 
                                  levels = c("Control", "StageC", "Foxo6"))
Control_idx <- which(dbObj$samples$Treatment == "Control")
StageC_idx <- which(dbObj$samples$Treatment == "StageC")
Foxo6_idx <- which(dbObj$samples$Treatment == "Foxo6")


dbObj <- dba.contrast(dbObj,group1 = StageC_idx,group2 =Control_idx,name1 = "StageC",name2 = "Control")
dbObj <- dba.contrast(dbObj, group1 = Foxo6_idx,  group2 = StageC_idx, name1="Foxo6",  name2="StageC")
dbObj <- dba.contrast(dbObj, group1 = Foxo6_idx,  group2 = Control_idx, name1="Foxo6",  name2="Control")
dba.show(dbObj, bContrasts = TRUE)
dbObj <- dba.analyze(dbObj, method = DBA_ALL_METHODS)
save(dbObj,file="ATACseq/dbObj_diffbind.RData")
load("ATACseq/dbObj_diffbind.RData")
StageC_Control <- dba.report(dbObj,  method=DBA_EDGER, contrast = 1, th=1)
Foxo6_StageC <- dba.report(dbObj,  method=DBA_EDGER, contrast = 2, th=1)
Foxo6_Control <- dba.report(dbObj,  method=DBA_EDGER, contrast = 3, th=1)

# ATACseq diffbind peak ---------------------------------------------------
StageC_Control.peak <- as.data.frame(StageC_Control)
colnames(StageC_Control.peak)[which(colnames(StageC_Control.peak)=="Fold")]<-"log2FC"
StageC_Control.peak$id<-paste0("chr",StageC_Control.peak$seqnames, ":",
                               StageC_Control.peak$start, "-",
                               StageC_Control.peak$end)


StageC_Control.peak <- StageC_Control.peak %>%
  mutate(type = case_when(
    log2FC > 0.5 & FDR < 0.05 ~ "up",
    log2FC < -0.5 & FDR < 0.05 ~ "down",
    TRUE ~ "nodiff"
  ))
table(StageC_Control.peak$type)
StageC_Control.peak$seqnames<-paste("chr",StageC_Control.peak$seqnames,sep="")
write.table(StageC_Control.peak, file="ATACseq/result/StageC_Control_ATAC.bed", sep="\t", quote=F, row.names=F, col.names=F)
#整理差异ATAC-peak:Foxo6_StageC
Foxo6_StageC.peak <- as.data.frame(Foxo6_StageC)
colnames(Foxo6_StageC.peak)[which(colnames(Foxo6_StageC.peak)=="Fold")]<-"log2FC"
Foxo6_StageC.peak$id<-paste0("chr",Foxo6_StageC.peak$seqnames, ":",
                             Foxo6_StageC.peak$start, "-",
                             Foxo6_StageC.peak$end)


Foxo6_StageC.peak <- Foxo6_StageC.peak %>%
  mutate(type = case_when(
    log2FC > 0.5 & FDR < 0.05 ~ "up",
    log2FC < -0.5 & FDR < 0.05 ~ "down",
    TRUE ~ "nodiff"
  ))
table(Foxo6_StageC.peak$type)
Foxo6_StageC.peak$seqnames<-paste("chr",Foxo6_StageC.peak$seqnames,sep="")


write.table(Foxo6_StageC.peak, file="ATACseq/result/Foxo6_StageC_ATAC.bed", sep="\t", quote=F, row.names=F, col.names=F)

# annotatePeak -------------------------------------------------------------------



StageC_Control.peak <- readPeakFile('ATACseq/result/StageC_Control_ATAC.bed')
StageC_Control_peakAnno <- annotatePeak(StageC_Control.peak, tssRegion = c(-4000, 4000), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                                        annoDb = "org.Mm.eg.db")

StageC_Control_peakAnno_gene<-StageC_Control_peakAnno@anno
StageC_Control_peakAnno_gene<-as.data.frame(StageC_Control_peakAnno_gene)
colnames(StageC_Control_peakAnno_gene)[6:15]<-c("width1","strand1","Conc","Conc_StageC","Conc_Control","log2FC","p.value","FDR","id","type")
StageC_Control_peakAnno_gene$id<-paste0(StageC_Control_peakAnno_gene$seqnames, ":",
                                        StageC_Control_peakAnno_gene$start, "-",
                                        StageC_Control_peakAnno_gene$end)
StageC_Control_peakAnno_gene$SYMBOL<-gsub("Atp5if1","Atpif1",StageC_Control_peakAnno_gene$SYMBOL)
StageC_Control_peakAnno_gene_4kb<-StageC_Control_peakAnno_gene%>%filter(grepl("Promoter", StageC_Control_peakAnno_gene$annotation))
sheets = list("4kb" = StageC_Control_peakAnno_gene_4kb,"ALL" = StageC_Control_peakAnno_gene)
write.xlsx(sheets,"ATACseq/result/StageC_Control_peakAnnotation.xlsx")


#差异ATAC-peak注释到基因:Foxo6_StageC


Foxo6_StageC.peak <- readPeakFile('ATACseq/result/Foxo6_StageC_ATAC.bed')
Foxo6_StageC_peakAnno <- annotatePeak(Foxo6_StageC.peak, tssRegion = c(-4000, 4000), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                                      annoDb = "org.Mm.eg.db")

Foxo6_StageC_peakAnno_gene<-Foxo6_StageC_peakAnno@anno
Foxo6_StageC_peakAnno_gene<-as.data.frame(Foxo6_StageC_peakAnno_gene)
colnames(Foxo6_StageC_peakAnno_gene)[6:15]<-c("width1","strand1","Conc","Conc_Foxo6","Conc_StageC","log2FC","p.value","FDR","id","type")
Foxo6_StageC_peakAnno_gene$id<-paste0(Foxo6_StageC_peakAnno_gene$seqnames, ":",
                                      Foxo6_StageC_peakAnno_gene$start, "-",
                                      Foxo6_StageC_peakAnno_gene$end)
Foxo6_StageC_peakAnno_gene$SYMBOL<-gsub("Atp5if1","Atpif1",Foxo6_StageC_peakAnno_gene$SYMBOL)
Foxo6_StageC_peakAnno_gene_4kb<-Foxo6_StageC_peakAnno_gene%>%filter(grepl("Promoter", Foxo6_StageC_peakAnno_gene$annotation))
sheets = list("4kb" = Foxo6_StageC_peakAnno_gene_4kb,"ALL" = Foxo6_StageC_peakAnno_gene)
write.xlsx(sheets,"ATACseq/result/Foxo6_StageC_peakAnnotation.xlsx")
# Table_S4-------------------------------------------------------------------

dds_LNAFg_StageC_Control<-read.xlsx("RNAseq/result/LNAFg/dds_LNAFg_StageC_Control.xlsx")
dds_LNAFg_StageC_Control_up<-dds_LNAFg_StageC_Control%>%filter(dds_LNAFg_StageC_Control$DifferentialExpression=="Up-regulated") %>%pull(gene_name) 
dds_LNAFg_StageC_Control_down<-dds_LNAFg_StageC_Control%>%filter(dds_LNAFg_StageC_Control$DifferentialExpression=="Down-regulated") %>%pull(gene_name) 

H3K4me3_peakAnno_gene<-read.xlsx("Chipseq/result/H3K4me3_peakAnno_gene.xlsx")
H3K27me3_peakAnno_gene<-read.xlsx("Chipseq/result/H3K27me3_peakAnno_gene.xlsx")

H3K27me3_positive <- H3K27me3_peakAnno_gene %>%filter(grepl("Promoter", annotation)) %>%pull(SYMBOL) %>%unique()
H3K4me3_high <- H3K4me3_peakAnno_gene %>%filter(grepl("Promoter", annotation)) %>%filter(grepl("High", strength_group))%>%pull(SYMBOL) %>%unique()
H3K4me3_low <- H3K4me3_peakAnno_gene %>%filter(grepl("Promoter", annotation)) %>%filter(grepl("Low", strength_group))%>%pull(SYMBOL) %>%unique()

StageC_Control_ATAC_peakAnno_gene_4kb<-read.xlsx("ATACseq/result/StageC_Control_peakAnnotation.xlsx",sheet="4kb")
up_ATAC<-StageC_Control_ATAC_peakAnno_gene_4kb%>%filter(type %in% c("up"))%>%pull(SYMBOL) %>%unique()
down_ATAC<-StageC_Control_ATAC_peakAnno_gene_4kb%>%filter(type %in% c("down"))%>%pull(SYMBOL) %>%unique()

#K4 low & LNAFg-StageC up & K27 positive & ATAC up
up_summary<-sort(intersect(intersect(intersect(H3K4me3_low,dds_LNAFg_StageC_Control_up),H3K27me3_positive),up_ATAC))
#K4 low & LNAFg-StageC down & K27 positive & ATAC down
down_summary<-sort(intersect(intersect(intersect(H3K4me3_low,dds_LNAFg_StageC_Control_down),H3K27me3_positive),down_ATAC))
##summary genes
gene_summary<-c(up_summary,down_summary)
RNA<-dds_LNAFg_StageC_Control[,c("gene_name","log2FoldChange","padj","DifferentialExpression")]
K4<-H3K4me3_peakAnno_gene[,c("SYMBOL","strength_group","annotation")]
K4<-K4%>%filter(grepl("Promoter", annotation)) 
K27<-H3K27me3_peakAnno_gene[,c("SYMBOL","annotation")]
K27<-K27%>%filter(grepl("Promoter", annotation)) 
ATAC<-StageC_Control_ATAC_peakAnno_gene_4kb[,c("SYMBOL","log2FC","FDR","type","annotation")]
ATAC<-ATAC%>%filter(grepl("Promoter", annotation))%>%filter(!grepl("nodiff", type)) 
colnames(RNA)<-paste("RNA",colnames(RNA),sep = "_")
colnames(K4)<-paste("K4",colnames(K4),sep = "_")
colnames(K27)<-paste("K27",colnames(K27),sep = "_")
colnames(ATAC)<-paste("ATAC",colnames(ATAC),sep = "_")

RNA_up<-RNA%>%filter(RNA_gene_name %in% up_summary)
K4_up<-K4%>%filter(K4_SYMBOL %in% up_summary)
K27_up<-K27%>%filter(K27_SYMBOL %in% up_summary)
ATAC_up<-ATAC%>%filter(ATAC_type %in% c("up"))%>%filter(ATAC_SYMBOL %in% up_summary)
ATAC_up <- ATAC_up %>%group_by(ATAC_SYMBOL) %>%slice_min(ATAC_FDR, n = 1, with_ties = FALSE) %>%ungroup()
ATAC_up<-as.data.frame(ATAC_up)


up_summary<-merge(ATAC_up,RNA_up,by.x="ATAC_SYMBOL",by.y="RNA_gene_name")
up_summary<-up_summary[,c("ATAC_SYMBOL","ATAC_log2FC","ATAC_FDR","RNA_log2FoldChange","RNA_padj","RNA_DifferentialExpression")]
colnames(up_summary)[c(1,6)]<-c("SYMBOL","Regulated")
up_summary <- up_summary %>% arrange(ATAC_FDR)



RNA_down<-RNA%>%filter(RNA_gene_name %in% down_summary)
K4_down<-K4%>%filter(K4_SYMBOL %in% down_summary)
K27_down<-K27%>%filter(K27_SYMBOL %in% down_summary)
ATAC_down<-ATAC%>%filter(ATAC_type %in% c("down"))%>%filter(ATAC_SYMBOL %in% down_summary)
ATAC_down <- ATAC_down %>%group_by(ATAC_SYMBOL) %>%slice_min(ATAC_FDR, n = 1, with_ties = FALSE) %>%ungroup()
ATAC_down<-as.data.frame(ATAC_down)


down_summary<-merge(ATAC_down,RNA_down,by.x="ATAC_SYMBOL",by.y="RNA_gene_name")
down_summary<-down_summary[,c("ATAC_SYMBOL","ATAC_log2FC","ATAC_FDR","RNA_log2FoldChange","RNA_padj","RNA_DifferentialExpression")]
colnames(down_summary)[c(1,6)]<-c("SYMBOL","Regulated")
down_summary <- down_summary %>% arrange(ATAC_FDR)


Table_S4<-rbind(down_summary,up_summary)




# -------------------------4.metabolomics ---------------------------------------------------------------- -------------------------------------------------------------------
metabolomics <- read.xlsx("metabolomics/metabolomics_data.xlsx")
metabolomics<-na.omit(metabolomics)
rownames(metabolomics) <-metabolomics$Metabolite
metabolomics <- metabolomics[, !colnames(metabolomics) %in% c("Metabolite")]
metabolomics<-log2(metabolomics)
all_col <- colnames(metabolomics)
sample_cols <- all_col[str_detect(all_col, "^Norm|^HFSC")]
met_anno <- rownames(metabolomics)
metabolomics <- metabolomics[, sample_cols]
sample_info <- data.frame(
  sample = sample_cols,
  group = case_when(
    str_detect(sample_cols, "^Norm") ~ "Control",
    str_detect(sample_cols, "^HFSC") ~ "StageC"
  ))
rownames(sample_info)<-sample_info$sample
pca_obj <- opls(t(metabolomics), predI = 2, algoC = "nipals",scaleC = "standard")
pca_score <- getScoreMN(pca_obj) %>% as.data.frame()
pca_score$sample <- rownames(pca_score)
pca_score <- left_join(pca_score, sample_info, by = "sample")
PCA_plot<-ggplot(pca_score, aes(x = p1, y = p2, color = group)) +
  geom_point(size = 5) +
  labs(x = paste0("PC1 (", round(pca_obj@modelDF$R2X[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(pca_obj@modelDF$R2X[2]*100, 1), "%)")) +
  theme_bw() +
  theme(panel.grid   = element_blank())+
  scale_color_manual(values = c("Control" = "black", "StageC" = "red"))+
  geom_text_repel(aes(label = sample),
                  size = 3, show.legend = FALSE,
                  box.padding = 0.4, point.padding = 0.3,
                  max.overlaps = Inf,
                  segment.color = "grey50")
metabolomics <- metabolomics[, !colnames(metabolomics) %in% c("HFSC4")]
sample_info<-sample_info[!rownames(sample_info) %in% c("HFSC4"),]
oplsda_obj <- opls(t(metabolomics), 
                   sample_info$group, 
                   predI = 1, 
                   algoC = "nipals",
                   permI = 200,
                   scaleC = "pareto")
print(oplsda_obj@summaryDF)
vip_df <- getVipVn(oplsda_obj) %>% as.data.frame()
colnames(vip_df) <- "VIP"
vip_df$Metabolite <- rownames(vip_df)
samples <- colnames(metabolomics)
info    <- sample_info[match(samples, sample_info$sample), ] 
design <- model.matrix(~ 0 + group, data = info)
colnames(design) <- c("Control", "StageC")
fit  <- lmFit(metabolomics, design)
contr <- makeContrasts(StageC - Control, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
limma_res <- topTable(fit2, number = Inf, sort.by = "none")
limma_res$Metabolite<-rownames(limma_res)
res_all <- limma_res %>%
  left_join(vip_df, by = "Metabolite") 
res_all <- res_all %>%
  mutate(
    sig_Metabolite = case_when(
      adj.P.Val < 0.05 & VIP > 1 & logFC >  0.58 ~ "up",
      adj.P.Val < 0.05 & VIP > 1 & logFC < -0.58 ~ "down",
      TRUE ~ "not_sig"
    )
  )
table(res_all$sig_Metabolite)
metabolomics$Metabolite<-rownames(metabolomics)
res_limma <- merge(metabolomics,res_all,by.x="Metabolite",by.y="Metabolite")
write.xlsx(res_limma, "metabolomics/result/metabolomics_diff_result.xlsx", rowNames = FALSE)
# Figure2A -----------------------------------------------------
res_limma<-read.xlsx("metabolomics/result/metabolomics_diff_result.xlsx")
rownames(res_limma)<-res_limma$Metabolite
res_limma_diff <- res_limma[,sample_info$sample]
res_heat_met<-res_limma_diff[c("2-Hydroxyglutaric acid","9E-tetradecenoic acid","5Z-Dodecenoic acid","Butyric acid","Octanoic acid",
                               "Decanoic acid","Glyceric acid","Glutaconic acid","2-Methylbutyroylcarnitine","Isoleucine","Phenylalanine",
                               "Homocitrulline","Hippuric acid","1-Methylhistidine","Oxoglutaric acid","Pyroglutamic acid","Methylcysteine","Carnosine"),]
tab_col <- table(sample_info$group)
tab_col <- tab_col[tab_col > 0]
gaps_col <- cumsum(tab_col)
gaps_col <- gaps_col[-length(gaps_col)]
if (length(gaps_col) == 0) gaps_col <- NULL
ann_col_col <- list(group = c( Control = "black",StageC = "red"))
F2_A<-ComplexHeatmap::pheatmap(res_heat_met,
                         scale = 'row',       
                         cluster_rows = FALSE, 
                         cluster_cols = FALSE,
                         annotation_col = sample_info["group"],
                         show_rownames = TRUE,
                         annotation_colors = ann_col_col,
                         color = colorRampPalette(c("#004488", "#4488BB", "#F4F6F5", "#CC5533", "#942003"))(100),
                         border_color = NA,
                         show_colnames = FALSE,
                         heatmap_legend_param = list(title = "Zscore"),fontsize = 12)
F2_A



# -------------------------5.metabolome --------------------------------------------------------------------
metabolome <- read.xlsx("metabolomics/FINAL_metabolite_data_updated 5 group with L2hg D2hg update.xlsx",sheet="Converted data")
sample_info_all<-data.frame(sample=metabolome$Sample,group=metabolome$Group)
all_sample<-sample_info_all$sample
sample_info_all$group<-gsub("F6","LNAFg-FoxO6",sample_info_all$group) 
sample_info_all$group <- factor(
  sample_info_all$group,
  levels = c("Control", "LNAFg",  "LNAFg-FoxO6","LNAFg-L2hgdh", "LNAFg-D2hgdh")
)
sample_info_all <- sample_info_all %>% arrange(group)
rownames(sample_info_all) <- sample_info_all$sample   # 可选：重置行号
metabolome<-metabolome[,-c(1:2)]
metabolome_num <- apply(metabolome, 2, as.numeric)
rownames(metabolome_num)<-all_sample
metabolome_num<-as.data.frame(metabolome_num)

metabolome<-t(metabolome_num)
metabolome_ratio<-metabolome[grep("/", rownames(metabolome)),]
metabolome_ratio<-log2(metabolome_ratio)
metabolome<-metabolome[-grep("/", rownames(metabolome)),]
metabolome2<-rbind(metabolome,metabolome_ratio)
metabolome2_sub <- metabolome2[,sample_info_all$sample]

# Figure4--------------------------------------------------------------------
group_colors <- c("Control"="black", "LNAFg"="red", "LNAFg-FoxO6"="#F1AA1B" , 
                  "LNAFg-L2hgdh"="#4858A6","LNAFg-D2hgdh"="#A30C5E")

annot_col <- data.frame(
  group = sample_info_all$group[match(colnames(metabolome2_sub), sample_info_all$sample)],
  row.names = colnames(metabolome2_sub)
)
annot_col$group <- factor(annot_col$group, levels =c("Control", "LNAFg","LNAFg-FoxO6","LNAFg-L2hgdh",  "LNAFg-D2hgdh"))

ann_colors <- list(group = group_colors)

cnt_ctrl  <- sum(sample_info_all$group=="Control")
cnt_HTAF     <- sum(sample_info_all$group=="LNAFg")
cnt_FoxO6     <- sum(sample_info_all$group=="LNAFg-FoxO6")
cnt_L2hgdh     <- sum(sample_info_all$group=="LNAFg-L2hgdh")
cnt_D2hgdh    <- sum(sample_info_all$group=="LNAFg-D2hgdh")

gaps_col  <- c(cnt_ctrl, cnt_ctrl+cnt_HTAF, cnt_ctrl+cnt_HTAF+cnt_FoxO6, cnt_ctrl+cnt_HTAF+cnt_FoxO6+cnt_L2hgdh)
gaps_row <- sum(!grepl("/", rownames(metabolome2_sub)))
F4<-pheatmap(
  mat = metabolome2_sub,
  annotation_col = annot_col,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("#004488", "#4488BB", "#F4F6F5", "#CC5533", "#942003"))(50),
  gaps_col = gaps_col,
  gaps_row = gaps_row,
  border_color = NA,
  fontsize = 10,  
  fontsize_row = 12,  
  fontsize_col = 10, 
  heatmap_legend_param = list(title = "Z-score", title_position = "topcenter")
)
F4
