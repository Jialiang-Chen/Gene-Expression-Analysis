# https://bioconductor.org/packages/3.11/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data

# 0. preparation 
# load library
library(tidyr)
library(dplyr)
library(stringr)
library(DESeq2)
library(clusterProfiler)
library(ggplot2)
library(tibble)
library(ggrepel)
library(wesanderson)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(survival)
library(survminer)
library(contsurvplot)
library(pammtools)
# library(pathviewr)

citation("DESeq2")
citation("clusterProfiler")
citation("survival")
citation("survminer")

# set directory
path <- "/Users/chenjialiang/Desktop/MSc AI for Medicine&Medical Research/ANAT40040-Bio Principles & Cellular Org/Assignment 2"
setwd(path)

# untar folder.
folder_name <- "brca_tcga_pan_can_atlas_2018.tar.gz"
folder <- paste(path, folder_name, sep = "/")
untar(folder)

# go to new path
new_dir <- paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )
setwd(new_dir)

# set cut off
alpha <- 0.05
padj.cutoff <- 0.05
lfc.cutoff <- 0 # 0.58 has also been tested, considering relatively small gene count, the threshold is removed


# 1. read files and match patient ids
# read the RNASeq file
data_Rnaseq <- read.delim("data_mrna_seq_v2_rsem.txt")

# read the Patient Data file
data_patient <- read.delim("data_clinical_patient.txt", skip = 4, header = TRUE)

# read the Copy Number Aberrations Data file
data_cna <- read.delim("data_cna.txt")

# plot histogram to visualize explore the data
hist(as.numeric(data_cna[data_cna$Hugo_Symbol == "ERBB2",-c(1,2)]))

# match the RNASeq patient ids with the CNA ids and the Patient Data ids
for (i in 3:dim(data_Rnaseq)[2]){
  pat_barcode <- colnames(data_Rnaseq)[i] 
  pat_barcode = substr(pat_barcode, 1, 12)
  pat_barcode = gsub("\\.", "-",pat_barcode)
  colnames(data_Rnaseq)[i] <- pat_barcode
}

for (i in 3:dim(data_cna)[2]){
  pat_barcode <- colnames(data_cna)[i] 
  pat_barcode <- substr(pat_barcode, 1, 12)
  pat_barcode <- gsub("\\.", "-",pat_barcode)
  colnames(data_cna)[i] <- pat_barcode
}

# match patients in rnaseq to patients in cna
rna_cna_id <- which(is.element(colnames(data_Rnaseq[, ]), colnames(data_cna[, ])))

# select only the rna cases which have cna data
data_Rnaseq <- data_Rnaseq[, rna_cna_id]

# remove duplicated genes
keep <- !duplicated(data_Rnaseq[, 1])
data_Rnaseq <- data_Rnaseq[keep, ]

# remove unknown genes
data_Rnaseq <- filter(data_Rnaseq, Hugo_Symbol != "")

# 2. create metadata
# CNA level of ERBB2
erbb2 <- data_cna |>
  filter(Hugo_Symbol == "ERBB2") |>
  pivot_longer(cols = colnames(data_cna)[3]:colnames(data_cna)[dim(data_cna)[2]]
               , names_to = "PATIENT_ID", values_to = "ERBB2_CNA_LEVEL")

# build assay
# assay <- as.matrix(data_Rnaseq[,-c(2)]) |>
#   data.frame() |> 
#   filter(Hugo_Symbol != "") |>
#   column_to_rownames(var = 'Hugo_Symbol')
assay <- tibble(data_Rnaseq[,-2])
colnames(assay)[1] <- "Gene"

assay <- assay |>
  column_to_rownames(var = 'Gene')

# build metadata
metadata <- matrix(NA, dim(assay)[2], 2)
pat_ids <- erbb2$PATIENT_ID
col_cna <- which(colnames(erbb2) == "ERBB2_CNA_LEVEL")

for (i in 1:dim(assay)[2]){
  idx <- which(colnames(assay)[i] == pat_ids)
  metadata[i, 1] <- colnames(assay)[i]
  metadata[i, 2] <- as.numeric(erbb2[idx, 4]) > 0 # greater then 0 means amplified, i.e., ERBB2+
}
metadata[is.na(metadata)] <- FALSE

colnames(metadata) <- c("PATIENT_ID","Subtype")

metadata[, 2] <- str_replace_all(metadata[, 2], "TRUE", "ERBB2 Amplified")
metadata[, 2] <- str_replace_all(metadata[, 2], "FALSE", "Others")

metadata <- metadata |>
  data.frame() |>
  column_to_rownames(var = 'PATIENT_ID')
all(colnames(data_Rnaseq[, -c(1, 2)]) %in% rownames(metadata)) # check if PATIENT_ID matches

# check number of ERBB2 amplified vs others
table(metadata)

# build DESeq object
assay[is.na(assay)] <- 0  # impute with 0 the NA
assay[assay < 0] <- 0 # correct negative value to 0

dds <- DESeqDataSetFromMatrix(countData = round(assay),
                              colData = metadata,
                              design = ~ Subtype)
dds
dds_size_factor <- estimateSizeFactors(dds)
dds_size_factor

# filter out low count genes
# smallestGroupSize <- 3
# keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize # keep genes with at least "smallestGroupSize" reads exceeding 10
# dds <- dds[keep,]

# 3. differential expression analysis
# set factor level
dds$Subtype <- factor(dds$Subtype, levels = c("ERBB2 Amplified","Others"))

# normalize data
dds <- DESeq(dds)

normalized_counts <- counts(dds_size_factor, normalized = TRUE)
normalized_counts <- normalized_counts |>
  data.frame() |>
  rownames_to_column(var = "Gene") |>
  as_tibble()

# plot normalized count of ERBB2
erbb2_plt <- plotCounts(dds, gene = "ERBB2", intgroup = "Subtype", returnData = TRUE)
ggplot(erbb2_plt, aes(x = Subtype, y = count, color = Subtype)) +
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  # geom_text_repel(aes(label = rownames(erbb2_plt))) +
  # theme_bw() +
  theme(legend.position="bottom") +
  labs(x = "", y = "Normalized Count", color = "Subtype"
       , title = "ERBB2 Normalized Count", subtitle = "ERBB2 Amplified vs Others") +
  scale_color_manual(values = c("#F8AFA8", "#74A089"))

# get result
res <- results(dds, contrast = c("Subtype","ERBB2 Amplified","Others"), alpha = alpha)
summary(res)
res
res_shrunken <- lfcShrink(dds, res = res, contrast = c("Subtype","ERBB2 Amplified","Others"), type = "normal")
summary(res_shrunken)
res_shrunken

# check meaning of each column
class(res)
mcols(res, use.names = TRUE)
mcols(res)$description

# MA plot
plotMA(res, ylim = c(-2, 2), xlab = "Mean of Normalized Counts", ylab = "Log Fold Change"
       , colNonSig = "grey", colSig = "#85D4E3") 
abline(h=c(-1,1), col="grey", lwd=2)
title(main = "A. MA Plot with Unshrunken Results", font.main = 1, adj = 0)

plotMA(res_shrunken, ylim = c(-2, 2), xlab = "Mean of Normalized Counts", ylab = "Log Fold Change"
       , colNonSig = "grey", colSig = "#85D4E3")
abline(h=c(-1,1), col="grey", lwd=2)
title(main = "B. MA Plot with Shrunken Results", font.main = 1, adj = 0)

# obtain significantly differentially expressed gene and sort by log2FoldChange
signif <- which(res_shrunken$padj < padj.cutoff & abs(res_shrunken$log2FoldChange) >= lfc.cutoff)

deg <- res_shrunken[signif, ] |>
  data.frame() |>
  rownames_to_column(var = "Gene") |>
  as_tibble() |>
  arrange(desc(abs(log2FoldChange)))

deg_ <- res[which(res$padj < padj.cutoff), ] |>
  data.frame() |>
  rownames_to_column(var = "Gene") |>
  as_tibble() |>
  arrange(desc(abs(log2FoldChange)))

# filter top 10 by absolute log2FoldChange
deg_top10 <- deg[1:10,]
deg_top10_gene <- pull(deg_top10, Gene)
write.xlsx(deg_top10, "Top 10 Differentially Expressed Genes.xlsx")

# get normalized counts for top 10 differentially expressed gene
deg_top10_norm <- normalized_counts |>
  filter(Gene %in% deg_top10$Gene)

# gather the columns to have normalized counts to a single column
gathered_deg_top10 <- deg_top10_norm |>
  gather(colnames(deg_top10_norm)[2:length(colnames(deg_top10_norm))], key = "PATIENT_ID", value = "Normalized_Counts")
metadata_ <- rownames_to_column(metadata, var = "PATIENT_ID") 
gathered_deg_top10$PATIENT_ID <- gsub("\\.", "-", gathered_deg_top10$PATIENT_ID)

gathered_deg_top10 <- inner_join(metadata_, gathered_deg_top10) # merge with metadata to get subtype group

ggplot(gathered_deg_top10) +
  geom_boxplot(aes(x = Gene, y = Normalized_Counts, color = Subtype)) +
  scale_y_log10() +
  labs(x = "Genes", y = "log10 Normalized Counts", color = "Subtype"
       , title = "Top 10 Differentially Expressed Gene Normalized Count", subtitle = "ERBB2 Amplified vs Others") +
  theme(legend.position="bottom") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("#F8AFA8", "#74A089"))


# volcano plot
res_shrunken_vol <- res_shrunken |>
  data.frame() |>
  rownames_to_column(var = "Gene") |>
  as_tibble() |>
  mutate(threshold = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff) |>
  mutate(expression = ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated")) |>
  drop_na() |>
  mutate(genelabel = "") |>
  #arrange(padj)
  arrange(desc(abs(log2FoldChange)))

res_shrunken_vol$genelabel[1: 10] <- res_shrunken_vol$Gene[1: 10]

res_shrunken_vol <- res_shrunken_vol |>
  arrange(padj)

res_shrunken_vol$genelabel[1: 10] <- res_shrunken_vol$Gene[1: 10]

ggplot(res_shrunken_vol, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = expression)) +
  geom_text_repel(aes(label = genelabel), size = 3) +
  ggtitle("Gene Expression Level") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  labs(color = "Differential Expression", subtitle = "Up-regulated vs Down-regulated") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.25))) +
  scale_color_manual(values = c("#9A8822", "#FAD77B"))

# res_vol <- res |>
#   data.frame() |>
#   rownames_to_column(var = "Gene") |>
#   as_tibble() |>
#   mutate(threshold = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff) |>
#   drop_na() |>
#   mutate(genelabel = "") |>
#   arrange(padj)
# 
# res_vol$genelabel[1: 10] <- res_vol$Gene[1: 10]
# 
# ggplot(res_vol, aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(aes(colour = threshold)) +
#   geom_text_repel(aes(label = genelabel), size = 3) +
#   ggtitle("Gene Expression Level") +
#   xlab("log2 fold change") + 
#   ylab("-log10 adjusted p-value") +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.5)),
#         axis.title = element_text(size = rel(1.25))) +
#   scale_color_manual(values = c("#9A8822", "#FAD77B"))
  

# 4. pathway enrichment analysis
# separate up/down regulated genes
dup <- deg[deg$log2FoldChange > 0., ]
ddown <- deg[deg$log2FoldChange < 0., ]

# get Entrez gene id for pathway enrichment
entrez_all <- data_Rnaseq |>
  filter(Hugo_Symbol %in% deg$Gene) |>
  select(Entrez_Gene_Id)
entrez_all <- entrez_all[, 1]

entrez_up <- data_Rnaseq |>
  filter(Hugo_Symbol %in% dup$Gene) |>
  select(Entrez_Gene_Id)
entrez_up <- entrez_up[, 1]

entrez_down <- data_Rnaseq |>
  filter(Hugo_Symbol %in% ddown$Gene) |>
  select(Entrez_Gene_Id)
entrez_down <- entrez_down[, 1]

# KEGG pathway over-representation analysis
all_paths <- enrichKEGG(gene = entrez_all, organism = 'hsa', pvalueCutoff = 0.05)
head(all_paths)

up_paths <- enrichKEGG(gene = entrez_up, organism = 'hsa', pvalueCutoff = 0.05)
head(up_paths)

down_paths <- enrichKEGG(gene = entrez_down, organism = 'hsa', pvalueCutoff = 0.05)
head(down_paths)

dotplot(all_paths, showCategory = 10, title = "Enriched Pathways")
dotplot(up_paths, showCategory = 10, title = "Enriched Pathways")
dotplot(down_paths, showCategory = 10, title = "Enriched Pathways")

cnetplot(all_paths, categorySize="pvalue", foldChange=entrez_all) 
cnetplot(up_paths, categorySize="pvalue", foldChange=entrez_up)
cnetplot(down_paths, categorySize="pvalue", foldChange=entrez_down)
# browseKEGG(all_paths, 'hsa04080') 

# pathview(gene.data  = gene,
#          pathway.id = "hsa04151",
#          species    = "hsa",
#          limit      = list(gene=max(abs(gene)), cpd=1))



# 5. principal component analysis

# calculate the variance stabilised transformed expression values
vsd <- vst(dds, blind = FALSE)
colData(vsd)

pca <- plotPCA(vsd, intgroup = "Subtype", returnData = TRUE)
ggplot(pca) +
  geom_point(aes(x = PC1, y = PC2, color = Subtype)) +
  scale_color_manual(values = c("#9A8822", "#85D4E3")) +
  labs(title = "PCA Plot", subtitle = "ERBB2 Amplified vs Others") +
  theme(legend.position="bottom")
  

vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat))
df <- as.data.frame(cbind(metadata, pca$x))
df <- rownames_to_column(df)
df <- mutate_at(df, c("PC3", "PC4"), as.numeric)
ggplot(df, aes(x = PC3, y = PC4, color = Subtype)) + geom_point()


# 6. gene expression cluster

a <- as.data.frame(head(deg))
a <- column_to_rownames(a, var = "Gene")
deg.dist1 <- dist(a, method = "euclidean")
deg.dist1
deg.clust1 <- hclust1(deg.dist)
plot(deg.clust1)

deg.dist2 <- 1-cor(t(a))
deg.dist2
deg.clust2 <- hclust(as.dist(deg.dist2))
plot(deg.clust2)

deg_top100 <- deg[1:100,]
b <- data.frame(t(deg_top100[-1]))
colnames(b) <- deg_top100$Gene
b1 <- as.matrix(b)

deg.cor <- cor(b1, method = "pearson")
deg.hc <- hclust(as.dist(1-deg.cor), method = "complete")

heat.colors <- brewer.pal(11,"RdBu")
pheatmap(deg.cor, col = heat.colors, cluster_rows = TRUE, cluster_cols = TRUE, scale = "none", fontsize = 2)

pheatmap(deg.cor
         , color = heat.colors
         , cluster_rows = TRUE
         , show_rownames = TRUE
         , border_color = NA
         , font_size = 10
         , scale = "row"
         , fontsize_row = 10
         , height = 20
         , show_colnames = TRUE
)



# heatmap 
norm_deg <- normalized_counts |>
  filter(Gene %in% deg$Gene, Gene != "") |> 
  data.frame() |>
  column_to_rownames(var = "Gene")

deg_top10_heat <- deg_top10_norm |>
  data.frame() |>
  select(-Gene) |>
  as.matrix()

rownames(deg_top10_heat) <- deg_top10_norm$Gene

deg_top10_heat <- deg_top10_heat |>
  t() |>
  scale() |>
  t()

#heat.colors <- wes_palette("Royal2", 5)
heat.colors <- brewer.pal(11,"RdBu")


metadata_$PATIENT_ID <- gsub("-", "\\.", metadata_$PATIENT_ID)
annotation <- metadata_ |>
  #mutate_at(metadata_$PATIENT_ID, gsub("-", "\\.")) |>
  data.frame(row.names = "PATIENT_ID")

pheatmap(deg_top10_heat
         , color = heat.colors
         , cluster_rows = TRUE
         , show_rownames = TRUE
         , border_color = NA
         , font_size = 10
         , scale = "row"
         , fontsize_row = 10
         , height = 20
         , show_colnames = TRUE
         , annotation = annotation
)

# 7. survival analysis
# prepare data
vsd_deg <- vsd |>
  assay() |>
  data.frame() |>
  rownames_to_column(var = "Gene") |>
  filter(Gene %in% deg$Gene) |>
  column_to_rownames(var = "Gene") |>
  as.matrix()

metadata_$PATIENT_ID <- gsub("\\.", "-", metadata_$PATIENT_ID)

sur_data <- data_patient |>
  mutate(OS_STATUS_ = as.numeric(substr(OS_STATUS,1,1))) |>
  filter(PATIENT_ID %in% metadata_$PATIENT_ID) 

sur_data <- merge(sur_data, metadata_, by = "PATIENT_ID")

data_combined <- data.frame(t(vsd_deg), sur_data)
data_combined$Subtype <- factor(data_combined$Subtype, levels = c("Others", "ERBB2 Amplified"))

# check categorical data
unique(data_combined$Subtype)
class(data_combined$Subtype)

# run cox regression
cox_res <- coxph(Surv(OS_MONTHS, OS_STATUS_) ~ Subtype + KRT12 + PZP + BTN1A1 + NEUROD2 + MYOC + PNMT + VCX3A + LYPD4 + ZPBP2 + KRT1
                   , data = data_combined, x = TRUE)
cox_res
summary(cox_res)

# check statistics
cox.zph(cox_res)

# check data pattern overtime
ggcoxzph(cox.zph(cox_res), point.size = 1)

# hazard ratio
ggforest(cox_res, data = data_combined)

# survival curve by subtype
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ Subtype, data = data_combined)
cox_fit
#plot(cox_fit)
ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.title="Subtype",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Breast Cancer Survival Curve ", 
           subtitle = "ERBB2 Amplified vs Others",
           risk.table.height=.3,
           xlab = "Time (Month)")

# survival curve by gene
plot_surv_area(time = "OS_MONTHS"
               , status = "OS_STATUS_"
               , variable = "KRT12"
               , data = data_combined
               , model = cox_res
               , title = "Brest Cancer Survival Curve - KRT12"
               , xlab = "Time (Month)"
               , start_color = "lightblue"
               , end_color = "orchid2")

plot_surv_area(time = "OS_MONTHS"
               , status = "OS_STATUS_"
               , variable = "ZPBP2"
               , data = data_combined
               , model = cox_res
               , title = "Brest Cancer Survival Curve - ZPBP2"
               , xlab = "Time (Month)"
               , start_color = "lightblue"
               , end_color = "orchid2")

plot_surv_area(time = "OS_MONTHS"
               , status = "OS_STATUS_"
               , variable = "LYPD4"
               , data = data_combined
               , model = cox_res
               , title = "Brest Cancer Survival Curve - LYPD4"
               , xlab = "Time (Month)"
               , start_color = "lightblue"
               , end_color = "orchid2")

plot_surv_area(time = "OS_MONTHS"
               , status = "OS_STATUS_"
               , variable = "VCX3A"
               , data = data_combined
               , model = cox_res
               , title = "Brest Cancer Survival Curve - VCX3A"
               , xlab = "Time (Month)"
               , start_color = "lightblue"
               , end_color = "orchid2")