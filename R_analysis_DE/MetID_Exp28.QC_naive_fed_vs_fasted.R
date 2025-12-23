
#================================ Loading libraries ================================= #
#updateR() 
library(installr)
library("DESeq2")
library(ggplot2)
library(ggrepel)
library(FactoMineR) ## for PCA
library(factoextra) ## for visualization functions
library(pheatmap) # for heatmap
library(RColorBrewer)
library("corrplot")

#================================ Spot checking================================= #
# Example for KO check expression : 
fc <- read.delim("/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/featurecounts/featurecounts.txt", comment.char = "#", check.names = FALSE)
head(fc)

#================================ Loading & description of DS ================================= #
# Loading featurecounts table
rownames(fc) <- fc$Geneid
head(fc)
col <- colnames(fc)
writeLines(col, "exp_design.txt")

# Keep only the counts from featurecounts
rawcounts <- fc[, -(1:6)]
rawcounts <- as.matrix(rawcounts)
storage.mode(rawcounts) <- "numeric"
bam_names <- colnames(fc)[-(1:6)] # gives bam names
num <- c(sub(pattern = ".Aligned.out.sorted.bam","\\1",bam_names))
bam_names <- paste0("mouse_",num)

# Remove NA values from rawcounts table
rawcounts <- na.omit(rawcounts)
is.na(rawcounts)

# How many reads do I have per sample? 
colSums(rawcounts)

# Description of the experimental design
expdsgn <- read.delim("/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/MetID_Exp28.experimentaldesign.txt",row.names = NULL, nrows = 2, check.names = FALSE)
expdsgn<- (t(expdsgn))

Diet <- rownames(expdsgn)
LearningStep <- expdsgn[,1]

colData <- data.frame(
  food = Diet,
  exposure = LearningStep,
  row.names = bam_names)

# Explicit all the different conditions : Diet x Learning i.e. 6
colData$condition <- interaction(colData$food,colData$exposure)

cond_levels <- levels(colData$condition)

condition_colors <- setNames(
  brewer.pal(n = length(cond_levels), name = "Set2"),
  cond_levels) # create a color palette for each condition

colData$colors <- condition_colors[colData$condition]
colData
colData$food <- factor(colData$food)
colData$exposure <- factor(colData$exposure)
#colData[order(colData$condition), ]

# Control colData !
# Make sure that rownames in colData are matching with column names in rawcounts 
colnames(rawcounts) <- bam_names
all((colnames(rawcounts)) %in% rownames(colData))
# Are they in the same order ?
all((colnames(rawcounts)) == rownames(colData))

# Construct DESeq object
# - counts : matrix obtained from featurecounts, (rows = gene ID x col = sample with different conditions)
# - colData : a data frame containing the experimental design, (row names = colnames(counts), colnames = experimental design)
# - design factor: the formula expresses how the counts for each gene depend on the variables in colData. Many R formula are valid, including designs with multiple variables, e.g., ~ group + condition, and designs with interactions, e.g., ~ genotype + treatment + genotype:treatment
# N.B. : here we want to compare the food Fasting vs Feeding condition

# 2. Hyp : learning step and food are independant - FALSE
dds <- DESeqDataSetFromMatrix(countData = rawcounts, colData = colData, design = ~ food*exposure)

dds
# Pre-filtering (optionnal)
# keeping rows that have at least 10 reads total ?
keep <- rowSums(counts(dds)) >= 0
dds <- dds[keep,]

#================================ Exploring DS ================================= # 
# First quality check
# Make bar plot of library sizes 
libsize <- colSums(rawcounts)
colnames(libsize)


png("/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/barplot_rawcounts.png", width = 1200, height = 800, res = 150)

par(mar=c(5,5,4,10))
par(xpd=NA)

barplot(libsize, 
        names.arg = rownames(colData),
        col = colData$colors,
        fill = colData$condition,
        las = 2, 
        main = "Barplot of the library size")
legend("topright",
       inset = c(-0.35, 0),
       legend = levels(colData$condition),
       fill = condition_colors)

dev.off()


# Draw distribution of raw counts (log2)
png("/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/boxplot_rawcounts.png", width = 1200, height = 800, res = 150)

par(mar=c(5,5,4,10))
par(xpd=NA)

boxplot(log2(1+rawcounts),        
        col = colData$colors,
        names = bam_names,
        las=2, 
        ylab="raw counts (log2)",main='Distribution of the raw counts - log2 transformation')

legend("topright",
       inset = c(-0.35, 0),
       legend = levels(colData$condition),
       fill = condition_colors)
dev.off()

# log2 transform
log_counts <- log2(1 + rawcounts)

# distance entre samples (colonnes)
dist_samples <- dist(t(log_counts))

# hierarchical clustering
hc_samples <- hclust(dist_samples, method = "complete")

# plot
png("/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/clustering_samples.png", width = 1200, height = 800, res = 150)

plot(hc_samples, main = "Hierarchical clustering of samples")
dev.off()

png("/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/corrplot.png",
    width = 1200, height = 1000, res = 150)

M<-cor(log_counts)
colnames(M)<-bam_names
rownames(M)<-bam_names
head(round(M,2))
corrplot(M,
         method="color", 
         order = "hclust", 
         main="Sample hierarchical clustering on log2 transformed rawcounts",
         mar=c(5,4,4,2))
dev.off()


#================================ Comparing Fed vs Fasted naive ================================= # 
# Set the factor level : L’ordre des niveaux détermine la comparaison statistique qui sera faite (R classe automatiquement les niveaux par ordre alphabétique)
# Le modèle est toujours construit à partir de la référence : sinon classé par ordre alphabétique

dds_naive <- dds[, dds$exposure == "naive"]
dds_naive$exposure <- droplevels(dds_naive$exposure) 
colData(dds_naive) <- droplevels(colData(dds_naive))
table(dds_naive$exposure)
levels(dds_naive$exposure)

dds_naive$food <- relevel(dds_naive$food, ref = "fed")
design(dds_naive) <- ~ food
dds_naive <- DESeq(dds_naive)
res_naive <- results(dds_naive)

dds_naive <- estimateSizeFactors(dds_naive) #estimate size factor
sf_naive <- sizeFactors(dds_naive)

png("/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/sizefactors_naive.png",
    width = 1200, height = 1000, res = 150)
barplot(sf_naive, las = 2, ylab = "Size factor", main = "DESeq2 size factors")
abline(h = 1, lty = 2)
dev.off()

# Plot PCA on vst
#Pourquoi VST avant la PCA ? Les comptes RNA-seq ont :
# - une variance dépendante de la moyenne
# - une distribution très asymétrique
# VST (Variance Stabilizing Transformation) : blind = TRUE : ignore le design => exploration pure
# - normalise (size factors)
# - stabilise la variance
# - rend les distances euclidiennes pertinentes
vsd_naive <- vst(dds_naive, blind = TRUE)

ntop.pca <- DESeq2::plotPCA(vsd_naive,  ntop=1000)
ntop.pca

# Plot PCA using FactoMine and factoextra
rawcounts.vsd_naive <- assay(vsd_naive) #extract vst matrix 
rawcounts.vsd_naive <- data.frame(rawcounts.vsd_naive)
gvar_naive <- apply(rawcounts.vsd_naive, 1, var)    # variance par gène
mostvargenes_naive <- head(order(gvar_naive, decreasing = TRUE), 1000) #select most variable genes
res_pca_naive <- PCA(t(rawcounts.vsd_naive[mostvargenes_naive, ]), ncp=3, graph = FALSE)

fviz_pca_ind(res_pca_naive,
             geom = "point",
             label = "ind",                    # show sample names
              habillage = colData(vsd_naive)$food,    # color by condition
             repel = TRUE,                     # avoid label overlap
             addEllipses = TRUE,               # ellipse per group
             ellipse.level = 0.95,             # 95% confidence
             title = "RNAseq PCA: Top 1000 Variable Genes",
             ggtheme = theme_minimal(base_size = 14))

fviz_pca_ind(res_pca_naive,
             geom = "point",
             repel = TRUE)  # labels intelligents

fviz_pca_ind(res_pca_naive,
             geom = "point",
             habillage = colData(vsd_naive)$food,    # color by condition
             repel = TRUE)

# Plot PCA using ggplot
coord <- as.data.frame(res_pca_naive$ind$coord)
coord$sample <- rownames(coord)
coord$condition <- colData(vsd_naive)$food

eig <- res_pca_naive$eig
pc1_var <- round(eig[1, 2], 1)  # % variance expliquée par PC1
pc2_var <- round(eig[2, 2], 1)  # % variance expliquée par PC2

png("/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/Fed_vs_Fasted_firstQC/PCA_naive.png", width = 1200, height = 800, res = 150)

ggplot(coord, aes(x = Dim.1, y = Dim.2, label = sample)) +
  geom_point(colour = "black",   # bordure noire
             fill = "brown1",  ) +
  geom_text_repel() +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_line(size = 0.1),
    text = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
  ) +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  ggtitle("PCA RNAseq - Top 1000 Variable Genes")

# Définir palette personnalisée
my_palette <- c("fed" = "brown2", "fasted" = "darkseagreen")  # adapter à tes conditions

ggplot(coord, aes(x = Dim.1, y = Dim.2, label = sample)) +
  # Points avec bordure noire et remplissage coloré selon condition
  geom_point(aes(fill = condition),
             shape = 21,   # cercle remplissable
             colour = "black",
             size = 3,
             stroke = 0.3) +
  
  # Labels des samples
  geom_text_repel(size = 4, box.padding = 0.5,segment.color = NA) +
  
  # Palette manuelle
  scale_fill_manual(
    name = "Food",   # nom de la légende
    values = my_palette
  ) +
  
  # Thème et axes
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_line(size = 0.1),
    text = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
  ) +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  ggtitle("PCA RNAseq - Top 1000 Variable Genes, naive stage")
dev.off()

# Fit your model : 
dds_naive <- DESeq(dds_naive)
# using pre-existing size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

res_naive <- results(dds_naive)
res_naive # display the results
summary(res_naive)

resdf_naive <- as.data.frame(res_naive)                 # convert to data.frame
resdf_naive$gene <- rownames(resdf_naive)    

## update colnames of count table to your sample names
#colnames(rawcounts) <- as.character(colData[colnames(rawcounts), "sname"])

res0.01_naive <- results(dds_naive, alpha = 0.01)
summary(res0.01_naive)

# Contrasts
resultsNames(dds_naive)

#================================ Visualization ================================= #
#================================ MA plot ================================= #
# Example : annotation of the top 5 genes with smallest pvalue
top_genes <- head(resdf_naive[order(resdf_naive$padj), ], 5) # Pas vraiment informatif : show plot with L2FC and significative pval

# Make sure gene column exists
resdf_naive$gene <- rownames(resdf_naive)

# Top 5 genes by smallest adjusted p-value
top_genes <- head(resdf_naive[order(resdf_naive$padj), ], 5)

png("/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/Fed_vs_Fasted_firstQC/MA_naive.png", width = 1200, height = 800, res = 150)
ggplot(resdf_naive, aes(x = baseMean, y = log2FoldChange)) +
  
  # Tous les gènes (non significatifs)
  geom_point(
    shape = 21,
    colour = "black",   # bordure
    fill = "grey80",    # remplissage
    size = 1.8,
    stroke = 0.1,
    alpha = 0.8
  ) +
  
  # Gènes significatifs
  geom_point(
    data = subset(res_naive, padj < 0.01),
    aes(x = baseMean, y = log2FoldChange),
    shape = 21,
    colour = "black",   # bordure noire
    fill = "brown1",       # remplissage rouge
    size = 1.8,
    stroke = 0.1
  ) +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "MA Plot - padj < 0.01",
    x = "Mean of normalized counts",
    y = "Log2 fold change"
  )
dev.off()

#================================ Volcano plot ================================= #
resdf_naive$negLog10Padj <- -log10(res_naive$padj)
top_genes <- head(res_naive[order(resdf_naive$padj), ], 10) #imo pas important
resdf_naive$gene <- rownames(resdf_naive)
top_genes <- head(
  resdf_naive[order(resdf_naive$padj), ],
  5
)
resdf_naive <- subset(resdf_naive, !is.na(padj))

# Préparation
resdf_naive$gene <- rownames(resdf_naive)
resdf_naive$negLog10Padj <- -log10(resdf_naive$padj)
resdf_naive <- subset(resdf_naive, !is.na(padj))

top_genes <- head(
  resdf_naive[order(resdf_naive$padj), ],
  5
)

# Plot
png("/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/Fed_vs_Fasted_firstQC/Volcano_naive.l2fc_.padj_.01.png", width = 1200, height = 800, res = 150)
ggplot(resdf_naive, aes(x = log2FoldChange, y = negLog10Padj)) +
  
  geom_point(
    shape = 21,
    colour = "black",
    fill = "grey80",
    size = 1.8,
    stroke = 0.1,
    alpha = 0.8
  ) +
  
  geom_point(
    data = subset(resdf_naive, padj < 0.01& abs(log2FoldChange) > 1),
    shape = 21,
    colour = "black",
    fill = "red",
    size = 1.8,
    stroke = 0.1
  ) +
  
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             color = "grey50") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "grey50") +
  
  theme_bw() +
  labs(
    title = "Volcano plot – Naive (fasted vs fed)",
    x = "Log2 fold change",
    y = "-log10(adjusted p-value)"
  )
dev.off()

#================================ Applying Threshold ================================= #
genes_pval0.01_LFC1.2 <- rownames(res_naive)[ !is.na(res_naive$padj) &
                                                !is.na(res_naive$log2FoldChange) &
                                                res_naive$padj < 0.01 & abs(res_naive$log2FoldChange) > 0.263]
genes_pval0.01_LFC1.2
lgenes_pval0.01_LFC1.2 <- list(genes_pval0.01_LFC1.2)
genes_pval0.01_LFC2 <- rownames(res_naive)[ !is.na(res_naive$padj) &
                                              !is.na(res_naive$log2FoldChange) &
                                              res_naive$padj < 0.01 & abs(res_naive$log2FoldChange) > 1]
length(genes_pval0.01_LFC1.2)
lgenes_pval0.01_LFC2 <- list(genes_pval0.01_LFC2)

writeLines(
  genes_pval0.01_LFC1.2,
  con = "/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/Fed_vs_Fasted_firstQC/Jiaqi_fed_vs_fasted.naive.DEgenes_pval0.01_LFC1.2.txt"
)

writeLines(
  genes_pval0.01_LFC2,
  con = "/media/user/LINUX/linux/kenza_work/Jiaqi/RNAseq/MetID_Exp28_feedfast_memory_wt_06112025/R_analysis/Fed_vs_Fasted_firstQC/Jiaqi_fed_vs_fasted.naive.DEgenes_pval0.01_LFC2.txt"
)

install.packages("readxl") # CRAN version

#Nate results
library(readxl)
res_pval0.01_f2 <- read_excel("/media/user/LINUX/linux/kenza_work/Nate/FedFastWT-RNAseq.Feb-June2020/RNAseqAnalysis/Nate_results/DESeq_FedvsFast_padj_0.01_f2.xlsx")
Nate0.01_2 <- list(res_pval0.01_f2$Gene)
res_pval0.01_f1.2 <- read_excel("/media/user/LINUX/linux/kenza_work/Nate/FedFastWT-RNAseq.Feb-June2020/RNAseqAnalysis/Nate_results/DESeq_FedvsFast_padj_0.01_f1.2.xlsx")
Nate0.01_1.2 <- list(res_pval0.01_f1.2$Gene)

library(ggVennDiagram)
all.equal(res_pval0.01_f2$Gene, genes_pval0.01_LFC2)
mapply(identical,res_pval0.01_f2$Gene, genes_pval0.01_LFC2)

library(eulerr)
x <-list(Jiaqi = genes_pval0.01_LFC2,nate = res_pval0.01_f2$Gene)
fit <- euler(x)  # fits proportional areas
plot(fit, main="pval 0.01 & |LFC|>2", quantities = TRUE)  # basic plot with number

plot(fit, quantities = TRUE)  # basic plot with numbers
x <-list(Jiaqi = genes_pval0.01_LFC1.2,nate = res_pval0.01_f1.2$Gene)
fit <- euler(x)  # fits proportional areas
plot(fit, main="pval 0.01 & |LFC|>1.2", quantities = TRUE)  # basic plot with number

length(genes_pval0.01_LFC1.2)
length(genes_pval0.01_LFC2)
