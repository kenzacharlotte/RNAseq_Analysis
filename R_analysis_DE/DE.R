setwd("//172.23.2.37/panasyuk lab/DATA%20PROJECTS/Nuclear%20Transcription_class%203%20PI3K/NT_Exp85_PI3P_probe_MEFs_EBSS_Methionine_Starvation_Stimulation/EBSS%20n1/R")

setwd("smb://172.23.2.37/panasyuk lab/DATA PROJECTS/Nuclear Transcription_class 3 PI3K/NT_Exp66_Nuclear_PI3P_vps34inh1 & OICR/PI3P probe_quantifications")
install.packages("installr")

#install.packages("installr")

#================================ Loading libraries ================================= #
#updateR() 
library(installr)
library("DESeq2")
library(ggplot2)
library(ggrepel)
library(FactoMineR) ## for PCA
library(factoextra) ## for visualization functions
library(pheatmap) # for heatmap

#================================ Spot checking================================= #
# Example for KO check expression : 
fc <- read.delim("/media/user/LINUX/linux/kenza_work/Nate_FedFastWT-RNAseq.Feb-June2020/RNAseqAnalysis/featurecounts/featurecounts.txt", comment.char = "#", check.names = FALSE)
head(fc)

#================================ DESeq analysis - My code ================================= #
#================================ Build your object ================================= #

# Loading featurecounts table
fc <- read.delim("/media/user/LINUX/linux/kenza_work/Nate_FedFastWT-RNAseq.Feb-June2020/RNAseqAnalysis/featurecounts/featurecounts.txt", comment.char = "#", check.names = FALSE)
rownames(fc) <- fc$Geneid
head(fc)

# Keep only the counts from featurecounts
rawcounts <- fc[, -(1:6)]
rawcounts <- as.matrix(rawcounts)
storage.mode(rawcounts) <- "numeric"
bam_names <- colnames(fc)[-(1:6)] # gives bam names

# Remove NA values from rawcounts table
rawcounts <- na.omit(rawcounts)
is.na(rawcounts)

# How many reads do I have per sample? 
colSums(rawcounts)

# First quality check
# Make bar plot of library sizes 
libsize <- colSums(rawcounts)
barplot(libsize, las = 2, main = "Barplot of lib size")

# Draw distribution of raw counts (log2)
boxplot(log2(1+rawcounts), las=2, ylab="raw counts (log2)",col="gray50", pch=16)

# Sample information : describe your experimental design
colData <- data.frame(
  diet = c("fed", "fed", "fed", "fast", "fast", "fast"),
  row.names = bam_names)

colData

# Control colData !
# Make sure that rownames in colData are matching with column names in rawcounts 
all((colnames(rawcounts)) %in% rownames(colData))
# Are they in the same order ?
all((colnames(rawcounts)) == rownames(colData))

# Construct DESeq object
# - counts : matrix obtained from featurecounts, (rows = gene ID x col = sample with different conditions)
# - colData : a data frame containing the experimental design, (row names = colnames(counts), colnames = experimental design)
# - design factor: the formula expresses how the counts for each gene depend on the variables in colData. Many R formula are valid, including designs with multiple variables, e.g., ~ group + condition, and designs with interactions, e.g., ~ genotype + treatment + genotype:treatment
# N.B. : here we want to compare the diet Fasting vs Feeding condition
dds <- DESeqDataSetFromMatrix(countData = rawcounts, colData = colData, design = ~ diet)

dds
# Pre-filtering (optionnal)
# keeping rows that have at least 10 reads total ?
keep <- rowSums(counts(dds)) >= 0
dds <- dds[keep,]


# Set the factor level : L’ordre des niveaux détermine la comparaison statistique qui sera faite (R classe automatiquement les niveaux par ordre alphabétique)
# Le modèle est toujours construit à partir de la référence : sinon classé par ordre alphabétique
dds$diet <- relevel(dds$diet, ref = "fed") # fast compared to fed
dds$diet # our reference level : better to be explicit

#================================ Explore dataset - Visualisation ================================= #
# Size factor :  Un size factor est un coefficient de normalisation par échantillon
# => corrige les différences de profondeur de séquençage et de composition.
# Deux échantillons peuvent avoir :
# - un nombre total de reads différent
# - un taux de mapping différent
# - une bibliothèque plus “riche” que l’autre
dds <- estimateSizeFactors(dds) #estimate size factor
sf <- sizeFactors(dds)

barplot(sf, las = 2, ylab = "Size factor", main = "DESeq2 size factors")
abline(h = 1, lty = 2)

# Plot PCA on vst
#Pourquoi VST avant la PCA ? Les comptes RNA-seq ont :
# - une variance dépendante de la moyenne
# - une distribution très asymétrique
# VST (Variance Stabilizing Transformation) : blind = TRUE : ignore le design => exploration pure
# - normalise (size factors)
# - stabilise la variance
# - rend les distances euclidiennes pertinentes
vsd <- vst(dds, blind = TRUE)

## Plot raw counts, log2counts, and vst normalized data to see effect
par(mfrow=c(1,3))

plot(counts(dds, normalized=FALSE)[,1],counts(dds, normalized=FALSE)[,4],
     pch=16, cex=0.3, xlim=c(0,20e3), ylim=c(0,20e3), main = "raw counts")

plot(log2(counts(dds, normalized=FALSE)[,3:4]+1),
     pch=16, cex=0.3, main = "log2 counts")

#QUESTION NATE : PK ce plot ? 
plot(assay(vsd)[,3:4], 
     pch=16, cex=0.3, main="vst normalized counts")

# Plot PCA with DESEQ (ugly)
ntop.pca <- DESeq2::plotPCA(vsd, intgroup = c("diet"), ntop=1000)
ntop.pca

# Plot PCA using FactoMine and factoextra
rawcounts.vsd <- assay(vsd) #extract vst matrix 
rawcounts.vsd <- data.frame(rawcounts.vsd)
gvar <- apply(rawcounts.vsd, 1, var)    # variance par gène
mostvargenes <- head(order(gvar, decreasing = TRUE), 1000) #select most variable genes
res_pca <- PCA(t(rawcounts.vsd[mostvargenes, ]), ncp=3, graph = FALSE)

fviz_pca_ind(res_pca,
             geom = "point",
             label = "ind",                    # show sample names
             habillage = colData(vsd)$diet,    # color by condition
             repel = TRUE,                     # avoid label overlap
             addEllipses = TRUE,               # ellipse per group
             ellipse.level = 0.95,             # 95% confidence
             title = "RNAseq PCA: Top 1000 Variable Genes",
             ggtheme = theme_minimal(base_size = 14))

fviz_pca_ind(res_pca,
             geom = "point",
             repel = TRUE)  # labels intelligents

fviz_pca_ind(res_pca,
             geom = "point",
             habillage = colData(vsd)$diet,
             repel = TRUE)

# Plot PCA using ggplot
coord <- as.data.frame(res_pca$ind$coord)
coord$sample <- rownames(coord)
coord$condition <- colData(vsd)$diet

eig <- res_pca$eig
pc1_var <- round(eig[1, 2], 1)  # % variance expliquée par PC1
pc2_var <- round(eig[2, 2], 1)  # % variance expliquée par PC2

ggplot(coord, aes(x = Dim.1, y = Dim.2, color = condition, label = sample)) +
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
my_palette <- c("fed" = "brown2", "fast" = "darkseagreen")  # adapter à tes conditions
coord$sample_label <- c("Fed_1","Fed_2","Fed_3","Fast_4","Fast_5","Fast_6")
ggplot(coord, aes(x = Dim.1, y = Dim.2, label = sample_label)) +
  
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
    name = "Diet condition",   # nom de la légende
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
  ggtitle("PCA RNAseq - Top 1000 Variable Genes")

#================================ Create your DESeq object ================================= #
# Fit your model : 
dds <- DESeq(dds)
# using pre-existing size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

res <- results(dds)
res # display the results
summary(res)

resdf <- as.data.frame(res)                 # convert to data.frame
resdf$gene <- rownames(resdf)    

## update colnames of count table to your sample names
#colnames(rawcounts) <- as.character(colData[colnames(rawcounts), "sname"])

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# Contrasts
resultsNames(dds)

#================================ Applying Threshold ================================= #
genes_pval0.01_LFC1.2 <- rownames(res)[ !is.na(res$padj) &
                                          !is.na(res$log2FoldChange) &
                                          res$padj < 0.01 & abs(res$log2FoldChange) > 0.263]
genes_pval0.01_LFC1.2
lgenes_pval0.01_LFC1.2 <- list(genes_pval0.01_LFC1.2)
genes_pval0.01_LFC2 <- rownames(res)[ !is.na(res$padj) &
                                           !is.na(res$log2FoldChange) &
                                           res$padj < 0.01 & abs(res$log2FoldChange) > 1]
genes_pval0.01_LFC2
lgenes_pval0.01_LFC2 <- list(genes_pval0.01_LFC2)

install.packages("readxl") # CRAN version

#Nate results

library(readxl)
res_pval0.01_f2 <- read_excel("/media/user/LINUX/linux/kenza_work/Nate_FedFastWT-RNAseq.Feb-June2020/RNAseqAnalysis/Nate_results/res_DEseq2_(Fast_WT)_vs_(Fed_WT)_P50_p(SC)0.01_f2_(1439).xlsx")
Nate0.01_2 <- list(res_pval0.01_f2$Gene)
res_pval0.01_f1.2 <- read_excel("/media/user/LINUX/linux/kenza_work/Nate_FedFastWT-RNAseq.Feb-June2020/RNAseqAnalysis/Nate_results/res_DEseq2_(Fast_WT)_vs_(Fed_WT)_P50_p(SC)0.01_f1.2_(3657).xlsx")
Nate0.01_1.2 <- list(res_pval0.01_f1.2$Gene)


#Comparison
install.packages("ggVennDiagram")
library(ggVennDiagram)
all.equal(res_pval0.01_f2$Gene, genes_pval0.01_LFC2)
mapply(identical,res_pval0.01_f2$Gene, genes_pval0.01_LFC2)

install.packages("eulerr")
library(eulerr)
x <-list(kenza = genes_pval0.01_LFC2,nate = res_pval0.01_f2$Gene)
fit <- euler(x)  # fits proportional areas
plot(fit, main="pval 0.01 & |LFC|>2", quantities = TRUE)  # basic plot with number

plot(fit, quantities = TRUE)  # basic plot with numbers
x <-list(kenza = genes_pval0.01_LFC1.2,nate = res_pval0.01_f1.2$Gene)
fit <- euler(x)  # fits proportional areas
plot(fit, main="pval 0.01 & |LFC|>1.2", quantities = TRUE)  # basic plot with number

ggVennDiagram(c(lgenes_pval0.01_LFC2,Nate0.01_2))+scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
ggVennDiagram(c(lgenes_pval0.01_LFC1.2,Nate0.01_1.2))+scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

length(genes_pval0.01_LFC1.2)
#Union
A <- c(genes_pval0.01_LFC1.2)
B <- c(res_pval0.01_f1.2$Gene)
length(union(A,B))
extA <- setdiff(A, B) #genes in A and not in B = kenza's genes
extB <- setdiff(B, A) #genes in B and not in A = nate's genes
length(extA)
length(extB)

writeLines(extB, "natesgenes_nonoverlapping.txt")
writeLines(extA, "kenzasgenes_nonoverlapping.txt")

#================================ Visualization ================================= #
#================================ MA plot ================================= #
# Example : annotation of the top 5 genes with smallest pvalue
top_genes <- head(res[order(resdf$padj), ], 5) # Pas vraiment informatif : show plot with L2FC and significative pval

# Make sure gene column exists
resdf$gene <- rownames(resdf)

# Top 5 genes by smallest adjusted p-value
top_genes <- head(resdf[order(resdf$padj), ], 5)

ggplot(resdf, aes(x = baseMean, y = log2FoldChange)) +
  
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
    data = subset(resdf, padj < 0.05),
    aes(x = baseMean, y = log2FoldChange),
    shape = 21,
    colour = "black",   # bordure noire
    fill = "brown1",       # remplissage rouge
    size = 1.8,
    stroke = 0.1
  ) +
  
  # Top 5 genes
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 4,
    box.padding = 0.5
  ) +
  
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "MA Plot",
    x = "Mean of normalized counts",
    y = "Log2 fold change"
  )

plotMA(res, ylim=c(-5,5), alpha=0.05,)

#================================ Volcano plot ================================= #
resdf$negLog10Padj <- -log10(resdf$padj)
top_genes <- head(resdf[order(resdf$padj), ], 10) #imo pas important

ggplot(resdf, aes(x = log2FoldChange, y = negLog10Padj)) +

#=============================  
  # Tous les gènes
  geom_point(
    shape = 21,
    colour = "black",
    fill = "grey80",
    size = 1.8,
    stroke = 0.1,
    alpha = 0.8
  ) +
  
  # Gènes significatifs
  geom_point(
    data = subset(resdf, padj < 0.05 & abs(log2FoldChange) > 1),
    aes(x = log2FoldChange, y = negLog10Padj),
    shape = 21,
    colour = "black",
    fill = "red",
    size = 1.8,
    stroke = 0.1
  ) +
  
  # Annotation top 5 gènes
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 4,
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  
  # Seuils
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             color = "grey50") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "grey50") +
  
  theme_bw() +
  labs(
    title = "Volcano plot",
    x = "Log2 fold change",
    y = "-log10(adjusted p-value)"
  )

#=============================
max_fc <- max(abs(resdf$log2FoldChange), na.rm = TRUE) # Pour centrer la figure

ggplot(resdf, aes(x = log2FoldChange, y = negLog10Padj)) +
  
  geom_point(
    shape = 21,
    colour = "black",
    fill = "grey80",
    size = 1.8,
    stroke = 0.3
  ) +
  
  geom_point(
    data = subset(resdf, padj < 0.05 & abs(log2FoldChange) > 1),
    shape = 21,
    colour = "black",
    fill = "red",
    size = 1.8,
    stroke = 0.3
  ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  
  coord_cartesian(xlim = c(-max_fc, max_fc)) +
  
  theme_bw() +
  labs(
    x = "Log2 fold change",
    y = "-log10(adjusted p-value)",
    title = "Volcano plot"
  )
#=============================
# Préparation des threshold
resdf$negLog10Padj <- -log10(resdf$padj)

# Définition des seuils
lfc_cutoff <- 1
padj_cutoff <- 0.05

# Classification des gènes
resdf$regulation <- "NS"
resdf$regulation[resdf$padj < padj_cutoff & resdf$log2FoldChange >  lfc_cutoff] <- "Up"
resdf$regulation[resdf$padj < padj_cutoff & resdf$log2FoldChange < -lfc_cutoff] <- "Down"

# Comptage
n_up   <- sum(resdf$regulation == "Up")
n_down <- sum(resdf$regulation == "Down")

# Top gènes à annoter (optionnel)
top_genes <- head(resdf[order(resdf$padj), ], 5)

max_fc <- max(abs(resdf$log2FoldChange), na.rm = TRUE)
ymax_data <- max(resdf$negLog10Padj, na.rm = TRUE)
ymax_plot <- ymax_data * 1.3   # espace garanti pour les annotations
xmax <- max(abs(resdf$log2FoldChange), na.rm = TRUE)

ggplot(resdf, aes(x = log2FoldChange, y = negLog10Padj)) +
  
  # Tous les points
  geom_point(
    aes(fill = regulation),
    shape = 21,
    colour = "black",
    size = 1.8,
    stroke = 0.1,
    alpha = 0.8
  ) +
  
  # Couleurs
  scale_fill_manual(
    values = c(
      "Up"   = "brown1",
      "Down" = "cornflowerblue",
      "NS"   = "grey80"
    )
  ) +
  
  # Seuils
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed", color = "grey50") +
  
  # Annotation top gènes
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 4,
    box.padding = 0.5
  ) +
  coord_cartesian(xlim = c(-max_fc, max_fc),
                  ylim = c(0, 200)) +
  
  theme_bw() +
  theme(
    legend.position = "none"
  ) +
  
  labs(
    title = "Volcano plot",
    x = "Log2 fold change",
    y = "-log10(adjusted p-value)"
  )


topGenes <- head(order(res$padj), 20)
plotMA(res, ylim=c(-5,5))
points(res$baseMean[topGenes], res$log2FoldChange[topGenes], col="red", cex=1.5)

#================================ Heatmap ================================= #
# Matrice VST
mat <- assay(vsd)

# With pheatmap 
# Sélection des gènes les plus variables
gvar <- apply(mat, 1, var)
top_genes <- head(order(gvar, decreasing = TRUE), 1000)

mat_top <- mat[top_genes, ]

# Z-score par gène (ligne)
mat_scaled <- t(scale(t(mat_top)))

# Annotation des colonnes
annotation_col <- data.frame(
  diet = colData(vsd)$diet
)
rownames(annotation_col) <- colnames(mat_scaled)

# All genes pheatmap
# Ask Nate !
pheatmap(
  mat = (degs_all),
  annotation_col = df_met,
  annotation_colors = list(genotype=c(WT="black",Vps15KO="purple"),
                           methionine=c(all="violetred1",starve_4h = "violetred2", stim_1h = "violetred3", stim_4h = "violetred4")),
  scale = "row",
  show_colnames = FALSE,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = mako(10),
  cutree_cols = 2,
  cutree_rows = 6,
  main = "transcriptional response to methionine starvation and stimulation")

# With complex heatmap
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "RNA-seq heatmap (top 50 variable genes)"
)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

mat <- assay(vsd)

# Gènes les plus variables
gvar <- apply(mat, 1, var)
top_genes <- head(order(gvar, decreasing = TRUE), 50)
mat_top <- mat[top_genes, ]

# Z-score
mat_scaled <- t(scale(t(mat_top)))

# Annotation colonnes
ha <- HeatmapAnnotation(
  diet = colData(vsd)$diet,
  col = list(diet = c("fed" = "red", "fast" = "blue"))
)

Heatmap(
  mat_scaled,
  name = "Z-score",
  top_annotation = ha,
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  show_row_names = FALSE,
  show_column_names = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  column_title = "RNA-seq heatmap (top 50 variable genes)"
)

#================================ Annotation ================================= #
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Mm.eg.db)   # ou org.Hs.eg.db selon ton organisme
library(DOSE)           # pour visualisation
library(enrichplot)

# Ajouter les IDs de gènes si nécessaire
resdf$gene <- rownames(resdf)

# Filtrer DEGs significatifs
deg <- subset(resdf, padj < 0.05 & abs(log2FoldChange) > 1)
gene_list <- deg$gene

# Si tes gènes sont SYMBOL, convertir en ENTREZID pour clusterProfiler
library(biomaRt)

# Exemple pour mouse SYMBOL -> ENTREZID
entrez_ids <- mapIds(
  org.Mm.eg.db,
  keys = gene_list,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

entrez_ids <- na.omit(entrez_ids)

#================================ Annotation NATE================================= #

### if KEGG enrichment does not work, do the following to run KEGG locally ###
#KEGG (download kegg.db from  [https://www.bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz])
#then : install.packages("~/Desktop/path where KEGG.db downloaded/KEGG.db_3.2.4.tar.gz", repos = NULL, type = "source")
#then use the "use_internal_data = TRUE" command on enrichKEGG
################################################
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Mm.eg.db)   # ou org.Hs.eg.db selon ton organisme
library(DOSE)           # pour visualisation
library(enrichplot)

# Ajouter les IDs de gènes si nécessaire
resdf$gene <- rownames(resdf)

# Filtrer DEGs significatifs
deg <- subset(resdf, padj < 0.05 & abs(log2FoldChange) > 1)
gene_list <- deg$gene

# Si tes gènes sont SYMBOL, convertir en ENTREZID pour clusterProfiler
library(biomaRt)

# Exemple pour mouse SYMBOL -> ENTREZID
entrez_ids <- mapIds(
  org.Mm.eg.db,
  keys = gene_list,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

entrez_ids <- na.omit(entrez_ids)

# NATE
## load geneList into DOSE analysis ## 
gene <- names(geneList) [abs(geneList) > 1]

###### GO PATHWAY ANALYSIS ######
genes_GO <- enrichGO(gene,
                     keyType = "ENTREZID",
                     universe = names(geneList),
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     readable = TRUE)
write.csv(genes_GO, file = "DOSE_gsea/DOSE_output/limma_DEG_upreg_gfp_stim4h_GOBP.csv")

goGSEA <- gseGO(geneList = geneList,
                keyType = "ENTREZID",
                OrgDb = org.Mm.eg.db,
                ont = "BP", # Biological Processes
                #pvalueCutoff = 0.05,
                by = "fgsea",
                verbose = TRUE)
write.csv(goGSEA, file = "DOSE_gsea/DOSE_output/limma_DEG_gfp_vs_cre_stim4h_goGSEA.csv")

## plots for GO BP ##
gseaplot2(goGSEA, geneSetID = c("GO:0023056","GO:0048584", "GO:0009893"), ES_geom = "line",base_size = 8, rel_heights = c(1.3,0.3,0.5),
          color = c("dodgerblue4", "dodgerblue3","cornflowerblue"),
          title = "GO - limma_DEG_upreg_gfp_stim1h, DEG limma")+
  theme_minimal()
ggsave("DOSE_output_plots/limma_DEG_upreg_gfp_starve_gseaplot_lines.png")


###### KEGG PATHWAY ANALYSIS ######
genes_KEGG <- enrichKEGG(gene,
                         organism = 'mmu',
                         minGSSize = 1,
                         use_internal_data = F)
write.csv(genes_KEGG, file = "DOSE_gsea/DOSE_output/limma_DEG_upreg_gfp_stim4h_KEGG.csv")

genes_gseaKEGG <- gseKEGG(geneList = geneList,
                          organism = 'mmu',
                          minGSSize = 1,
                          verbose = TRUE,
                          use_internal_data = F)
write.csv(genes_gseaKEGG, file = "DOSE_gsea/DOSE_output/limma_DEG_upreg_gfp_stim4h_gseKEGG.csv")

## plots for KEGG  ##
gseaplot2(genes_gseaKEGG, geneSetID = c(1), ES_geom = "line",base_size = 8, rel_heights = c(1.3,0.3,0.5),
          color = c("purple"),
          title = "KEGG - limma_DEG_upreg_gfp_stim1h, DEG limma")
ggsave("DOSE_output_plots/limma_DEG_upreg_gfp_starve_gseaplot_lines.png")



###### REACOME ANALYSIS ######
genes_reactome <- gsePathway(geneList = geneList,
                             pvalueCutoff = 0.05,
                             organism = "mouse",
                             verbose = TRUE)
write.csv(genes_reactome, file = "DOSE_gsea/DOSE_output/limma_DEG_upreg_gfp_stim4h_gseReactome.csv")
## plots for Reactome ##
gseaplot2(genes_reactome, geneSetID = c(1:5,7), ES_geom = "dot", base_size = 8, rel_heights = c(1.5,0.5,0.5),
          title = "Reactome - upreg in feno-Vps15LKO vs Vps15LKO, DEG limma")



#########===== Additional GSEA enrichment plot visualization using GseaVis =====#########
# https://junjunlab.github.io/gseavis-manual/basic-usage.html#marking-gene-names #
# GseaVis can work with gseKEGG, gseGO, gsePathway (reactome) objects
library(GseaVis)

gseaNb(object = goGSEA,
       geneSetID = c("GO:0043066","GO:0051246", "GO:0010604"),
       curveCol = c("dodgerblue4","dodgerblue3","cornflowerblue"),
       subPlot = 2)

gseaNb(object = goGSEA,
       geneSetID = c("GO:0032502","GO:0019538", "GO:0033554"),
       curveCol = c("dodgerblue4", "dodgerblue3","cornflowerblue","cadetblue2", "skyblue1","lightblue"),
       subPlot = 2)

# to add gene name labels to plot
#  mygenes <- c("gene1", "gene2","gene3")
#  gseaNb(object = goGSEA,
#      geneSetID = c("GO:0023056","GO:0048584", "GO:0009893"),
#      subPlot = 2,
#       addgene = mygenes)


#================================ NATE'S CODEEEEE ================================= #

## remove NA values from rawcounts table
rawcounts <- na.omit(rawcounts)
is.na(rawcounts)

## How many reads do I have per sample? ##
colSums(rawcounts)

## make bar plot of library sizes ##
libsize <- colSums(rawcounts)
barplot(libsize, las = 2, main = "Barplot of lib size")

## Draw distribution of raw counts (log2)
boxplot(log2(1+rawcounts),
        las=2, ylab="raw counts (log2)",col="gray50", pch=16)

# =================================  Variance stabilization (vst) (only use vst for clustering and PCA plots) ================================= #
#install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)
## Load Data
dds <- DESeqDataSetFromMatrix(rawcounts, DataFrame(condition=splan$subtype),
                              ~condition)
## Estmate size factors
dds <- dds[rowSums(counts(dds)) > 0]

## Run vst normalization
rld <-vst(dds, blind=TRUE)
vst_matrix <- assay(rld)
#save vst matrix as a csv file
dir.create("vst_normalized_counts", showWarnings = FALSE)

write.csv(vst_matrix, file = "vst_normalized_counts/vst_normalized_counts.csv", row.names = TRUE)

file.exists("vst_normalized_counts")
# TRUE

## r vst plotting
par(mfrow=c(1,3))
plot(counts(dds, normalized=FALSE)[,3:4],
     pch=16, cex=0.3, xlim=c(0,20e3), ylim=c(0,20e3), main = "raw counts")

plot(log2(counts(dds, normalized=FALSE)[,3:4]+1),
     pch=16, cex=0.3, main = "log2 counts")

plot(assay(rld)[,3:4], 
     pch=16, cex=0.3, main="vst normalized counts")




