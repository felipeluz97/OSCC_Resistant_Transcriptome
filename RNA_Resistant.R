
# if (!requireNamespace("edgeR", quietly = TRUE)) {
#   install.packages("BiocManager")
#   BiocManager::install("edgeR")
# }
BiocManager::install("clusterProfiler")
BiocManager::install("singscore")

library(edgeR)
library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(singscore)


# Load count matrix
counts <- read.csv("~/projetos/felipeluz/dieila/RNAseq/counts.csv", row.names = 1)
# 6. Map ENSG tp HGNC with biomaRt
entrez_id = mapIds(
  org.Hs.eg.db,
  keys = rownames(counts),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

counts$SYMBOL = entrez_id

counts = counts %>%
  dplyr::filter(SYMBOL %in% counts$SYMBOL[!is.na(counts$SYMBOL)]) # filter out NAs

counts$GENE = rownames(counts)
counts$GENE[duplicated(counts$SYMBOL)] # to know duplicated genes

counts = counts %>%
  dplyr::filter(!SYMBOL %in% counts$SYMBOL[duplicated(counts$SYMBOL)]) # filter out duplicated ENSEMBL

rownames(counts) = counts$SYMBOL
counts$SYMBOL = NULL
counts$GENE = NULL

# 2. Define groups
group <- factor(c("parental", "parental", "intermed", "intermed", "resistant", "resistant"),
                levels = c("parental", "intermed", "resistant"))

# 3. Make DGEList
dge <- DGEList(counts = counts, group = group)

# 4. Filter out low expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# 5. Normalize
dge <- calcNormFactors(dge)

# 6. Define design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 7. Estimate dispersion
dge <- estimateDisp(dge, design)

# 8. Fit model
fit <- glmFit(dge, design)

# 9. Define contrasts
contrasts <- list(
  Intermed_vs_Parental   = c(-1, 1, 0),
  Resistant_vs_Parental  = c(-1, 0, 1),
  Resistant_vs_Intermed  = c( 0, -1, 1)
)


# 7. Loop
for (name in names(contrasts)) {
  contrast_vector <- contrasts[[name]]
  
  # Try
  lrt <- glmLRT(fit, contrast = contrast_vector)
  results <- topTags(lrt, n = Inf)$table
  
  # Save
  file_name <- paste0("edgeR_results_", name, "_HGNC.csv")
  write.csv(results, file = file_name, row.names = TRUE)
  
  message("Arquivo salvo: ", file_name)
}


#####################################################################
############### FUNCTIONAL ENRICGMENT ###############################
#####################################################################

res_par <- read.csv("~/projetos/felipeluz/dieila/RNAseq/edgeR_results_Resistant_vs_Parental_HGNC.csv")
int_par <- read.csv("~/projetos/felipeluz/dieila/RNAseq/edgeR_results_Intermed_vs_Parental_HGNC.csv")
res_int <- read.csv("~/projetos/felipeluz/dieila/RNAseq/edgeR_results_Resistant_vs_Intermed_HGNC.csv")

results = list(res_par,int_par,res_int)
names  = c("resXpar","intXpar","resXint")


sig_genes <- results[[i]][results[[i]]$FDR < 0.05, ]
gene_list <- sig_genes$X

ego_BP <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",         # Biological Process 
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

ego_MF <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "MF",         # Molecular function 
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

png(paste("~/projetos/felipeluz/dieila/RNAseq/dotplot_GO_top5",names[i], ".png"), height = 800, width = 1500)
dotplot(ego_BP, showCategory=5) + ggtitle("GO - BP") | dotplot(ego_MF, showCategory=5) + ggtitle("GO - MF")
dev.off()


#####################################################################
###############           VOLCANO           #########################
#####################################################################


results[[i]]$significant <- with(results[[i]], ifelse(FDR < 0.01 & abs(logFC) > 5, "yes", "no"))

label_genes <- results[[i]][results[[i]]$significant == "yes", ]
label_genes <- label_genes[order(-log10(label_genes$PValue)), ]
label_genes <- label_genes[1:30,]

png(paste("~/projetos/felipeluz/dieila/RNAseq/volcano",names[i], ".png"), height = 500, width = 1000)
ggplot(results[[i]], aes(x = logFC, y = -log10(PValue), color = significant)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(
    data = label_genes,
    aes(label = X), 
    size = 5,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano plot: Intermed vs Parental",
    x = "log2 Fold Change",
    y = "-log10 P-value"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")
dev.off()

#####################################################################
###############         HEATMAP      ###############################
#####################################################################

top_genes <- head(sig_genes[order(sig_genes$FDR), ], 50)$X
logCPM <- cpm(dge, log=TRUE)
heatmap_data <- logCPM[rownames(logCPM) %in% top_genes, ]
colnames(heatmap_data) <- paste(group, 1:2, sep="_")
set.seed(123)

png("~/projetos/felipeluz/dieila/RNAseq/heatmap.png", height = 800, width = 800)
pheatmap(heatmap_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "row"
         )
dev.off()

#####################################################################
###############        PCA            ###############################
#####################################################################

logCPM_all <- cpm(dge, log=TRUE)
pca <- prcomp(t(logCPM_all), scale. = TRUE)


pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     group = group)
png("~/projetos/felipeluz/dieila/RNAseq/PCA.png", height = 500, width = 500)
ggplot(pca_df, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=4) +
  theme_minimal() +
  labs(title = "PCA dos samples (logCPM)",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)")) +
  scale_color_brewer(palette = "Set2")
dev.off()


#####################################################################
###############        MARKERS            ###########################
#####################################################################


# Unificar por hgnc_symbol
merged <- merge(res_par[, c("X", "logFC", "FDR")],
                int_par[, c("X", "logFC", "FDR")],
                by = "X", suffixes = c("_ResPar", "_IntPar"))

merged <- merge(merged,
                res_int[, c("X", "logFC", "FDR")],
                by = "X")

colnames(merged)[c(5,6)] <- c("logFC_ResInt", "FDR_ResInt")

#   Resistance markers
resistant_markers <- merged[
  merged$logFC_ResPar > 5  &
    merged$FDR_ResPar < 0.01,  ]

resistant_markers <- resistant_markers[order(resistant_markers$FDR_ResPar), ]

# Intermed  markers
intermed_markers <- merged[
  merged$logFC_IntPar > 5 &
    merged$FDR_IntPar < 0.01 , ] # NAo tem 

# Visualize
head(resistant_markers)


# Export
write.csv(resistant_markers, "~/projetos/felipeluz/dieila/RNAseq/markers_resistant.csv", row.names = FALSE)
resistant_markers = read.csv("~/projetos/felipeluz/dieila/RNAseq/markers_resistant.csv")

genes_resistant_Di = resistant_markers$X[1:30]


genes_resistant_paper = read.csv("~/projetos/felipeluz/dieila/RNAseq/resistance_genes_paper.csv")
genes_resistant_paper = genes_resistant_paper$HUGO_Gene



#######################
# sing Score  #########
#######################


rankData <- rankGenes(logCPM)

################################
### Plot rank desity ###########
################################
plotRankDensity(rankData[,2,drop = FALSE], upSet = genes_resistant_paper,downSet = genes_resistant_Di, isInteractive = FALSE)


## -- Score Giomo resistance
scoredf_Di <- simpleScore(rankData, 
                          upSet = genes_resistant_Di,
                          centerScore=FALSE)


scoredf_Di$sample = rownames(scoredf_Di)






## -- Score paper resistance
scoredf_paper <- simpleScore(rankData, 
                             upSet = genes_resistant_paper,
                             centerScore=FALSE)


scoredf_paper$sample = rownames(scoredf_paper)

########################################################
### Dieila Resistance related genes heatmap ############
########################################################

# Genes resistance Dieila

Genes_Di_log_cpm_data  = logCPM[genes_resistant_Di,]



png("~/projetos/felipeluz/dieila/RNAseq/heatmap_50_genes_Di.png", height = 800, width = 800)
pheatmap(Genes_Di_log_cpm_data,
         cluster_cols = TRUE,
         scale = "row",
         show_colnames = TRUE,
         cutree_cols = 3)

dev.off()

######################################
## Plot Landscape ####################
######################################
plot = plotScoreLandscape(scoredf_paper, scoredf_Di, 
                          scorenames = c('paper','Di'),hexMin = 1)




projectScoreLandscape(plotObj = plot,scoredf_paper, scoredf_Di,
                      subSamples = rownames(scoredf_Di), 
                      sampleLabels = rownames(scoredf_Di),
                      isInteractive = FALSE)

Genes_Paper_log_cpm_data  = logCPM[rownames(logCPM)[rownames(logCPM) %in% genes_resistant_paper],]
png("~/projetos/felipeluz/dieila/RNAseq/heatmap_812_genes_paper.png", height = 800, width = 800)
pheatmap(Genes_Paper_log_cpm_data,
         cluster_cols = TRUE,
         scale = "row",
         show_colnames = TRUE,
         cutree_cols = 3)

dev.off()

