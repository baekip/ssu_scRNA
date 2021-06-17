#!/usr/local/bin/R
#
suppressMessages({
	library(Seurat)
	library(dplyr)
	library(patchwork)
	library(Matrix)
	library(Cairo)
})

set.seed(123)
args = commandArgs(TRUE)
options (bitmapType = 'cairo')
sample.id = args[1]
input.path = args[2]
output.path = args[3]
basic.path = paste0 (output.path, "/Basic/")

###load QC.Rda
sample.filt <- readRDS (file = paste0 (input.path, "/Rdata/", sample.id, ".Filt.Rda"))

#-------------------------------------------------------------------------------
##Normalizaing the data
sample.norm <- NormalizeData(sample.filt, 
														 normalization.method = "LogNormalize", 
														 scale.factor = 10000)
sample.fvf <- FindVariableFeatures(sample.norm,
																	 selection.method = "vst",
																	 nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sample.fvf), 10)
top10
# plot variable features with and without labels
CairoPNG (file = paste0 (basic.path, "VariableFeaturePlot.png"), 
          width = 800, height = 600)
plot1 <- VariableFeaturePlot(sample.fvf)
plot2 <- LabelPoints(plot = plot1, 
										 points = top10, 
										 repel = TRUE)
plot1 + plot2
dev.off()
##Scaling the data
all.genes <- rownames(sample.fvf)
sample.scale <- ScaleData(sample.fvf,
													features = all.genes)

#-------------------------------------------------------------------------------
sample.pca <- RunPCA(sample.scale, 
										 features = VariableFeatures(object = sample.scale))
# Examine and visualize PCA results a few different ways
print(sample.pca[["pca"]], 
			dims = 1:5, nfeatures = 5)
VizDimLoadings(sample.pca, 
							 dims = 1:2, reduction = "pca")
DimPlot(sample.pca, reduction = "pca")
DimHeatmap(sample.pca, dims = 1, cells = 500, balanced = TRUE)

CairoPNG (file = paste0(basic.path, "DimHeatmap.png"), width = 800, height = 600)
DimHeatmap(sample.pca, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

## Determine the 'dimensionality; of the dataset
sample.jack <- JackStraw(sample.pca, num.replicate = 100)
sample.score <- ScoreJackStraw(sample.jack, dims = 1:20)

CairoPNG (file = paste0(basic.path, "JackStrawPlot.png"), width = 800, height = 600)
JackStrawPlot(sample.score, dims = 1:15)
dev.off()

CairoPNG (file = paste0(basic.path, "ElbowPlot.png"), width = 800, height = 600)
ElbowPlot(sample.score)
dev.off()

#-------------------------------------------------------------------------------
##Cluster the cells
sample.find <- FindNeighbors(sample.score, dims = 1:10)
sample.cluster <- FindClusters(sample.find, 
															 resolution = 0.5)
head(Idents(sample.cluster), 5)

##Run non-linear dimensional reduction (UMAP)
sample.umap <- RunUMAP(sample.cluster, dims = 1:10)

CairoPNG (file = paste0(basic.path, "merge_umap.png"), width = 800, height = 600)
p1 <- DimPlot(sample.umap, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sample.umap, reduction = "umap", label = TRUE)
p1 + p2
dev.off()

##Run non-linear dimensional reduction (tSNE)
sample.tsne <- RunTSNE(sample.cluster, dims = 1:10)

CairoPNG (file = paste0 (basic.path, "cluster_plot.png"), width=800, height=600)
DimPlot(object = sample.tsne, reduction = "tsne", group.by = "orig.ident")
dev.off()

CairoPNG (file = paste0 (basic.path, "merge_tsne.png"), width = 1500, height = 600)
p1 <- DimPlot(sample.tsne, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(sample.tsne, reduction = "tsne", label = TRUE)
p3 <- DimPlot(sample.tsne, reduction = "tsne", split.by = "orig.ident")
p1 + p2 + p3
dev.off()

CairoPNG (file = paste0 (basic.path, "merge_tsne.split.png"), width = 800, height = 600)
DimPlot(sample.tsne, reduction = "tsne", split.by = "orig.ident")
dev.off()

##write FindAllMarkers
sample.markers <- FindAllMarkers (sample.tsne,
																	logfc.threshold = 0, 
																	test.use = "roc")
write.table (sample.markers, 
						 file = paste0 (basic.path, "expression.cluster.xls"),
						 col.names = TRUE, row.names = FALSE,
						 quote = FALSE, sep = "\t")

#-------------------------------------------------
saveRDS (sample.tsne, file = paste0 (output.path, "/Rdata/", sample.id, ".tSNE.Rda"))

