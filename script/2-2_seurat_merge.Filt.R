#!/usr/local/bin/R
suppressMessages({
	library(Seurat)
	library(dplyr)
	library(patchwork)
	library(Matrix)
	library(Cairo)
})

args = commandArgs(TRUE)
options (bitmapType = 'cairo')
sample.id = args[1]
project.path = args[2]
input.path = paste0 (project.path, "/result/2-1_Seurat_QC/", sample.id)
output.path = paste0 (project.path, "/result/2-2_Seurat_Filt/", sample.id)
filt.path = paste0 (output.path, "/Filt/")
rdata.path = paste0 (output.path, "/Rdata/")

dir.create (output.path, showWarnings = FALSE)
dir.create (filt.path, showWarnings = FALSE)
dir.create (rdata.path, showWarnings = FALSE)
#------------------------------------------------
min.nGene = 200 
max.nGene = 7000 
min.nUMI = -Inf
max.nUMI = Inf
min.mito = -Inf
max.mito = 5 

###load QC.Rda
sample.merge <- readRDS (file = paste0 (input.path, "/Rdata/", sample.id, ".Merge.Rda"))

#-------------------------------------------------
sample.filt <- subset(sample.merge,
											subset = min.nGene < nFeature_RNA & nFeature_RNA < max.nGene & min.nUMI < nCount_RNA & nCount_RNA < max.nUMI & min.mito < percent.mt & percent.mt < max.mito)


CairoPNG (file = paste0 (filt.path, "Filterted_Distribution_of_detected_genes.png"), width = 800, height = 600)
hist(sample.filt@meta.data$nFeature_RNA, breaks=100, 
		 main = paste0 ("Distribution of detected genes - ", sample.id),
		 xlab = "Genes with at least on tag")
dev.off()

CairoPNG (file = paste0 (filt.path, "Filtered_Density_distribution_of_detected_genes.png"), width = 800, height = 600)
d<-density(sample.filt@meta.data$nFeature_RNA)
plot (d, main = paste0 ("Density distribution of detected genes - ", sample.id))
polygon (d, col="red", border="blue")
dev.off()

CairoPNG (file = paste0 (filt.path, "Filtered_Expression_sum_per_cell.png"), width = 800, height = 600)
hist(sample.filt@meta.data$nCount_RNA, breaks=100,
		 main = paste0 ("Expression sum per cell - ", sample.id),
		 xlab = "Sum of expression")
dev.off()

CairoPNG (file = paste0 (filt.path, "Filtered_Density_expression_sum_per_cell.png"), width = 800, height = 600)
d<-density(sample.filt@meta.data$nCount_RNA)
plot(d, main= paste0 ("Density expression sum per cell - ", sample.id ))
polygon (d, col="red", border="blue")
dev.off()

CairoPNG (file = paste0 (filt.path, "Filtered_vlnplot.png"), width=800, height=600)
VlnPlot(sample.filt,
				features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
				ncol = 3)
dev.off()

CairoPNG (file = paste0 (filt.path, "Filtered_FeatureScatter.png"), width = 800, height = 600)
plot1 <- FeatureScatter(sample.filt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sample.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

#-------------------------------------------------
saveRDS (sample.filt, file = paste0 (output.path, "/Rdata/", sample.id, ".Filt.Rda"))



