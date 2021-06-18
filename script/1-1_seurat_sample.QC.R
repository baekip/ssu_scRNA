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
input.path = paste0(project.path, "/data/", sample.id)
output.path = paste0(project.path, "/result/1-1_Seurat_QC/", sample.id)
qc.path = paste0 (output.path, "/QC/")
rdata.path = paste0 (output.path, "/Rdata/")
sample.id

dir.create(output.path, showWarnings = FALSE)
dir.create(qc.path, showWarnings = FALSE)
dir.create(rdata.path, showWarnings = FALSE)
################################################################################
work.path = input.path
sample.data <- Read10X(data.dir = work.path)
sample.orig <- CreateSeuratObject (counts = sample.data,
																	 min.cells = 3,
																	 min.features = 200,
																	 project = paste0 ("10X_", sample.id))

sample.orig[["percent.mt"]] <- PercentageFeatureSet(sample.orig, pattern = "^MT-|^mt-")
#-------------------------------------------------
CairoPNG (file = paste0 (qc.path, "Distribution_of_detected_genes.png"), width = 800, height = 600)
hist(sample.orig@meta.data$nFeature_RNA, breaks=100, 
		 main = paste0 ("Distribution of detected genes - ", sample.id),
		 xlab = "Genes with at least on tag")
dev.off()

CairoPNG (file = paste0 (qc.path, "Density_distribution_of_detected_genes.png"), width = 800, height = 600)
d<-density(sample.orig@meta.data$nFeature_RNA)
plot (d, main = paste0 ("Density distribution of detected genes - ", sample.id))
polygon (d, col="red", border="blue")
dev.off()

CairoPNG (file = paste0 (qc.path, "Expression_sum_per_cell.png"), width = 800, height = 600)
hist(sample.orig@meta.data$nCount_RNA, breaks=100,
		 main = paste0 ("Expression sum per cell - ", sample.id),
		 xlab = "Sum of expression")
dev.off()

CairoPNG (file = paste0 (qc.path, "Density_expression_sum_per_cell.png"), width = 800, height = 600)
d<-density(sample.orig@meta.data$nCount_RNA)
plot(d, main= paste0 ("Density expression sum per cell - ", sample.id ))
polygon (d, col="red", border="blue")
dev.off()

CairoPNG (file = paste0 (qc.path, "vlnplot.png"), width=800, height=600)
VlnPlot(sample.orig,
				features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
				ncol = 3)
dev.off()

CairoPNG (file = paste0 (qc.path, "FeatureScatter.png"), width = 800, height = 600)
plot1 <- FeatureScatter(sample.orig, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sample.orig, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

capture.output (median(sample.orig@meta.data$nCount_RNA), file = paste0 (qc.path, "median_ExpressionSum.txt"))
capture.output (median(sample.orig@meta.data$nFeature_RNA), file = paste0 (qc.path, "median_DetectedGenes.txt"))

saveRDS (sample.orig, file = paste0 (output.path, "/Rdata/", sample.id, ".QC.Rda"))
