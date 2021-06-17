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
input.path = args[2]
output.path = args[3]

qc.path = paste0 (output.path, "/QC/")

###load QC.Rda
hiPSC = readRDS (file = paste0 (input.path, "/hiPSC/Rdata/hiPSC.QC.Rda"))
Scl = readRDS (file = paste0 (input.path, "/Scl/Rdata/Scl.QC.Rda"))
Cp = readRDS (file = paste0 (input.path, "/Cp/Rdata/Cp.QC.Rda"))
D1 = readRDS (file = paste0 (input.path, "/D1/Rdata/D1.QC.Rda"))
D7 = readRDS (file = paste0 (input.path, "/D7/Rdata/D7.QC.Rda"))
D14 = readRDS (file = paste0 (input.path, "/D14/Rdata/D14.QC.Rda"))
D28 = readRDS (file = paste0 (input.path, "/D28/Rdata/D28.QC.Rda"))
D42 = readRDS (file = paste0 (input.path, "/D42/Rdata/D42.QC.Rda"))

sample.merge <- merge (hiPSC, y = c(Scl, Cp, D1, D7, D14, D28, D42),
											 add.cell.ids = c("hiPSC", "Scl", "Cp", "D1", "D7", "D14", "D28", "D42"),
											 project = "NCom")

###make QC Figure
CairoPNG (file = paste0 (qc.path, "Distribution_of_detected_genes.png"), width = 800, height = 600)
hist(sample.merge@meta.data$nFeature_RNA, breaks=100, 
		 main = paste0 ("Distribution of detected genes - ", sample.id),
		 xlab = "Genes with at least on tag")
dev.off()

CairoPNG (file = paste0 (qc.path, "Density_distribution_of_detected_genes.png"), width = 800, height = 600)
d<-density(sample.merge@meta.data$nFeature_RNA)
plot (d, main = paste0 ("Density distribution of detected genes - ", sample.id))
polygon (d, col="red", border="blue")
dev.off()

CairoPNG (file = paste0 (qc.path, "Expression_sum_per_cell.png"), width = 800, height = 600)
hist(sample.merge@meta.data$nCount_RNA, breaks=100,
		 main = paste0 ("Expression sum per cell - ", sample.id),
		 xlab = "Sum of expression")
dev.off()

CairoPNG (file = paste0 (qc.path, "Density_expression_sum_per_cell.png"), width = 800, height = 600)
d<-density(sample.merge@meta.data$nCount_RNA)
plot(d, main= paste0 ("Density expression sum per cell - ", sample.id ))
polygon (d, col="red", border="blue")
dev.off()

CairoPNG (file = paste0 (qc.path, "vlnplot.png"), width=800, height=600)
VlnPlot(sample.merge,
				features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
				ncol = 3)
dev.off()

CairoPNG (file = paste0 (qc.path, "FeatureScatter.png"), width = 800, height = 600)
plot1 <- FeatureScatter(sample.merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sample.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
#-------------------------------------------------
saveRDS (sample.merge, file = paste0 (output.path, "/Rdata/Merge.Rda"))

