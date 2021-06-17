#!/usr/local/bin/R
suppressMessages({
  library(monocle)
  library(Seurat)
  library(Matrix)
  library(DDRTree)
  library(Cairo)
  library(scales)
  library(gridExtra)
})
#Requirement--------------------------------------------------------------------
set.seed(123)
options (bitmapType='cairo')
sample.id = args[1]
in.path = args[2]
out.path = args[3]
ncom.tsne <- readRDS (file = paste0(in.path, "/Rdata/", sample.id, ".tSNE.Rda"))
#1.Make monocle data set--------------------------------------------------------
data <- as(as.matrix(ncom.tsne@assays$RNA@data), 'sparseMatrix')
pd <- new ('AnnotatedDataFrame', data = ncom.tsne@meta.data)
fData <- data.frame (gene_short_name = row.names (data),
                     row.names = row.names(data))
fd <- new ('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds.SF <- estimateSizeFactors(monocle_cds) 
monocle_cds.D <- estimateDispersions(monocle_cds.SF)

saveRDS (monocle_cds.D, 
         file = paste0(out.path, "/Rdata/", sample.id, ".monocle_cds.D.rds"))

#2.Step-------------------------------------------------------------------------
#Genes that aren't highly expressed enought will not be used for clustering, since 
#they may not give meaningful signal and would only add noise

disp_table <- dispersionTable(monocle_cds.D)
clustering_genes <- subset (disp_table, mean_expression >= 0.1)
monocle_cds.SOF <- setOrderingFilter (monocle_cds.D, clustering_genes$gene_id)

#check
tiff(file = paste0(out.path, "/plot_ordering_genes.tiff"),
     units = "in", width = 5, height = 5, res = 300)
plot_ordering_genes(monocle_cds.SOF)
dev.off()

tiff(file = paste0(out.path, "/plot_pc_variance_explained.tiff"),
     units = "in", width = 5, height = 5, res = 300)
plot_pc_variance_explained(monocle_cds.SOF, return_all = F)
dev.off()


monocle_cds.rD <- reduceDimension (monocle_cds.SOF, num_dim = 40, 
                                   reduction_method = 'tSNE',
                                   norm_method = "log")
monocle_cds.CC <- clusterCells (monocle_cds.rD, method = 'louvain')

plot_cell_clusters (monocle_cds.CC, color_by = "orig.ident")
plot_cell_clusters (monocle_cds.CC, color_by = "Cluster")

saveRDS (monocle_cds.CC, 
         file = paste0(out.path, "/Rdata/", sample.id, ".monocle_cds.CC.rds"))

#Genes likely to be informative of ordering of cells along the pseudotime trajectory
#will be selected for pseudotime inference
#3.Step-------------------------------------------------------------------------
str(monocle_cds.CC)

diff_genes <- differentialGeneTest(monocle_cds.CC,
                                   fullModelFormulaStr = "~ RNA_snn_res.0.5",
                                   cores = 18)
saveRDS (diff_genes, 
         file = paste0(out.path, "/Rdata/", sample.id, ".diff_genes_RNA_snn_res.0.5_orig.ident.rds"))
ordering_genes <- row.names(subset(diff_genes, qval < 0.01))[order(diff_genes$qval)][1:1000]
monocle_cds.CC.SOF <- setOrderingFilter(monocle_cds.CC, ordering_genes)
monocle_cds.CC.rD <- reduceDimension(monocle_cds.CC.SOF, max_components = 2, 
                                     method = 'DDRTree')
monocle_cds.CC.OC <- orderCells(monocle_cds.CC.rD)



###Rearrange Levels-------------------------------------------------------------
#monocle_cds.CC.OC <- readRDS (file = paste0(out.path, "/Rdata/", sample.id, ".monocle_cds.check.rds"))
monocle_cds.CC.OC$orig.ident <- factor (x=monocle_cds.CC.OC$orig.ident,
                                        levels = c("10X_hiPSC", "10X_Scl", "10X_Cp", 
                                                   "10X_D1", "10X_D7", "10X_D14", "10X_D28", "10X_D42"))
monocle_cds.CC.OC$orig.ident
#Idents(monocle_cds.CC.OC)

tiff(file = paste0(out.path, "/total_trjactory_sample.tiff"),
     units = "in", width = 5, height = 5, res = 300)
plot_cell_trajectory(monocle_cds.CC.OC,
                     color_by="orig.ident",
                     #color_by="Pseudotime",
                     #color_by = "seurat_clusters",
                     show_branch_points=FALSE,
                     show_backbone = FALSE,
                     show_tree =F,
                     cell_size = 0.5) +
  guides(x = "none", y = "none") +
  labs(x = NULL, y = NULL) +
  theme(#legend.position="bottom",
    legend.title = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size = unit(1,'cm'), 
    legend.key.width = unit(0.7,"cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.text = element_text (face = "bold", size = 8, family="Arial")) 
dev.off()

tiff(file = paste0(out.path, "/total_trjactory_pseudotime.tiff"),
     units = "in", width = 5, height = 5, res = 300)
plot_cell_trajectory(monocle_cds.CC.OC,
                     color_by = "Pseudotime",
                     show_branch_points=FALSE,
                     show_backbone = FALSE,
                     show_tree =F,
                     cell_size = 1) +
  guides(x = "none", y = "none") +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position="bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 12, family="Arial"),
    legend.key.width = unit(0.8,"cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.text = element_text (face = "bold", size = 8, family="Arial")) +
  guides (color = guide_colorbar (title.position = "top"))
#Here Monocle2 will first project the data to 
dev.off()
saveRDS(monocle_cds.CC.OC, file = paste0(out.path, "/Rdata/", sample.id, ".monocle_cds.CC.OC.rds"))


#-------------------------------------------------------------------------------
monocle_cds.CC.OC <- readRDS(file = paste0(out.path, "/Rdata/", sample.id, ".monocle_cds.CC.OC.rds"))


tiff(file = paste0(out.path, "/total_trjactory_chodrogenic_markers.tiff"),
     units = "in", width = 7, height = 7, res = 300)
plot_cell_trajectory(monocle_cds.CC.OC,
                     markers = c("ACAN", "COL2A1", "COMP", "SOX9"),
                     cell_size = 1,
                     use_color_gradient = TRUE,
                     show_cell_names =FALSE,
                     show_branch_points=FALSE) + 
  guides(fill=guide_legend(title="New Legend Title"), x = "none", y = "none") +
  scale_colour_gradient2 (low = muted("blue"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0,
                         guide = "colourbar",
                         na.value = "grey50") +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position="bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 12, family="Arial"),
    legend.key.width = unit(0.8,"cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.text = element_text (face = "bold", size = 8, family="Arial"))
dev.off()



tiff(file = paste0(out.path, "/total_trjactory_neurogenic_markers.tiff"),
     units = "in", width = 7, height = 7, res = 300)
plot_cell_trajectory(monocle_cds.CC.OC,
                     markers = c("NES", "OTX2", "SOX2", "WNT3A"),
                     cell_size = 1,
                     use_color_gradient = TRUE,
                     show_branch_points=FALSE) + 
  guides(fill=guide_legend(title="New Legend Title"), x = "none", y = "none") +
  scale_color_gradient2 (low = muted("blue"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0,
                         guide = "colourbar",
                         na.value = "grey50") +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position="bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 12, family="Arial"),
    legend.key.width = unit(0.8,"cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.text = element_text (face = "bold", size = 8, family="Arial"))
dev.off()

#-------------------------------------------------------------------------------
monocle_cds.CC.OC <- readRDS(file = paste0(out.path, "/Rdata/", sample.id, ".monocle_cds.CC.OC.rds"))



plot_cell_trajectory(monocle_cds.CC.OC,
                     markers = "OTX2",
                     cell_size = 1,
                     use_color_gradient = TRUE,
                     show_branch_points=FALSE) + 
  guides(fill=guide_legend(title="New Legend Title"), x = "none", y = "none") +
  scale_color_gradient2 (low = muted("blue"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0,
                         guide = "colourbar",
                         na.value = "grey50") +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position="none",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 12, family="Arial"),
    legend.key.width = unit(0.8,"cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.text = element_text (face = "bold", size = 8, family="Arial"))

out.path
for (i in Fig5_D_gene_list){
  print (i)
  plot_customized(monocle_cds.CC.OC, i, out.path)
}

#Figure------------------------------------------------------------------------
#
ncom.tsne$orig.ident <- factor (x= ncom.tsne$orig.ident,
                                levels = c("10X_hiPSC", "10X_Scl", "10X_Cp", 
                                  "10X_D1", "10X_D7", "10X_D14", "10X_D28", "10X_D42"))
ncom.tsne$seurat_clusters  <- factor (x = ncom.tsne$seurat_clusters,
                                      levels = as.numeric(levels(ncom.tsne$seurat_clusters)))

DimPlot(ncom.tsne)


p1 <- DimPlot(ncom.tsne) +
  labs (x=NULL, y= NULL) + 
  guides(x = "none", y = "none") +
  theme (axis.line.x = element_blank(),
         axis.line.y = element_blank(),
         legend.text = element_text (face = "bold", size = 9, family = "Arial"),
         legend.position = "left") 

p2 <- dittoDimPlot(ncom.tsne, "orig.ident",
             main = NULL,
             theme = theme_void()) + 
  labs (fill = "", x=NULL, y= NULL) + 
  theme (legend.position = "none") 


p3 <- dittoBarPlot(ncom.tsne, 
             "orig.ident", group.by = "orig.ident",
             scale = "count",
             var.labels.reorder = c(7,8,1,2,6,3,4,5),
             x.reorder = c(7,8,1,2,6,3,4,5),
             main = NULL) +
  ylab ("Count of each samples") +
  scale_y_continuous (
    expand = expansion(mult = c(0,0.1))) +
  scale_x_discrete (expand = expansion(mult = c(0, 0.1))) +
  theme (axis.title.x = element_blank(),
         axis.line.x = element_blank(),
         axis.text.x = element_text(face = "bold", size = 9, family = "Arial"),
         
         axis.ticks.y = element_line(),
         axis.line.y = element_blank(),
         axis.text.y = element_text (face = "bold", size = 9, family = "Arial"),
         axis.title.y = element_text(face = "bold", size = 12, family = "Arial"),
         
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line (colour = "black"),
         
         legend.text = element_text (face = "bold", size = 9, family = "Arial"),
         legend.position = "none") 

p4 <- dittoBarPlot(ncom.tsne, 
             var.labels.reorder = c(7, 8, 1, 2, 6, 3, 4, 5), 
             x.reorder = c(1,2,7,8,9,10,11,12,13,14,3,4,5,6),
             "orig.ident", group.by = "seurat_clusters", scale = "percent",
             main = NULL) +
  ylab ("Proportion of each clsuters") +
  labs (fill = NULL) +
  
  scale_y_continuous (labels = scales::percent,
                      expand = expansion(mult = c(0,0.1))) +
  scale_x_discrete (expand = expansion(mult = c(0, 0.1))) +
  theme (axis.title.x = element_blank(),
         axis.line.x = element_blank(),
         axis.text.x = element_text(face = "bold", size = 9, family = "Arial"),
         
         axis.ticks.y = element_line(),
         axis.line.y = element_blank(),
         axis.text.y = element_text (face = "bold", size = 9, family = "Arial"),
         axis.title.y = element_text(face = "bold", size = 12, family = "Arial"),
         
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_line (colour = "black"),
         
         legend.text = element_text (face = "bold", size = 6, family = "Arial"),
         legend.key.size = unit(0.3,'cm'), 
         legend.key.width = unit(0.3,"cm"),
         legend.key.height = unit(0.3, "cm"),
         legend.position = "bottom") 


tiff(file = paste0(out.path, "/total_summary_figure.tiff"),
     units = "in", width = 14, height = 6, res = 300)
grid.arrange(p1, p2, p3, p4, ncol=4)
dev.off()
#function-----------------------------------------------------------------------
plot_customized <- function (seurat, gene, output_path){
  tiff(file = paste0(output_path, "/total_trjactory_", gene, ".tiff"), 
       units = "in", width = 7, height = 7, res = 300)
  print(plot_cell_trajectory(seurat,
                       markers = gene,
                       cell_size = 1,
                       use_color_gradient = TRUE,
                       show_branch_points=FALSE) + 
    guides(fill=guide_legend(title="New Legend Title"), x = "none", y = "none") +
    scale_color_gradient2 (low = muted("blue"),
                           mid = "white",
                           high = muted("red"),
                           midpoint = 0,
                           guide = "colourbar",
                           na.value = "grey50") +
    labs(x = NULL, y = NULL) +
    theme(
      legend.position="none",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold", size = 12, family="Arial"),
      legend.key.width = unit(0.8,"cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.text = element_text (face = "bold", size = 8, family="Arial"))
  )
  dev.off()
}


