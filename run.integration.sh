#!/bin/bash
#############################
###Requirement Information###
#############################
project_path=$(pwd)
script_path="${project_path}/script/"
out_path="${project_path}/result/"

###Run Seurat Integration TGF-b Multi-Sample###
#Step1: seurat multi-sample QC
mkdir -p ${out_path}/2-1_Seurat_QC/
Rscript ${script_path}/2-1_seurat_merge.QC.R ALL ${project_path}
#Step2: seurat multi-sample Filt 
mkdir -p ${out_path}/2-2_Seurat_Filt/
Rscript ${script_path}/2-2_seurat_merge.Filt.R ALL ${project_path}
#Step3: seurat multi-sample Basic
mkdir -p ${out_path}/2-3_Seurat_Basic/
Rscript ${script_path}/2-3_seurat_merge.Basic.R ALL ${project_path}
#Step4: monocle multi-sample Basic
mkdir -p ${out_path}/3-1_Monocle_Result/
Rscript ${script_path}/3-1_monocle.Basic.R ALL ${project_path}

###Run Seurat Integration TGF-b Multi-Sample###
#Step1: seurat multi-sample QC
mkdir -p ${out_path}/2-1_Seurat_QC/
Rscript ${script_path}/2-1_seurat_merge.QC.R ALL_C59 ${project_path}
#Step2: seurat multi-sample Filt 
mkdir -p ${out_path}/2-2_Seurat_Filt/
Rscript ${script_path}/2-2_seurat_merge.Filt.R ALL_C59 ${project_path}
#Step3: seurat multi-sample Basic
mkdir -p ${out_path}/2-3_Seurat_Basic/
Rscript ${script_path}/2-3_seurat_merge.Basic.R ALL_C59 ${project_path}
#Step4: monocle multi-sample Basic
mkdir -p ${out_path}/3-1_Monocle_Result/
Rscript ${script_path}/3-1_monocle.Basic.R ALL_C59 ${project_path}

