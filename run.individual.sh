#!/bin/bash
#############################
###Requirement Information###
#############################
project_path=$(pwd)
script_path="${project_path}/script/"
out_path="${project_path}/result/"
sample_list=("hiPSC Scl Cp D1 D7 D14 D28 D42")
C59_sample_list=("C59_D14 C59_D28 C59_D42 C59_D7")


###Run Seurat TGF-b Sample### 
for sample in ${sample_list};
do
				echo "${sample}"
				###Step1: seurat QC
				echo "###Step1: Seurat QC###"
				mkdir -p ${out_path}/1-1_Seurat_QC/
				Rscript ${script_path}/1-1_seurat_sample.QC.R ${sample} ${project_path}
				###Step2: seurat Filt
				echo "###Step2: Seurat Filt###"
				mkdir -p ${out_path}/1-2_Seurat_Filt/
				Rscript ${script_path}/1-2_seurat_sample.Filt.R ${sample} ${project_path}
				###Step3: seurat Basic
				echo "###Step3: Seurat Basic###"
				mkdir -p ${out_path}/1-3_Seurat_Basic/
				Rscript ${script_path}/1-3_seurat_sample.Basic.R ${sample} ${project_path}
done

###Run Seurat Wnt inhibitor###
for sample in ${C59_sample_list};
do
				echo "${sample}"
				###Step1: seurat QC
				echo "###Step1: Seurat QC###"
				mkdir -p ${out_path}/1-1_Seurat_QC/
				Rscript ${script_path}/1-1_seurat_sample.QC.R ${sample} ${project_path}
				###Step2: seurat Filt
				echo "###Step2: Seurat Filt###"
				mkdir -p ${out_path}/1-2_Seurat_Filt/
				Rscript ${script_path}/1-2_seurat_sample.Filt.R ${sample} ${project_path}
				###Step3: seurat Basic
				echo "###Step3: Seurat Basic###"
				mkdir -p ${out_path}/1-3_Seurat_Basic/
				Rscript ${script_path}/1-3_seurat_sample.Basic.R ${sample} ${project_path}
done
