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
				echo "Processing Sample:${sample}"
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
###Run Seurat Integration TGF-b Multi-Sample###
###Step1: seurat multi-sample QC
mkdir -p ${out_path}/2-1_Seurat_QC/
Rscript ${script_path}/2-1_seurat_merge.QC.R ALL ${project_path}
###Step2: seurat multi-sample Filt 
mkdir -p ${out_path}/2-2_Seurat_Filt/
Rscript ${script_path}/2-2_seurat_merge.Filt.R ALL ${project_path}
###Step3: seurat multi-sample Basic
mkdir -p ${out_path}/2-3_Seurat_Basic/
Rscript ${script_path}/2-3_seurat_merge.Basic.R ALL ${project_path}
END
###Step4: monocle multi-sample Basic
mkdir -p ${out_path}/3-1_Monocle_Result/
Rscript ${script_path}/3-1_monocle.Basic.R ALL ${project_path}

###Run Seurat Wnt inhibitor###
for sample in ${C59_sample_list};
do
				echo "Processing Sample:${sample}"
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
###Run Seurat Integration TGF-b Multi-Sample###
###Step1: seurat multi-sample QC
mkdir -p ${out_path}/2-1_Seurat_QC/
Rscript ${script_path}/2-1_seurat_merge.QC.R ALL_C59 ${project_path}
###Step2: seurat multi-sample Filt 
mkdir -p ${out_path}/2-2_Seurat_Filt/
Rscript ${script_path}/2-2_seurat_merge.Filt.R ALL_C59 ${project_path}
###Step3: seurat multi-sample Basic
mkdir -p ${out_path}/2-3_Seurat_Basic/
Rscript ${script_path}/2-3_seurat_merge.Basic.R ALL_C59 ${project_path}

###Step4: monocle multi-sample Basic
mkdir -p ${out_path}/3-1_Monocle_Result/
Rscript ${script_path}/3-1_monocle.Basic.R ALL_C59 ${project_path}

