# ssu_scRNA
## 숭실대 단일세포유전체분석 강의 
### 1. Contents
```bash
./ssu_scRNA
├── README.md
├── data #cellranger mtx 파일 
├── result #결과들이 쌓일 폴더
├── run.integration.sh #integration 실행 bash 파일
├── run.individual.sh #individual 실행 bash 파일
└── script #실행 R script 파일
```

### 2. Sample Information
 Sample 은 총 12샘플이고, 카테고리는 Progenitor, Original Protocol, C59(Wnt inhibitor)처리한 Protocol로 나눌 수 있다.
```data
./data
#Progenitor
├── hiPSC 
├── Scl #Sclerotome
├── Cp #Chondroprogenitor
#Original Protocol
├── D1 #time point:day1
├── D7 #time point:day7
├── D14 #time point:day14
├── D28 #time point:day28
├── D42 #time point:day42
#C59 (Wnt inhibitor)
├── C59_D7 #time point:day7
├── C59_D14 #time point:day14
├── C59_D28 #time point:day28
└── C59_D42 #time point:day42
```

### 3. Script Information
 프로세싱은 2가지 방법으로 진행된다. 첫번째는 individual processing 두번째는 integration processing이다.
integration processing 은 TGF-beta (Original Protocol) 처리만 한 샘플과 TGF-beta + C59 (Wnt inhibition) 처리한 샘플들을 integration하여 분석 진행했다. 
```bash
./script
#1.individual processing
├── 1-1_seurat_sample.QC.R
├── 1-2_seurat_sample.Filt.R
├── 1-3_seurat_sample.Basic.R
#2.integration processing
├── 2-1_seurat_merge.C59.QC.R
├── 2-1_seurat_merge.QC.R
├── 2-2_seurat_merge.Filt.R
├── 2-3_seurat_merge.Basic.R
└── 3-1_monocle.Basic.R
```

### 4. Run
 실행 파일은 2 프로세스로 진행되고, result 폴더에는 figure와 txt 파일이 저장된다. run.integration.sh 를 실행할때는 샘플의 individual QC 결과(bash run.individual.sh) 가 있어야 한다. 
1. run.individual.sh: 12 sample 의 individual sample 의 seurat process
2. run.integration.sh: TGF-b, Wnt inhibitor samples 의 seurat, monocle process
```bash
$ git clone https://github.com/baekip/ssu_scRNA.git
$ cd ssu_scRNA
$ bash run.individual.sh #individual 12 sample 실행
$ bash run.integration.sh #integration 실행
```

### 5. Result
 Result도 individual 과 integration 으로 나눠서 저장된다. Seurat 결과는 3 step으로 나눠 QC, Filt, Basic 의 폴더에 결과와 Rda 파일이 각각 저장이 되며, monocle 결과는 한 폴더에 저장된다.
```
./result
├── 1-1_Seurat_QC #individual 12 sample Seurat QC result
├── 1-2_Seurat_Filt #individual 12 sample Seurat Filt result
├── 1-3_Seurat_Basic #individual 12 sample Seurat Basic result
├── 2-1_Seurat_QC #integration Seurat QC result
├── 2-2_Seurat_Filt #integration Seurat Filt result
├── 2-3_Seurat_Basic #integration Seurat Basic result
└── 3-1_Monocle_Result #integration Monocle result
```

