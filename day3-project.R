#：step 1：准备相关的数据
#：1）参考基因组，下载完成后，解压
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
#：mouse
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-and-mm10-2020-A.tar.gz
#：2）下载相应的SRA数据，其中前者为N0，后者为N1，
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR16796879/SRR16796879
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR16796880/SRR16796880
#：3）SRA转fastq
fastq-dump --split-files --origfmt --gzip SRR16796879
fastq-dump --split-files --origfmt --gzip SRR16796880
#：查看CB和UMI序列
les /data2/kkwan/ESCC_WGBS_RNA/ESCC_LNM/scRNA_data/00.fastq/SRR16796879/SRR16796879_1.fastq.gz |awk '/TTTTTTTT/{print substr($0,1,27)}' |les
#：4）准备白名单序列文件 10X v2和V3不太一样
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz

#：step 2：STAR solo计算表达矩阵
#：1）准备STAR solo相关的文件，包括参考基因组、参考注释文件等
/data2/kkwan/software/STAR  --runMode genomeGenerate --runThreadN 24  --genomeDir /data2/kkwan/software/cellranger_data/star_solo_data --genomeFastaFiles /data2/kkwan/software/cellranger_data/refdata-gex-GRCh38-2020-A/fasta/genome.fa  --sjdbGTFfile /data2/kkwan/software/cellranger_data/refdata-gex-GRCh38-2020-A/genes/genes.gtf --genomeSAsparseD 3
#：2）使用STAR进行序列比对，要求输入参考基因组序列、参考的注释基因，reads文件、白名单文件
/data2/kkwan/software/STAR --genomeDir /data2/kkwan/software/cellranger_data/star_solo_data/ --readFilesIn /data2/kkwan/ESCC_WGBS_RNA/ESCC_LNM/scRNA_data/00.fastq/SRR16796879/SRR16796879_1.fastq.gz /data2/kkwan/ESCC_WGBS_RNA/ESCC_LNM/scRNA_data/00.fastq/SRR16796879/SRR16796879_2.fastq.gz --soloType CB_UMI_Simple --soloBarcodeMate 1 --clip3pNbases 26 0 --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --readFilesCommand zcat --soloCBwhitelist /data2/kkwan/software/cellranger_data/cellbarcode_whitelist/737K-august-2016.txt
#-
/data2/kkwan/software/STAR --genomeDir /data2/kkwan/software/cellranger_data/star_solo_data/ --readFilesIn /data2/kkwan/ESCC_WGBS_RNA/ESCC_LNM/scRNA_data/00.fastq/SRR16796880/SRR16796880_1.fastq.gz /data2/kkwan/ESCC_WGBS_RNA/ESCC_LNM/scRNA_data/00.fastq/SRR16796880/SRR16796880_2.fastq.gz --soloType CB_UMI_Simple --soloBarcodeMate 1 --clip3pNbases 26 0 --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --readFilesCommand zcat --soloCBwhitelist /data2/kkwan/software/cellranger_data/cellbarcode_whitelist/737K-august-2016.txt
#：step 3：整理STAR solo输出的结果，可以直接被读入到 seurat 对象
library(Seurat) #BiocManager::install("Seurat")
library(ggplot2) # install.packages("ggplot2")
library(SingleR) # BiocManager::install("SingleR")
library(dplyr) #: install.packages("dplyr")
library(celldex) #:BiocManager::install("celldex")
library(RColorBrewer) #: install.packages("RColorBrewer")
library(SingleCellExperiment) #:BiocManager::install("SingleCellExperiment")

#-------------
setwd("I:/单细胞培训/相关数据/练习题/day3");

#:------------------------------
#:------------------------------
#：1、读入数据
# Load the PBMC dataset
Read10X(data.dir = "scRNA-data/SRR16796863-matrix")->SRR16796863_data;
Read10X(data.dir = "scRNA-data/SRR16796864-matrix")->SRR16796864_data;
Read10X(data.dir = "scRNA-data/SRR16796883-matrix")->SRR16796883_data;
Read10X(data.dir = "scRNA-data/SRR16796884-matrix")->SRR16796884_data;

# ：创建 seurat 对象
CreateSeuratObject(counts = SRR16796863_data, project = "ESCC_N0", min.cells = 3, min.features = 200)->SRR16796863_seurat_obj;
CreateSeuratObject(counts = SRR16796864_data, project = "ESCC_N0", min.cells = 3, min.features = 200)->SRR16796864_seurat_obj;
CreateSeuratObject(counts = SRR16796883_data, project = "ESCC_N1", min.cells = 3, min.features = 200)->SRR16796883_seurat_obj;
CreateSeuratObject(counts = SRR16796884_data, project = "ESCC_N1", min.cells = 3, min.features = 200)->SRR16796884_seurat_obj;
#：添加批次信息
"N0_r1"->SRR16796863_seurat_obj$batch;
"N0_r2"->SRR16796864_seurat_obj$batch;
"N1_r1"->SRR16796883_seurat_obj$batch;
"N1_r2"->SRR16796884_seurat_obj$batch;

#------------------------------------
#:------------------------------
#：2、基本QC
#：a. 计算线粒体基因比例
PercentageFeatureSet(SRR16796863_seurat_obj, pattern = "^MT-")->SRR16796863_seurat_obj[["percent.mt"]];
PercentageFeatureSet(SRR16796883_seurat_obj, pattern = "^MT-")->SRR16796883_seurat_obj[["percent.mt"]];
#：b. 计算核糖体基因比例
PercentageFeatureSet(SRR16796863_seurat_obj, pattern = "^RP[SL]")->SRR16796863_seurat_obj[["percent.rb"]]
PercentageFeatureSet(SRR16796883_seurat_obj, pattern = "^RP[SL]")->SRR16796883_seurat_obj[["percent.rb"]];
#：以上得到的两个数值均存放在meta.data的切片中
"Failed"->SRR16796863_seurat_obj$QC_type;
which(SRR16796863_seurat_obj$nFeature_RNA > 200 & SRR16796863_seurat_obj$nFeature_RNA < 2500 & SRR16796863_seurat_obj$percent.mt < 5)->qc_pass_index;
"Pass"->SRR16796863_seurat_obj$QC_type[qc_pass_index];
#
"Failed"->SRR16796883_seurat_obj$QC_type;
which(SRR16796883_seurat_obj$nFeature_RNA > 200 & SRR16796883_seurat_obj$nFeature_RNA < 2500 & SRR16796883_seurat_obj$percent.mt < 5)->qc_pass_index;
"Pass"->SRR16796883_seurat_obj$QC_type[qc_pass_index];
#:---过滤低质量的基因
rowMeans(SRR16796863_seurat_obj)->x;
SRR16796863_seurat_obj[which(x>=0.1),]->SRR16796863_seurat_obj;
rowMeans(SRR16796883_seurat_obj)->x;
SRR16796883_seurat_obj[which(x>=0.1),]->SRR16796883_seurat_obj;

##################################################################################
#:合并seurat对象
merge(subset(SRR16796863_seurat_obj,QC_type=="Pass"),subset(SRR16796883_seurat_obj,QC_type=="Pass"))->ESCC_N0N1_seurat_obj;
#：未做整合之前的结果
NormalizeData(ESCC_N0N1_seurat_obj)->ESCC_N0N1_seurat_obj;
FindVariableFeatures(ESCC_N0N1_seurat_obj, selection.method = "vst", nfeatures = 2000)->ESCC_N0N1_seurat_obj
ScaleData(ESCC_N0N1_seurat_obj, verbose = FALSE)->ESCC_N0N1_seurat_obj
RunPCA(ESCC_N0N1_seurat_obj, npcs = 30, verbose = FALSE)->ESCC_N0N1_seurat_obj;
RunUMAP(ESCC_N0N1_seurat_obj, reduction = "pca", dims = 1:30)->ESCC_N0N1_seurat_obj
FindNeighbors(ESCC_N0N1_seurat_obj, reduction = "pca", dims = 1:30)->ESCC_N0N1_seurat_obj
FindClusters(ESCC_N0N1_seurat_obj, resolution = 0.5)->ESCC_N0N1_seurat_obj
#：展示
DimPlot(ESCC_N0N1_seurat_obj, reduction = "umap", group.by = "batch")
#--------------------------------------------------------------
#--------------------------------------------------------------
#：去除condition效应
#：按照条件来分割seurat对象
SplitObject(ESCC_N0N1_seurat_obj, split.by = "batch")->ESCC_N0N1_list;

#step 1. 对每个数据集进行均一化，并识别显著变化的前2000个基因
lapply(X = ESCC_N0N1_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})->ESCC_N0N1_list;
#step 2. 挑选数据集之间共有的高度变化的基因
SelectIntegrationFeatures(object.list = ESCC_N0N1_list)->ESCC_N0N1_seleted_features;
#step 3. 然后，我们使用FindIntegrationAnchors()函数识别锚点，该函数将Seurat对象的列表作为输入，
FindIntegrationAnchors(object.list = ESCC_N0N1_list, anchor.features = ESCC_N0N1_seleted_features)->ESCC_N0N1_anchors;
#step 4. 使用这些锚点将两个数据集整合在一起。
IntegrateData(anchorset = ESCC_N0N1_anchors)->ESCC_N0N1_combined
#-----
#：对合并后的数据集进行标准分析
ESCC_N0N1_combined@assays # 查看有哪些assay
"integrated"->DefaultAssay(ESCC_N0N1_combined);
# Run the standard workflow for visualization and clustering
ScaleData(ESCC_N0N1_combined, verbose = FALSE)->ESCC_N0N1_combined
RunPCA(ESCC_N0N1_combined, npcs = 30, verbose = FALSE)->ESCC_N0N1_combined;
RunUMAP(ESCC_N0N1_combined, reduction = "pca", dims = 1:30)->ESCC_N0N1_combined
FindNeighbors(ESCC_N0N1_combined, reduction = "pca", dims = 1:30)->ESCC_N0N1_combined
FindClusters(ESCC_N0N1_combined, resolution = 0.5)->ESCC_N0N1_combined
#：展示
DimPlot(ESCC_N0N1_combined, reduction = "umap", group.by = "batch")->p1_combined;
DimPlot(ESCC_N0N1_combined, reduction = "umap", label = TRUE, repel = TRUE)->p2_combined;
p1_combined + p2_combined
#------------------------------------------------------
#-------细胞类型注释：
celldex::MonacoImmuneData()->celldex_ref;
as.SingleCellExperiment(DietSeurat(ESCC_N0N1_combined))->ESCC_N0N1_combined_sc_obj;
#：
SingleR(test = ESCC_N0N1_combined_sc_obj,assay.type.test = 1,ref = celldex_ref,labels = celldex_ref$label.main)->ESCC_N0N1_sc_celldex_main
SingleR(test = ESCC_N0N1_combined_sc_obj,assay.type.test = 1,ref = celldex_ref,labels = celldex_ref$label.fine)->ESCC_N0N1_sc_celldex_fine;
# 将注释的细胞类型添加到 seurat 对象中
ESCC_N0N1_sc_celldex_main$pruned.labels->pruned_labels;
rownames(ESCC_N0N1_sc_celldex_main)->names(pruned_labels);
AddMetaData(ESCC_N0N1_combined,pruned_labels,col.name="Monaco_main_pruned_labels")->ESCC_N0N1_combined;
#：另一种添加注释的方法
ESCC_N0N1_sc_celldex_fine$pruned.labels->ESCC_N0N1_combined$Monaco_fine_pruned_labels
#:-----------------------
#：这里同时使用HumanPrimaryCellAtlasData来注释
HumanPrimaryCellAtlasData()->celldex_ref2;
#：
SingleR(test = ESCC_N0N1_combined_sc_obj,assay.type.test = 1,ref = celldex_ref2,labels = celldex_ref2$label.main)->ESCC_N0N1_sc_celldex_main2
SingleR(test = ESCC_N0N1_combined_sc_obj,assay.type.test = 1,ref = celldex_ref2,labels = celldex_ref2$label.fine)->ESCC_N0N1_sc_celldex_fine2;
ESCC_N0N1_sc_celldex_main2$pruned.labels->ESCC_N0N1_combined$HPC_main_pruned_labels;
ESCC_N0N1_sc_celldex_fine2$pruned.labels->ESCC_N0N1_combined$HPC_fine_pruned_labels;

#-------------------
#：画图展示
#：查看细胞的identify信息，DimPlot画图需要指定
Idents(ESCC_N0N1_combined);
SetIdent(ESCC_N0N1_combined, value = "HPC_main_pruned_labels")->ESCC_N0N1_combined;
DimPlot(ESCC_N0N1_combined, label = T , repel = T, label.size = 3) + NoLegend()->p1_annotated;
#
SetIdent(ESCC_N0N1_combined, value = "HPC_fine_pruned_labels")->ESCC_N0N1_combined;
DimPlot(ESCC_N0N1_combined, label = T , repel = T, label.size = 3) + NoLegend()->p2_annotated;
#
SetIdent(ESCC_N0N1_combined, value = "Monaco_main_pruned_labels")->ESCC_N0N1_combined;
DimPlot(ESCC_N0N1_combined, label = T , repel = T, label.size = 3) + NoLegend()->p3_annotated;
#
SetIdent(ESCC_N0N1_combined, value = "Monaco_fine_pruned_labels")->ESCC_N0N1_combined;
DimPlot(ESCC_N0N1_combined, label = T , repel = T, label.size = 3) + NoLegend()->p4_annotated;
#---比较N0和N1之间细胞群体差异
table(ESCC_N0N1_combined$batch,ESCC_N0N1_combined$integrated_snn_res.0.5)
table(ESCC_N0N1_combined$Monaco_main_pruned_labels,ESCC_N0N1_combined$batch)
table(ESCC_N0N1_combined$Monaco_fine_pruned_labels,ESCC_N0N1_combined$batch)

#：比较注释的结果
table(ESCC_N0N1_combined$integrated_snn_res.0.5,ESCC_N0N1_combined$Monaco_main_pruned_labels)
table(ESCC_N0N1_combined$integrated_snn_res.0.5,ESCC_N0N1_combined$Monaco_fine_pruned_labels)
#----热图展示不同类群之间的注释情况
pheatmap(table(ESCC_N0N1_combined$integrated_snn_res.0.5,ESCC_N0N1_combined$batch),cluster_rows=F,cluster_cols=F,scale="row")
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(250),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(250))->fill_colors;
pheatmap(table(ESCC_N0N1_combined$integrated_snn_res.0.5,ESCC_N0N1_combined$HPC_fine_pruned_labels),cluster_rows=F,cluster_cols=F,scale="row",col=fill_colors)


#：展示特定的基因
as.data.frame(top10_markers);
FeaturePlot(ESCC_N0N1_combined, features = c("GP9", "SPARC"))+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

#:------------------------------------------------------------
#：8、差异表达分析
# 首先设置差异表达分析的组别：以聚类的结果来定
SetIdent(ESCC_N0N1_combined,value="integrated_snn_res.0.5")->ESCC_N0N1_combined;
FindAllMarkers(ESCC_N0N1_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)->ESCC_N0N1_all_markers;
ESCC_N0N1_all_markers %>%group_by(cluster) %>%slice_max(n = 2, order_by = avg_log2FC)
#：也可以指定使用哪个slot的矩阵来展示
#VlnPlot(pbmc_seurat_obj, features = c("RPS12", "FCGR3A"), slot = "counts", log = TRUE)# data,scale.data
#：聚类图上展示
FeaturePlot(pbmc_seurat_obj, features = c("RPS12", "FCGR3A"))
#：也可以通过热图展示
ESCC_N0N1_all_markers%>%group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC) -> top10_markers;
DoHeatmap(ESCC_N0N1_combined, features = top10_markers$gene) + NoLegend()



























