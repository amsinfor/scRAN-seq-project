#:
library(Seurat)
library(SeuratData)# devtools::install_github('satijalab/seurat-data')
library(patchwork)


#：设置工作路径
setwd("F:/projects-work/单细胞培训/相关数据/练习题/day2")


#;--------准备相关的数据
#：查看seurat提供的数据集有哪些
AvailableData()
#：安装 ifnb 数据
InstallData("ifnb")
#：载入数据
LoadData("ifnb")
#：查看ifnb数据集（实际上为seurat对象）的信息
table(ifnb$stim)
#：按照条件来分割seurat对象
SplitObject(ifnb, split.by = "stim")->ifnb_list;

#step 1. 对每个数据集进行均一化，并识别显著变化的前2000个基因
lapply(X = ifnb_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})->ifnb_list;
#step 2. 挑选数据集之间共有的高度变化的基因
SelectIntegrationFeatures(object.list = ifnb_list)->ifnb_seleted_features;
#step 3. 然后，我们使用FindIntegrationAnchors()函数识别锚点，该函数将Seurat对象的列表作为输入，
FindIntegrationAnchors(object.list = ifnb_list, anchor.features = ifnb_seleted_features)->ifnb_anchors;
#step 4. 使用这些锚点将两个数据集整合在一起。
IntegrateData(anchorset = ifnb_anchors)->ifnb_combined
#-----
#：对合并后的数据集进行标准分析
ifnb_combined@assays # 查看有哪些assay
"integrated"->DefaultAssay(ifnb_combined);
# Run the standard workflow for visualization and clustering
ScaleData(ifnb_combined, verbose = FALSE)->ifnb_combined
RunPCA(ifnb_combined, npcs = 30, verbose = FALSE)->ifnb_combined;
RunUMAP(ifnb_combined, reduction = "pca", dims = 1:30)->ifnb_combined
FindNeighbors(ifnb_combined, reduction = "pca", dims = 1:30)->ifnb_combined
FindClusters(ifnb_combined, resolution = 0.5)->ifnb_combined
#：展示
DimPlot(ifnb_combined, reduction = "umap", group.by = "stim")->p1_combined;
DimPlot(ifnb_combined, reduction = "umap", label = TRUE, repel = TRUE)->p2_combined;
p1_combined + p2_combined


#------------
#：作为对照，对未合并的数据集使用相同的标准分析流程并展示
#：设置默认的assay
"RNA"->DefaultAssay(ifnb_combined);
#：均一化
NormalizeData(ifnb_combined, verbose = FALSE)->ifnb_combined
#：识别高度变化的基因
FindVariableFeatures(ifnb_combined, selection.method = "vst", nfeatures = 2000, verbose = F)->ifnb_combined;
#：降维之前的数据归一化
ScaleData(ifnb_combined, verbose = F)->ifnb_combined;
#：PCA降维，该结果用于聚类使用
RunPCA(ifnb_combined, npcs = 30, verbose = FALSE)->ifnb_combined;
#：UMAP降维，该结果用于展示
RunUMAP(ifnb_combined, reduction = "pca", dims = 1:30)->ifnb_combined
#：执行聚类，seurat基于图的聚类算法
FindNeighbors(ifnb_combined, reduction = "pca", dims = 1:30)->ifnb_combined
FindClusters(ifnb_combined, resolution = 0.5)->ifnb_combined
#：展示
DimPlot(ifnb_combined, reduction = "umap", group.by = "stim")->p1_uncombined;
DimPlot(ifnb_combined, reduction = "umap", label = TRUE, repel = TRUE)->p2_uncombined;
p1_uncombined + p2_uncombined




#----------------------------------------------
#--使用 hamony
library(harmony)#install.packages("harmony")
#：标准分析
ifnb->ifnb_harmony
NormalizeData(ifnb_harmony, verbose = F)->ifnb_harmony;
FindVariableFeatures(ifnb_harmony, selection.method = "vst", nfeatures = 2000, verbose = F)->ifnb_harmony;
ScaleData(ifnb_harmony, verbose = F)->ifnb_harmony;
RunPCA(ifnb_harmony, npcs = 30, verbose = F)->ifnb_harmony;
RunUMAP(ifnb_harmony, reduction = "pca", dims = 1:30, verbose = F)->ifnb_harmony;
#：执行harmony，harmony的运行需要PCA的结果
#：运行该过程，需要指定哪一变量为待消除的因素，比如 condition、batch、source等
RunHarmony(ifnb_harmony, "stim")->ifnb_harmony;
#: 查看harmony的结果
Embeddings(ifnb_harmony, 'harmony')->ifnb_harmony_embeddings;
#：使用harmony结果进行展示
RunUMAP(ifnb_harmony, reduction = "harmony",dims = 1:30, verbose = F)->ifnb_harmony;
#: 执行聚类
FindNeighbors(ifnb_harmony,reduction = "harmony", k.param = 10, dims = 1:30)->ifnb_harmony;
FindClusters(ifnb_harmony,resolution=0.5)->ifnb_harmony; 
DimPlot(ifnb_harmony, reduction = "umap", group.by = "stim")->p1_harmony;
#: 展示聚类的结果
SetIdent(ifnb_harmony,value = "RNA_snn_res.0.5")->ifnb_harmony;
DimPlot(ifnb_harmony, reduction = "umap", label = TRUE, repel = TRUE)->p2_harmony;
#：

#-------------------------------------------
#-：使用LIGER
library(rliger)#install.packages('rliger')
library(SeuratWrappers)#remotes::install_github('satijalab/seurat-wrappers') # 需要在4.2以上版本来做 
#：1、系统性比较人类PBMC的细胞，这些细胞来自不同方法的数据集，将它们进行整合后可以获得完整的PBMC细胞类型景观
#：step 1：准备相关数据
InstallData("pbmcsca")
data("pbmcsca")
# 注意要求`liger` 版本高于 0.5.0
pbmcsca->pbmcsca_rliger
#：查看哪些方法
table(pbmcsca_rliger$Method)

#：step 2：标准分析
NormalizeData(pbmcsca_rliger)->pbmcsca_rliger
FindVariableFeatures(pbmcsca_rliger)->pbmcsca_rliger;
#：归一化用于降维，这里使用的是非中心化方法（rliger使用NMF方法，需要非负值输入）
ScaleData(pbmcsca_rliger, split.by = "Method", do.center = FALSE)->pbmcsca_rliger;
#：step 3：执行rliger，指定消除的变量为 Method
RunOptimizeALS(pbmcsca_rliger, k = 20, lambda = 5, split.by = "Method")->pbmcsca_rliger;
RunQuantileNorm(pbmcsca_rliger, split.by = "Method")->pbmcsca_rliger;
# 
#：step 4：聚类
FindNeighbors(pbmcsca_rliger, reduction = "iNMF", dims = 1:20)->pbmcsca_rliger;
FindClusters(pbmcsca, resolution = 0.3)->pbmcsca_rliger;
# 展示
RunUMAP(pbmcsca_rliger, dims = 1:ncol(pbmcsca_rliger[["iNMF"]]), reduction = "iNMF")->pbmcsca_rliger;
DimPlot(pbmcsca_rliger, group.by = c("Method", "ident", "CellType"), ncol = 3)

#---
#：2、整合不同条件下的单细胞数据，使用stimulated组和control组的PBMC数据集
# 这里使用ifnb数据，
#：InstallData("ifnb")
LoadData("ifnb")->ifnb_rliger
# step 1：标准分析
NormalizeData(ifnb_rliger)->ifnb_rliger;
FindVariableFeatures(ifnb_rliger)->ifnb_rliger;
#：归一化，
ScaleData(ifnb_rliger, split.by = "stim", do.center = FALSE)->ifnb_rliger;
#：step 2：执行liger，指定消除的变量为stim
RunOptimizeALS(ifnb_rliger, k = 20, lambda = 5, split.by = "stim")->ifnb_rliger;
RunQuantileNorm(ifnb_rliger, split.by = "stim")->ifnb_rliger;
# step 3：聚类
FindNeighbors(ifnb_rliger, reduction = "iNMF", dims = 1:20)->ifnb_rliger;
FindClusters(ifnb_rliger, resolution = 0.55)->ifnb_rliger;
# 展示
RunUMAP(ifnb_rliger, dims = 1:ncol(ifnb_rliger[["iNMF"]]), reduction = "iNMF")->ifnb_rliger
DimPlot(ifnb_rliger, group.by = c("stim", "ident", "seurat_annotations"), ncol = 3)

#：3、整合不同来源的数据集
#InstallData("panc8")
LoadData("panc8")->panc8_rliger;
# step1：标准分析
NormalizeData(panc8_rliger)->panc8_rliger;
FindVariableFeatures(panc8_rliger)->panc8_rliger;
#：归一化
ScaleData(panc8_rliger, split.by = "replicate", do.center = FALSE)->panc8_rliger;
#：step 2：执行liger，指定消除的变量为 replicate
RunOptimizeALS(panc8_rliger, k = 20, lambda = 5, split.by = "replicate")->panc8_rliger;
RunQuantileNorm(panc8_rliger, split.by = "replicate")->panc8_rliger;
# step 3：聚类
FindNeighbors(panc8_rliger, reduction = "iNMF", dims = 1:20)->panc8_rliger;
FindClusters(panc8_rliger, resolution = 0.4)->panc8_rliger;
# 展示
RunUMAP(panc8_rliger, dims = 1:ncol(panc8_rliger[["iNMF"]]), reduction = "iNMF")->panc8_rliger;
DimPlot(panc8_rliger, group.by = c("replicate", "ident", "celltype"), ncol = 3)
















