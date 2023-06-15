#：使用 seurat 包对单细胞数据进行分析
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
#--

setwd("F:/projects-work/单细胞培训/相关数据/练习题/day2");
#:------------------------------
#:------------------------------
#：1、读入数据
# Load the PBMC dataset
Read10X(data.dir = "pbmc3k_matrices/hg19")->pbmc_data;
# 查看pbmc数据的形式
pbmc_data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
#：矩阵中的.值代表0（没有检测到）。由于scRNA-seq矩阵中的大多数值都是0，Seurat尽可能地使用稀疏矩阵表示。这使得Drop-seq/inDrop/10x数据的内存和速度大大节省。
#：可以用object.size(as.matrix(pbmc_data))来查看对象占用的内存大小

# ：创建 seurat 对象
CreateSeuratObject(counts = pbmc_data, project = "pbmc3k", min.cells = 3, min.features = 200)->pbmc_seurat_obj;
#：seurat对象的结果
str(pbmc_seurat_obj)
#：行对应基因，也称为特征（feature），列对应的是细胞，也就是样本（meta data）
#：查看细胞信息
head(pbmc_seurat_obj@meta.data)
#:------------------------------
#:------------------------------
#：2、基本QC
#：a. 计算线粒体基因比例
PercentageFeatureSet(pbmc_seurat_obj, pattern = "^MT-")->pbmc_seurat_obj[["percent.mt"]];
#：b. 计算核糖体基因比例
PercentageFeatureSet(pbmc_seurat_obj, pattern = "^RP[SL]")->pbmc_seurat_obj[["percent.rb"]]
#：以上得到的两个数值均存放在meta.data的切片中
head(pbmc_seurat_obj@meta.data);
#：如果有doublets注释的结果，可以通过以下命令将此信息添加到seurat对象的meta data中，AddMetaData(pbmc_seurat_obj,doublets)
#：要保证doublets的细胞名称与seurat对象的meta data名称一致
#：在这个例子中，过滤那些UMI count>2500或<200的细胞 线粒体count比例>5%的细胞
# 可视化 QC metrics as a violin plot
VlnPlot(pbmc_seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)
#：FeatureScatter通常用于可视化特征-特征关系，但也可用于对象的任何计算，即对象元数据中的列，PC分数等。
FeatureScatter(pbmc_seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")->f_plot1
FeatureScatter(pbmc_seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")->f_plot2;
f_plot1 + f_plot2
#：练习
FeatureScatter(pbmc_seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(pbmc_seurat_obj, feature1 = "percent.rb", feature2 = "percent.mt")
#：----
#：过滤低质量细胞
"Failed"->pbmc_seurat_obj$QC_type;
which(pbmc_seurat_obj$nFeature_RNA > 200 & pbmc_seurat_obj$nFeature_RNA < 2500 & pbmc_seurat_obj$percent.mt < 5)->qc_pass_index;
"Pass"->pbmc_seurat_obj$QC_type[qc_pass_index];
#：展示
VlnPlot(subset(pbmc_seurat_obj, QC_type == "Pass"), features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1)

#:------------------------------
#:------------------------------
#：3、均一化QC之后的表达值
subset(pbmc_seurat_obj,QC_type=="Pass")->pbmc_seurat_obj;
NormalizeData(pbmc_seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)->pbmc_seurat_obj;
#：默认情况下，许多参数可以省略，因此上述代码通过以下方式实现,
NormalizeData(pbmc_seurat_obj)->pbmc_seurat_obj;
#：4、鉴定细胞间差异最大的基因
FindVariableFeatures(pbmc_seurat_obj, selection.method = "vst", nfeatures = 2000)->pbmc_seurat_obj;
# Identify the 10 most highly variable genes
head(VariableFeatures(pbmc_seurat_obj), 10)->top10_vst_genes;

# plot variable features with and without labels
VariableFeaturePlot(pbmc_seurat_obj)->plot1;
LabelPoints(plot = plot1, points = top10_vst_genes, repel = TRUE)->plot2;
plot1 + plot2

#:------------------------------
#:------------------------------
#：4、对均一化的表达值进行缩放，这是所有降维分析的标准步骤
rownames(pbmc_seurat_obj)->all.genes
ScaleData(pbmc_seurat_obj, features = all.genes)->pbmc_seurat_obj;
#：说明：上述缩放命令也可以不加上features参数，默认使用前2000个差异最大的基因来做
#：如果要去掉某些细胞，也可以指定vars.to.regress参数
ScaleData(pbmc_seurat_obj, vars.to.regress = "percent.mt")->pbmc_seurat_obj;
#：新版seurat建议使用SCTransform()，SCTransform()函数也包括一个vars.to.regress参数。

#:------------------------------
#:------------------------------
#：5、PCA做降维分析，PCA是一种线性降维算法
RunPCA(pbmc_seurat_obj, features = VariableFeatures(object = pbmc_seurat_obj))->pbmc_seurat_obj;
#：可视化：展示每个主成分的重要性组成（这里用loading表示） 
VizDimLoadings(pbmc_seurat_obj, dims = 1:2, reduction = "pca")
#：二维散点图
DimPlot(pbmc_seurat_obj, reduction = "pca")
#：热图展示每个主成分的重要性组成
DimHeatmap(pbmc_seurat_obj, dims = 1:3, cells = 500, balanced = TRUE)

#:------------------------------
#:------------------------------
#：6、如何挑选恰当的主成分进行细胞聚类
#：对于大数据集来说，JackStraw 过程可能需要很长的时间
JackStraw(pbmc_seurat_obj, num.replicate = 100)->pbmc_seurat_obj
ScoreJackStraw(pbmc_seurat_obj, dims = 1:20)->pbmc_seurat_obj
#：展示 JackStraw 结果
JackStrawPlot(pbmc_seurat_obj, dims = 1:15)
#：在这个例子中，我们可以观察到PC9-10周围的 "肘部"，这表明大部分的真实信号是在前10个PC中捕获的。
ElbowPlot(pbmc_seurat_obj)

#:------------------------------
#:------------------------------
#：7、细胞聚类
FindNeighbors(pbmc_seurat_obj, dims = 1:10)->pbmc_seurat_obj;
FindClusters(pbmc_seurat_obj, resolution = 0.5)->pbmc_seurat_obj;
#：查看前10个细胞的聚类结果
head(Idents(pbmc_seurat_obj), 10)
#：使用非线性降维技术来展示聚类的细胞，推荐使用与细胞聚类相同的PCs
RunUMAP(pbmc_seurat_obj, dims = 1:10)->pbmc_seurat_obj;
#：
DimPlot(pbmc_seurat_obj, reduction = "umap",repel=T,label=T)
#：特定基因在细胞亚群上表达水平
FeaturePlot(pbmc_seurat_obj, features = c("CST3", "MS4A1", "CD79A", "GZMB"))
FeaturePlot(pbmc_seurat_obj, features = "percent.mt") & theme(plot.title = element_text(size=10))
#：小提琴图
VlnPlot(pbmc_seurat_obj,features = "percent.mt") & theme(plot.title = element_text(size=10))
VlnPlot(pbmc_seurat_obj,features = c("CST3", "MS4A1")) & theme(plot.title = element_text(size=10))

#:------------------------------
#:------------------------------
#：8、差异表达分析
# cluster 2相对于其他亚类的差异表达基因
FindMarkers(pbmc_seurat_obj, ident.1 = 2, min.pct = 0.25)->cluster2.markers
head(cluster2.markers)
# cluster 5 相对于 clusters 0 和 3的差异表达基因
FindMarkers(pbmc_seurat_obj, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)->cluster5.markers
head(cluster5.markers)
# find markers for every cluster compared to all remaining cells, report only the positive ones
FindAllMarkers(pbmc_seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)->pbmc_all_markers;
pbmc_all_markers %>%group_by(cluster) %>%slice_max(n = 2, order_by = avg_log2FC)
#：Seurat有几种差异表达的test，可以用test.use参数来设置，比如ROC test返回一个标记物的 "分类能力"（范围从0-随机，到1-完美）。
FindMarkers(pbmc_seurat_obj, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)->cluster0.markers
#：展示marker基因在不同亚类细胞上表达情况
VlnPlot(pbmc_seurat_obj, features = c("RPS12", "FCGR3A"))
#：也可以指定使用哪个slot的矩阵来展示
VlnPlot(pbmc_seurat_obj, features = c("RPS12", "FCGR3A"), slot = "counts", log = TRUE)# data,scale.data
#：聚类图上展示
FeaturePlot(pbmc_seurat_obj, features = c("RPS12", "FCGR3A"))
#：也可以通过热图展示
pbmc_all_markers%>%group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC) -> top10_markers;
DoHeatmap(pbmc_seurat_obj, features = top10_markers$gene) + NoLegend()

#:------------------------------
#:------------------------------
#：9、细胞类型注释
#：将 seurat 对象转换为 SingleCellExperiment 对象
celldex::MonacoImmuneData()->celldex_ref;
as.SingleCellExperiment(DietSeurat(pbmc_seurat_obj))->pbmc_sc_obj;
#：
SingleR(test = pbmc_sc_obj,assay.type.test = 1,ref = celldex_ref,labels = celldex_ref$label.main)->pbmc_sc_celldex_main
SingleR(test = pbmc_sc_obj,assay.type.test = 1,ref = celldex_ref,labels = celldex_ref$label.fine)->pbmc_sc_celldex_fine;
# 将注释的细胞类型添加到 seurat 对象中
pbmc_sc_celldex_main$pruned.labels->pruned_labels;
rownames(pbmc_sc_celldex_main)->names(pruned_labels);
AddMetaData(pbmc_seurat_obj,pruned_labels,col.name="Monaco_main_pruned_labels")->pbmc_seurat_obj;
#：另一种添加注释的方法
pbmc_sc_celldex_fine$pruned.labels->pbmc_seurat_obj@meta.data$Monaco_fine_pruned_labels


#：画图展示
#：查看细胞的identify信息，DimPlot画图需要指定
Idents(pbmc_seurat_obj);
SetIdent(pbmc_seurat_obj, value = "Monaco_main_pruned_labels")->pbmc_seurat_obj;
DimPlot(pbmc_seurat_obj, label = T , repel = T, label.size = 3) + NoLegend()
#
SetIdent(pbmc_seurat_obj, value = "Monaco_fine_pruned_labels")->pbmc_seurat_obj;
DimPlot(pbmc_seurat_obj, label = T , repel = T, label.size = 3) + NoLegend()
#：比较注释的结果
table(pbmc_seurat_obj@meta.data$RNA_snn_res.0.5,pbmc_seurat_obj@meta.data$Monaco_main_pruned_labels)
table(pbmc_seurat_obj@meta.data$RNA_snn_res.0.5,pbmc_seurat_obj@meta.data$Monaco_fine_pruned_labels)
#：展示特定的基因
as.data.frame(top10_markers);
FeaturePlot(pbmc_seurat_obj, features = c("GP9", "SPARC"))+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))


























