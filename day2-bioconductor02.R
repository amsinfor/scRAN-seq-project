#---单细胞数据的标准分析：聚类和差异表达基因分析
library(pcaMethods)
library(SC3) # BiocManager::install("SC3")
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
set.seed(642237);
#--使用小鼠胚胎细胞数据集，因为它提供了每个细胞的类型信息，方便我们比较不同聚类方法的性能
setwd("F:/projects-work/单细胞培训/相关数据/练习题/day2")
#：载入：
readRDS("F:/projects-work/单细胞培训/相关数据/data/deng/deng-reads.rds")->SCE_mouse_embryo
#:查看该SCE对象的结构
SCE_mouse_embryo
#：查看细胞类型
head(colData(SCE_mouse_embryo))
table(colData(SCE_mouse_embryo)$cell_type2)
#：PCA展示结果，默认使用 logcounts assay数据做降维
runPCA(SCE_mouse_embryo)->SCE_mouse_embryo;
plotPCA(SCE_mouse_embryo, colour_by = "cell_type2",size_by="log10_total_counts")
#: plotReducedDim(object,dimred,ncomponents = 2,percentVar = NULL,colour_by = NULL,shape_by = NULL,size_by = NULL,by_exprs_values = "logcounts",text_by = NULL,text_size = 5,text_colour = "black",   label_format = c("%s %i", " (%i%%)"),other_fields = list(),swap_rownames = NULL,   ... )
#-------------------------------
#--- 1. 聚类分析
#-- 1）SC3方法
#： SC3的优势在于它可以直接摄取SingleCellExperiment对象，首先可以用SC3估计一个聚类的数量：
sc3_estimate_k(SCE_mouse_embryo)->SCE_mouse_embryo;
#： 查看估计的cluster数量
metadata(SCE_mouse_embryo)$sc3$k_estimation
#：SC3预测的细胞类型的数量比原始数据注释中的要少。如果将不同细胞类型的早期、中期和晚期合并在一起，这个结果正好对应6种细胞类型
#：测试 plotPCA(SCE_mouse_embryo, colour_by = "cell_type1",size_by="log10_total_counts")
#：运行
sc3(SCE_mouse_embryo, ks = 10, biology = TRUE, n_cores = 4)->SCE_mouse_embryo;
#：查看SC3的结果，注意：参考k要与ks的值对应
#：??SC3 查看帮助
#：a. 基于shiny的交互式结果展示
sc3_interactive(SCE_mouse_embryo) 
#：b：单独查看各项结果，一致性矩阵
sc3_plot_consensus(SCE_mouse_embryo, k = 10,show_pdata = c("cell_type1","log10_total_features","sc3_10_clusters","sc3_10_log2_outlier_score")) 
#：轮廓系数是对共识矩阵的对角线性量化测量。平均轮廓宽度（显示在轮廓图的左下方）从0到1变化，其中1表示完全分块对角共识矩阵，而0表示没有分块对角结构的情况。当平均轮廓宽度接近于1时，最佳聚类效果得以实现。
sc3_plot_silhouette(SCE_mouse_embryo, k = 10)
#：稳定性指数显示了每个聚类在选定的ks范围内的稳定程度。稳定性指数在0和1之间变化，其中1意味着在不同的k中，每个解决方案都出现相同的集群。
sc3_plot_cluster_stability(SCE_mouse_embryo, k = 10)
#：差异表达是用非参数Kruskal-Wallis检验计算的。显著的P值表明，至少有一个簇的基因表达随机地支配着另一个簇。SC3提供了一个调整后P值<0.01的所有差异表达基因的列表，并绘制了P值最低的50个基因的基因表达谱。请注意，聚类后差异表达的计算会在p值的分布中引入偏差，因此我们建议只使用p值对基因进行排序。
sc3_plot_de_genes(SCE_mouse_embryo, k = 10,show_pdata = c("cell_type1","log10_total_features","sc3_10_clusters","sc3_10_log2_outlier_score")) 
#：为了找到标记基因，对于每个基因，根据平均群组表达值构建一个二元分类器。然后用基因表达量的等级计算分类器的预测值。接受者操作特征（ROC）曲线下的面积被用来量化预测的准确性。通过使用Wilcoxon签名等级测试为每个基因分配一个P值。默认情况下，ROC曲线下面积（AUROC）>0.85且P值<0.01的基因被选中，每个集群的前10个标记基因在此热图中被可视化。
sc3_plot_markers(SCE_mouse_embryo, k = 10,p.val=1e-15)
#：PCA展示 SC3 聚类的结果
plotPCA(SCE_mouse_embryo, colour_by = "sc3_10_clusters",size_by="log10_total_counts")
#--- 2）K-means和TSNE方法
library(factoextra)
library(cluster)
#-确定最优K值：
logcounts(SCE_mouse_embryo)->data_f;
#：注意：k-means聚类输入数据要求行为样本，列为基因
fviz_nbclust(t(data_f), kmeans, method = "wss")
#-另一种方式：非常耗时，约30min
clusGap(t(data_f),FUN = kmeans,nstart = 25,K.max = 10,B = 50)->df_gap_stat
fviz_gap_stat(df_gap_stat)
#- 运行
kmeans(t(data_f), centers = 10, nstart = 25)->data_df_km;
data_df_km
#：展示
fviz_cluster(data_df_km, data = t(data_f))
#：tSNE展示
paste("KC",data_df_km$cluster,sep="")->colData(SCE_mouse_embryo)$KCluster;
runTSNE(SCE_mouse_embryo)->SCE_mouse_embryo;
plotPCA(SCE_mouse_embryo, colour_by = "KCluster",size_by="log10_total_counts")
#：对比
table(SCE_mouse_embryo$cell_type1,SCE_mouse_embryo$KCluster)

#--- 3）层次聚类方法
#：首先进行Z-score转换
apply(data_f, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))->data_f_zscore;
# hierarchical clustering
as.dist((1 - cor(t(data_f_zscore), method = "pearson"))/2)->data_f_zscore_dd;
hclust(data_f_zscore_dd, method = "average")->data_f_zscore_dd_hc;
#：确定聚类数量
0->num_singleton
1->kk;
for (i in 2:dim(data_f_zscore)[2]) {
    clusters <- cutree(data_f_zscore_dd_hc, k = i)
    clustersizes <- as.data.frame(table(clusters))
    singleton_clusters <- which(clustersizes$Freq < 2)
    if (length(singleton_clusters) <= num_singleton) {
        i->kk;
    } else {
        break;
    }
}
#：PCA展示
paste("HC",clusters,sep="")->colData(SCE_mouse_embryo)$HCluster;
plotPCA(SCE_mouse_embryo, colour_by = "HCluster",size_by="log10_total_counts")
#-------------------------------
#--- 2. 差异表达分析
#：说明：为了测试不同的单细胞差异表达方法，这里使用Blischak数据集。在这个实验中，除了单细胞数据外，还有每个细胞系的bulk RNA-seq数据。
#：通过比较bulk RNA-seq的差异表达基因与单细胞方法的差异表达基因来评估几种scRNA算法的准确性
library(scRNA.seq.funcs)
library(edgeR)
#library(monocle)
library(MAST) # BiocManager::install("MAST")
library(ROCR)
set.seed(2232)
#--
#----------
#：1）读入表达谱数据
read.table("F:/projects-work/单细胞培训/相关数据/data/tung/molecules.txt", sep = "\t")->expr_matrix;
read.table("F:/projects-work/单细胞培训/相关数据/data/tung/annotation.txt", sep = "\t", header = TRUE)->cell_pheno;
#：只保留两个样本的单细数据：NA19101 和 NA19239
c(grep("NA19101",colnames(expr_matrix)),grep("NA19239",colnames(expr_matrix)))->keep_samples;
expr_matrix[,keep_samples]->expr_matrix_testd;
cell_pheno[which(cell_pheno$individual%in%c("NA19239","NA19101")),]->cell_pheno_testd;
#：分组和重复
cell_pheno_testd$individual->test_group;
cell_pheno_testd$batch->test_batch;

# remove genes that aren't expressed in at least 6 cells
(rowSums(expr_matrix_testd > 0) > 5)->keep_genes;
expr_matrix_testd[keep_genes,]->expr_matrix_testd;
# Library size normalization
apply(expr_matrix_testd,2,function(cx){cx/sum(cx)*1000000})->expr_matrix_testd_norm;

# Variant of CPM for datasets with library sizes of fewer than 1 mil molecules
#：2）差异表达分析
#：a) 非参数显著性检验的Kolmogorov-Smirnov Test，也称为KS-test
apply(expr_matrix_testd_norm, 1, function(gx) {
		gx[test_group == "NA19101"]->gx_group1_values;
		gx[test_group == "NA19239"]->gx_group2_values;
        ks.test(gx_group1_values, gx_group2_values)->gx_test;
		c(mean(gx_group1_values,na.rm=T),mean(gx_group2_values,na.rm=T),gx_test$p.value)
})->expr_matrix_testd_norm_ks_test;
t(expr_matrix_testd_norm_ks_test)->expr_matrix_testd_norm_ks_test;
c("group1_mean","group2_mean","P.value")->colnames(expr_matrix_testd_norm_ks_test)
data.frame("gName"=rownames(expr_matrix_testd_norm),expr_matrix_testd_norm_ks_test)->expr_matrix_testd_norm_ks_test;
# multiple testing correction
p.adjust(expr_matrix_testd_norm_ks_test$P.value, method = "fdr")->expr_matrix_testd_norm_ks_test$Padj;
#：b) 非参数显著性检验的Wilcox/Mann-Whitney-U Test
apply(expr_matrix_testd_norm, 1, function(gx) {
		gx[test_group == "NA19101"]->gx_group1_values;
		gx[test_group == "NA19239"]->gx_group2_values;
        wilcox.test(gx_group1_values, gx_group2_values)->gx_test;
		c(mean(gx_group1_values,na.rm=T),mean(gx_group2_values,na.rm=T),gx_test$p.value)
})->expr_matrix_testd_norm_ms_test;
t(expr_matrix_testd_norm_ms_test)->expr_matrix_testd_norm_ms_test;
c("group1_mean","group2_mean","P.value")->colnames(expr_matrix_testd_norm_ms_test)
data.frame("gName"=rownames(expr_matrix_testd_norm),expr_matrix_testd_norm_ms_test)->expr_matrix_testd_norm_ms_test;
# multiple testing correction
p.adjust(expr_matrix_testd_norm_ms_test$P.value, method = "fdr")->expr_matrix_testd_norm_ms_test$Padj;
#---
#:3）edgeR：采用负二项分布模型
#：要求输入的是原始的read count类型数据
DGEList(counts = expr_matrix_testd, norm.factors = rep(1, length(expr_matrix_testd[1,])), group = test_group)->edgeR_dge;
factor(test_group)->edgeR_group;
model.matrix(~ edgeR_group)->edgeR_design;
estimateDisp(edgeR_dge, design = edgeR_design, trend.method = "none")->edgeR_dge;
glmFit(edgeR_dge, edgeR_design)->edgeR_fit;
glmLRT(edgeR_fit)->edgeR_fit_res;
#-
data.frame("gName"=rownames(edgeR_fit_res$table),edgeR_fit_res$table)->expr_matrix_testd_edgeR;
p.adjust(expr_matrix_testd_edgeR$PValue, method = "fdr")->expr_matrix_testd_edgeR$Padj;
#---
#:4）MAST：采用0膨胀模型；需要指出的是，该算法耗时比较长，如果结合R语言的并行计算，时间会大大缩短
#：log转换
log2(expr_matrix_testd + 1)->expr_matrix_testd_log;
#：构建特征信息
data.frame(names = rownames(expr_matrix_testd_log))->MAST_fData;
rownames(expr_matrix_testd_log)->rownames(MAST_fData);
#：构建细胞信息
data.frame(cond = test_group)->MAST_cData;
colnames(expr_matrix_testd_log)->rownames(MAST_cData) 
#：构建MAST对象
FromMatrix(as.matrix(expr_matrix_testd_log), MAST_cData, MAST_fData)->MAST_obj;
scale(colSums(assay(MAST_obj) > 0))->colData(MAST_obj)$cngeneson
factor(colData(MAST_obj)$cond)->MAST_cond;
# 构建模型
zlm(~ MAST_cond + cngeneson, MAST_obj)->MAST_zlmCond;
summary(MAST_zlmCond, doLRT = "MAST_condNA19239")->MAST_summaryCond
as.data.frame(MAST_summaryCond$datatable)->MAST_summaryDt
#--
MAST_summaryDt[!is.na(MAST_summaryDt[,4]),]->MAST_summaryDt;
p.adjust(MAST_summaryDt[,4], method = "fdr")->MAST_summaryDt$Padj;
#------------
#：其他方法：BPSC包：https://github.com/nghiavtr/BPSC
#install.packages("devtools")
library("devtools")
install_github("nghiavtr/BPSC")
library("BPSC")
#：查看帮助
vignette("BPSC")
#--准备输入数据：这里只取1个replicate的单细胞
c(grep("NA19101.r1",colnames(expr_matrix_testd_norm)),grep("NA19239.r1",colnames(expr_matrix_testd_norm)))->keep_samples;
expr_matrix_testd_norm[,keep_samples]->bpsc_data;
cell_pheno_testd$batch[which(cell_pheno_testd$batch%in%c("NA19101.r1","NA19239.r1"))]->bpsc_group;
cell_pheno_testd$sample_id[which(cell_pheno_testd$batch%in%c("NA19101.r1","NA19239.r1"))]->names(bpsc_group)
#：构建design模型
which(bpsc_group == "NA19101.r1")->bpsc_control_cells
model.matrix(~bpsc_group)->bpsc_design;
#：执行计算：这里 control ID是 bpsc_data 中的样本列下标，需要保证顺序一致
BPglm(data=bpsc_data, controlIds=bpsc_control_cells, design=bpsc_design, coef=2, estIntPar=FALSE, useParallel = FALSE)->bpsc_res;
#p.adjust(bpsc_res$PVAL, method = "fdr")->bpsc_res$Padj;

#：SCDE包：http://hms-dbmi.github.io/scde/
#：scde软件包实现了一套用于分析单细胞RNA-seq数据的统计方法。scde为单细胞RNA-seq测量值拟合了单个误差模型。
#：这些模型随后可用于评估细胞组之间的差异性表达，以及其他类型的分析。
#：scde软件包还包含pagoda框架，它应用路径和基因组过度分散分析来确定单细胞之间转录异质性的各个方面。
library(devtools)
devtools::install_version('flexmix', '2.3-13')
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)
#：供练习





















