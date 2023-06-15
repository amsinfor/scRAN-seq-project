#---练习使用bioconductor包来操作单细胞数据：
#--引入相关的包
library(SingleCellExperiment);#BiocManager::install("SingleCellExperiment")






#------------------------------------------
#——-1. 读入数据
setwd("F:/projects-work/单细胞培训/相关数据/练习题/day2");
read.table("sce-data/molecules.txt", sep = "\t",header=T)->sce_counts;
read.table("sce-data/annotation.txt", sep = "\t", header = TRUE)->sce_annotation;
#------------------------------------------
#--创建 SingleCellExperiment 对象
#-- 读入的count变量要求是matrix类型
#-查看count类型：
class(sce_counts);
as.matrix(sce_counts[,-1])->sce_counts_matrix;
sce_counts$gName->rownames(sce_counts_matrix);
#--2. 创建：
SingleCellExperiment(assays = list(counts = sce_counts_matrix),colData = sce_annotation)->sce_obj;
#———删除原始的count数据，释放内存
rm(sce_counts,sce_counts_matrix,sce_annotation)

#———————---如果是标准的10X数据，还可以使用以下方法创建 SingleCellExperiment 对象
library(DropletUtils)
# importing the raw count data
read10xCounts("data/pbmc_1k_raw")->sce_obj;
#————————————
#------------------------------------------
#------2. 查看SCE对象的属性
class(colData(sce_obj))
class(counts(sce_obj))
#------------------------------------------
#-----SCE对象中表达矩阵的操作
#-- new_matrix->assay(sce, "name_of_new_assay")
log2(counts(sce_obj) + 1)->assay(sce_obj, "logcounts");
#-----展示logcounts数据
logcounts(sce_obj)[1:10, 1:4]
#------------------------------------------
#----3. 表达数据的统计
#--1. 每个细胞（所有基因）的平均counts
colMeans(counts(sce_obj))
colMeans(counts(sce_obj))->colData(sce_obj)$mean_count
#--练习：
#1）在colData中添加一个名为“total_counts”的新列，其中每个单元格的计数总和。
colSums(counts(sce_obj))->colData(sce_obj)$total_counts;

#2）创建一个名为“cpm”（Counts-Per-Million）的新测量值，其中包含将counts矩阵除以百万分之一的总counts所得到的结果。
#:：计算方式是针对每个细胞而言，首先得到每个细胞的total count值，该值除以100万进行归一化，然后每个基因的count再除以该值
apply(counts(sce_obj),2,function(cx){cx/sum(cx)*1000000})->assay(sce_obj,"cpm")

#3）如何访问这个新测量值？
assay(sce_obj,"cpm")
colSums(cpm(sce_obj))[1:10]
#------------------------------------------
#----4. 取SCE对象的子集：可以是表达矩阵取子集也可以是细胞注释信息取子集，如果有基因注释信息当然也可以取子集
# 对于表达矩阵而言，通过下标方式：
sce_obj[1:3, ] # the first 3 genes, keep all cells
sce_obj[, 1:3] # the first 3 cells, keep all genes
sce_obj[1:3, 1:2] # the first 3 genes and first 2 cells
# 也可以是名称的方式
sce_obj[c("ENSG00000069712", "ENSG00000237763"), ]
sce_obj[, c("NA19098.r1.A01", "NA19098.r1.A03")]
sce_obj[c("ENSG00000069712", "ENSG00000237763"), c("NA19098.r1.A01", "NA19098.r1.A03")]
#- 条件筛选的方式取子集
# 计算每个基因在所有细胞上平均表达值
rowMeans(counts(sce_obj))->gene_means
sce_obj[gene_means > 0.01, ] #取平均表达值大于0.01以上的基因
# 问题：如果对细胞进行筛选呢？通常我们希望每个细胞检测到的基因数量大于某个阈值，比如至少有2000个基因，
colSums(counts(sce_obj)>0)->total_detected_per_cell;
sce_obj[,total_detected_per_cell>5000]
#--练习：
#创建一个名为sce_filtered的新对象，要求
#1）保留那些至少具有25000个total count的细胞，
(colSums(counts(sce_obj))>25000)->filtered_cells;
#2）保留那些在至少一半的细胞中count>=5的基因
(rowSums(counts(sce_obj)>=5)/ncol(sce_obj)>=0.5)->filtered_genes;
#3）取子集
sce_obj[filtered_genes,filtered_cells]->sce_filtered;
#------------------------------------------
#----5. 数据的展示：有多种展示方法，这里我们借助bioconductor的scater包来实现，该包针对SCE对象而设计的
library(scater)
#--展示每个批次的total counts
ggcells(sce_obj, aes(x = batch, y = total_counts)) + geom_boxplot(fill = 'green') + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#--展示给定基因在每个批次上的count
ggcells(sce_obj,aes(x=batch,y=ENSG00000242485))+geom_boxplot(fill="blue")+ theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#-练习1）：展示mean count与count sd的关系
colSds(counts(sce_obj))->sce_obj$sd_count;
colVars(counts(sce_obj))->sce_obj$var_count;
ggcells(sce_obj,aes(x=mean_count,y=sd_count))+geom_point(aes(color=batch))
ggcells(sce_obj,aes(x=mean_count,y=var_count))+geom_point(aes(color=batch))
#-练习1）：展示两个基因的相关性
ggcells(sce_obj,aes(x=ENSG00000242485,y=ENSG00000175756))+geom_point(aes(color=batch))


##########################################################################
#-
#----5. 数据质量控制
library(scater)
library(SingleCellExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
#————去掉ERCC基因，该基因为加入的标准参考片段
sce_obj[grep("^ERCC-",rownames(sce_obj)), ]->altExp(sce_obj,"ERCC") 
sce_obj[grep("^ERCC-",rownames(sce_obj),invert = T), ]->sce_obj


mapIds(org.Hs.eg.db, keys=rownames(sce_obj), keytype="ENSEMBL", columns="SYMBOL",column="SYMBOL")->gene_names
gene_names->rowData(sce_obj)$SYMBOL;
#--去掉为NA的基因
sce_obj[!is.na(rowData(sce_obj)$SYMBOL),]->sce_obj;


#--借助 EnsDb.Hsapiens.v86 去掉线粒体基因
genes(EnsDb.Hsapiens.v86)->ensdb_genes
#: 获得线粒体上的基因
table(seqnames(ensdb_genes))
ensdb_genes[seqnames(ensdb_genes) == "MT"]$gene_id->MT_names
rownames(sce_obj) %in% MT_names->is_mito;
#————————基本QC：主要信息统计
perCellQCMetrics(sce_obj,subsets=list(Mito=is_mito))->sce_obj_cell_summary
perFeatureQCMetrics(sce_obj)->sce_obj_feature_summary;
#--将基本信息添加到 sce_obj 对象
addPerCellQC(sce_obj, subsets=list(Mito=is_mito))->sce_obj
addPerFeatureQC(sce_obj)->sce_obj
#--1. 对于细胞的过滤
#-查看分布情况，便于手动过滤
hist(sce_obj$total,breaks = 100,xlab="Total reads per cell")
abline(v = 25000, col = "red")
#-
hist(sce_obj_cell_summary$detected, breaks = 100)
abline(v = 7000, col = "red")
#--通过程序自动判断过滤的截断值
#：以下示例表示通过程序判断每个细胞总reads数量的截断值，该值可以用于去除那些测序量不足的单细胞
isOutlier(sce_obj_cell_summary$sum, log=TRUE, type="lower")->qc.lib2;
attr(qc.lib2, "thresholds")
#：当然可以使用以下函数一步操作，实现对细胞、基因过滤
quickPerCellQC(sce_obj_cell_summary, sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))->cell_reasons
colSums(as.matrix(cell_reasons))
cell_reasons$discard->sce_obj$discard;
#--画图展示过滤指标与排除的细胞情况
#：线粒体基因reads数量与总reads数量
plotColData(sce_obj, x="sum", y="subsets_Mito_percent", colour_by="discard")
#：检测到的基因数量与总reads数量
plotColData(sce_obj, x="sum", y="detected", colour_by="discard")
#：ERCC检测的reads比例与线粒体检测的reads比例
plotColData(sce_obj, x="altexps_ERCC_percent", y="subsets_Mito_percent",colour_by="discard")
#：当然也可以分批次（每个个体）来展示某个指标
library(scales)
plotColData(sce_obj, x="sum", y="detected", colour_by="discard", other_fields = "individual") + facet_wrap(~individual) + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))
#：练习：按重复次数来展示检测的基因数量
#plotColData(sce_obj, x="sum", y="detected", colour_by="discard", other_fields = "replicate") + facet_wrap(~replicate)  + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))
#--2. 对于基因的过滤
plotHighestExprs(sce_obj, exprs_values = "counts", feature_names_to_plot = "SYMBOL", colour_cells_by="detected")
#-过滤: 保留在2个或更多细胞中检测到的基因表达值>1
(nexprs(sce_obj,byrow = TRUE,detection_limit = 1) >= 2)->keep_feature
!keep_feature->rowData(sce_obj)$discard;
table(rowData(sce_obj)$discard)

#--3. 排除低质量基因和细胞
sce_obj[!rowData(sce_obj)$discard,! colData(sce_obj)$discard]->sce_obj_qc;
#--: 对基因的read count进行log转换
log2(counts(sce_obj)+1)->assay(sce_obj,"logcounts_raw")
#--4. 展示QC前和QC后的区别
# 未做QC
set.seed(2390213)
runTSNE(sce_obj, exprs_values = "logcounts_raw", perplexity = 130)->sce_obj;
plotTSNE(sce_obj, colour_by = "batch", size_by = "detected", shape_by = "individual")
# QC之后
set.seed(2390213)
sce_obj_qc
runTSNE(sce_obj_qc, exprs_values = "logcounts_raw", perplexity = 130)->sce_obj_qc;
sce_obj_qc
plotTSNE(sce_obj_qc, colour_by = "batch", size_by = "detected", shape_by = "individual")

#---------------------------
#————————：干扰因素的识别和排除
#-1. 检测的基因数与主成分关系
#：首先是使用 logcounts_raw 数据进行PCA分析
runPCA(sce_obj_qc, exprs_values = "logcounts_raw")->sce_obj_qc
dim(reducedDim(sce_obj_qc, "PCA"))
plotPCA(sce_obj_qc, colour_by = "batch", size_by = "sum", shape_by = "individual")
#：获得前10个PC值
getExplanatoryPCs(sce_obj_qc,variables = "sum")
#：展示所有PC的贡献
plotExplanatoryPCs(sce_obj_qc,variables = "sum")+xlab("PC correlation with the number of detected genes") 
#：结果说明：可以看到PC1几乎可以完全（86%）由UMI总计数（测序深度）解释，这是scRNA-seq中一个众所周知的问题。
#assay(sce_obj_qc, "logcounts_raw")->logcounts(sce_obj_qc) 
#NULL->logcounts(sce_obj_qc)
#-2. 多个变量的相关性评估
plotExplanatoryVariables(sce_obj_qc,exprs_values = "logcounts_raw",variables = c("detected","sum","batch","individual","altexps_ERCC_percent","subsets_Mito_percent"))
#；检测到的基因数（再次）和测序深度（计数数）对许多基因有很大的解释力，因此这些变量是在规范化步骤中调理出来的很好的候选变量，
#：或包括在下游统计模型中。ERCCs的表达似乎也是一个重要的解释变量。上述展示的结果显示，一个明显特征是批次解释多于个体。
#：这告诉我们数据的技术和生物变异性是什么
#---------------------------------------
#————————：两个主要干扰因素：文库大小和批次效应
#----1. 文库大小因素的消除：均一化
#：首先用PCA展示未均一化的样本
plotPCA(sce_obj_qc, colour_by = "batch", size_by = "detected", shape_by = "individual")
#：1）CPM方法均一化
log2(calculateCPM(sce_obj_qc) + 1)->logcounts(sce_obj_qc)
runPCA(sce_obj_qc)->sce_obj_qc;
plotPCA(sce_obj_qc, colour_by = "batch", size_by = "detected", shape_by = "individual")
#   说明：SCE对象中名为logcounts的assay是大多数绘图和降维函数的默认值。
#   后续分析中，将使用不同均一化方法来替换这个assay。每次重新进行均一化或运行PCA时，名为logcounts的assay和PCA reducedDim对象都会被替换。
#使用相对对数表达（RLE）画图可以非常有用地评估归一化程序是否成功。https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0191629
plotRLE(sce_obj_qc, exprs_values = "logcounts_raw",colour_by = "batch") + ggtitle("RLE plot for logcounts_raw")
plotRLE(sce_obj_qc, exprs_values = "logcounts",colour_by = "batch") + ggtitle("RLE plot for log2(CPM) counts")
#：2）使用scran包来对做均一化
library(scran) # BiocManager::install("scran")
quickCluster(sce_obj_qc, min.size = 30)->sce_obj_qc.qclust;
#:查看聚类数量
table(sce_obj_qc.qclust);
#：使用上述聚类结果来计算size factor。下述操作为colData增加了一列，名为sizeFactor，这些值被logNormCounts使用。
computeSumFactors(sce_obj_qc, clusters = sce_obj_qc.qclust)->sce_obj_qc;
#:
logNormCounts(sce_obj_qc)->sce_obj_qc;
#:展示均一化后的效果
runPCA(sce_obj_qc)->sce_obj_qc;
plotPCA(sce_obj_qc, colour_by = "batch",size_by = "detected", shape_by = "individual")
#：RLE展示
plotRLE(sce_obj_qc, exprs_values = "logcounts",colour_by = "batch")
#  说明：有时scran会产生负的或零的size factor，这会破坏均一化表达矩阵原来数值结构。可以通过以下方式查看scran计算的size factor：
summary(sizeFactors(sce_obj_qc))
#  对于这个数据集，所有的size factor都是正值；因此这个均一化结果满足进一步分析要求。如果scran均一化出了负的size factor，可以尝试增加聚类数量（min.size参数），直到它们都变成正值。
#——————————————
library(scater)
library(scran)
library(sva)
library(batchelor) #BiocManager::install("batchelor")
library(kBET)# library(devtools);install_github('theislab/kBET')


#----2. 去除批次效应，多个重复之间的干扰因素去除
set.seed(1234567)
#: 1) 使用 combat 方法
ComBat(logcounts(sce_obj_qc),batch = sce_obj_qc$replicate)->assay(sce_obj_qc, "combat") 
#：展示批次校正后的效果
runPCA(sce_obj_qc, exprs_values = "combat", ncomponents = 20)->tmp_
plotPCA(tmp_,colour_by = "batch",size_by = "detected",shape_by = "individual") +ggtitle("ComBat")

#：练习：将检测到的基因数量作为一个变量，对数据进行校正，保存到在SCE对象的assay slot中，combat_tf 中
# ComBat(logcounts(sce_obj_qc),batch = sce_obj_qc$detected)->assay(sce_obj_qc, "combat_tf")
#：2）使用 MNN 方法；batchelor包提供
fastMNN(sce_obj_qc,batch = sce_obj_qc$replicate)->mnn_out
assay(mnn_out,'reconstructed')->assay(sce_obj_qc, "mnn")
#：展示
runPCA(sce_obj_qc, exprs_values = "mnn", ncomponents = 20)->tmp_
plotPCA(tmp_,colour_by = "batch",size_by = "detected",shape_by = "individual") +ggtitle("MNN")
#---------------------------------------
#————————：评估干扰因素的排除效果
#：1）PCA画图展示，可以用不同颜色对应技术重复，形状对应不同的生物样本（个体）。生物样本和交错批次的分离表明已经消除了技术变异
#：注意：这里用来做PCA的数据，都是log2-cpm均一化以后的。
runPCA(sce_obj_qc, exprs_values = "combat_tf", ncomponents = 20)->tmp_
plotPCA(tmp_,colour_by = "batch",size_by = "detected",shape_by = "individual") +ggtitle("Normalization by combat_tf")
#: 2）RLE方法
#：计算每个细胞中每个基因的RLE值，根据RLE定义，可以写出相应的计算代码如下：
get_cell_RLE_values<-function(expd_matrix){
	apply(expd_matrix,1,function(ex){
		if (median(unlist(ex)) > 0) {
            log( (ex + 1) / ( median(unlist(ex)) + 1 ) ) / log(2)
        } else {
            rep(NA, times = length(ex))
        }
	})->RLE_matrix;
	t(RLE_matrix)->RLE_matrix;
	#--
	apply(RLE_matrix, 2, median, na.rm = T)->cell_RLE;
    return(cell_RLE)
}
#--
lapply(assayNames(sce_obj_qc),function(nx){
	suppressWarnings(get_cell_RLE_values(assay(sce_obj_qc, nx)))
})->res;
#-画图
par(mar=c(6,4,1,1))
boxplot(res, las=2)
#：3）kBET方法
calculate_kBET_results <- function(SCE_obj){
    unique(as.character(SCE_obj$individual))->indiv;
    assayNames(SCE_obj)->norm_names # Get all normalization names
    list()->results;
    for (i in indiv){ 
        for (j in norm_names){
            tmp <- kBET(
                df = t(assay(SCE_obj[,SCE_obj$individual== i], j)), #第i个个体的第j个均一化的矩阵
                batch = SCE_obj$batch[SCE_obj$individual==i], #批次，这里对应的是三个重复
                heuristic = TRUE, 
                verbose = FALSE, 
                addTest = FALSE, 
                plot = FALSE)
            tmp$summary$kBET.observed[1]->results[[i]][[j]]
        }
    }
    return(do.call(rbind.data.frame, results))
}
#-
calculate_kBET_results(sce_obj_qc)->eff_debatching
eff_debatching
#-展示：高拒绝率意味着更强的批次效应，值越大批次效应越强
library(pheatmap)
pheatmap(eff_debatching);



















