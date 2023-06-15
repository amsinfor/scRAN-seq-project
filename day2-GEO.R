#-R语言分析GEO表达数据
#-----------------------------------------
#--设置工作路径：这一步可不操作
setwd("F:/projects-work/单细胞培训/相关数据/练习题/day2")

#-----------------------------------------
#--导入相关的package
library(affy)# BiocManager::install("affy")
library(hgu133a.db);# BiocManager::install("hgu133a.db")


#-----------------------------------------
#-----读取CEL文件
ReadAffy(celfile.path="GSE46517/GSE46517_RAW/")->celdat;
#--均一化：mas5和rma两种方法
#mas5(celdat)->celdat_normalized;
rma(celdat)->celdat_normalized;
exprs(celdat_normalized)->celdat_normalized_exp
colnames(celdat_normalized_exp)->expd_names;
unlist(lapply(expd_names,function(x){unlist(strsplit(x,split="\\."))->res;res[1]}))->expd_names;
expd_names->colnames(celdat_normalized_exp);
data.frame("DATA"=rownames(celdat_normalized_exp),celdat_normalized_exp)->celdat_normalized_exp;


#-----------------------------------------
#----将探针注释到基因上：
as.character(celdat_normalized_exp$DATA)->expd_probes;	
annotate::getSYMBOL(expd_probes,"hgu133a.db")->expd_probes_substr_SYMBOL;
which(is.na(expd_probes_substr_SYMBOL))->na_index;
names(expd_probes_substr_SYMBOL)[na_index]->expd_probes_substr_SYMBOL[na_index]
data.frame("Probe"=celdat_normalized_exp[,1],"gName"=expd_probes_substr_SYMBOL,celdat_normalized_exp[,-1])->GSE46517_normalized_exp;

#-----------------------------------------
#----样本临床特征信息处理：
c("Sample_title","Sample_geo_accession","Sample_last_update_date","Sample_source_name_ch1","Sample_characteristics_ch1","Sample_supplementary_file")->infors;
readLines(gzfile("GSE46517/GSE46517_series_matrix.txt.gz",'r'))->myd_info;
c()->res;
c()->res_colnames;
for(s in infors){
	myd_info[grep(s,myd_info)]->s_lines;
	for(sl in s_lines){
		unlist(strsplit(sl,split="\t"))->s_lines.split;
		gsub("\"","",s_lines.split)->s_lines.split;
		if(grepl(": ",sl)){	
			unlist(strsplit(s_lines.split[2],split=": "))[1]->s;
			unlist(lapply(s_lines.split,function(x){
				unlist(strsplit(x,split=": "))[2]
			}))->x_res;
			c(res,x_res[-1])->res;
		}else if(grepl("ftp",sl)){
			unlist(lapply(s_lines.split,function(x){
				gsub("\\.gz","",basename(x));
			}))->x_basenames;
			c(res,x_basenames[-1])->res;
		}else{
			c(res,s_lines.split[-1])->res;
		}
		c(res_colnames,s)->res_colnames;
	}	
}
matrix(res,ncol=length(res_colnames),byrow=F)->res.matrix;
res_colnames->colnames(res.matrix);
as.data.frame(res.matrix)->GSE46517_infor;
#---样本信息处理
GSE46517_infor[,c(1,2,5)]->GSE46517_infor;
c("SampleID","GSM_ID","SampleType")->colnames(GSE46517_infor);
#--修改样本分组的信息
barplot(table(GSE46517_infor$SampleType))
unlist(lapply(GSE46517_infor$SampleType,function(sx){unlist(strsplit(sx,split=" "))[1]}))->GSE46517_infor$SampleType;
#--将临床信息与表达谱整合到一起，这一步稍微复杂
#-首先：对表达矩阵的列名称进行处理，以保证与样本信息数据框中样本名称一致
unlist(lapply(colnames(GSE46517_normalized_exp),function(cx){unlist(strsplit(cx,split="_"))[1]}))->colnames(GSE46517_normalized_exp)
#-表达谱矩阵进行转换
GSE46517_normalized_exp[,-c(1,2)]->GSE46517_normalized_exp_matrix;
paste(GSE46517_normalized_exp$gName,GSE46517_normalized_exp$Probe,sep="|")->rownames(GSE46517_normalized_exp_matrix);
data.frame("GSM_ID"=colnames(GSE46517_normalized_exp_matrix),t(GSE46517_normalized_exp_matrix))->tmp_;
merge(GSE46517_infor,tmp_,by.x="GSM_ID",by.y="GSM_ID")->GSE46517_exp_factor;

#-----------------------------------------
#----样本结构展示：用于判断有无异常样本，若有最好是去掉异常值后再进行后续分析
#---方法1：PCA
library(ggbiplot)
#-首先调整样本顺序：
GSE46517_normalized_exp_matrix[,GSE46517_exp_factor$GSM_ID]->GSE46517_normalized_exp_matrix;
prcomp(t(GSE46517_normalized_exp_matrix),center = TRUE,scale = TRUE)->GSE46517_normalized_exp_matrix.pca#row->samples, column->features;
ggbiplot(GSE46517_normalized_exp_matrix.pca,ellipse=TRUE,circle=TRUE,groups=GSE46517_exp_factor$SampleType,varname.size=0,var.axes=FALSE)->pca_p;
#-设置主题
pca_p+theme_bw()
#-设置颜色
pca_p+theme_bw()+scale_color_manual(values=c("red","green","blue","black"));
#-展示样本名称：设置labels参数
ggbiplot(GSE46517_normalized_exp_matrix.pca,ellipse=TRUE,circle=TRUE,labels=GSE46517_exp_factor$GSM_ID, groups=GSE46517_exp_factor$SampleType,varname.size=0,var.axes=FALSE)+theme_bw()+scale_color_manual(values=c("red","green","blue","black"))

#---方法2：tSNE方法
library(Rtsne)
t(GSE46517_normalized_exp_matrix[,GSE46517_exp_factor$GSM_ID])->Rtsne_data;
## do tsne: Executing the algorithm on curated data	
Rtsne(Rtsne_data, dims = 2, perplexity=round((nrow(Rtsne_data)-1)/3.1), verbose=TRUE, max_iter = 500)->tsne
## Plotting 
c("tSNE1","tSNE2")->colnames(tsne$Y);
cbind(GSE46517_exp_factor,tsne$Y)->GSE46517_exp_factor;
#-
library(ggpubr)
ggscatter(GSE46517_exp_factor,x="tSNE1",y="tSNE2",color="SampleType");
#-设置颜色和形状
ggscatter(GSE46517_exp_factor,x="tSNE1",y="tSNE2",color="SampleType",palette="npg",shape="SampleType")
#-显示样本名称：
ggscatter(GSE46517_exp_factor,x="tSNE1",y="tSNE2",color="SampleType",palette="npg",shape="SampleType",label=GSE46517_exp_factor$GSM_ID,font.label = c(8, "plain"))

#-----------------------------------------
#----差异表达分析：limma包
library(limma);
#1）不去掉异常值：
#-只对转移和非转移两组样本进行比较：
GSE46517_exp_factor$GSM_ID[which(GSE46517_exp_factor$SampleType%in%c("Primary","Metastatic"))]->keep_samples;
GSE46517_normalized_exp_matrix[,keep_samples]->expd;
#-构建表型对象：
data.frame("Condition"=GSE46517_exp_factor$SampleType[which(GSE46517_exp_factor$SampleType%in%c("Primary","Metastatic"))])->phenoData;
keep_samples->rownames(phenoData);
new("AnnotatedDataFrame",data = phenoData)->phenoData;
#- 构建 ExpressionSet 对象必须是matrix
ExpressionSet(assayData=as.matrix(expd),phenoData = phenoData)->experimentData;
pData(experimentData)$Condition->condition;
factor(condition,levels=c("Primary","Metastatic"))->condition;
model.matrix(~condition)->design;
lmFit(experimentData,design = design)->fit;
eBayes(fit)->fit_eBayes;
topTable(fit_eBayes,number=nrow(expd))->GSE46517_exp_limma_res_no;
#-----------------------
#2）去掉异常值：
#c("GSM1131667")->removed_samples;
#setdiff(keep_samples,removed_samples)->keep_samples_filter;
#-
#GSE46517_normalized_exp_matrix[,keep_samples_filter]->expd;
#-构建表型对象：
#data.frame("Condition"=GSE46517_exp_factor$SampleType[which(GSE46517_exp_factor$GSM_ID%in%keep_samples_filter)])->phenoData;
#GSE46517_exp_factor$GSM_ID[which(GSE46517_exp_factor$GSM_ID%in%keep_samples_filter)]->rownames(phenoData);
#new("AnnotatedDataFrame",data = phenoData)->phenoData;
#- 构建 ExpressionSet 对象必须是matrix
#ExpressionSet(assayData=as.matrix(expd),phenoData = phenoData)->experimentData;
#pData(experimentData)$Condition->condition;
#model.matrix(~condition)->design;
#c("Intercept","Condition_Compare")->colnames(design)
#lmFit(experimentData,design = design)->fit;
#eBayes(fit)->fit_eBayes;
#topTable(fit_eBayes,number=nrow(expd))->GSE46517_exp_limma_res_yes;

#-----------------------------------------
#----差异表达基因的功能富集分析
#：对limma的结果进行进一步处理，添加基因和探针
rownames(GSE46517_exp_limma_res_no)->tmp_;
tmp_->GSE46517_exp_limma_res_no$nameProbe;
unlist(lapply(tmp_,function(tx){unlist(strsplit(tx,split="\\|"))[1]}))->GSE46517_exp_limma_res_no$gName;
unlist(lapply(tmp_,function(tx){unlist(strsplit(tx,split="\\|"))[2]}))->GSE46517_exp_limma_res_no$Probe;
#--画火山图：
#-标记基因的表达上调、下调
"Not"->GSE46517_exp_limma_res_no$DEG_type;
"Up"->GSE46517_exp_limma_res_no$DEG_type[GSE46517_exp_limma_res_no$adj.P.Val<0.05&GSE46517_exp_limma_res_no$logFC>0]
"Down"->GSE46517_exp_limma_res_no$DEG_type[GSE46517_exp_limma_res_no$adj.P.Val<0.05&GSE46517_exp_limma_res_no$logFC<0]
table(GSE46517_exp_limma_res_no$DEG_type)
-log10(GSE46517_exp_limma_res_no$adj.P.Val)->GSE46517_exp_limma_res_no$nlog_adj_P;
#-
ggscatter(GSE46517_exp_limma_res_no,x="logFC",y="adj.P.Val",yscale="-log10");
ggscatter(GSE46517_exp_limma_res_no,x="logFC",y="nlog_adj_P",color="DEG_type")+geom_vline(xintercept=c(0),lty=2);
#-----------
#--取padj<0.05的基因，或者是|logFC|>2的
GSE46517_exp_limma_res_no[GSE46517_exp_limma_res_no$adj.P.Val<0.05,]->GSE46517_exp_limma_res_no_filter;
#- up-DEGs 和 download-DEGs
GSE46517_exp_limma_res_no_filter$gName[GSE46517_exp_limma_res_no_filter$logFC>0]->up_DEGs;
GSE46517_exp_limma_res_no_filter$gName[GSE46517_exp_limma_res_no_filter$logFC<0]->down_DEGs;
#- 去掉重复的基因：因为每个基因有多个探针对应，因此会出现重复基因。
unique(up_DEGs)->up_DEGs;
unique(down_DEGs)->down_DEGs;
length(up_DEGs);
length(down_DEGs);
#- 去掉那些没有注释到基因的探针
up_DEGs[-grep("_",up_DEGs)]->up_DEGs;
down_DEGs[-grep("_",down_DEGs)]->down_DEGs;
length(up_DEGs);
length(down_DEGs);
#----------------------
#--- GO富集
library(R.utils);#解决 fail to download KEGG data 问题
#R.utils::setOption("clusterProfiler.download.method",'auto')#或者
options(clusterProfiler.download.method = "wininet")
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)#combine multiple plots
#-- up_DEGs
enrichGO(gene=up_DEGs,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05)->GO_KEGG.selected_genes.GO;
#----------------------
#--- KEGG富集
bitr(up_DEGs,fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb=org.Hs.eg.db)->myd.gName.df;
enrichKEGG(gene=myd.gName.df[,3],organism="hsa",keyType="ncbi-geneid",pvalueCutoff=0.05)->myd.gName.df.KEGG
#--
dotplot(GO_KEGG.selected_genes.GO, showCategory=20,label_format = 100) + ggtitle("BP for up-DEGs")->p1_GO
dotplot(myd.gName.df.KEGG, showCategory=20,label_format = 100) + ggtitle("KEGG for up-DEGs")->p1_KEGG;
#-GO+KEGG: merged plot 
plot_grid(p1_GO,p1_KEGG,ncol=2,rel_widths = c(4,4.5);
#-- gsea富集分析
bitr(up_DEGs,fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb=org.Hs.eg.db)->myd.gName.df;
merge(GSE46517_exp_limma_res_no_filter,myd.gName.df,by.y="SYMBOL",by.x="gName",all=F)->easy_df
# 只需要基因名称、fold change和 ENTREZID 三列信息
easy_df[,c(1,2,11)]->easy_df;
c("SYMBOL","foldChange","ENTREZID")->colnames(easy_df);
easy_df[order(easy_df$foldChange, decreasing = T),]->easy_df;#降序排序
easy_df$foldChange->gene_expr#把foldchange按照从大到小提取出来
easy_df$ENTREZID->names(gene_expr)#给上面提取的foldchange对应上ENTREZID
gseKEGG(gene_expr,organism = "hsa",pvalueCutoff=0.05) ->gene_expr_kk;
#--按照enrichment score从高到低排序，便于查看富集通路
library(enrichplot)
#gene_expr_kk[order(gene_expr_kk$enrichmentScore, decreasing = T),]->gene_expr_kk;
gseaplot2(gene_expr_kk, c("hsa03010","hsa05202","hsa04970","hsa00190"), color = colorspace::rainbow_hcl(2), pvalue_table = TRUE)




##########################################################
#---------------------------------------------------------
#----无CEL格式的表达数据处理：以GSE79691为例
library(illuminaHumanv4.db)#for Illumina HumanHT-12 V4.0 expression beadchip
#： 将上述处理series matrix的代码封装，用函数形式呈现，便于其他数据集的使用
prepare_series_matrix_data<-function(infile){
	readLines(gzfile(infile,'r'))->infile_lines;
	grep("series_matrix_table_begin",infile_lines)->skip_lines;
	read.table(infile,skip=skip_lines,sep="\t",header=T,stringsAsFactors=F,nrows=length(infile_lines)-skip_lines-2)->myd;
	"DATA"->names(myd)[1];
	return(myd);
}
process_series_matrix<-function(series_file,infors){
	readLines(gzfile(series_file,'r'))->myd_info;
	c()->res;
	c()->res_colnames;
	for(s in infors){
		myd_info[grep(s,myd_info)]->s_lines;
		for(sl in s_lines){
			unlist(strsplit(sl,split="\t"))->s_lines.split;
			gsub("\"","",s_lines.split)->s_lines.split;
			if(grepl(": ",sl)){	
				unlist(strsplit(s_lines.split[2],split=": "))[1]->s;
				unlist(lapply(s_lines.split,function(x){
					unlist(strsplit(x,split=": "))[2]
				}))->x_res;
				c(res,x_res[-1])->res;
			}else if(grepl("ftp",sl)){
				unlist(lapply(s_lines.split,function(x){
					gsub("\\.gz","",basename(x));
				}))->x_basenames;
				c(res,x_basenames[-1])->res;
			}else{
				c(res,s_lines.split[-1])->res;
			}
			c(res_colnames,s)->res_colnames;
		}
	}
	matrix(res,ncol=length(res_colnames),byrow=F)->res.matrix;
	res_colnames->colnames(res.matrix);
	as.data.frame(res.matrix)->res.matrix
	return(res.matrix);
}
#--
map_illuminaHuman2SYMBOL<-function(expd,db){
	as.character(expd[,1])->expd_probes;
	if(db=="illuminaHumanv3.db"){
		unlist(as.list(illuminaHumanv3SYMBOL))->illuminaHumanv4SYMBOL_list;
	}else if(db=="illuminaHumanv4.db"){
		unlist(as.list(illuminaHumanv4SYMBOL))->illuminaHumanv4SYMBOL_list;
	}else if(db=="illuminaHumanv2.db"){
		unlist(as.list(illuminaHumanv2SYMBOL))->illuminaHumanv4SYMBOL_list;
	}else if(db=="IlluminaHumanMethylation27k"){
		unlist(as.list(IlluminaHumanMethylation27kSYMBOL))->illuminaHumanv4SYMBOL_list;
	}else if(db=="IlluminaHumanMethylation450k"){
		#unlist(as.list(IlluminaHumanMethylation450kSYMBOL))->illuminaHumanv4SYMBOL_list;
		print("not work for 450k!");flush.console();return(NULL);
	}
	illuminaHumanv4SYMBOL_list[as.character(expd_probes)]->expd_probes_symbol;
	#select(illuminaHumanv4.db,keys = expd_probes,columns=c("SYMBOL"),keytype="PROBEID")->expd_probes_symbol;
	which(is.na(expd_probes_symbol))->na_index;
	expd_probes[na_index]->expd_probes_symbol[na_index]
	data.frame("Probe"=expd_probes,"gName"=expd_probes_symbol,expd[,-1])->expd_t;
	return(expd_t);
}
#---
process_series_matrix("GSE79691/GSE79691_series_matrix.txt.gz",c("Sample_title","Sample_geo_accession","Sample_last_update_date","Sample_source_name_ch1","Sample_characteristics_ch1","Sample_supplementary_file"))->GSE79691_targets;
c("Sample_title","A0_Samples","Sample_last_update_date","Sample_source_name_ch1","Response","Sample_supplementary_file")->colnames(GSE79691_targets);
#--
prepare_series_matrix_data("GSE79691/GSE79691_series_matrix.txt.gz")->GSE79691_expd;
map_illuminaHuman2SYMBOL(GSE79691_expd,"illuminaHumanv4.db")->GSE79691_expd;
#---------------------------------------------------
#---以下几项自行练习
#1）数据集的结构
#2）分组信息
#3）差异表达分析
#4）功能富集分析



























