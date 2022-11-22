rm(list = ls())  ## 魔幻操作，一键清空~

exp_file = 'D:/Download/GEO analysis/20221101/download/gene_count_matrix.csv'
exp_file = 'D:/Download/GEO analysis/20221101/DEG/gene_trans.csv'
Compare = "D:/Download/GEO analysis/20221101/DEG/group_file.txt"
workdir = 'D:/Download/GEO analysis/20221101/DEG'
setwd(workdir)

#去批次分析 
count <- read.csv(exp_file, header = T, row.names = 1,check.names = F)
#过滤掉低表达的基因
Expr <- count[rowSums(count)>1,]
Expr = data.frame(t(count))
Expr = data.frame(t(data))
##载入样品信息，含批次
batch <- read.table(Compare, row.names=1, header=T, comment.char = "", check.names=F)
batch = group
batch[1:10,2] = 1
batch[11:20,2] = 2
batch[21:25,2] = 1
batch[26:30,2] = 2
colnames(batch)[2] = 'batch'
batch[,1] <- as.factor(batch$group)##type表示的是生物学上的不同处理
batch[,2] <- as.factor(batch$batch)##batch表示不同的批次信息
#row.names(batch) <- batch$sample

library(FactoMineR)## 不做校正的PCA分析,没有请先安装
library(factoextra)
pre.pca <- PCA(t(Expr),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = batch$batch,
             habillage=batch$batch,
             palette = "Dark2",
             addEllipses = TRUE,
             legend.title="Group")
#fviz_pca_ind(pre.pca)
##limma包自带的removeBatchEffect函数 消除批次效应
library(limma)
design <- model.matrix(~group,data = batch)
rB_Expr <- removeBatchEffect(Expr,batch = batch$batch,design = design)
##差异分析
fit <- lmFit( rB_Expr, design )
fit2 <- eBayes(fit)
allDiff3 = topTable(fit2,adjust = 'fdr',coef = 2,number = Inf)
diffLab3 <- subset(allDiff3,abs(logFC)>1 & adj.P.Val < 0.05)
dim(diffLab3)###查看差异基因的个数
###################
##PCA分析
af2.pca <- PCA(t(rB_Expr),graph = FALSE)
fviz_pca_ind(af2.pca,
             geom= "point",
             col.ind = batch$batch,
             habillage=batch$batch,
             palette = "Dark2",
             addEllipses = TRUE,
             legend.title="Group"  )


##消除批次效应  sva包中Combat函数
design <- model.matrix(~group,data = batch)

library(sva)
combat_Expr <- ComBat_seq(as.matrix(Expr),batch = batch$batch)
##combat_Expr就是校正后的数据
######################################
#差异分析
fit <- lmFit(combat_Expr,design)
fit2 <- eBayes(fit)
allDiff2 = topTable(fit2,adjust = 'fdr',coef = 2,number = Inf)
diffLab2 <- subset(allDiff2,abs(logFC)>1 & adj.P.Val < 0.05)
dim(diffLab2)###查看差异基因的个数
###PCA分析
library(FactoMineR)
af1.pca <- PCA(t(combat_Expr),graph = FALSE)
fviz_pca_ind(af1.pca,
             geom= "point",
             col.ind = batch$batch,
             habillage=batch$batch,
             palette = "Dark2",
             addEllipses = TRUE,
             legend.title="Group" )

write.csv(combat_Expr, file = 'gene_sva.csv', quote=F, row.names = T)
