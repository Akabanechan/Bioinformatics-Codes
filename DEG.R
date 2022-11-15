library(data.table)
library(DESeq2)

exp_file = "D:/Download/GEO analysis/20221101/download/gene_count_matrix.csv" #表达矩阵
Compare = "D:/Download/GEO analysis/20221101/DEG/group_file.txt"

workdir = "D:/Download/GEO analysis/20221101/DEG"

setwd(workdir)

###设定阈值
fdr = 0.05
logFC = 1

datExpr = read.csv(exp_file, header = T, row.names = 1,check.names = F)
datExpr <- datExpr[rowSums(datExpr)>0,]
datExpr <- datExpr[apply(datExpr, 1, var)!=0,]

m.mad <- apply(datExpr,1,mad)
# dataVar <- data[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
# dataVar <- data[which(m.mad > 0),]
# data <- as.data.frame(t(dataVar))

datExpr2 = datExpr[rev(order(m.mad)),]
datExpr = datExpr2[1:500,]
data = data.frame(t(datExpr))

for (i in 1:nrow(data)){
  rownames(data)[i] = substr(rownames(data)[i],1,nchar(rownames(data)[i])-6)
}

data = data.frame(t(data))

#载入分组
group = read.table(Compare, row.names=1, header=T, comment.char = "", check.names=F)
#group$group = as.character(group$group)

for (i in 1:nrow(group)){
  rownames(group)[i] = paste(substr(rownames(group)[i],1,2),'_',substr(rownames(group)[i],4,nchar(rownames(group)[i])),sep = '')
}

#筛选
#group[,2] = grepl(pattern = '6h',rownames(group))
#group = group[group[,2]==TRUE,1]

group[1:10,1] = substr(group[1:10,1],1,5)
group[11:20,1] = substr(group[11:20,1],1,4)
#group[21:30,1] = substr(group[21:30,1],1,7)
#group = group[1:20,]


#
group[1:10,2] = 12
group[11:20,2] = 6
group = group[1:20,]
colnames(group)[2] = 'type'
group$group <- factor(group$group)

data = data[,1:20]

#差异分析
dds = DESeqDataSetFromMatrix(countData = data, colData = group, design = ~ group)
dds = DESeq(dds)

res = results(dds, contrast = c('group', 'mv12h', 'mv6h'))
apply(res,2, function(x) sum(is.na(x)))
#NA值是因为使用results()提取差异分析结果时，大于alpha值（这里是默认的0.1）的矫正后p-value都会被当做是NA。因此，我们将这些padj值都设为1
res$padj[is.na(res$padj)] <- 1
table(is.na(res$padj))
summary(res)

deg.ord1 <- res[order(res$padj, decreasing = F),]
deg.sig1 <- subset(deg.ord1,abs(deg.ord1$log2FoldChange)>logFC & deg.ord1$padj<fdr)

write.csv(res,"all_DEGs.csv",quote = F)
write.csv(deg.sig1,"DEGs_sig.csv",quote = F)

plotMA(res, alpha = 0.05, colSig = 'red', colLine = 'skyblue')

filter_up = subset(res, pvalue < 0.05 & log2FoldChange > logFC) #过滤上调基因
filter_down = subset(res, pvalue < 0.05 & log2FoldChange < -logFC) #过滤下调基因

print(paste('差异基因上调数量', nrow(filter_up)))
print(paste('差异基因下调数量', nrow(filter_down)))

write.csv(filter_up,'DEGs_up.csv',quote = F)
write.csv(filter_down,'DEGs_down.csv',quote = F)



