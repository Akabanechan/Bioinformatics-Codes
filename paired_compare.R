library(ggplot2)
library(ggpubr)
library(tidyverse)
library(rstatix)

#exp_file = "D:/Download/GEO analysis/20221101/download/gene_count_matrix.csv"
expr_file = "D:/Download/GEO analysis/20221101/download/genes.TPM_Cross.matrix"
#trait_file = "D:/Download/GEO analysis/20221101/WGCNA/trait.txt"
Compare = "D:/Download/GEO analysis/20221101/DEG/group_file.txt"
workdir = 'D:/Download/GEO analysis/20221101/DEG'
setwd(workdir)

#读入表型文件

datExpr = read.table(expr_file,header=T,row.names=1,stringsAsFactors = F,comment.char = "",check.names=F)
#datExpr = read.csv(exp_file, header = T, row.names = 1,check.names = F)
datExpr <- datExpr[rowSums(datExpr)>0,]
datExpr <- datExpr[apply(datExpr, 1, var)!=0,]

m.mad <- apply(datExpr,1,mad)
# dataVar <- data[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
# dataVar <- data[which(m.mad > 0),]
# data <- as.data.frame(t(dataVar))

datExpr2 = datExpr[rev(order(m.mad)),]
datExpr = datExpr2[1:500,]
data = data.frame(t(datExpr))

data = log2(data+1)

for (i in 1:nrow(data)){
  rownames(data)[i] = substr(rownames(data)[i],1,nchar(rownames(data)[i])-6)
}

#载入分组
group = read.table(Compare, row.names=1, header=T, comment.char = "", check.names=F)
group$group = as.character(group$group)

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

#组合
expdata = data.frame(group[,1],data[match(rownames(group),rownames(data)),])
colnames(expdata)[1] = 'group'

expdata = expdata[1:20,1:2]
colnames(expdata)[2] = 'expr'

for (i in 1:10){
  expdata[i,3] = i
  expdata[i+10,3] = i
}
colnames(expdata)[3] = 'pair'

#expdata$group <- factor(expdata$group)

#使用 ggplot2 包绘制箱线图
p <- ggplot(expdata, aes(x = group, y = expr)) +
  geom_boxplot(aes(fill = group), show.legend = FALSE, width = 0.6) +  #绘制箱线图
  scale_fill_manual(values = c('#FE7280', '#AC88FF')) +  #箱线图的填充色
  geom_point(size = 2) +  #绘制样本点
  geom_line(aes(group = pair), color = 'gray', lwd = 0.5) +  #绘制配对样本间连线
  ##以下是ggplot2的主题设置，修改边框、背景、标题、字体等
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black', linewidth = 1), panel.background = element_blank(), 
        plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 20, color = 'black'), axis.title = element_text(size = 20, color = 'black')) +
  labs(x = '', y = 'log2(count+1)', title = '', subtitle = 'mv_12h vs mv_6h')

my_compare = list(c('mv12h','mv6h'))
p + stat_compare_means(comparisons = my_compare,
                       label = 'p.signif',
                       method = 't.test')

p
