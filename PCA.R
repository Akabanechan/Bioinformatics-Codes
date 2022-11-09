
workdir = "D:/Download/GEO analysis/20221101/DEG"
exp_file = "D:/Download/GEO analysis/20221101/download/gene_count_matrix.csv"
Compare = "D:/Download/GEO analysis/20221101/DEG/group_file.txt"

setwd(workdir)


#载入数据
data = read.csv(exp_file, header = T, row.names = 1,check.names = F)
m.mad <- apply(data,1,mad)
dataVar <- data[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataVar <- data[which(m.mad > 0),]
data2 = data[rev(order(m.mad)),]
data = data2[1:500,]

data <- as.data.frame(t(dataVar))
data = data.frame(t(data))

for (i in 1:nrow(data)){
  rownames(data)[i] = substr(rownames(data)[i],1,nchar(rownames(data)[i])-6)
}

#载入分组
group = read.table(Compare, row.names=1, header=T, comment.char = "", check.names=F)

for (i in 1:nrow(group)){
  rownames(group)[i] = paste(substr(rownames(group)[i],1,2),'_',substr(rownames(group)[i],4,nchar(rownames(group)[i])),sep = '')
}

group[,2] = rownames(group)
group[,3] = grepl(pattern = '12h',group[,1])
group = group[group[,3]==TRUE,1:2]

expdata = data.frame(group[,1],data[match(rownames(group),rownames(data)),])
colnames(expdata)[1] = 'group'

pca1 <- prcomp(expdata[,2:ncol(expdata)],center = TRUE,scale. = TRUE)


#iris_input <- iris # 使用R自带iris数据集（150*5，行为样本，列为特征）
#rownames(iris_input) <- paste("sample",1:nrow(iris_input),sep = "") 

#pca1 <- prcomp(iris_input[,-ncol(iris_input)],center = TRUE,scale. = TRUE)

df1 <- pca1$x # 提取PC score
df1 <- as.data.frame(df1) # 注意：如果不转成数据框形式后续绘图时会报错

summ1 <- summary(pca1)
summ1


xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")


library(ggplot2)
p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color = expdata$group))+
  stat_ellipse(aes(fill = expdata$group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  #scale_fill_manual(values = c("purple","orange","pink"))+
  #scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
ggsave(p.pca1,filename = "PCA3000.pdf")


