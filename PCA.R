
workdir = "D:/Download/GEO analysis/20221101/DEG"
exp_file = "D:/Download/GEO analysis/20221101/download/gene_count_matrix.csv"
Compare = "D:/Download/GEO analysis/20221101/DEG/group_file.txt"

setwd(workdir)


#载入数据
data = read.csv(exp_file, header = T, row.names = 1,check.names = F)

#去除总和为0的行
data <- data[rowSums(data)>0,]
#去掉方差为0的行
data <- data[apply(data, 1, var)!=0,]



m.mad <- apply(data,1,mad)
# dataVar <- data[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
# dataVar <- data[which(m.mad > 0),]
# data <- as.data.frame(t(dataVar))

data2 = data[rev(order(m.mad)),]
#data = data2[1:500,]
data = data.frame(t(data))

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



#组合
expdata = data.frame(group[,1],data[match(rownames(group),rownames(data)),])
colnames(expdata)[1] = 'group'

pca1 <- prcomp(expdata[,2:ncol(expdata)],center = TRUE,scale. = TRUE)

print(str(pca1))
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
p.pca1
ggsave(p.pca1,filename = "PCA_all.pdf")


library(factoextra)
# 碎石图展示每个主成分的贡献
fviz_eig(pca1, addlabels = TRUE)
# #PCA样品聚类信息展示
# fviz_pca_ind(pca1, repel=T) 
# #根据分组上色并绘制

# fviz_pca_ind(pca1, col.ind=data$group, mean.point=F, addEllipses = T, legend.title="Groups")

#根据分组上色并绘制95%置信区间,最后需要AI对标注文字位置进行修改
ind.p = fviz_pca_ind(pca1, 
                     title="",
                     col.ind=expdata$group, 
                     geom.ind = c("point"),
                     mean.point=F, 
                     #palette = c("#00AFBB", "#FC4E07"),
                     addEllipses = T, 
                     legend.title="Groups", 
                     ellipse.type="confidence",
                     ellipse.level=0.95,
                     pointsize = 10)+
  xlab(xlab1)+
  ylab(ylab1)+
  theme_bw()+
  xlim(-150,120)+
  #geom_text(aes(label=rownames(df1)),size=8,nudge_x = 25)+#调整样本标签的位置和大小
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 20,colour = "black"),
        axis.text.y = element_text(size = 20,colour = "black"),
        legend.key = element_blank(),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 25,face = "bold"),
        axis.title.y = element_text(size = 25,face = "bold"),
        legend.text = element_text(size = 15),#调整图例文字大小
        legend.position = c(0.01,0.99),legend.justification = c(0.01,0.99),#调整图例位置
        text = element_text(family = "serif"))


ind.p




####ggplot的方法画图
pca_df$Group = c(rep("CMV",5),rep("Control",5))
ggplot(data = pca_df,aes(x=PC1,y=PC2,color=Group))+
  geom_point(size=8)+
  xlab(paste("PC1: ","(", x_per,"%"," variance",")",sep=""))+
  ylab(paste("PC2: ","(", y_per,"%"," variance",")",sep=""))+
  theme_bw()+
  xlim(-200,120)+
  geom_text(aes(label=rownames(pca_df)),size=8,nudge_x = 25)+#调整样本标签的位置和大小
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 20,colour = "black"),
        axis.text.y = element_text(size = 20,colour = "black"),
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 25,face = "bold"),
        axis.title.y = element_text(size = 25,face = "bold"),
        legend.text = element_text(size = 15),#调整图例文字大小
        legend.position = c(0.01,0.99),legend.justification = c(0.01,0.99),#调整图例位置
        text = element_text(family = "serif"))+
  geom_vline(xintercept=0,lty=3,col="black",lwd=0.5)+
  geom_hline(yintercept=0,lty=3,col="black",lwd=0.5)+
  stat_ellipse(level = 0.95)





pca1$x<-data.frame(pca1$x)
# 设置颜色，有几个分组就写几个颜色
colors <- c("red","blue","green","black")
colors <- colors[as.numeric(as.factor(group[,1]))]
# 设置点形状，仅供参考
# shape<-16:18
# shape<-shape[as.numeric(as.factor(dfGroup[,1]))]

#计算PC值，并替换列名，用来替换坐标轴上的标签
pVar <- pca1$sdev^2/sum(pca1$sdev^2)
pVar = round(pVar,digits = 3)
colnames(pca1$x) = c(
  paste0("PC1 (",as.character(pVar[1] * 100 ),"%)"),
  paste0("PC2 (",as.character(pVar[2] * 100 ),"%)"),
  paste0("PC3 (",as.character(pVar[3] * 100 ),"%)"),
  "PC4",
  "PC5",
  "PC6",
  "PC7",
  "PC8",
  "PC9"
)

#BiocManager::install("scatterplot3d")
library(scatterplot3d)
# 绘图
s3d <- scatterplot3d(pca1$x[,1:3],
                     pch = 16,       # 点形状
                     color=colors,   # 点颜色
                     cex.symbols = 2,# 点大小
                     angle = 60,
                     type = "h"
)


# 设置图例
legend(s3d$xyz.convert(7.5, 3, 4.5),
       legend = unique(group[,1]),
       col =  c("red","blue","green","black"),
       pch = 16,
       inset = -0.1,
       bty = "n" ,
       bg = "transparent",
       xpd = T
)


# # 设置文字标注
# text(s3d$xyz.convert(pca1$x[,c(1,2,3)] + 2),
#      labels = row.names(pca1$x),
#      cex = 0.8,col = "black")
