

exp_file = 'D:/Download/GEO analysis/新建文件夹/exp_matrix.csv'

data = read.csv(exp_file, header = T, row.names = 1,check.names = F)
data[1:4,1:4]
data = data.frame(t(data))
data[1:4,1:4]
data <- data[rowSums(data)>0,]
data <- data[apply(data, 1, var)!=0,]

m.mad <- apply(data,1,mad)
data2 = data[rev(order(m.mad)),]
data = data2[1:5,]
colnames(data) = c('a','b','c','d','e')
#转换为列为基因
data = data.frame(t(data))
#标准化
data = scale(data)
#pearson
cor_pearson = cor(data, method = 'pearson')

library(corrplot)

corrplot(cor_pearson, method = 'number', number.cex = 0.8, diag = FALSE, tl.cex = 0.8)
corrplot(cor_pearson, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n')

cor.test( ~ a + b, 
          data=data,
          method = "spearman",
          continuity = FALSE,
          conf.level = 0.95)
plot(a ~ b, 
     data=data, 
     pch=16)
