library(biomaRt)

expr_file = "D:/Download/GEO analysis/20221101/download/genes.TPM_Cross.matrix"
exp_file = 'D:/Download/GEO analysis/20221101/download/gene_count_matrix.csv'
Compare = "D:/Download/GEO analysis/20221101/DEG/group_file.txt"
workdir = 'D:/Download/GEO analysis/20221101/DEG'
setwd(workdir)

#读取基因数据
datExpr = read.table(expr_file,header=T,row.names=1,stringsAsFactors = F,comment.char = "",check.names=F)
datExpr = read.csv(exp_file, header = T, row.names = 1,check.names = F)

datExpr <- datExpr[rowSums(datExpr)>0,]
datExpr <- datExpr[apply(datExpr, 1, var)!=0,]

data = data.frame(t(datExpr))

for (i in 1:nrow(data)){
  rownames(data)[i] = substr(rownames(data)[i],1,nchar(rownames(data)[i])-6)
}

#读取转换数据
human = useMart("ensembl",dataset="hsapiens_gene_ensembl")
mouse = useMart("ensembl",dataset="mmusculus_gene_ensembl")

#转换
ensembl_ids = data.frame(colnames(data))

gene_symbols<- getBM(attributes=c('ensembl_gene_id','mgi_symbol'), 
                   filters= 'ensembl_gene_id', values = ensembl_ids, mart = mouse)

data2 = data[,match(gene_symbols[,1],colnames(data))]
colnames(data2) = gene_symbols[,2]

data = data2
