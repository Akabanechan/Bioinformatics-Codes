library(data.table)
library('ggplot2')
workdir = 'D:/Download/GEO analysis/新建文件夹'
setwd(workdir)

#读入表型文件
pd <- fread("TcgaTargetGTEX_phenotype.txt.gz")
pd <- as.data.frame.matrix(pd)
#pd_N <- pd[pd[,4] == "Brain",]
pd_N = pd
pd_n <- as.character(pd_N$sample)
for (i in 1:nrow(pd_N)){
  if (substr(pd_N[i,1],14,14) == 0 & substr(pd_N[i,1],1,4) == 'TCGA'){
    pd_N[i,8] = 'Tumor'}
  else{
    pd_N[i,8] = 'Normal'}
}
colnames(pd_N)[8] = 'group'
colnames(pd_N)[4] = 'organ'




#读入表达矩阵文件
exp <- fread("TcgaTargetGtex_gene_expected_count.gz")
exp <- as.data.frame.matrix(exp)
exp_n <- exp[,colnames(exp) %in% pd_n]
exp_p = exp_n[4,]
gene = rownames(exp_p)
exp_p = data.frame(t(exp_p))
exp_p[,2] = rownames(exp_p)
colnames(exp_p) = c('expr','sample')



#组合
melted = data.frame(pd_N[match(exp_p[,2],pd_N[,1]),],exp_p[,1])
melter = data.frame(melted[,1],melted[,4],melted[,8:9])
colnames(melter) = c('sample','organ','group','expr')
for (i in 1:nrow(melter)){
  if (melter[i,2]=='Sympathetic\xcaNervous System'){
    melter[i,2] = 'Sympathetic Nervous System'}
  else if(melter[i,2]=='') {
    melter[i,2] = 'Blank'}
}



pdf(file=paste(workdir,"/","表达谱.pdf",sep=""), width=24, height=18, onefile = FALSE)
ggplot(melter)+
  geom_boxplot(aes(organ,expr,fill = group))+
  scale_fill_manual(values = c('Tumor' = '#d6503a','Normal' = '#5488ef'))+
  theme_bw()+
  theme(axis.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()
