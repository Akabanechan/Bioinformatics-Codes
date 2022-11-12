library(data.table)
library(ggplot2)
library(tidyverse)
library(rstatix)
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


#筛选健康与癌症都有的部位
a = unique(melter[,'organ'])
for (i in unique(melter[,'organ'])){
  if (length(unique(melter[melter[,'organ']==i,'group'])) == 1){
    a = a[-which(a==i)]
  }
}

melter = melter[which(melter[,'organ'] %in% a),]


#计算组内p值
melter$organ <- factor(melter$organ)

melt_p_val1 <- melter %>% 
  group_by(organ) %>% 
  wilcox_test(formula  = expr ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='organ')


#绘图
pdf(file=paste(workdir,"/","表达谱.pdf",sep=""), width=24, height=18, onefile = FALSE)
ggplot(melter)+
  geom_boxplot(aes(organ,expr,fill = group))+
  scale_fill_manual(values = c('Tumor' = '#d6503a','Normal' = '#5488ef'))+
  stat_pvalue_manual(melt_p_val1,label = '{p.signif}',tip.length = 0)+
  labs(x='Organ',y='Expr',
       caption = "Visualization by <span style='color:#DD6449'>gtrobot</span>")+
  guides(fill=guide_legend(title = 'Cancer'))+
  theme_bw()+
  theme(axis.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = c(0.7,0.1),
        legend.direction = 'horizontal')
  #scale_y_continuous(limits = c(0,20),breaks = seq(0,20,5),expand = c(0,0))
dev.off()



