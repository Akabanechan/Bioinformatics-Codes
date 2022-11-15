library(ggplot2)

degs = 'D:/Download/GEO analysis/20221101/DEG/all_DEGs.csv'
workdir = 'D:/Download/GEO analysis/20221101/DEG'

setwd(workdir)

###设定阈值
fdr = 0.05
logFC = 1

#读入数据
degs = read.csv(degs, header = T, row.names = 1,check.names = F)
degs = data.frame(degs[,2],degs[,6])
colnames(degs) = c('logFC','FDR')

#增加上下调
degs[fdr > degs[,"FDR"]  &  degs[,"logFC"] >= logFC, ncol(degs)+1] = "Up"
degs[fdr > degs[,"FDR"]  & -logFC >= degs[,"logFC"], ncol(degs)] = "Down"
degs[degs[,"FDR"] >= fdr | logFC > abs(degs[,"logFC"]) , ncol(degs)] = "Normal"
colnames(degs)[ncol(degs)] = "Regulate"


temp1 = degs[,c("FDR","logFC","Regulate")]
temp1[,"FDR"] = -log10(temp1$FDR)
colnames(temp1)=c("-log10FDR","logFC","Regulate")
temp1$Regulate=factor(temp1$Regulate, levels=c("Up","Down","Normal"), order=T)

#绘图
P_volcano=ggplot(temp1,aes(x=temp1$logFC,y=temp1[,"-log10FDR"]))+
  geom_point(aes(color=temp1$Regulate))+
  scale_color_manual(values =c("Up" = "red", "Down" = "blue", "Normal" = "grey"))+
  labs(x="log2FC",y="-log10FDR")+
  geom_hline(yintercept=-log10(fdr),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  xlim(-5,5)+
  theme(plot.title = element_text(size = 25,face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, face = "bold"),
        legend.position = 'right',
        legend.key.size=unit(0.8,'cm'),
        
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),
        
        axis.title.x = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
        
        panel.background = element_rect(fill = "transparent",colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = "black"))


pdf(file=paste(workdir,"/","volcano.pdf",sep=""), width=18, height=15)
print(P_volcano)
dev.off()

png(file=paste(workdir,"/","volcano.png",sep=""), width=5400,height=4500,res=300)
print(P_volcano)
dev.off()
