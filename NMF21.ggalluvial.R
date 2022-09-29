


#���ð�
library(ggalluvial)
library(ggplot2)
library(dplyr)

clusterFile="cluster.txt"         #���ͽ���ļ�
subtypeFile="immune type.txt"       #���߷����ļ�
setwd("~/desktop/zpx/2")     #���ù���Ŀ¼

#��ȡ���ͽ���ļ�
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#��ȡ���߷����ļ�
subtype=read.table(subtypeFile, header=T, sep="\t", check.names=F, row.names=1)
#rownames(subtype)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(subtype))

#��Ʒȡ����
sameSample=intersect(row.names(subtype), row.names(cluster))
subtype=subtype[sameSample,,drop=F]
subtype[,1]=gsub(".+\\((.+?)\\)", "\\1", subtype[,1])
colnames(subtype)=c("Immune subtype")
cluster=cluster[sameSample,,drop=F]
data=cbind(cluster, subtype)

#׼��ɣ��ͼ�����ļ�
corLodes=to_lodes_form(data, axes=1:ncol(data), id = "Cohort")

#�õ�����ļ�
pdf(file="ggalluvial.pdf", width=5, height=10)
mycol=rep(c("#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 #��aes.flow����������ɫ��forward˵����ɫ��ǰ�����״ͼһ�£�backward˵���ͺ������״ͼһ�¡�
  	 geom_flow(width = 2/10,aes.flow = "backward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 #size=3���������С
	 geom_text(stat = "stratum", size = 3,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ��������
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()

