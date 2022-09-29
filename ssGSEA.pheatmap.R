

library(pheatmap)
setwd("C:\\Users\\biowolf\\Desktop\\ssGSEA\\09.heatmap")      #设置工作目录
rt=read.table("ssgseaOut.txt",sep="\t",header=T,row.names=1,check.names=F)    #读取文件

Type=read.table("cluster.txt",sep="\t",check.names=F,header=F)
Type[,2]=paste0("Cluster",Type[,2])
Type=Type[order(Type[,2]),]
rt=rt[,as.vector(Type[,1])]

cluster=as.data.frame(Type[,2])
row.names(cluster)=Type[,1]
colnames(cluster)="Cluster"

pdf("heatmap.pdf",height=5,width=9)
pheatmap(rt, annotation=cluster, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3)
dev.off()



#定义每个cluster分组信息，需要修改
Cluster1="Immunity_H"
Cluster2="Immunity_M"
Cluster3="Immunity_L"

#输出分组结果
a=c()
a[Type[,2]=="Cluster1"]=Cluster1
a[Type[,2]=="Cluster2"]=Cluster2
a[Type[,2]=="Cluster3"]=Cluster3
clusterOut=cbind(Type,a)
write.table(clusterOut,file="cluster.Immunity.txt",sep="\t",quote=F,col.names=F,row.names=F)

