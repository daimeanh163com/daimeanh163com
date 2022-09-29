
library(limma)
workDir="C:\\Users\\lexb4\\Desktop\\m6A\\10.ConsensusClusterPlus"        #����Ŀ¼
setwd(workDir)
rt=read.table("m6Aexp.txt",sep="\t",header=T,check.names=F)              #��ȡ�����ļ�

#һ��������ֶ��У�ȡ��ֵ
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]


#����
maxK=9
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123456,
              plot="png")

#������
clusterNum=2                  #�ּ��࣬�����жϱ�׼�ж�
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster,file="cluster.txt",sep="\t",quote=F,col.names=F)