

setwd("C:\\Users\\biowolf\\Desktop\\ssGSEA\\07.ssGSEA")          #���ù���Ŀ¼
inputFile="symbol.txt"                                         #�����ļ�
gmtFile="immune.gmt"                                           #GMT�ļ�

#���ð�
library(GSVA)
library(limma)
library(GSEABase)

#��ȡ�����ļ������������ļ�����
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())

#ssgsea����
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#����ssGSEA score��������
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#��ssGSEA score���н���
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut.txt",sep="\t",quote=F,col.names=F)
