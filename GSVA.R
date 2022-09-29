
setwd("C:\\Users\\lexb4\\Desktop\\GSVA\\03.GSVA")
inputFile="input.txt"
gmtFile="c2.cp.kegg.v6.2.symbols.gmt"

library(GSVA)
library(limma)
library(GSEABase)

rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=normalizeBetweenArrays(mat)

c3gsc2=getGmt( gmtFile, 
               collectionType=BroadCollection(category="c3"), 
               geneIdType=SymbolIdentifier())
gsvaOut=gsva(mat, 
             c3gsc2, 
             min.sz=10, 
             max.sz=500, 
             verbose=TRUE,
             parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaOut),gsvaOut)
write.table(gsvaOut,file="gsvaOut.txt",sep="\t",quote=F,col.names=F)




