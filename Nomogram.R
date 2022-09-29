

library(rms)
setwd("~/desktop/8")                         #���ù���Ŀ¼
rt=read.table("input1.txt",sep="\t",header=T,row.names=1,check.names=F)           #��ȡ�����ļ�

#���ݴ��
dd <- datadist(rt)
options(datadist="dd")

#���ɺ���
f <- cph(Surv(futime, fustat) ~ CDE_therapy+GeneExp_Subtype+age_at_initial_pathologic_diagnosis+gender+cancer_status+riskScore, x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)

#����nomogram
nom <- nomogram(f, fun=list(function(x) surv(365, x), function(x) surv(730, x), function(x) surv(1095, x)), 
    lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  

#nomogram���ӻ�
pdf(file="nomogram.pdf",height=6,width=25)
plot(nom)
dev.off()
