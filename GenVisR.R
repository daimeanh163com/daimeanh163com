

setwd("~/desktop/m")

library(GenVisR)
rt=read.table("lowmRNAsi.txt",header=T,sep="\t",check.names=F,quote="")
pdf(file="waterfall2.pdf",height=6,width=10)
waterfall(rt,fileType = 'Custom',mainRecurCutoff = 0.06, variant_class_order = levels(factor(rt$variant_class)))
dev.off()


