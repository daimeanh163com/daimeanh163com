
#install.packages("corrplot")

setwd("C:\\Users\\lexb4\\Desktop\\GEOimmune\\10.corHeatmap")      #设置工作目录
rt=read.table("CIBERSORT.filter.txt",sep="\t",header=T,row.names=1,check.names=F)

library(corrplot)
pdf("corHeatmap.pdf",height=13,width=13)              #保存图片的文件名称
corrplot(corr=cor(rt),
         method = "color",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         )
dev.off()

