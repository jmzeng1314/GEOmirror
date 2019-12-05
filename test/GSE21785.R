rm(list = ls())
options(stringsAsFactors = F)
library(GEOmirror)
geoChina('GSE21785')
load('GSE21785_eSet.Rdata')
exp <- exprs(gset[[1]])
exp[1:4,1:4]
pd <- pData(gset[[1]])
anno = gset[[1]]@annotation
group_list =c(rep("Tubulus",6),rep("Glomerulus",6))
group_list=factor(group_list,levels = c("Tubulus","Glomerulus"))
boxplot(exp,las=2,col=group_list)
