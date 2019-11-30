# GEOmirror

A package for researchers in mainland China to downlaod the GEO dataset,  which is just a **replacement** of getGEO function from GEOquery package.

### just one function: geoChina

Install the **development version** from Github:

```r
library(devtools)
install_github("jmzeng1314/GEOmirror")
library(GEOmirror)
```

Then use it to download GEO dataset, as below :

```r
eSet=geoChina('GSE1009') 
eSet=geoChina('GSE27533') 
eSet=geoChina('GSE95166') 
```

### easy download GEO dataset and go through downstream analysis 

Once you download the ExpressionSet of GEO dataset, you can access the expression matrix and phenotype data: 

```r
## download GSE95166 data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95166
#eSet=getGEO('GSE95166', destdir=".", AnnotGPL = F, getGPL = F)[[1]]
library(GEOmirror)
eSet=geoChina('GSE95166')
eSet
eSet=eSet[[1]]


probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)

## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])
# https://www.ncbi.nlm.nih.gov/pubmed/31430288

groupList=factor(c(rep('npc',4),rep('normal',4)))
table(groupList)
eSet@annotation
# GPL15314	Arraystar Human LncRNA microarray V2.0 (Agilent_033010 Probe Name version)

```

The simple example for downstream analysis of expression profiling by array:

```r
genes_expr=probes_expr
library("FactoMineR")
library("factoextra")
dat.pca <- PCA(t(genes_expr) , graph = FALSE)
dat.pca
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = groupList,
             addEllipses = TRUE,
             legend.title = "Groups"
)
library(limma)
design=model.matrix(~factor(groupList))
design
fit=lmFit(genes_expr,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,n=Inf)
head(DEG)
# We observed that 2107 lncRNAs were upregulated
# while 2090 lncRNAs were downregulated by more than 2-fold,
# NKILA among these downregulated lncRNAs (Fig 1A, GSE95166).

## for volcano plot
df=DEG
attach(df)
df$v= -log10(P.Value)
df$g=ifelse(df$P.Value>0.05,'stable',
            ifelse( df$logFC >1,'up',
                    ifelse( df$logFC < -1,'down','stable') )
)
table(df$g)
df$name=rownames(df)
head(df)
ggpubr::ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
          label = "name", repel = T,
          label.select =head(rownames(df)),
          palette = c("#00AFBB", "#E7B800", "#FC4E07") )
detach(df)


x=DEG$logFC
names(x)=rownames(DEG)
cg=c(names(head(sort(x),100)),
     names(tail(sort(x),100)))
cg
library(pheatmap)
n=t(scale(t(genes_expr[cg,])))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
ac=data.frame(groupList=groupList)
rownames(ac)=colnames(n)  
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac)
```

Actually, I should change the probe ID to gene names, but it's not this package should do, so I just omit it.

### warnings

- 1.only the expression profiling by array datasets will be offered by our package.
- 2.it's free for all of us, so we cann't make sure the internet connection will always be better than GEOquery. 

