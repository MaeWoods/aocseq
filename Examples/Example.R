setwd("~/aocseqExample/aocseq-main/")
roxygen2::roxygenise()
library(devtools)
build()
options(repos="https://CRAN.R-project.org")

##Install BiocParrallel
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocParallel")

library(BiocParallel)

install.packages("~/aocseqExample/aocseq_0.1.0.tar.gz", repos = NULL, type = "source")
library(aocseq)
library(Seurat)
library(SeuratObject)
BiocManager::install('glmGamPoi')
library('glmGamPoi')

######Citeseq example
setwd("~/aocseqExample")
gex.path=c("AOCSEQ_Example/SeqData/CD5")
marker.gene=c("IFNG","CRTAM","GZMB","CCR7","CCL1","SELL","CD28","CD27")
save.dir="AOCSEQ_Example/RDS/CD5.rds"
my.sample.name=c("CD5")
my.demultiplex.index = c(8, 9, 10)
my.nameshashtags = c('donor_1', 'donor_2','donor_3')
my.n.samples.ht = 3
my.hashtags = 3
my.sample.name=c("CD5")

CellData=CombineData(gex.path,
                     marker.gene,
                     file.saved=save.dir,
                     sample.name=my.sample.name,
                     demultiplex=TRUE,
                     demultiplex.index=my.demultiplex.index,
                     nameshashtags=my.nameshashtags,
                     hashtags=my.hashtags,
                     n.samples.ht = my.n.samples.ht,
                     nvariable_features=500,
                     QC_plots = TRUE)

my.ident.list=c("orig.ident","Threshold_CCL1")
my.save.dir="AOCSEQ_Example"

UMAPReduce(
  Clonal_Obs[[2]],
  save.dir=my.save.dir,
  ident.list=my.ident.list
)
