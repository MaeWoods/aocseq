#########################################################
####Load dependencies                               #####
#########################################################
library(ggplot2)
library(gplots)
library(dplyr)
library(cowplot)
library(devtools)
#install.packages("viridis")
library(viridis)
library(devtools)
library(purrr)
library(RVenn)
library(rJava)
library(UpSetR)
library(tidyverse)
library(venneuler)
library(grid)
library(crone)

markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
upsetnames <- c("IFNGCMVCD4","CRTAMCMVCD4","CD69CMVCD4","TNFCMVCD4","CD70CMVCD4","TNFRSF9CMVCD4",
                "IFNGCMVCD8","CRTAMCMVCD8","CD69CMVCD8","TNFCMVCD8","CD70CMVCD8","TNFRSF9CMVCD8",
                "IFNGEBVCD8","CRTAMEBVCD8","CD69EBVCD8","TNFEBVCD8","CD70EBVCD8","TNFRSF9EBVCD8")

Nsetsin=3
Nmarkersin=6

upsetpath="../Figures/Genes_lower_3hour.pdf"
combination_path = "../SupplementaryTables/SummaryGenes/Genes_lower_3hour"
path_list = list("../SupplementaryTables/Genes_lower_3hour/CMV_IFNG/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/CMV_CRTAM/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/CMV_CD69/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/CMV_TNF/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/CMV_CD70/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/CMV_TNFRSF9/CD4genes_upregulated.csv",
                 
                 "../SupplementaryTables/Genes_lower_3hour/CMV_IFNG/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/CMV_CRTAM/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/CMV_CD69/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/CMV_TNF/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/CMV_CD70/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/CMV_TNFRSF9/CD8genes_upregulated.csv",
                 
                 "../SupplementaryTables/Genes_lower_3hour/EBV_IFNG/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/EBV_CRTAM/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/EBV_CD69/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/EBV_TNF/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/EBV_CD70/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_lower_3hour/EBV_TNFRSF9/CD8genes_upregulated.csv")

SummarySignature(Nsetsin,Nmarkersin,path_list,upsetnames,upsetpath,combination_path)



upsetpath="../Figures/Genes_upper_3hour.pdf"
combination_path = "../SupplementaryTables/SummaryGenes/Genes_upper_3hour"
path_list = list("../SupplementaryTables/Genes_upper_3hour/CMV_IFNG/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/CMV_CRTAM/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/CMV_CD69/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/CMV_TNF/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/CMV_CD70/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/CMV_TNFRSF9/CD4genes_upregulated.csv",
                 
                 "../SupplementaryTables/Genes_upper_3hour/CMV_IFNG/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/CMV_CRTAM/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/CMV_CD69/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/CMV_TNF/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/CMV_CD70/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/CMV_TNFRSF9/CD8genes_upregulated.csv",
                 
                 "../SupplementaryTables/Genes_upper_3hour/EBV_IFNG/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/EBV_CRTAM/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/EBV_CD69/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/EBV_TNF/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/EBV_CD70/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_upper_3hour/EBV_TNFRSF9/CD8genes_upregulated.csv")

SummarySignature(Nsetsin,Nmarkersin,path_list,upsetnames,upsetpath,combination_path)




upsetpath="../Figures/Genes_all6hour_new4.pdf"
combination_path = "../SupplementaryTables/SummaryGenes/Genes_all6hour"
path_list = list("../SupplementaryTables/Genes_all6hour/CMV_IFNG/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/CMV_CRTAM/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/CMV_CD69/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/CMV_TNF/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/CMV_CD70/CD4genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/CMV_TNFRSF9/CD4genes_upregulated.csv",
                 
                 "../SupplementaryTables/Genes_all6hour/CMV_IFNG/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/CMV_CRTAM/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/CMV_CD69/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/CMV_TNF/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/CMV_CD70/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/CMV_TNFRSF9/CD8genes_upregulated.csv",
                 
                 "../SupplementaryTables/Genes_all6hour/EBV_IFNG/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/EBV_CRTAM/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/EBV_CD69/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/EBV_TNF/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/EBV_CD70/CD8genes_upregulated.csv",
                 "../SupplementaryTables/Genes_all6hour/EBV_TNFRSF9/CD8genes_upregulated.csv")

SummarySignature(Nsetsin,Nmarkersin,path_list,upsetnames,upsetpath,combination_path)


SummarySignature <- function(Nsetsin,Nmarkersin,path_list,upsetnames,upset_path,combination_path){
  
StoreGeneLists = list(rep(0,Nsetsin*Nmarkersin))
StoreGeneLists <- vector(mode = "list", length = Nsetsin*Nmarkersin)
StoreGeneLists[1] = read_csv(path_list[[1]])$x
for(s in 2:length(path_list)){
  StoreGeneLists[s] = list(read_csv(path_list[[s]])$x)
}

unionGeneSet=StoreGeneLists[1][[1]]
for(h in 2:length(StoreGeneLists)){
  unionGeneSet = union(unionGeneSet,StoreGeneLists[h][[1]])
}

allgenes=unionGeneSet
upsetdata <- data.frame(matrix(ncol = Nsetsin*Nmarkersin, nrow = length(allgenes)))

for(q in 1:(Nsetsin*Nmarkersin)){
for(k in 1:length(allgenes)){
  if(is.na(match(allgenes[k],StoreGeneLists[[q]]))){
    upsetdata[k,q]=0
  }
  else{
    upsetdata[k,q]=1
  }

}
}

UpsetDataForC = subset(upsetdata,
(rowSums(upsetdata[,1:Nmarkersin])>0)&
(rowSums(upsetdata[,(Nmarkersin+1):(2*Nmarkersin)])>0)&
(rowSums(upsetdata[,(2*Nmarkersin+1):(3*Nmarkersin)])>0))

CreateMatrix <- rep(0,(dim(UpsetDataForC)[1]*dim(UpsetDataForC)[2]))
for(j in 1:(dim(UpsetDataForC)[1])){
  for(k in 1:(dim(UpsetDataForC)[2])){
    CreateMatrix[dim(UpsetDataForC)[2]*(j-1) + k] = UpsetDataForC[j,k]
  }
}

UniqueGenes=dim(UpsetDataForC)[1]
Genedatain=CreateMatrix
myinteger=GetOptimalMarkers(UniqueGenes,Nmarkersin,Nsetsin,Genedatain)
s1=t(unname(UpsetDataForC[(myinteger+1),]))[,]
maximalGeneSetindexes=which(s1==1,s1)
maximalGeneSet=StoreGeneLists[maximalGeneSetindexes[1]][[1]]
for(h in 2:length(maximalGeneSetindexes)){
  maximalGeneSet = intersect(maximalGeneSet,StoreGeneLists[maximalGeneSetindexes[h]][[1]])
}
names(upsetdata) <- upsetnames

pdf(upsetpath,width=25,height=10)
p=upset(upsetdata,sets=names(upsetdata),
      nsets = (Nsetsin*Nmarkersin), nintersects=NA, number.angles = 45, keep.order=TRUE, point.size = 3.5, line.size = 1,
      mainbar.y.label = "Number of common genes", sets.x.label = "Total number of genes", text.scale = c(1.2, 1.3, 1, 1, 1.2, 1.4)
)
print(p)
dev.off()

write.csv(maximalGeneSet,paste(combination_path,".csv",sep=""))

}
