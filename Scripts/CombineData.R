#install.packages("Seurat")
library(Seurat)
library(patchwork)
library(ggplot2)
library(gplots)
library(cowplot)
library(devtools)
library(plyr)
library(matrixStats)
library("ggplot2")
library("reshape2")
library(magrittr)
library(dplyr)    
library(readr)
library("readxl")
library(UpSetR)
library(tidyverse)
library(grid)
library(crone)
library(ggvenn)
setwd("/Users/maewoodsphd/mVSTManuscript/Rscripts")

###################################################
####Unstimulated ##################################
###################################################

markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
mask=c("TRBV.","TRAV.")
save_dir="../RDS/ThreeHourStim.rds"
dir_rna1="../outs_t1/UMO"
dir_vdj1="../vdj_t1/UMO/all_contig_annotations.csv"

dir_rna2="../outs_t1/EBV"
dir_vdj2="../vdj_t1/EBV/all_contig_annotations.csv"

dir_rna3="../outs_t1/CMV"
dir_vdj3="../vdj_t1/CMV/all_contig_annotations.csv"

dir_rna4="../outs_t1/BKV"
dir_vdj4="../vdj_t1/BKV/all_contig_annotations.csv"

ThreeHourStim=CombineData(dir_rna1,dir_rna2,dir_rna3,dir_rna4,
dir_vdj1,dir_vdj2,dir_vdj3,dir_vdj4,markerGene,mask,save_dir)

markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
mask=c("TRBV.","TRAV.")
save_dir="../RDS/SixHourStim.rds"

dir_rna1="../outs_t2/UMO"
dir_vdj1="../vdj_t2/UMO/all_contig_annotations.csv"

dir_rna2="../outs_t2/EBV"
dir_vdj2="../vdj_t2/EBV/all_contig_annotations.csv"

dir_rna3="../outs_t2/CMV"
dir_vdj3="../vdj_t2/CMV/all_contig_annotations.csv"

dir_rna4="../outs_t2/BKV"
dir_vdj4="../vdj_t2/BKV/all_contig_annotations.csv"
save_dir="../RDS/SixHourStim.rds"

SixHourStim=CombineData(dir_rna1,dir_rna2,dir_rna3,dir_rna4,
                          dir_vdj1,dir_vdj2,dir_vdj3,dir_vdj4,markerGene,mask,save_dir)

CombineData <- function(dir_rna1,dir_rna2,dir_rna3,dir_rna4,
                        dir_vdj1,dir_vdj2,dir_vdj3,dir_vdj4,markerGene,mask,save_dir){
clonal.data <- Read10X(data.dir = dir_rna1)
Total_TCellsUMO=length(clonal.data@p)
Clonal_TUMO <- CreateSeuratObject(counts = clonal.data,project = "UMO")
Clonal_TUMO[["percent.mt"]] <- PercentageFeatureSet(Clonal_TUMO, pattern = "^MT-")
Clonal_TUMO<- subset(Clonal_TUMO, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 5)
RnaStoreUMO=Clonal_TUMO@assays$RNA@counts

#VDJ list with barcodes
UMOTCRarray <- read.csv(dir_vdj1)
tcrUMO_cell=subset(UMOTCRarray,productive=="true" & is_cell=="true")
TCRlistUMO=tcrUMO_cell$barcode
joint.bcsUMO <- intersect(colnames(RnaStoreUMO), TCRlistUMO)
#Remove barcodes from list with abberant GE
tcrhashUMO=subset(UMOTCRarray,UMOTCRarray$barcode %in% joint.bcsUMO)
mvsts.UMO <- RnaStoreUMO[, joint.bcsUMO]
Clonal_TUMO=Clonal_TUMO[,joint.bcsUMO]

###remove genes that are not of interest
for(c in 1:length(mask)){
Allgenes=row.names(Clonal_TUMO)
Gset=Allgenes[grep(mask[c],Allgenes)]
NewGeneSet=setdiff(Allgenes,Gset)
IndexesGTR=match(NewGeneSet,row.names(Clonal_TUMO))
Clonal_TUMO <- Clonal_TUMO[IndexesGTR,]
}

####Genes of interest
Gene_indUMO=match(markerGene,row.names(mvsts.UMO))
CD8_UMO=match("CD8A",row.names(mvsts.UMO))
CD4_UMO=match("CD4",row.names(mvsts.UMO))
CD8cells=rep(0,length(mvsts.UMO[Gene_indUMO[1],]))
CD4cells=rep(0,length(mvsts.UMO[Gene_indUMO[1],]))

vec1=mvsts.UMO[CD4_UMO,]
vec2=mvsts.UMO[CD8_UMO,]

for(g in 1:length(CD8cells)){
  if((vec2[g]==0)&(vec1[g]>0)){
    CD4cells[g]=1
  }
  else if((vec2[g]>0)&(vec1[g]==0)){
    CD8cells[g]=1
  }
  else{}
}
rm(vec1)
rm(vec2)

Clonal_TUMO=AddMetaData(Clonal_TUMO, CD8cells, col.name = 'CD8cells')
Clonal_TUMO=AddMetaData(Clonal_TUMO, CD4cells, col.name = 'CD4cells')

rm(CD8cells)
rm(CD4cells)
#Check for errors
plot(log10(mvsts.UMO[Gene_indUMO[3],]),log10(mvsts.UMO[Gene_indUMO[2],]),xlab=paste(Gene_indUMO[3]," (log10 UMIs)",sep=""),ylab=paste(Gene_indUMO[2]," (log10 UMIs)",sep=""),col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))
plot(log10(mvsts.UMO[CD8_UMO,]),log10(mvsts.UMO[CD4_UMO,]),xlab="CD8",ylab="CD4",col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))
print("Starting scTRansform sample 1...")
Clonal_TUMO <- SCTransform(Clonal_TUMO, vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n = 7000)
Gene_indUMO=match(markerGene,row.names(Clonal_TUMO))
d1=as.matrix(Clonal_TUMO[["SCT"]]@data)[Gene_indUMO,]
alist <- 1:length(Gene_indUMO)
cutoff=lapply(alist, function(alist) quantile(d1[alist,],.975)[[1]])

Thresholds=matrix("unassigned",nrow=dim(Clonal_TUMO)[2],ncol=length(Gene_indUMO))
for(k in 1:length(Gene_indUMO)){
  vec1=Clonal_TUMO[Gene_indUMO[k],][["SCT"]]@data
for(j in 1:dim(Clonal_TUMO)[2]){
  if(vec1[j]>cutoff[k]){
    Thresholds[j,k]="UMO_high"
  }
}
  Clonal_TUMO=AddMetaData(Clonal_TUMO, Thresholds[,k], col.name = paste("Threshold_",markerGene[k],sep=""))
}

rm(mvsts.UMO)
rm(RnaStoreUMO)
rm(clonal.data)

##########################################
####EBV ##################################
##########################################

clonal.data <- Read10X(data.dir = dir_rna2)
Total_TCellsEBV=length(clonal.data@p)
Clonal_TEBV <- CreateSeuratObject(counts = clonal.data,project = "EBV")
Clonal_TEBV[["percent.mt"]] <- PercentageFeatureSet(Clonal_TEBV, pattern = "^MT-")
Clonal_TEBV<- subset(Clonal_TEBV, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 5)
RnaStoreEBV=Clonal_TEBV@assays$RNA@counts

#VDJ list with barcodes
EBVTCRarray <- read.csv(dir_vdj2)
tcrEBV_cell=subset(EBVTCRarray,productive=="true" & is_cell=="true")
TCRlistEBV=tcrEBV_cell$barcode
joint.bcsEBV <- intersect(colnames(RnaStoreEBV), TCRlistEBV)
#Remove barcodes from list with abberant GE
tcrhashEBV=subset(EBVTCRarray,EBVTCRarray$barcode %in% joint.bcsEBV)
mvsts.EBV <- RnaStoreEBV[, joint.bcsEBV]
Clonal_TEBV=Clonal_TEBV[,joint.bcsEBV]

###remove genes that are not of interest
for(c in 1:length(mask)){
  Allgenes=row.names(Clonal_TEBV)
"  Gset=Allgenes[grep(mask[c],Allgenes)]"
"  NewGeneSet=setdiff(Allgenes,Gset)"
"  IndexesGTR=match(NewGeneSet,row.names(Clonal_TEBV))"
"  Clonal_TEBV <- Clonal_TEBV[IndexesGTR,]"
}

####Genes of interest
"Gene_indEBV=match(markerGene,row.names(mvsts.EBV))"
"CD8_EBV=match("CD8A",row.names(mvsts.EBV))"
"CD4_EBV=match("CD4",row.names(mvsts.EBV))"
"CD8cells=rep(0,length(mvsts.EBV[Gene_indEBV[1],]))"
"CD4cells=rep(0,length(mvsts.EBV[Gene_indEBV[1],]))"

"vec1=mvsts.EBV[CD4_EBV,]"
"vec2=mvsts.EBV[CD8_EBV,]"

for(g in 1:length(CD8cells)){
  if((vec2[g]==0)&(vec1[g]>0)){
    CD4cells[g]=1
  }
  else if((vec2[g]>0)&(vec1[g]==0)){
    CD8cells[g]=1
  }
  else{}
}
rm(vec1)
rm(vec2)

"Clonal_TEBV=AddMetaData(Clonal_TEBV, CD8cells, col.name = 'CD8cells')"
"Clonal_TEBV=AddMetaData(Clonal_TEBV, CD4cells, col.name = 'CD4cells')"

rm(CD8cells)
rm(CD4cells)
#Check for errors
"plot(log10(mvsts.EBV[Gene_indEBV[3],]),log10(mvsts.EBV[Gene_indEBV[2],]),xlab=paste(Gene_indEBV[3]," (log10 UMIs)",sep=""),ylab=paste(Gene_indEBV[2]," (log10 UMIs)",sep=""),col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))"
"plot(log10(mvsts.EBV[CD8_EBV,]),log10(mvsts.EBV[CD4_EBV,]),xlab="CD8",ylab="CD4",col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))"
"print("Starting scTRansform sample 2...")"
"Clonal_TEBV <- SCTransform(Clonal_TEBV, vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n = 7000)"
"Gene_indEBV=match(markerGene,row.names(Clonal_TEBV))"
"d1=as.matrix(Clonal_TEBV[["SCT"]]@data)[Gene_indEBV,]"
alist <- 1:length(Gene_indEBV)

"Thresholds=matrix("unassigned",nrow=dim(Clonal_TEBV)[2],ncol=length(Gene_indEBV))"
for(k in 1:length(Gene_indEBV)){
"  vec1=Clonal_TEBV[Gene_indEBV[k],][["SCT"]]@data"
  for(j in 1:dim(Clonal_TEBV)[2]){
    if(vec1[j]>cutoff[k]){
"      Thresholds[j,k]="EBV_high""
    }
  }
"  Clonal_TEBV=AddMetaData(Clonal_TEBV, Thresholds[,k], col.name = paste("Threshold_",markerGene[k],sep=""))"
}

rm(mvsts.EBV)
rm(RnaStoreEBV)
rm(clonal.data)

###################################################
####CMV ###########################################
###################################################

clonal.data <- Read10X(data.dir = dir_rna3)
Total_TCellsCMV=length(clonal.data@p)
"Clonal_TCMV <- CreateSeuratObject(counts = clonal.data,project = "CMV")"
"Clonal_TCMV[["percent.mt"]] <- PercentageFeatureSet(Clonal_TCMV, pattern = "^MT-")"
"Clonal_TCMV<- subset(Clonal_TCMV, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 5)"
RnaStoreCMV=Clonal_TCMV@assays$RNA@counts

#VDJ list with barcodes
CMVTCRarray <- read.csv(dir_vdj3)
"tcrCMV_cell=subset(CMVTCRarray,productive=="true" & is_cell=="true")"
TCRlistCMV=tcrCMV_cell$barcode
"joint.bcsCMV <- intersect(colnames(RnaStoreCMV), TCRlistCMV)"
#Remove barcodes from list with abberant GE
"tcrhashCMV=subset(CMVTCRarray,CMVTCRarray$barcode %in% joint.bcsCMV)"
"mvsts.CMV <- RnaStoreCMV[, joint.bcsCMV]"
"Clonal_TCMV=Clonal_TCMV[,joint.bcsCMV]"

###remove genes that are not of interest
for(c in 1:length(mask)){
  Allgenes=row.names(Clonal_TCMV)
"  Gset=Allgenes[grep(mask[c],Allgenes)]"
"  NewGeneSet=setdiff(Allgenes,Gset)"
"  IndexesGTR=match(NewGeneSet,row.names(Clonal_TCMV))"
"  Clonal_TCMV <- Clonal_TCMV[IndexesGTR,]"
}

####Genes of interest
"Gene_indCMV=match(markerGene,row.names(mvsts.CMV))"
"CD8_CMV=match("CD8A",row.names(mvsts.CMV))"
"CD4_CMV=match("CD4",row.names(mvsts.CMV))"
"CD8cells=rep(0,length(mvsts.CMV[Gene_indCMV[1],]))"
"CD4cells=rep(0,length(mvsts.CMV[Gene_indCMV[1],]))"

"vec1=mvsts.CMV[CD4_CMV,]"
"vec2=mvsts.CMV[CD8_CMV,]"

for(g in 1:length(CD8cells)){
  if((vec2[g]==0)&(vec1[g]>0)){
    CD4cells[g]=1
  }
  else if((vec2[g]>0)&(vec1[g]==0)){
    CD8cells[g]=1
  }
  else{}
}
rm(vec1)
rm(vec2)

Clonal_TCMV=AddMetaData(Clonal_TCMV, CD8cells, col.name = 'CD8cells')
Clonal_TCMV=AddMetaData(Clonal_TCMV, CD4cells, col.name = 'CD4cells')

rm(CD8cells)
rm(CD4cells)
#Check for errors
plot(log10(mvsts.CMV[Gene_indCMV[3],]),log10(mvsts.CMV[Gene_indCMV[2],]),xlab=paste(Gene_indCMV[3],"(log10 UMIs)",sep=""),ylab=paste(Gene_indCMV[2],"(log10 UMIs)",sep=""),col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))
plot(log10(mvsts.CMV[CD8_CMV,]),log10(mvsts.CMV[CD4_CMV,]),xlab="CD8",ylab="CD4",col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))
print("Starting scTRansform sample 3...")
Clonal_TCMV <- SCTransform(Clonal_TCMV, vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n = 7000)
Gene_indCMV=match(markerGene,row.names(Clonal_TCMV))
d1=as.matrix(Clonal_TCMV[["SCT"]]@data)[Gene_indCMV,]
alist <- 1:length(Gene_indCMV)

Thresholds=matrix("unassigned",nrow=dim(Clonal_TCMV)[2],ncol=length(Gene_indCMV))
for(k in 1:length(Gene_indCMV)){
  vec1=Clonal_TCMV[Gene_indCMV[k],][["SCT"]]@data
  for(j in 1:dim(Clonal_TCMV)[2]){
    if(vec1[j]>cutoff[k]){
      Thresholds[j,k]="CMV_high"
    }
  }
  Clonal_TCMV=AddMetaData(Clonal_TCMV, Thresholds[,k], col.name = paste("Threshold_",markerGene[k],sep=""))
}

rm(mvsts.CMV)
rm(RnaStoreCMV)
rm(clonal.data)
###################################################
####BkV ###########################################
###################################################

clonal.data <- Read10X(data.dir = dir_rna4)
Total_TCellsBKV=length(clonal.data@p)
"Clonal_TBKV <- CreateSeuratObject(counts = clonal.data,project = "BKV")"
"Clonal_TBKV[["percent.mt"]] <- PercentageFeatureSet(Clonal_TBKV, pattern = "^MT-")"
"Clonal_TBKV<- subset(Clonal_TBKV, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & percent.mt < 5)"
RnaStoreBKV=Clonal_TBKV@assays$RNA@counts

#VDJ list with barcodes
BKVTCRarray <- read.csv(dir_vdj4)
"tcrBKV_cell=subset(BKVTCRarray,productive=="true" & is_cell=="true")"
TCRlistBKV=tcrBKV_cell$barcode
"joint.bcsBKV <- intersect(colnames(RnaStoreBKV), TCRlistBKV)"
#Remove barcodes from list with abberant GE
"tcrhashBKV=subset(BKVTCRarray,BKVTCRarray$barcode %in% joint.bcsBKV)"
"mvsts.BKV <- RnaStoreBKV[, joint.bcsBKV]"
"Clonal_TBKV=Clonal_TBKV[,joint.bcsBKV]"

###remove genes that are not of interest
for(c in 1:length(mask)){
  Allgenes=row.names(Clonal_TBKV)
"  Gset=Allgenes[grep(mask[c],Allgenes)]"
"  NewGeneSet=setdiff(Allgenes,Gset)"
"  IndexesGTR=match(NewGeneSet,row.names(Clonal_TBKV))"
"  Clonal_TBKV <- Clonal_TBKV[IndexesGTR,]"
}

####Genes of interest
"Gene_indBKV=match(markerGene,row.names(mvsts.BKV))"
"CD8_BKV=match("CD8A",row.names(mvsts.BKV))"
"CD4_BKV=match("CD4",row.names(mvsts.BKV))"
"CD8cells=rep(0,length(mvsts.BKV[Gene_indBKV[1],]))"
"CD4cells=rep(0,length(mvsts.BKV[Gene_indBKV[1],]))"

"vec1=mvsts.BKV[CD4_BKV,]"
"vec2=mvsts.BKV[CD8_BKV,]"

for(g in 1:length(CD8cells)){
  if((vec2[g]==0)&(vec1[g]>0)){
    CD4cells[g]=1
  }
  else if((vec2[g]>0)&(vec1[g]==0)){
    CD8cells[g]=1
  }
  else{}
}
rm(vec1)
rm(vec2)

"Clonal_TBKV=AddMetaData(Clonal_TBKV, CD8cells, col.name = 'CD8cells')"
"Clonal_TBKV=AddMetaData(Clonal_TBKV, CD4cells, col.name = 'CD4cells')"

rm(CD8cells)
rm(CD4cells)
#Check for errors
"plot(log10(mvsts.BKV[Gene_indBKV[3],]),log10(mvsts.BKV[Gene_indBKV[2],]),xlab=paste(Gene_indBKV[3]," (log10 UMIs)",sep=""),ylab=paste(Gene_indBKV[2]," (log10 UMIs)",sep=""),col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))"
"plot(log10(mvsts.BKV[CD8_BKV,]),log10(mvsts.BKV[CD4_BKV,]),xlab="CD8",ylab="CD4",col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))"
"print("Starting scTRansform sample 4...")"
"Clonal_TBKV <- SCTransform(Clonal_TBKV, vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n = 7000)"
"Gene_indBKV=match(markerGene,row.names(Clonal_TBKV))"
"d1=as.matrix(Clonal_TBKV[["SCT"]]@data)[Gene_indBKV,]"
alist <- 1:length(Gene_indBKV)

"Thresholds=matrix("unassigned",nrow=dim(Clonal_TBKV)[2],ncol=length(Gene_indBKV))"
for(k in 1:length(Gene_indBKV)){
"  vec1=Clonal_TBKV[Gene_indBKV[k],][["SCT"]]@data"
  for(j in 1:dim(Clonal_TBKV)[2]){
    if(vec1[j]>cutoff[k]){
"      Thresholds[j,k]="BKV_high""
    }
  }
"  Clonal_TBKV=AddMetaData(Clonal_TBKV, Thresholds[,k], col.name = paste("Threshold_",markerGene[k],sep=""))"
}

rm(mvsts.BKV)
rm(RnaStoreBKV)
rm(clonal.data)

#####################################################
##         Generate clonotypes by TCRB (na)        ##
#####################################################

"monoB=rep(0,dim(Clonal_TUMO)[2])"
"monoTRBNA=rep("unassigned",dim(Clonal_TUMO)[2])"
"monoTRBAA=rep("unassigned",dim(Clonal_TUMO)[2])"
"monoclono=rep("unassigned",dim(Clonal_TUMO)[2])"

##TRB loop
for(j in 1:dim(Clonal_TUMO)[2]){
  
  monoB[j]=colnames(Clonal_TUMO)[j]
"  set_umo=subset(tcrhashUMO,barcode==colnames(Clonal_TUMO)[j] & productive=="true" & chain=="TRB")"
  TCRs=set_umo
  len=length(set_umo$reads)
  if(len==0){
    
  }
  else if(len==1){
    monoTRBNA[j]=set_umo$cdr3_nt
    monoTRBAA[j]=set_umo$cdr3
  }
  else{
    position= which.max(set_umo$reads)
    monoTRBNA[j]=set_umo$cdr3_nt[position]
    monoTRBAA[j]=set_umo$cdr3[position]
  }
  
}

"cmvB=rep(0,dim(Clonal_TUMO)[2])"
"cmvTRBNA=rep("unassigned",dim(Clonal_TCMV)[2])"
"cmvTRBAA=rep("unassigned",dim(Clonal_TCMV)[2])"
"cmvclono=rep("unassigned",dim(Clonal_TCMV)[2])"

##TRB loop
for(j in 1:dim(Clonal_TCMV)[2]){
  
  cmvB[j]=colnames(Clonal_TCMV)[j]
"  set_CMV=subset(tcrhashCMV,barcode==colnames(Clonal_TCMV)[j] & productive=="true" & chain=="TRB")"
  TCRs=set_CMV
  len=length(set_CMV$reads)
  if(len==0){
    
  }
  else if(len==1){
    cmvTRBNA[j]=set_CMV$cdr3_nt
    cmvTRBAA[j]=set_CMV$cdr3
  }
  else{
    position= which.max(set_CMV$reads)
    cmvTRBNA[j]=set_CMV$cdr3_nt[position]
    cmvTRBAA[j]=set_CMV$cdr3[position]
  }
  
}

"ebvB=rep(0,dim(Clonal_TUMO)[2])"
"ebvTRBNA=rep("unassigned",dim(Clonal_TEBV)[2])"
"ebvTRBAA=rep("unassigned",dim(Clonal_TEBV)[2])"
"ebvclono=rep("unassigned",dim(Clonal_TEBV)[2])"

##TRB loop
for(j in 1:dim(Clonal_TEBV)[2]){
  
  ebvB[j]=colnames(Clonal_TEBV)[j]
"  set_EBV=subset(tcrhashEBV,barcode==colnames(Clonal_TEBV)[j] & productive=="true" & chain=="TRB")"
  TCRs=set_EBV
  len=length(set_EBV$reads)
  if(len==0){
    
  }
  else if(len==1){
    ebvTRBNA[j]=set_EBV$cdr3_nt
    ebvTRBAA[j]=set_EBV$cdr3
  }
  else{
    position= which.max(set_EBV$reads)
    ebvTRBNA[j]=set_EBV$cdr3_nt[position]
    ebvTRBAA[j]=set_EBV$cdr3[position]
  }
  
}

"bkvB=rep(0,dim(Clonal_TUMO)[2])"
"bkvTRBNA=rep("unassigned",dim(Clonal_TBKV)[2])"
"bkvTRBAA=rep("unassigned",dim(Clonal_TBKV)[2])"
"bkvclono=rep("unassigned",dim(Clonal_TBKV)[2])"

##TRB loop
for(j in 1:dim(Clonal_TBKV)[2]){
  
  bkvB[j]=colnames(Clonal_TBKV)[j]
"  set_BKV=subset(tcrhashBKV,barcode==colnames(Clonal_TBKV)[j] & productive=="true" & chain=="TRB")"
  TCRs=set_BKV
  len=length(set_BKV$reads)
  if(len==0){
    
  }
  else if(len==1){
    bkvTRBNA[j]=set_BKV$cdr3_nt
    bkvTRBAA[j]=set_BKV$cdr3
  }
  else{
    position= which.max(set_BKV$reads)
    bkvTRBNA[j]=set_BKV$cdr3_nt[position]
    bkvTRBAA[j]=set_BKV$cdr3[position]
  }
  
}


"IntersectTCRs=intersect(monoTRBNA,intersect(ebvTRBNA,intersect(bkvTRBNA,cmvTRBNA)))"

"Tailmono4=setdiff(monoTRBNA,IntersectTCRs)"
"IMTCRs4=c(IntersectTCRs,Tailmono4)"

"Tailmono3=setdiff(ebvTRBNA,IMTCRs4)"
"IMTCRs3=c(IMTCRs4,Tailmono3)"

"Tailmono2=setdiff(bkvTRBNA,IMTCRs3)"
"IMTCRs2=c(IMTCRs3,Tailmono2)"

"Tailmono1=setdiff(cmvTRBNA,IMTCRs2)"
"IMTCRs1=c(IMTCRs2,Tailmono1)"

"SizesMono=rep(0,length(IMTCRs1))"
for(j in 1:length(IMTCRs1)){
"  SizesMono[j]=length(subset(tcrhashUMO,cdr3_nt==IMTCRs1[j] & productive=="true" & chain=="TRB")$barcode)"
}

"OSM=order(SizesMono,decreasing=TRUE)"
"Nclono=rep("blank",length(IMTCRs1))"
for(f in 1: length(IMTCRs1)){
  if(OSM[f]==7){
    
  }
  else{
"    Nclono[OSM[f]]=paste("clonotype",f,sep="")"
  }
}

"NA_clonotypes_mono=data.frame(TCRs=IMTCRs1,clonoT=Nclono,freqT=SizesMono)"

#original arrays
"#monoB=rep(0,dim(Clonal_TUMO)[2])"
"#monoTRBNA=rep("unassigned",dim(Clonal_TUMO)[2])"
"#monoTRBAA=rep("unassigned",dim(Clonal_TUMO)[2])"
"#monoclono=rep("unassigned",dim(Clonal_TUMO)[2])"

"umoCname=rep("unassigned",dim(Clonal_TUMO)[2])"
"umoSize=rep(0,dim(Clonal_TUMO)[2])"

##TRB loop
for(j in 1:dim(Clonal_TUMO)[2]){
  
"  if(monoTRBNA[j]=="unassigned"){"
    
  }
  else{
"    umoSize[j]=subset(NA_clonotypes_mono,TCRs==monoTRBNA[j])$freqT"
"    umoCname[j]=subset(NA_clonotypes_mono,TCRs==monoTRBNA[j])$clonoT"
  }
  
}



"SizesEBV=rep(0,length(IMTCRs1))"
for(j in 1:length(IMTCRs1)){
"  SizesEBV[j]=length(subset(tcrhashEBV,cdr3_nt==IMTCRs1[j] & productive=="true" & chain=="TRB")$barcode)"
}

"OSM=order(SizesEBV,decreasing=TRUE)"
"Nclono=rep("blank",length(IMTCRs1))"
for(f in 1: length(IMTCRs1)){
  if(OSM[f]==7){
    
  }
  else{
"    Nclono[OSM[f]]=paste("clonotype",f,sep="")"
  }
}

"NA_clonotypes_EBV=data.frame(TCRs=IMTCRs1,clonoT=Nclono,freqT=SizesEBV)"

#original arrays
"#EBVB=rep(0,dim(Clonal_TEBV)[2])"
"#EBVTRBNA=rep("unassigned",dim(Clonal_TEBV)[2])"
"#EBVTRBAA=rep("unassigned",dim(Clonal_TEBV)[2])"
"#EBVclono=rep("unassigned",dim(Clonal_TEBV)[2])"

"EBVCname=rep("unassigned",dim(Clonal_TEBV)[2])"
"EBVSize=rep(0,dim(Clonal_TEBV)[2])"

##TRB loop
for(j in 1:dim(Clonal_TEBV)[2]){
  
"  if(ebvTRBNA[j]=="unassigned"){"
    
  }
  else{
"    EBVSize[j]=subset(NA_clonotypes_EBV,TCRs==ebvTRBNA[j])$freqT"
"    EBVCname[j]=subset(NA_clonotypes_EBV,TCRs==ebvTRBNA[j])$clonoT"
  }
  
}



"SizesBKV=rep(0,length(IMTCRs1))"
for(j in 1:length(IMTCRs1)){
"  SizesBKV[j]=length(subset(tcrhashBKV,cdr3_nt==IMTCRs1[j] & productive=="true" & chain=="TRB")$barcode)"
}

"OSM=order(SizesBKV,decreasing=TRUE)"
"Nclono=rep("blank",length(IMTCRs1))"
for(f in 1: length(IMTCRs1)){
  if(OSM[f]==7){
    
  }
  else{
"    Nclono[OSM[f]]=paste("clonotype",f,sep="")"
  }
}

"NA_clonotypes_BKV=data.frame(TCRs=IMTCRs1,clonoT=Nclono,freqT=SizesBKV)"

#original arrays
"#BKVB=rep(0,dim(Clonal_TBKV)[2])"
"#BKVTRBNA=rep("unassigned",dim(Clonal_TBKV)[2])"
"#BKVTRBAA=rep("unassigned",dim(Clonal_TBKV)[2])"
"#BKVclono=rep("unassigned",dim(Clonal_TBKV)[2])"

"BKVCname=rep("unassigned",dim(Clonal_TBKV)[2])"
"BKVSize=rep(0,dim(Clonal_TBKV)[2])"

##TRB loop
for(j in 1:dim(Clonal_TBKV)[2]){
  
"  if(bkvTRBNA[j]=="unassigned"){"
    
  }
  else{
"    BKVSize[j]=subset(NA_clonotypes_BKV,TCRs==bkvTRBNA[j])$freqT"
"    BKVCname[j]=subset(NA_clonotypes_BKV,TCRs==bkvTRBNA[j])$clonoT"
  }
  
}

"SizesCMV=rep(0,length(IMTCRs1))"
for(j in 1:length(IMTCRs1)){
"  SizesCMV[j]=length(subset(tcrhashCMV,cdr3_nt==IMTCRs1[j] & productive=="true" & chain=="TRB")$barcode)"
}

"OSM=order(SizesCMV,decreasing=TRUE)"
"Nclono=rep("blank",length(IMTCRs1))"
for(f in 1: length(IMTCRs1)){
  if(OSM[f]==7){
    
  }
  else{
"    Nclono[OSM[f]]=paste("clonotype",f,sep="")"
  }
}

"NA_clonotypes_CMV=data.frame(TCRs=IMTCRs1,clonoT=Nclono,freqT=SizesCMV)"

#original arrays
"#CMVB=rep(0,dim(Clonal_TCMV)[2])"
"#CMVTRBNA=rep("unassigned",dim(Clonal_TCMV)[2])"
"#CMVTRBAA=rep("unassigned",dim(Clonal_TCMV)[2])"
"#CMVclono=rep("unassigned",dim(Clonal_TCMV)[2])"

"CMVCname=rep("unassigned",dim(Clonal_TCMV)[2])"
"CMVSize=rep(0,dim(Clonal_TCMV)[2])"

##TRB loop
for(j in 1:dim(Clonal_TCMV)[2]){
  
"  if(cmvTRBNA[j]=="unassigned"){"
    
  }
  else{
"    CMVSize[j]=subset(NA_clonotypes_CMV,TCRs==cmvTRBNA[j])$freqT"
"    CMVCname[j]=subset(NA_clonotypes_CMV,TCRs==cmvTRBNA[j])$clonoT"
  }
  
}

"ClonotypesAll=data.frame(Spec=c(rep(paste("Threshold_",markerGene[1],sep=""),dim(Clonal_TEBV)[2]),"
"                                rep(paste("Threshold_",markerGene[1],sep=""),dim(Clonal_TCMV)[2])"
"                                ,rep(paste("Threshold_",markerGene[1],sep=""),dim(Clonal_TBKV)[2]),"
"                                rep(paste("Threshold_",markerGene[1],sep=""),dim(Clonal_TUMO)[2])),"
"                         Name=c(rep("ebv",dim(Clonal_TEBV)[2]),rep("cmv",dim(Clonal_TCMV)[2])"
"                                ,rep("bkv",dim(Clonal_TBKV)[2]),rep("umo",dim(Clonal_TUMO)[2])),"
"                         Clono=c(EBVCname,CMVCname,BKVCname,"
"                                 umoCname),"
"                         CS=c(EBVSize,CMVSize,BKVSize,umoSize),"
"                         cd4=c(Clonal_TEBV@meta.data$CD4cells,Clonal_TCMV@meta.data$CD4cells,"
"                               Clonal_TBKV@meta.data$CD4cells,Clonal_TUMO@meta.data$CD4cells),"
"                         cd8=c(Clonal_TEBV@meta.data$CD8cells,Clonal_TCMV@meta.data$CD8cells,"
"                               Clonal_TBKV@meta.data$CD8cells,Clonal_TUMO@meta.data$CD8cells))"


"Clonal_TEBV=AddMetaData(Clonal_TEBV, EBVCname, col.name = 'clonotype')"
"Clonal_TEBV=AddMetaData(Clonal_TEBV, EBVSize, col.name = 'countcln')"
"Clonal_TUMO=AddMetaData(Clonal_TUMO, umoCname, col.name = 'clonotype')"
"Clonal_TUMO=AddMetaData(Clonal_TUMO, umoSize, col.name = 'countcln')"
"Clonal_TBKV=AddMetaData(Clonal_TBKV, BKVCname, col.name = 'clonotype')"
"Clonal_TBKV=AddMetaData(Clonal_TBKV, BKVSize, col.name = 'countcln')"
"Clonal_TCMV=AddMetaData(Clonal_TCMV, CMVCname, col.name = 'clonotype')"
"Clonal_TCMV=AddMetaData(Clonal_TCMV, CMVSize, col.name = 'countcln')"

#####################################################
# Now add clonotype # and size to metadata vectors ##
#####################################################

###############
## Monocytes ##
###############

"xmonoclonoNAseq=rep("unassigned",dim(Clonal_TUMO)[2])"
"xmonoclonoAAseq=rep("unassigned",dim(Clonal_TUMO)[2])"

##TRB loop
for(j in 1:dim(Clonal_TUMO)[2]){
  
  monoB[j]=colnames(Clonal_TUMO)[j]
"  set_umo=subset(tcrhashUMO,barcode==colnames(Clonal_TUMO)[j] & productive=="true" & chain=="TRB")"
  TCRs=set_umo
  len=length(set_umo$reads)
  if(len==0){
    
  }
  else if(len==1){
    mTRBNA=set_umo$cdr3_nt
    mTRBAA=set_umo$cdr3
    xmonoclonoNAseq[j]=mTRBNA
    xmonoclonoAAseq[j]=mTRBAA
  }
  else{
    position= which.max(set_umo$reads)
    mTRBNA=set_umo$cdr3_nt[position]
    mTRBAA=set_umo$cdr3[position]
    xmonoclonoNAseq[j]=mTRBNA
    xmonoclonoAAseq[j]=mTRBAA
  }
  
}

"Clonal_TUMO=AddMetaData(Clonal_TUMO, xmonoclonoNAseq, col.name = 'cdr3_na')"
"Clonal_TUMO=AddMetaData(Clonal_TUMO, xmonoclonoAAseq, col.name = 'cdr3')"


###############
## EBV ##
###############

"xEBVclonoNAseq=rep("unassigned",dim(Clonal_TEBV)[2])"
"xEBVclonoAAseq=rep("unassigned",dim(Clonal_TEBV)[2])"
##TRB loop
for(j in 1:dim(Clonal_TEBV)[2]){
  
"  set_ebv=subset(tcrhashEBV,barcode==colnames(Clonal_TEBV)[j] & productive=="true" & chain=="TRB")"
  len=length(set_ebv$reads)
  if(len==0){
  }
  else if(len==1){
    mTRBNA=set_ebv$cdr3_nt
    mTRBAA=set_ebv$cdr3
    xEBVclonoNAseq[j]=mTRBNA
    xEBVclonoAAseq[j]=mTRBAA
  }
  else{
    position= which.max(set_ebv$reads)
    mTRBNA=set_ebv$cdr3_nt[position]
    mTRBAA=set_ebv$cdr3[position]
    xEBVclonoNAseq[j]=mTRBNA
    xEBVclonoAAseq[j]=mTRBAA
  }
}
"Clonal_TEBV=AddMetaData(Clonal_TEBV, xEBVclonoNAseq, col.name = 'cdr3_na')"
"Clonal_TEBV=AddMetaData(Clonal_TEBV, xEBVclonoAAseq, col.name = 'cdr3')"



###############
## CMV ##
###############

"xCMVclonoNAseq=rep("unassigned",dim(Clonal_TCMV)[2])"
"xCMVclonoAAseq=rep("unassigned",dim(Clonal_TCMV)[2])"
##TRB loop
for(j in 1:dim(Clonal_TCMV)[2]){
  
"  set_CMV=subset(tcrhashCMV,barcode==colnames(Clonal_TCMV)[j] & productive=="true" & chain=="TRB")"
  len=length(set_CMV$reads)
  if(len==0){
  }
  else if(len==1){
    mTRBNA=set_CMV$cdr3_nt
    mTRBAA=set_CMV$cdr3
    xCMVclonoNAseq[j]=mTRBNA
    xCMVclonoAAseq[j]=mTRBAA
  }
  else{
    position= which.max(set_CMV$reads)
    mTRBNA=set_CMV$cdr3_nt[position]
    mTRBAA=set_CMV$cdr3[position]
    xCMVclonoNAseq[j]=mTRBNA
    xCMVclonoAAseq[j]=mTRBAA
  }
}
"Clonal_TCMV=AddMetaData(Clonal_TCMV, xCMVclonoNAseq, col.name = 'cdr3_na')"
"Clonal_TCMV=AddMetaData(Clonal_TCMV, xCMVclonoAAseq, col.name = 'cdr3')"


###############
## BKV ##
###############

"xBKVclonoNAseq=rep("unassigned",dim(Clonal_TBKV)[2])"
"xBKVclonoAAseq=rep("unassigned",dim(Clonal_TBKV)[2])"
##TRB loop
for(j in 1:dim(Clonal_TBKV)[2]){
  
"  set_BKV=subset(tcrhashBKV,barcode==colnames(Clonal_TBKV)[j] & productive=="true" & chain=="TRB")"
  len=length(set_BKV$reads)
  if(len==0){
  }
  else if(len==1){
    mTRBNA=set_BKV$cdr3_nt
    mTRBAA=set_BKV$cdr3
    xBKVclonoNAseq[j]=mTRBNA
    xBKVclonoAAseq[j]=mTRBAA
  }
  else{
    position= which.max(set_BKV$reads)
    mTRBNA=set_BKV$cdr3_nt[position]
    mTRBAA=set_BKV$cdr3[position]
    xBKVclonoNAseq[j]=mTRBNA
    xBKVclonoAAseq[j]=mTRBAA
  }
}
"Clonal_TBKV=AddMetaData(Clonal_TBKV, xBKVclonoNAseq, col.name = 'cdr3_na')"
"Clonal_TBKV=AddMetaData(Clonal_TBKV, xBKVclonoAAseq, col.name = 'cdr3')"

rm(AgspecBKV)
rm(AgspecCMV)
rm(AgspecEBV)
rm(NA_clonotypes_BKV)
rm(NA_clonotypes_EBV)
rm(NA_clonotypes_CMV)
rm(prodTCRBKV)
rm(prodTCRCMV)
rm(prodTCREBV)
rm(set_BKV)
rm(set_EBV)
rm(set_umo)
rm(tcrBKV_cell)
rm(tcrCMV_cell)
rm(tcrEBV_cell)
rm(tcrhashBKV)
rm(tcrhashCMV)
rm(tcrhashUMO)
rm(tcrhashEBV)
rm(UMOTCRarray)
rm(BKVTCRarray)
rm(ClonotypesAll)
rm(CMVTCRarray)
rm(EBVTCRarray)
rm(NA_clonotypes_CMV)
rm(set_CMV)
rm(set_ebv)
rm(tcrUMO_cell)
rm(NA_clonotypes_mono)
rm(TCRs)

#scalers
rm(bkvB)
rm(bkvclono)
rm(BKVCname)
rm(BKVSize)
rm(bkvTRBAA)
rm(bkvTRBNA)
rm(CD4_BKV)
rm(CD4_CMV)
rm(CD4_UMO)
rm(CD8_BKV)
rm(CD8_CMV)
rm(CD8_EBV)
rm(CD8_UMO)
rm(cmvB)
rm(cmvclono)
rm(CMVCname)
rm(CMVSize)
rm(cmvTRBAA)
rm(cmvTRBNA)

rm(CD4_EBV)
rm(ebvB)
rm(ebvclono)
rm(EBVCname)
rm(EBVSize)
rm(ebvTRBAA)
rm(ebvTRBNA)
rm(f)
rm(g)
rm(IMTCRs1)
rm(IMTCRs2)
rm(IMTCRs3)
rm(IMTCRs4)
rm(IntersectTCRs)
rm(j)
rm(joint.bcsBKV)
rm(joint.bcsCMV)
rm(joint.bcsEBV)
rm(joint.bcsUMO)
rm(Total_TCellsBKV)
rm(umoCname)
rm(umoSize)
rm(vec1)
rm(xBKVclonoAAseq)
rm(xCMVclonoAAseq)
rm(xEBVclonoAAseq)
rm(xBKVclonoNAseq)
rm(xCMVclonoNAseq)
rm(xEBVclonoNAseq)
rm(xmonoclonoAAseq)
rm(xmonoclonoNAseq)

rm(Total_TCellsBKV)
rm(Total_TCellsEBV)
rm(Total_TCellsCMV)

rm(TCRlistBKV)
rm(TCRlistCMV)
rm(TCRlistEBV)
rm(TCRlistUMO)
rm(TNF_indUMO)
rm(Total_TCellsUMO)

rm(Tailmono1)
rm(Tailmono2)
rm(Tailmono3)
rm(Tailmono4)

rm(SizesBKV)
rm(SizesCMV)
rm(SizesEBV)
rm(SizesMono)

rm(position)
rm(allcells)
rm(len)
rm(monoB)
rm(monoclono)
rm(monoTRBAA)
rm(monoTRBNA)
rm(mTRBAA)
rm(mTRBNA)
rm(Nclono)
rm(OSM)

Clonal_TBKV_t1=Clonal_TBKV
Clonal_TEBV_t1=Clonal_TEBV
Clonal_TCMV_t1=Clonal_TCMV
Clonal_TUMO_t1=Clonal_TUMO

"Dataset=list(Clonal_TCMV,Clonal_TEBV,Clonal_TBKV,Clonal_TUMO)"
"saveRDS(Dataset,save_dir)"
return(Dataset)
}