###Read in marker intensity scan table
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
ThreeHourStim=readRDS("../RDS/ThreeHourStim.rds")
SixHourStim=readRDS("../RDS/SixHourStim.rds")


GetSpecificCells <- function(Clonaldata,expression,phenotype,tcrseq,moipos){
  if(phenotype=="CD8"){
    pindx=5
  if(length(subset(subset(Clonaldata,(`cdr3_na` %in% tcrseq))@meta.data[[pindx]],subset(Clonaldata,(`cdr3_na` %in% tcrseq))@meta.data[[pindx]]=="1"))>0){
    if(length(subset(subset(Clonaldata,(`cdr3_na` %in% tcrseq) & Clonaldata@meta.data[[pindx]]=="1")@meta.data[[8+moipos]],subset(Clonaldata,(`cdr3_na` %in% tcrseq) & Clonaldata@meta.data[[pindx]]=="1")@meta.data[[8+moipos]]==expression))>0){
      return(subset(Clonaldata,Clonaldata@meta.data[[8+moipos]]==expression & Clonaldata@meta.data[[pindx]]=="1" & (`cdr3_na` %in% tcrseq)))
    }
    else{
      return(0)
    }
  }
  else{
    return(0)
  }
  }
  else if(phenotype=="CD4"){
    pindx=6
  if(length(subset(subset(Clonaldata,(`cdr3_na` %in% tcrseq))@meta.data[[pindx]],subset(Clonaldata,(`cdr3_na` %in% tcrseq))@meta.data[[pindx]]=="1"))>0){
    if(length(subset(subset(Clonaldata,(`cdr3_na` %in% tcrseq) & Clonaldata@meta.data[[pindx]]=="1")@meta.data[[8+moipos]],subset(Clonaldata,(`cdr3_na` %in% tcrseq) & Clonaldata@meta.data[[pindx]]=="1")@meta.data[[8+moipos]]==expression))>0){
      return(subset(Clonaldata,Clonaldata@meta.data[[8+moipos]]==expression & Clonaldata@meta.data[[pindx]]=="1" & (`cdr3_na` %in% tcrseq)))
    }
    else{
      return(0)
    }
  }
  else{
    return(0)
  }
  }
  else{
    if(length(subset(subset(Clonaldata,(`cdr3_na` %in% tcrseq))@meta.data[[8+moipos]],subset(Clonaldata,(`cdr3_na` %in% tcrseq))@meta.data[[8+moipos]]==expression))>0){
      return(subset(Clonaldata,Clonaldata@meta.data[[8+moipos]]==expression & (`cdr3_na` %in% tcrseq)))
    }
    else{
      return(0)
    }
  }
}

ClonotypeTable_path = "../SupplementaryTables/SummaryTNF_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,9,10,11,17)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                       "3hrCMV",
                       "N3hrCMV")

Clonal_T=ThreeHourStim[[3]]
Threshold="BKV_high"
moipos=4
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/BKV_TNF"

GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)

ClonotypeTable_path = "../SupplementaryTables/SummaryCD70_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,9,10,11,17)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")

Clonal_T=ThreeHourStim[[3]]
Threshold="BKV_high"
moipos=5
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/BKV_CD70"

GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


ClonotypeTable_path = "../SupplementaryTables/SummaryTNFRSF9_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,9,10,11,17)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")

Clonal_T=ThreeHourStim[[3]]
Threshold="BKV_high"
moipos=6
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/BKV_TNFRSF9"

GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)

ClonotypeTable_path = "../SupplementaryTables/SummaryIFNG_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,9,10,11,17)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")

Clonal_T=ThreeHourStim[[3]]
Threshold="BKV_high"
moipos=1
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/BKV_IFNG"

GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)

#######CMV cells
ClonotypeTable_path = "../SupplementaryTables/SummaryIFNG_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,3,4,5,15)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[1]]
Threshold="CMV_high"
moipos=1
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/CMV_IFNG"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


ClonotypeTable_path = "../SupplementaryTables/SummaryCRTAM_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,3,4,5,15)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[1]]
Threshold="CMV_high"
moipos=2
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/CMV_CRTAM"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


ClonotypeTable_path = "../SupplementaryTables/SummaryCD69_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,3,4,5,15)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[1]]
Threshold="CMV_high"
moipos=3
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/CMV_CD69"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


ClonotypeTable_path = "../SupplementaryTables/SummaryTNF_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,3,4,5,15)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[1]]
Threshold="CMV_high"
moipos=4
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/CMV_TNF"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)



ClonotypeTable_path = "../SupplementaryTables/SummaryCD70_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,3,4,5,15)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[1]]
Threshold="CMV_high"
moipos=5
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/CMV_CD70"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)



ClonotypeTable_path = "../SupplementaryTables/SummaryTNFRSF9_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,3,4,5,15)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[1]]
Threshold="CMV_high"
moipos=6
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/CMV_TNFRSF9"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


#######EBV cells
ClonotypeTable_path = "../SupplementaryTables/SummaryIFNG_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,6,7,8,16)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[2]]
Threshold="EBV_high"
moipos=1
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/EBV_IFNG"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


ClonotypeTable_path = "../SupplementaryTables/SummaryCRTAM_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,6,7,8,16)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[2]]
Threshold="EBV_high"
moipos=2
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/EBV_CRTAM"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


ClonotypeTable_path = "../SupplementaryTables/SummaryCD69_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,6,7,8,16)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[2]]
Threshold="EBV_high"
moipos=3
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/EBV_CD69"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


ClonotypeTable_path = "../SupplementaryTables/SummaryTNF_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,6,7,8,16)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[2]]
Threshold="EBV_high"
moipos=4
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/EBV_TNF"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


ClonotypeTable_path = "../SupplementaryTables/SummaryCD70_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,6,7,8,16)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[2]]
Threshold="EBV_high"
moipos=5
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/EBV_CD70"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


ClonotypeTable_path = "../SupplementaryTables/SummaryTNFRSF9_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,6,7,8,16)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")
Clonal_T=ThreeHourStim[[2]]
Threshold="EBV_high"
moipos=6
markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")
deg_path = "../SupplementaryTables/Genes"
gene="/EBV_TNFRSF9"
GetGeneSignature(Clonal_T,Summarydf,Threshold,moipos,gene)


GetGeneSignature <- function(Clonal_T,Summarydf,Threshold,moipos,gene,upper_lim,lower_lim,bystander_lim){
  
  CD4MaxClonotypeUpper=quantile(subset(Summarydf,(Summarydf[,2]/100)*Summarydf[,5]>10)[,2],upper_lim)[[1]]
  CD4MaxClonotype=quantile(subset(Summarydf,(Summarydf[,2]/100)*Summarydf[,5]>10)[,2],lower_lim)[[1]]
  CD4MinClonotype=quantile(subset(Summarydf,Summarydf[,2]>0 & Summarydf[,5]>10)[,2],bystander_lim)[[1]]
  
  CD4TCRActEBV=list(subset(Summarydf,(Summarydf[,2]>=CD4MaxClonotype) & (Summarydf[,2]<=CD4MaxClonotypeUpper))$cdr3_na)
  CD4TCRNonEBV=list(subset(Summarydf,Summarydf[,2]<=CD4MinClonotype & Summarydf[,5]>10)$cdr3_na)
  
  CD8MaxClonotypeUpper=quantile(subset(Summarydf,(Summarydf[,3]/100)*Summarydf[,5]>10)[,3],upper_lim)[[1]]
  CD8MaxClonotype=quantile(subset(Summarydf,(Summarydf[,3]/100)*Summarydf[,5]>10)[,3],lower_lim)[[1]]
  CD8MinClonotype=quantile(subset(Summarydf,Summarydf[,3]>0 & Summarydf[,5]>10)[,3],bystander_lim)[[1]]
  
  CD8TCRActEBV=list(subset(Summarydf,(Summarydf[,3]>=CD8MaxClonotype) & (Summarydf[,3]<=CD8MaxClonotypeUpper))$cdr3_na)
  CD8TCRNonEBV=list(subset(Summarydf,Summarydf[,3]<=CD8MinClonotype & Summarydf[,5]>10)$cdr3_na)
  
  MaxClonotypeUpper=quantile(subset(Summarydf,(Summarydf[,4]/100)*Summarydf[,5]>10)[,4],upper_lim)[[1]]
  MaxClonotype=quantile(subset(Summarydf,(Summarydf[,4]/100)*Summarydf[,5]>10)[,4],lower_lim)[[1]]
  MinClonotype=quantile(subset(Summarydf,Summarydf[,4]>0 & Summarydf[,5]>10)[,4],bystander_lim)[[1]]
  
  TCRActEBV=list(subset(Summarydf,(Summarydf[,4]>=MaxClonotype) & (Summarydf[,4]<=MaxClonotypeUpper))$cdr3_na)
  TCRNonEBV=list(subset(Summarydf,Summarydf[,4]<=MinClonotype & Summarydf[,5]>10)$cdr3_na)
  
  unassignedTCRs = setdiff(levels(factor(Clonal_T$cdr3_na)),c(CD4TCRActEBV[[1]],
                                                              CD4TCRNonEBV[[1]],
                                                              CD8TCRActEBV[[1]],
                                                              CD8TCRNonEBV[[1]],
                                                              TCRActEBV[[1]],
                                                              TCRNonEBV[[1]]))
  
  #------------------------------------#
  #-----Set all the specific cells-----#
  #------------------------------------#
  if(as.character(CD4TCRActEBV)=="character(0)"){
    Arrangemeta=data.frame(bc=colnames(Clonal_T),Idx=1:length(colnames(Clonal_T)),Viralgroup=rep(1,length(colnames(Clonal_T))))
    CD8_HighClones=GetSpecificCells(Clonal_T,Threshold,"CD8",CD8TCRActEBV[[1]],moipos)
    CD8_LowClones=GetSpecificCells(Clonal_T,"unassigned","CD8",CD8TCRActEBV[[1]],moipos)
    CD8_HighBystanderClones=GetSpecificCells(Clonal_T,Threshold,"CD8",CD8TCRNonEBV[[1]],moipos)
    CD8_LowBystanderClones=GetSpecificCells(Clonal_T,"unassigned","CD8",CD8TCRNonEBV[[1]],moipos)
    CD8_HighViralClones=GetSpecificCells(Clonal_T,Threshold,"CD8",unassignedTCRs,moipos)
    CD8_LowViralClones=GetSpecificCells(Clonal_T,"unassigned","CD8",unassignedTCRs,moipos)
    
    HighClones=GetSpecificCells(Clonal_T,Threshold,"all",TCRActEBV[[1]],moipos)
    LowClones=GetSpecificCells(Clonal_T,"unassigned","all",TCRActEBV[[1]],moipos)
    HighBystanderClones=GetSpecificCells(Clonal_T,Threshold,"all",TCRNonEBV[[1]],moipos)
    LowBystanderClones=GetSpecificCells(Clonal_T,"unassigned","all",TCRNonEBV[[1]],moipos)
    HighViralClones=GetSpecificCells(Clonal_T,Threshold,"all",unassignedTCRs,moipos)
    LowViralClones=GetSpecificCells(Clonal_T,"unassigned","all",unassignedTCRs,moipos)
    
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_HighClones))$Idx] = "CD8Specific_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_LowClones))$Idx] = "CD8Specific_Low"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% setdiff(colnames(CD8_HighBystanderClones),colnames(CD8_HighClones)))$Idx] = "CD8Bystander_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_LowBystanderClones))$Idx] = "CD8Bystander_Low"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_HighViralClones))$Idx] = "CD8Unspecified_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_LowViralClones))$Idx] = "CD8Unspecified_Low"
    
    Clonal_T=AddMetaData(Clonal_T, Arrangemeta$Viralgroup, col.name = 'SpecificityGroup')
    Idents(Clonal_T = Clonal_T) <- Clonal_T@meta.data$SpecificityGroup

    ##CD8 T cells
    DEGmarkers <- FindMarkers(Clonal_T, ident.1 = c("CD8Specific_High"), ident.2 = c("CD8Bystander_Low"),min.pct = 0)
    SigDEG=data.frame(gene=row.names(DEGmarkers),FC=DEGmarkers$avg_log2FC,pval=DEGmarkers$p_val_adj)
    SigDEG=subset(SigDEG, abs(FC)>0 & pval<0.01)
    Snew_df2Up=subset(SigDEG,FC>0)
    SpecificGLUP=Snew_df2Up$gene[order(Snew_df2Up$FC,decreasing=TRUE)]
    Snew_df2Down=subset(SigDEG,FC<0)
    SpecificGLDown=Snew_df2Down$gene[order(Snew_df2Down$FC,decreasing=TRUE)]
    allgenes=c(SpecificGLUP,SpecificGLDown)
    write.csv(SigDEG,paste(deg_path,paste(gene,"/CD8DEGs.csv",sep=""),sep=""))
    write.csv(allgenes,paste(deg_path,paste(gene,"/CD8genes.csv",sep=""),sep=""))
    write.csv(SpecificGLUP,paste(deg_path,paste(gene,"/CD8genes_upregulated.csv",sep=""),sep=""))
    write.csv(SpecificGLDown,paste(deg_path,paste(gene,"/CD8genes_downregulated.csv",sep=""),sep=""))
    ##All T cells
    
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(HighClones))$Idx] = "Specific_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(LowClones))$Idx] = "Specific_Low"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(HighBystanderClones))$Idx] = "Bystander_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(LowBystanderClones))$Idx] = "Bystander_Low"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(HighViralClones))$Idx] = "Unspecified_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(LowViralClones))$Idx] = "Unspecified_Low"
    
    Clonal_T=AddMetaData(Clonal_T, Arrangemeta$Viralgroup, col.name = 'SpecificityGroup')
    Idents(Clonal_T = Clonal_T) <- Clonal_T@meta.data$SpecificityGroup
    DEGmarkers <- FindMarkers(Clonal_T, ident.1 = c("Specific_High"), ident.2 = c("Bystander_Low"),min.pct = 0)
    SigDEG=data.frame(gene=row.names(DEGmarkers),FC=DEGmarkers$avg_log2FC,pval=DEGmarkers$p_val_adj)
    SigDEG=subset(SigDEG, abs(FC)>0 & pval<0.01)
    Snew_df2Up=subset(SigDEG,FC>0)
    SpecificGLUP=Snew_df2Up$gene[order(Snew_df2Up$FC,decreasing=TRUE)]
    Snew_df2Down=subset(SigDEG,FC<0)
    SpecificGLDown=Snew_df2Down$gene[order(Snew_df2Down$FC,decreasing=TRUE)]
    allgenes=c(SpecificGLUP,SpecificGLDown)
    write.csv(SigDEG,paste(deg_path,paste(gene,"/DEGs.csv",sep=""),sep=""))
    write.csv(allgenes,paste(deg_path,paste(gene,"/genes.csv",sep=""),sep=""))
    write.csv(SpecificGLUP,paste(deg_path,paste(gene,"/genes_upregulated.csv",sep=""),sep=""))
    write.csv(SpecificGLDown,paste(deg_path,paste(gene,"/genes_downregulated.csv",sep=""),sep=""))
  }
  else if(as.character(CD8TCRActEBV)=="character(0)"){
    CD4_HighClones=GetSpecificCells(Clonal_T,Threshold,"CD4",CD4TCRActEBV[[1]],moipos)
    CD4_LowClones=GetSpecificCells(Clonal_T,"unassigned","CD4",CD4TCRActEBV[[1]],moipos)
    CD4_HighBystanderClones=GetSpecificCells(Clonal_T,Threshold,"CD4",CD4TCRNonEBV[[1]],moipos)
    CD4_LowBystanderClones=GetSpecificCells(Clonal_T,"unassigned","CD4",CD4TCRNonEBV[[1]],moipos)
    CD4_HighViralClones=GetSpecificCells(Clonal_T,Threshold,"CD4",unassignedTCRs,moipos)
    CD4_LowViralClones=GetSpecificCells(Clonal_T,"unassigned","CD4",unassignedTCRs,moipos)
    
    HighClones=GetSpecificCells(Clonal_T,Threshold,"all",TCRActEBV[[1]],moipos)
    LowClones=GetSpecificCells(Clonal_T,"unassigned","all",TCRActEBV[[1]],moipos)
    HighBystanderClones=GetSpecificCells(Clonal_T,Threshold,"all",TCRNonEBV[[1]],moipos)
    LowBystanderClones=GetSpecificCells(Clonal_T,"unassigned","all",TCRNonEBV[[1]],moipos)
    HighViralClones=GetSpecificCells(Clonal_T,Threshold,"all",unassignedTCRs,moipos)
    LowViralClones=GetSpecificCells(Clonal_T,"unassigned","all",unassignedTCRs,moipos)
    
    Arrangemeta=data.frame(bc=colnames(Clonal_T),Idx=1:length(colnames(Clonal_T)),Viralgroup=rep(1,length(colnames(Clonal_T))))
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_HighClones))$Idx] = "CD4Specific_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_LowClones))$Idx] = "CD4Specific_Low"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% setdiff(colnames(CD4_HighBystanderClones),colnames(CD4_HighClones)))$Idx] = "CD4Bystander_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_LowBystanderClones))$Idx] = "CD4Bystander_Low"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_HighViralClones))$Idx] = "CD4Unspecified_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_LowViralClones))$Idx] = "CD4Unspecified_Low"
    
    Clonal_T=AddMetaData(Clonal_T, Arrangemeta$Viralgroup, col.name = 'SpecificityGroup')
    Idents(Clonal_T = Clonal_T) <- Clonal_T@meta.data$SpecificityGroup
    
    ##CD4 T cells
    DEGmarkers <- FindMarkers(Clonal_T, ident.1 = c("CD4Specific_High"), ident.2 = c("CD4Bystander_Low"),min.pct = 0)
    SigDEG=data.frame(gene=row.names(DEGmarkers),FC=DEGmarkers$avg_log2FC,pval=DEGmarkers$p_val_adj)
    SigDEG=subset(SigDEG, abs(FC)>0 & pval<0.01)
    Snew_df2Up=subset(SigDEG,FC>0)
    SpecificGLUP=Snew_df2Up$gene[order(Snew_df2Up$FC,decreasing=TRUE)]
    Snew_df2Down=subset(SigDEG,FC<0)
    SpecificGLDown=Snew_df2Down$gene[order(Snew_df2Down$FC,decreasing=TRUE)]
    allgenes=c(SpecificGLUP,SpecificGLDown)
    write.csv(SigDEG,paste(deg_path,paste(gene,"/CD4DEGs.csv",sep=""),sep=""))
    write.csv(allgenes,paste(deg_path,paste(gene,"/CD4genes.csv",sep=""),sep=""))
    write.csv(SpecificGLUP,paste(deg_path,paste(gene,"/CD4genes_upregulated.csv",sep=""),sep=""))
    write.csv(SpecificGLDown,paste(deg_path,paste(gene,"/CD4genes_downregulated.csv",sep=""),sep=""))
    ##All T cells
    
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(HighClones))$Idx] = "Specific_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(LowClones))$Idx] = "Specific_Low"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(HighBystanderClones))$Idx] = "Bystander_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(LowBystanderClones))$Idx] = "Bystander_Low"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(HighViralClones))$Idx] = "Unspecified_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(LowViralClones))$Idx] = "Unspecified_Low"
    
    Clonal_T=AddMetaData(Clonal_T, Arrangemeta$Viralgroup, col.name = 'SpecificityGroup')
    Idents(Clonal_T = Clonal_T) <- Clonal_T@meta.data$SpecificityGroup
    DEGmarkers <- FindMarkers(Clonal_T, ident.1 = c("Specific_High"), ident.2 = c("Bystander_Low"),min.pct = 0)
    SigDEG=data.frame(gene=row.names(DEGmarkers),FC=DEGmarkers$avg_log2FC,pval=DEGmarkers$p_val_adj)
    SigDEG=subset(SigDEG, abs(FC)>0 & pval<0.01)
    Snew_df2Up=subset(SigDEG,FC>0)
    SpecificGLUP=Snew_df2Up$gene[order(Snew_df2Up$FC,decreasing=TRUE)]
    Snew_df2Down=subset(SigDEG,FC<0)
    SpecificGLDown=Snew_df2Down$gene[order(Snew_df2Down$FC,decreasing=TRUE)]
    allgenes=c(SpecificGLUP,SpecificGLDown)
    write.csv(SigDEG,paste(deg_path,paste(gene,"/DEGs.csv",sep=""),sep=""))
    write.csv(allgenes,paste(deg_path,paste(gene,"/genes.csv",sep=""),sep=""))
    write.csv(SpecificGLUP,paste(deg_path,paste(gene,"/genes_upregulated.csv",sep=""),sep=""))
    write.csv(SpecificGLDown,paste(deg_path,paste(gene,"/genes_downregulated.csv",sep=""),sep=""))
  }
  else{
  CD4_HighClones=GetSpecificCells(Clonal_T,Threshold,"CD4",CD4TCRActEBV[[1]],moipos)
  CD4_LowClones=GetSpecificCells(Clonal_T,"unassigned","CD4",CD4TCRActEBV[[1]],moipos)
  CD4_HighBystanderClones=GetSpecificCells(Clonal_T,Threshold,"CD4",CD4TCRNonEBV[[1]],moipos)
  CD4_LowBystanderClones=GetSpecificCells(Clonal_T,"unassigned","CD4",CD4TCRNonEBV[[1]],moipos)
  CD4_HighViralClones=GetSpecificCells(Clonal_T,Threshold,"CD4",unassignedTCRs,moipos)
  CD4_LowViralClones=GetSpecificCells(Clonal_T,"unassigned","CD4",unassignedTCRs,moipos)
  
  CD8_HighClones=GetSpecificCells(Clonal_T,Threshold,"CD8",CD8TCRActEBV[[1]],moipos)
  CD8_LowClones=GetSpecificCells(Clonal_T,"unassigned","CD8",CD8TCRActEBV[[1]],moipos)
  CD8_HighBystanderClones=GetSpecificCells(Clonal_T,Threshold,"CD8",CD8TCRNonEBV[[1]],moipos)
  CD8_LowBystanderClones=GetSpecificCells(Clonal_T,"unassigned","CD8",CD8TCRNonEBV[[1]],moipos)
  CD8_HighViralClones=GetSpecificCells(Clonal_T,Threshold,"CD8",unassignedTCRs,moipos)
  CD8_LowViralClones=GetSpecificCells(Clonal_T,"unassigned","CD8",unassignedTCRs,moipos)
  
  HighClones=GetSpecificCells(Clonal_T,Threshold,"all",TCRActEBV[[1]],moipos)
  LowClones=GetSpecificCells(Clonal_T,"unassigned","all",TCRActEBV[[1]],moipos)
  HighBystanderClones=GetSpecificCells(Clonal_T,Threshold,"all",TCRNonEBV[[1]],moipos)
  LowBystanderClones=GetSpecificCells(Clonal_T,"unassigned","all",TCRNonEBV[[1]],moipos)
  HighViralClones=GetSpecificCells(Clonal_T,Threshold,"all",unassignedTCRs,moipos)
  LowViralClones=GetSpecificCells(Clonal_T,"unassigned","all",unassignedTCRs,moipos)
  
  ###--- Make a dataframe with barcodes and an empty vector that will
  ###--- be attached to the seurat matrix
  Arrangemeta=data.frame(bc=colnames(Clonal_T),Idx=1:length(colnames(Clonal_T)),Viralgroup=rep(1,length(colnames(Clonal_T))))
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_HighClones))$Idx] = "CD4Specific_High"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_LowClones))$Idx] = "CD4Specific_Low"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% setdiff(colnames(CD4_HighBystanderClones),colnames(CD4_HighClones)))$Idx] = "CD4Bystander_High"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_LowBystanderClones))$Idx] = "CD4Bystander_Low"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_HighViralClones))$Idx] = "CD4Unspecified_High"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_LowViralClones))$Idx] = "CD4Unspecified_Low"
  
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_HighClones))$Idx] = "CD8Specific_High"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_LowClones))$Idx] = "CD8Specific_Low"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% setdiff(colnames(CD8_HighBystanderClones),colnames(CD8_HighClones)))$Idx] = "CD8Bystander_High"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_LowBystanderClones))$Idx] = "CD8Bystander_Low"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_HighViralClones))$Idx] = "CD8Unspecified_High"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_LowViralClones))$Idx] = "CD8Unspecified_Low"

  Clonal_T=AddMetaData(Clonal_T, Arrangemeta$Viralgroup, col.name = 'SpecificityGroup')
  Idents(Clonal_T = Clonal_T) <- Clonal_T@meta.data$SpecificityGroup
  ##CD4 T cells
  DEGmarkers <- FindMarkers(Clonal_T, ident.1 = c("CD4Specific_High"), ident.2 = c("CD4Bystander_Low"),min.pct = 0)
  SigDEG=data.frame(gene=row.names(DEGmarkers),FC=DEGmarkers$avg_log2FC,pval=DEGmarkers$p_val_adj)
  SigDEG=subset(SigDEG, abs(FC)>0 & pval<0.01)
  Snew_df2Up=subset(SigDEG,FC>0)
  SpecificGLUP=Snew_df2Up$gene[order(Snew_df2Up$FC,decreasing=TRUE)]
  Snew_df2Down=subset(SigDEG,FC<0)
  SpecificGLDown=Snew_df2Down$gene[order(Snew_df2Down$FC,decreasing=TRUE)]
  allgenes=c(SpecificGLUP,SpecificGLDown)
  write.csv(SigDEG,paste(deg_path,paste(gene,"/CD4DEGs.csv",sep=""),sep=""))
  write.csv(allgenes,paste(deg_path,paste(gene,"/CD4genes.csv",sep=""),sep=""))
  write.csv(SpecificGLUP,paste(deg_path,paste(gene,"/CD4genes_upregulated.csv",sep=""),sep=""))
  write.csv(SpecificGLDown,paste(deg_path,paste(gene,"/CD4genes_downregulated.csv",sep=""),sep=""))
  ##CD8 T cells
  DEGmarkers <- FindMarkers(Clonal_T, ident.1 = c("CD8Specific_High"), ident.2 = c("CD8Bystander_Low"),min.pct = 0)
  SigDEG=data.frame(gene=row.names(DEGmarkers),FC=DEGmarkers$avg_log2FC,pval=DEGmarkers$p_val_adj)
  SigDEG=subset(SigDEG, abs(FC)>0 & pval<0.01)
  Snew_df2Up=subset(SigDEG,FC>0)
  SpecificGLUP=Snew_df2Up$gene[order(Snew_df2Up$FC,decreasing=TRUE)]
  Snew_df2Down=subset(SigDEG,FC<0)
  SpecificGLDown=Snew_df2Down$gene[order(Snew_df2Down$FC,decreasing=TRUE)]
  allgenes=c(SpecificGLUP,SpecificGLDown)
  write.csv(SigDEG,paste(deg_path,paste(gene,"/CD8DEGs.csv",sep=""),sep=""))
  write.csv(allgenes,paste(deg_path,paste(gene,"/CD8genes.csv",sep=""),sep=""))
  write.csv(SpecificGLUP,paste(deg_path,paste(gene,"/CD8genes_upregulated.csv",sep=""),sep=""))
  write.csv(SpecificGLDown,paste(deg_path,paste(gene,"/CD8genes_downregulated.csv",sep=""),sep=""))
  ##All T cells
  
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(HighClones))$Idx] = "Specific_High"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(LowClones))$Idx] = "Specific_Low"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(HighBystanderClones))$Idx] = "Bystander_High"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(LowBystanderClones))$Idx] = "Bystander_Low"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(HighViralClones))$Idx] = "Unspecified_High"
  Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(LowViralClones))$Idx] = "Unspecified_Low"
  
  Clonal_T=AddMetaData(Clonal_T, Arrangemeta$Viralgroup, col.name = 'SpecificityGroup')
  Idents(Clonal_T = Clonal_T) <- Clonal_T@meta.data$SpecificityGroup
  DEGmarkers <- FindMarkers(Clonal_T, ident.1 = c("Specific_High"), ident.2 = c("Bystander_Low"),min.pct = 0)
  SigDEG=data.frame(gene=row.names(DEGmarkers),FC=DEGmarkers$avg_log2FC,pval=DEGmarkers$p_val_adj)
  SigDEG=subset(SigDEG, abs(FC)>0 & pval<0.01)
  Snew_df2Up=subset(SigDEG,FC>0)
  SpecificGLUP=Snew_df2Up$gene[order(Snew_df2Up$FC,decreasing=TRUE)]
  Snew_df2Down=subset(SigDEG,FC<0)
  SpecificGLDown=Snew_df2Down$gene[order(Snew_df2Down$FC,decreasing=TRUE)]
  allgenes=c(SpecificGLUP,SpecificGLDown)
  write.csv(SigDEG,paste(deg_path,paste(gene,"/DEGs.csv",sep=""),sep=""))
  write.csv(allgenes,paste(deg_path,paste(gene,"/genes.csv",sep=""),sep=""))
  write.csv(SpecificGLUP,paste(deg_path,paste(gene,"/genes_upregulated.csv",sep=""),sep=""))
  write.csv(SpecificGLDown,paste(deg_path,paste(gene,"/genes_downregulated.csv",sep=""),sep=""))
  }
}



