#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Subsets a SEURAT object based on a list of meta data for the CD4 and CD8 subsets
#'
#' This function will read in Seurat objects processed by aocseq::AnnotateCellTypes
#' or any other software based on cell annotation and subset the array accordingly
#'
#' @param cell.data A Seurat object containing single cell RNA sequencing data.
#' @param expression Choice of "high" or "unassigned".
#' @param phenotype CD4 or CD8 classification.
#' @param cell.types List of cell types to subset single cells included in the analysis.
#' @param goi_v Gene of interest.
#' @return A subset of a Seurat object.
#' @concept Single cell analysis
#' @export
#'
GetSpecificCells <- function(
    cell.data,
    expression,
    phenotype,
    cell.types,
    goi_v
){
  if(phenotype=="CD8"){
    pindx=match("CD8cells",names(subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data))
    entry=match(paste("Threshold_",goi_v,sep=""),names(subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data))
    if(length(subset(subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data[[pindx]],subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data[[pindx]]=="1"))>0){
      if(length(subset(subset(cell.data,(`cdr3_na` %in% cell.types) & cell.data@meta.data[[pindx]]=="1")@meta.data[[entry]],subset(cell.data,(`cdr3_na` %in% cell.types) & cell.data@meta.data[[pindx]]=="1")@meta.data[[entry]]==expression))>0){
        return(subset(cell.data,(cell.data@meta.data[[entry]]==(expression)) & cell.data@meta.data[[pindx]]=="1" & (`cdr3_na` %in% cell.types)))
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
    pindx=match("CD4cells",names(subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data))
    entry=match(paste("Threshold_",goi_v,sep=""),names(subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data))
    if(length(subset(subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data[[pindx]],subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data[[pindx]]=="1"))>0){
      if(length(subset(subset(cell.data,(`cdr3_na` %in% cell.types) & cell.data@meta.data[[pindx]]=="1")@meta.data[[entry]],subset(cell.data,(`cdr3_na` %in% cell.types) & cell.data@meta.data[[pindx]]=="1")@meta.data[[entry]]==expression))>0){
        return(subset(cell.data,cell.data@meta.data[[entry]]==expression & cell.data@meta.data[[pindx]]=="1" & (`cdr3_na` %in% cell.types)))
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
    entry=match(paste("Threshold_",goi_v,sep=""),names(subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data))
    if(length(subset(subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data[[entry]],subset(cell.data,(`cdr3_na` %in% cell.types))@meta.data[[entry]]==expression))>0){

      return(subset(cell.data,cell.data@meta.data[[entry]]==expression & (`cdr3_na` %in% cell.types)))
    }
    else{
      return(0)
    }
  }
}


#' This function saves differential gene expression from multiple subsets from samples that have been processed and classified by aocseq.
#'
#' This function will take a set of annotation spreadsheets and export differential
#' expression between different cell types labelled by their marker genes
#' Currently the phenotypic separation is CD4 and CD8 T cell phenotypes, but the user is encouraged
#' to generate as many subsets as needed for downstream analysis.
#'
#' @param cell.data A Seurat object pre-processed with aocseq::CombineData.
#' @param annotation.path Character array. Directory of an aocseq clonotype annotation table.
#' @param save.dir Directory for storing differential genes.
#' @param goi Character array. Gene name, a marker of interest.
#' @param FClim Fold change cutoff for differential genes.
#' @param specific.pval If cells are classified by a zero inflated negative binomial (ZINB) model, sets the significance cutoff for specificity in each cell subset specified by the cell type metadata.
#' @param bystander.pval Same parameter as used for specific.pval, but sets the threshold for bystanders.
#' @param set.min.pct Same set.min.pct in Seurat FindMarkers function.
#' @param logfc.threshold Same logfc.threshold in Seurat FindMarkers function.
#'
#' @return return an assay containing predicted expression value in the data
#' slot
#' @concept Single cell analysis
#' @export
GetGeneSignature <- function(
    cell.data,
    annotation.path,
    save.dir,
    goi,
    FClim=0,
    specific.pval=10^(-5),
    bystander.pval=10^(-2),
    set.min.pct=0.25,
    logfc.threshold = 0){

  Classification=data.frame(read.csv(annotation.path))
  #-----------------------------------------#
  #-----Find specific TCRs based on goi-----#
  #-----------------------------------------#
  if(is.na(match(paste(levels(factor(cell.data@meta.data$orig.ident)),".abundance",sep=""),names(Classification)))==TRUE)
    stop("No matching column data in clonotype annotation table... To fix: Rerun GetGeneSignature using a Clonal Object with a matching orig.ident in the aocseq annotation table");gc();

  entry_idx=match(paste(levels(factor(cell.data@meta.data$orig.ident)),".abundance",sep=""),names(Classification))
  pheno_idx=match("phenotype",names(Classification))
  e_i=paste(levels(factor(cell.data@meta.data$orig.ident)),".abundance",sep="")
  e_i_s=paste(levels(factor(cell.data@meta.data$orig.ident)),".status",sep="")
  e_i_s_t=match(e_i_s,names(Classification))
  s=2
  idxs=s
  while(s<(length(grep(".status",names(Classification)))[[1]])){
    s=s+2
    idxs=append(idxs,s)

  }
  #idxs=idxs[-1]
  grep(".status",names(Classification))[idxs]

  indexofs=grep(e_i_s,names(Classification)[grep(".status",names(Classification))[idxs]])

  namesused=names(Classification)[grep(".status",names(Classification))[idxs]]

  abundanceused=names(Classification)[grep(".abundance",names(Classification))[idxs]]

  indexofa=grep(e_i,names(Classification)[grep(".abundance",names(Classification))[idxs]])

  CD4_Act=list(subset(
    subset(subset(Classification,
                  Classification[,e_i]!="-"),
          as.numeric(subset(Classification,
                 Classification[,e_i]!="-")[,e_i])<specific.pval),
    subset(subset(Classification,Classification[,e_i]!="-"),
           as.numeric(subset(Classification,Classification[,e_i]!="-")[,e_i])<specific.pval)[,"phenotype"]=="CD4")$cdr3_na)

  for(g in 1:length(namesused)){
    if(g!=indexofs){
  Act_restrict=subset(subset(Classification,(Classification[,namesused[g]]!="-" & Classification[,e_i_s]!="-")),
                          (as.numeric(subset(Classification,(Classification[,namesused[g]]!="-" & Classification[,e_i_s]!="-"))[,e_i_s])>as.numeric(subset(Classification,(Classification[,namesused[g]]!="-" & Classification[,e_i_s]!="-"))[,namesused[g]]))
  )$cdr3_na
    }
  }

  for(g in 1:length(abundanceused)){
    if(g!=indexofa){
      Act_restrict=append(Act_restrict,subset(subset(Classification,(Classification[,abundanceused[g]]!="-" & Classification[,e_i]!="-")),
                          (subset(Classification,(Classification[,abundanceused[g]]!="-" & Classification[,e_i]!="-"))[,e_i]>=subset(Classification,(Classification[,abundanceused[g]]!="-" & Classification[,e_i]!="-"))[,abundanceused[g]])
      )$cdr3_na)
    }
  }
  if(length(as.character(Act_restrict))==0){
  }else{
  CD4_Act=list(setdiff(CD4_Act[[1]],Act_restrict))
}
  CD4_Bys_1=list(subset(subset(subset(Classification,Classification[,e_i]!="-"),
                               as.numeric(subset(Classification,Classification[,e_i]!="-")[,e_i])>bystander.pval),
                        subset(subset(Classification,Classification[,e_i]!="-"),
                               as.numeric(subset(Classification,Classification[,e_i]!="-")[,e_i])>bystander.pval)[,"phenotype"]!=("CD8"))$cdr3_na)
  CD4_Bys_2=list(subset(subset(Classification,Classification[,e_i]=="-"),
              subset(Classification,Classification[,e_i]=="-")[,"phenotype"]!=("CD8"))$cdr3_na)
  CD4_Bys=c(CD4_Bys_1[[1]],CD4_Bys_2[[1]])

  CD8_Act=list(subset(
    subset(subset(Classification,
                  Classification[,e_i]!="-"),
           as.numeric(subset(Classification,
                  Classification[,e_i]!="-")[,e_i])<=specific.pval),
    subset(subset(Classification,Classification[,e_i]!="-"),
           as.numeric(subset(Classification,Classification[,e_i]!="-")[,e_i])<=specific.pval)[,"phenotype"]=="CD8")$cdr3_na)

  CD8_Bys_1=list(subset(subset(subset(Classification,Classification[,e_i]!="-"),as.numeric(subset(Classification,Classification[,e_i]!="-")[,e_i])>bystander.pval),subset(subset(Classification,Classification[,e_i]!="-"),as.numeric(subset(Classification,Classification[,e_i]!="-")[,e_i])>bystander.pval)[,"phenotype"]!=("CD4"))$cdr3_na)
  CD8_Bys_2=list(subset(subset(Classification,Classification[,e_i]=="-"),subset(Classification,Classification[,e_i]=="-")[,"phenotype"]!=("CD4"))$cdr3_na)
  #CD8_Bys_3=list(subset(subset(Classification,subset(Classification[,e_i],Classification[,e_i]=="-")),subset(Classification,subset(Classification[,e_i],Classification[,e_i]=="-"))[,"phenotype"]=="-")$cdr3_na)
  CD8_Bys=list(union(CD8_Bys_1[[1]],CD8_Bys_2[[1]]))#CD8_Bys_1#list(union(union(CD8_Bys_1[[1]],CD8_Bys_2[[1]]),CD8_Bys_3[[1]]))

    #list(subset(subset(Classification,is.na(Classification[,e_i])==TRUE),is.na(subset(Classification,is.na(Classification[,e_i])==TRUE)[,"phenotype"]))$cdr3_na)

    #list(union(CD8_Bys_1[[1]],CD8_Bys_2[[1]]))

  unassignedTCRs = setdiff(levels(factor(cell.data@meta.data$cdr3_na)),c(CD4_Act[[1]],
          CD4_Bys[[1]],CD8_Act[[1]],CD8_Bys[[1]]))

  #----------------------------#
  #-Set all the specific cells-#
  #----------------------------#
  if(as.character(CD8_Act)[[1]]=="character(0)"){
  }else{
    
    CD8_HighCells=GetSpecificCells(cell.data,"high","CD8",CD8_Act[[1]],goi)
    print(CD8_Act[[1]])
    CD8_LowBystanderCells=GetSpecificCells(cell.data,"unassigned","CD8",CD8_Bys[[1]],goi)
    CD8_HighViralCells=GetSpecificCells(cell.data,"high","CD8",unassignedTCRs,goi)
    CD8_LowViralCells=GetSpecificCells(subset(cell.data,CD8cells=="1"),"unassigned","CD8",setdiff(unassignedTCRs,"unassigned"),goi)
  }

  if(as.character(CD4_Act)[[1]]=="character(0)"){
    }else{
    CD4_HighCells=GetSpecificCells(cell.data,"high","CD4",CD4_Act[[1]],goi)
    #print(CD4_Bys)
    CD4_LowBystanderCells=GetSpecificCells(cell.data,"unassigned","CD4",CD4_Bys,goi)
    CD4_HighViralCells=GetSpecificCells(cell.data,"high","CD4",unassignedTCRs,goi)
    CD4_LowViralCells=GetSpecificCells(cell.data,"unassigned","CD4",unassignedTCRs,goi)
    }
  if(as.character(CD8_Act)[[1]]=="character(0)"){
  }else{
    Arrangemeta=data.frame(bc=colnames(cell.data),Idx=1:length(colnames(cell.data)),Viralgroup=rep(1,length(colnames(cell.data))))
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_HighCells))$Idx] = "CD8Specific_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_LowBystanderCells))$Idx] = "CD8Bystander_Low"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_HighViralCells))$Idx] = "CD8Unspecified_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_LowViralCells))$Idx] = "CD8Unspecified_Low"

    cell.data=AddMetaData(cell.data, Arrangemeta$Viralgroup, col.name = 'SpecificityGroup')
    Idents(cell.data = cell.data) <- cell.data@meta.data$SpecificityGroup


    ##CD8 T cells
    DEGmarkers <- FindMarkers(cell.data, ident.1 = c("CD8Specific_High"), ident.2 = c("CD8Bystander_Low"),min.pct=set.min.pct,test.use ="bimod",logfc.threshold = logfc.threshold)
    SigDEG=data.frame(gene=row.names(DEGmarkers),FC=DEGmarkers$avg_log2FC,pval=DEGmarkers$p_val_adj)
    SigDEG=subset(SigDEG, abs(FC)>FClim & pval<0.05)
    Snew_df2Up=subset(SigDEG,FC>FClim)
    SpecificGLUP=Snew_df2Up$gene[order(Snew_df2Up$FC,decreasing=TRUE)]
    Snew_df2Down=subset(SigDEG,FC<0)
    SpecificGLDown=Snew_df2Down$gene[order(Snew_df2Down$FC,decreasing=TRUE)]
    allgenes=c(SpecificGLUP,SpecificGLDown)
    write.csv(SigDEG,paste(save.dir,paste(goi,"_CD8DEGs.csv",sep=""),sep=""))
    write.csv(allgenes,paste(save.dir,paste(goi,"_CD8genes.csv",sep=""),sep=""))
    write.csv(SpecificGLUP,paste(save.dir,paste(goi,"_CD8genes_upregulated.csv",sep=""),sep=""))
    write.csv(SpecificGLDown,paste(save.dir,paste(goi,"_CD8genes_downregulated.csv",sep=""),sep=""))
}
    if(as.character(CD4_Act)[[1]]=="character(0)"){
    }else{
    Arrangemeta=data.frame(bc=colnames(cell.data),Idx=1:length(colnames(cell.data)),Viralgroup=rep(1,length(colnames(cell.data))))
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_HighCells))$Idx] = "CD4Specific_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_LowBystanderCells))$Idx] = "CD4Bystander_Low"
   # Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_HighViralCells))$Idx] = "CD4Unspecified_High"
   # Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_LowViralCells))$Idx] = "CD4Unspecified_Low"

    cell.data=AddMetaData(cell.data, Arrangemeta$Viralgroup, col.name = 'SpecificityGroup')
    Idents(cell.data = cell.data) <- cell.data@meta.data$SpecificityGroup


    ##CD4 T cells
    DEGmarkers <- FindMarkers(cell.data, ident.1 = c("CD4Specific_High"), ident.2 = c("CD4Bystander_Low"),min.pct=set.min.pct,test.use ="bimod",logfc.threshold = logfc.threshold)
    SigDEG=data.frame(gene=row.names(DEGmarkers),FC=DEGmarkers$avg_log2FC,pval=DEGmarkers$p_val_adj)
    SigDEG=subset(SigDEG, abs(FC)>FClim & pval<0.05)
    Snew_df2Up=subset(SigDEG,FC>FClim)
    SpecificGLUP=Snew_df2Up$gene[order(Snew_df2Up$FC,decreasing=TRUE)]
    Snew_df2Down=subset(SigDEG,FC<0)
    SpecificGLDown=Snew_df2Down$gene[order(Snew_df2Down$FC,decreasing=TRUE)]
    allgenes=c(SpecificGLUP,SpecificGLDown)
    write.csv(SigDEG,paste(save.dir,paste(goi,"_CD4DEGs.csv",sep=""),sep=""))
    write.csv(allgenes,paste(save.dir,paste(goi,"_CD4genes.csv",sep=""),sep=""))
    write.csv(SpecificGLUP,paste(save.dir,paste(goi,"_CD4genes_upregulated.csv",sep=""),sep=""))
    write.csv(SpecificGLDown,paste(save.dir,paste(goi,"_CD4genes_downregulated.csv",sep=""),sep=""))
}
}






