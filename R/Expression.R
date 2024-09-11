#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Subsets a SEURAT object based on a list of meta data for the CD4 and CD8 subsets
#'
#' This function will read in Seurat objects processed by traceseq::AnnotateClonotypes
#' or any other software based on cell annotation and subset the SEURAT array accordingly
#'
#' @param Clonal_Obs A Seurat object.
#' @param expression Choice of "high" or "unassigned".
#' @param phenotype CD4 or CD8 clasification.
#' @param tcrseq List of TCRs to mark clonotypes included in the analysis.
#' @param goi_v Gene of interest.
#' @return A subset of a Seurat object processed with aocseq.
#' @concept gene expression
#' @export
#'
GetSpecificCells <- function(
    Clonal_Obs,
    expression,
    phenotype,
    tcrseq,
    goi_v
){
  if(phenotype=="CD8"){
    pindx=match("CD8cells",names(subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data))
    entry=match(paste("Threshold_",goi_v,sep=""),names(subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data))
    if(length(subset(subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data[[pindx]],subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data[[pindx]]=="1"))>0){
      if(length(subset(subset(Clonal_Obs,(`cdr3_na` %in% tcrseq) & Clonal_Obs@meta.data[[pindx]]=="1")@meta.data[[entry]],subset(Clonal_Obs,(`cdr3_na` %in% tcrseq) & Clonal_Obs@meta.data[[pindx]]=="1")@meta.data[[entry]]==expression))>0){
        return(subset(Clonal_Obs,(Clonal_Obs@meta.data[[entry]]==(expression)) & Clonal_Obs@meta.data[[pindx]]=="1" & (`cdr3_na` %in% tcrseq)))
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
    pindx=match("CD4cells",names(subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data))
    entry=match(paste("Threshold_",goi_v,sep=""),names(subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data))
    if(length(subset(subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data[[pindx]],subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data[[pindx]]=="1"))>0){
      if(length(subset(subset(Clonal_Obs,(`cdr3_na` %in% tcrseq) & Clonal_Obs@meta.data[[pindx]]=="1")@meta.data[[entry]],subset(Clonal_Obs,(`cdr3_na` %in% tcrseq) & Clonal_Obs@meta.data[[pindx]]=="1")@meta.data[[entry]]==expression))>0){
        return(subset(Clonal_Obs,Clonal_Obs@meta.data[[entry]]==expression & Clonal_Obs@meta.data[[pindx]]=="1" & (`cdr3_na` %in% tcrseq)))
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
    entry=match(paste("Threshold_",goi_v,sep=""),names(subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data))
    if(length(subset(subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data[[entry]],subset(Clonal_Obs,(`cdr3_na` %in% tcrseq))@meta.data[[entry]]==expression))>0){

      return(subset(Clonal_Obs,Clonal_Obs@meta.data[[entry]]==expression & (`cdr3_na` %in% tcrseq)))
    }
    else{
      return(0)
    }
  }
}


#' Export differential gene expression from several different samples and phenotypes
#'
#' This function will take a set of annotation spreadsheets and export differential
#' expression between different cell types labelled by their marker genes
#' Currently the phenotypic separation is CD4 and CD8, but the reader is encouraged
#' to generate as many subsets as needed for downstream analysis
#'
#' @param Clonal_Obs A Seurat object pre-processed with aocseq::CombineData.
#' @param clonotype.path Character array. Directory of an aocseq clonotype annotation table.
#' @param save.dir Directory for storing differentially expressed genes.
#' @param goi Character array. Gene name, a marker of interest.
#' @param verbose Print progress
#'
#' @return return an assay containing predicted expression value in the data
#' slot
#' @concept integration
#' @export
GetGeneSignature <- function(
    Clonal_Obs,
    clonotype.path,
    save.dir,
    goi,
    FClim=10^(-3),
    specific.pval=10^(-5),
    bystander.pval=10^(-2),
    logfc.threshold = 0.01){

  Classification=data.frame(read.csv(clonotype.path))
  #-----------------------------------------#
  #-----Find specific TCRs based on goi-----#
  #-----------------------------------------#
  if(is.na(match(paste(levels(factor(Clonal_Obs@meta.data$orig.ident)),".abundance",sep=""),names(Classification)))==TRUE)
    stop("No matching column data in clonotype annotation table... To fix: Rerun GetGeneSignature using a Clonal Object with a matching orig.ident in the aocseq annotation table");gc();

  entry_idx=match(paste(levels(factor(Clonal_Obs@meta.data$orig.ident)),".abundance",sep=""),names(Classification))
  pheno_idx=match("phenotype",names(Classification))
  e_i=paste(levels(factor(Clonal_Obs@meta.data$orig.ident)),".abundance",sep="")
  e_i_s=paste(levels(factor(Clonal_Obs@meta.data$orig.ident)),".status",sep="")
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

  unassignedTCRs = setdiff(levels(factor(Clonal_Obs@meta.data$cdr3_na)),c(CD4_Act[[1]],
          CD4_Bys[[1]],CD8_Act[[1]],CD8_Bys[[1]]))

  #----------------------------#
  #-Set all the specific cells-#
  #----------------------------#
  if(as.character(CD8_Act)[[1]]=="character(0)"){
  }else{
    
    CD8_HighCells=GetSpecificCells(Clonal_Obs,"high","CD8",CD8_Act[[1]],goi)
    print(CD8_Act[[1]])
    CD8_LowBystanderCells=GetSpecificCells(Clonal_Obs,"unassigned","CD8",CD8_Bys[[1]],goi)
    CD8_HighViralCells=GetSpecificCells(Clonal_Obs,"high","CD8",unassignedTCRs,goi)
    CD8_LowViralCells=GetSpecificCells(subset(Clonal_Obs,CD8cells=="1"),"unassigned","CD8",setdiff(unassignedTCRs,"unassigned"),goi)
  }

  if(as.character(CD4_Act)[[1]]=="character(0)"){
    }else{
    CD4_HighCells=GetSpecificCells(Clonal_Obs,"high","CD4",CD4_Act[[1]],goi)
    #print(CD4_Bys)
    CD4_LowBystanderCells=GetSpecificCells(Clonal_Obs,"unassigned","CD4",CD4_Bys,goi)
    CD4_HighViralCells=GetSpecificCells(Clonal_Obs,"high","CD4",unassignedTCRs,goi)
    CD4_LowViralCells=GetSpecificCells(Clonal_Obs,"unassigned","CD4",unassignedTCRs,goi)
    }
  if(as.character(CD8_Act)[[1]]=="character(0)"){
  }else{
    Arrangemeta=data.frame(bc=colnames(Clonal_Obs),Idx=1:length(colnames(Clonal_Obs)),Viralgroup=rep(1,length(colnames(Clonal_Obs))))
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_HighCells))$Idx] = "CD8Specific_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_LowBystanderCells))$Idx] = "CD8Bystander_Low"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_HighViralCells))$Idx] = "CD8Unspecified_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD8_LowViralCells))$Idx] = "CD8Unspecified_Low"

    Clonal_Obs=AddMetaData(Clonal_Obs, Arrangemeta$Viralgroup, col.name = 'SpecificityGroup')
    Idents(Clonal_Obs = Clonal_Obs) <- Clonal_Obs@meta.data$SpecificityGroup


    ##CD8 T cells
    DEGmarkers <- FindMarkers(Clonal_Obs, ident.1 = c("CD8Specific_High"), ident.2 = c("CD8Bystander_Low"),min.pct = 0.0025,test.use ="bimod",min.cells.group=50,min.cells.feature=50,logfc.threshold = logfc.threshold)
    SigDEG=data.frame(gene=row.names(DEGmarkers),FC=DEGmarkers$avg_log2FC,pval=DEGmarkers$p_val_adj)
    SigDEG=subset(SigDEG, abs(FC)>FClim & pval<1)
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
    Arrangemeta=data.frame(bc=colnames(Clonal_Obs),Idx=1:length(colnames(Clonal_Obs)),Viralgroup=rep(1,length(colnames(Clonal_Obs))))
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_HighCells))$Idx] = "CD4Specific_High"
    Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_LowBystanderCells))$Idx] = "CD4Bystander_Low"
   # Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_HighViralCells))$Idx] = "CD4Unspecified_High"
   # Arrangemeta$Viralgroup[subset(Arrangemeta,bc %in% colnames(CD4_LowViralCells))$Idx] = "CD4Unspecified_Low"

    Clonal_Obs=AddMetaData(Clonal_Obs, Arrangemeta$Viralgroup, col.name = 'SpecificityGroup')
    Idents(Clonal_Obs = Clonal_Obs) <- Clonal_Obs@meta.data$SpecificityGroup


    ##CD4 T cells
    DEGmarkers <- FindMarkers(Clonal_Obs, ident.1 = c("CD4Specific_High"), ident.2 = c("CD4Bystander_Low"),min.pct = 0,test.use ="bimod",logfc.threshold = logfc.threshold)
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






