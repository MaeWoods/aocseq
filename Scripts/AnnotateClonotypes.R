#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Combine single cell data and annotate with cell labels based on functionality
#'
#' This function will read in Seurat objects processed by traceseq and generate a 
#' spreadsheet of clonotypes grouped by CDR3 beta sequence that displays the 
#' percentage of CD4 and CD8 cells with transcript expression that is greater than
#' the quantile cutoff of the control set in CombineData. A gene list used to estimate 
#' specificity can be chosen and is application specific. Other parameters are listed 
#' for debugging, but can be left as default values.
#'
#' @param gex.path path to gene expression data in the format of cellranger feature barcode matrices
#' @param vdj.path path to VDJ expression data in the format of cellranger csv files 
#' @param marker.gene list of marker genes used for specificity analysis
#' @param mask listing of genes to be dropped from analysis
#' @param save.dir directory for storing Seurat object RDS files
#' @param index.control index of the control sequencing sample
#' @param demultiplex Boolean value used to indicate if the data was hashtagged
#' @param hashtags Integer represents the number of hashtags used (if any)
#' @param verbose Print progress bars and output
#'
#' @return A Seurat object list containing metadata and VDJ annotations.
#' @concept annotation
#' @export
#'


ThreeHourStim=TetAssay
SixHourStim=readRDS("../RDS/SixHourStim.rds")

markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")

moi="TNFRSF9"
clonotype_path="../SupplementaryTables/SummaryTNFRSF9_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=6
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)


moi="CD70"
clonotype_path="../SupplementaryTables/SummaryCD70_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=5
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)

moi="TNF"
clonotype_path="../SupplementaryTables/SummaryTNF_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=4
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)


moi="CD69"
clonotype_path="../SupplementaryTables/SummaryCD69_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=3
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)


moi="CRTAM"
clonotype_path="../SupplementaryTables/SummaryCRTAM_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=2
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)


moi="IFNG"
clonotype_path="../SupplementaryTables/Tetramer/SummaryIFNG_2_5%.csv"
cell_path="../SupplementaryTables/Tetramer/Cells.csv"
moipos=1

index.control=1
threshold.cutoff=.975
n.batch=2
names.spreadsheet = c("Clone","3hr CD4 high CMV (%)","3hr CD8 high CMV (%)","3hr cells high CMV (%)"
                      ,"3hr CD4 high EBV (%)","3hr CD8 high EBV (%)","3hr cells high EBV (%)","3hr CD4 high BKV (%)",
                      "3hr CD8 high BKV (%)","3hr cells high BKV (%)", "3hr CD4 high US (%)",
                      "3hr CD8 high US (%)","3hr cells high US (%)","3hr Total cells CMV","3hr Total cells EBV",
                      "3hr Total cells BKV","3hr Total cells US",
                      "6hr CD4 high CMV (%)","6hr CD8 high CMV (%)","6hr cells high CMV (%)",
                      "6hr CD4 high EBV (%)","6hr CD8 high EBV (%)","6hr cells high EBV (%)","6hr CD4 high BKV (%)",
                      "6hr CD8 high BKV (%)","6hr cells high BKV (%)","6hr CD4 high US (%)",
                      "6hr CD8 high US (%)","6hr cells high US (%)","6hr Total cells CMV","6hr Total cells EBV",
                      "6hr Total cells BKV","6hr Total cells US")

AnnotateClonotypes(ThreeHourStim,ThreeHourStim,moi,clonotype_path,cell_path,moipos)


AnnotateClonotypes <- function(Clonal_Obs,moi,clonotype_path,moipos,n.batch,threshold.cutoff,index.control,conditions,names.spreadsheet){
  
  #match markers of interest
  
  unstim_t1 <- Clonal_Obs[[q*n.batch + k]]
  moi_samp[[q*n.batch + k]]=match(moi, row.names(Clonal_Obs[[q*n.batch + k]][["SCT"]]@data))
  cutoff[[q]]=quantile(Clonal_Obs[[index.control[[q]]]][["SCT"]]@data[moi_samp[[q*n.batch + k]],],threshold.cutoff)[[1]]
  
  #probably end loop here
  allclonotypes = union(Clonal_Obs[[1]]@meta.data$cdr3_na,Clonal_Obs[[2]]@meta.data$cdr3_na)
  intersectclonotypes=intersect(Clonal_Obs[[1]]@meta.data$cdr3_na,Clonal_Obs[[2]]@meta.data$cdr3_na)
  for(q in 1:conditions){
    for(k in 3:n.batch){
      allclonotypes = union(allclonotypes,Clonal_Obs[[q*n.batch + k]]@meta.data$cdr3_na)
      intersectclonotypes=intersect(intersectclonotypes,Clonal_Obs[[q*n.batch + k]]@meta.data$cdr3_na)
    }
  }
  
  allclonotypes = setdiff(levels(factor(allclonotypes)),"unassigned")
  intersectclonotypes = setdiff(levels(factor(intersectclonotypes)),"unassigned")
  
  sample.clonotypes=vector(mode = "list", length = conditions*n.batch)
  
  sample.clonotypes[[index.control[[1]]]]=setdiff(levels(factor(Clonal_Obs[[index.control[[q]]]]@meta.data$cdr3_na)),"unassigned")
  combine.clonotypes = sample.clonotypes[[index.control[[1]]]]
  for(q in 1:conditions){
    for(k in 3:n.batch){
      if(index.control[[1]]!=q*n.batch + k){
      f=q*n.batch + k
      sample.clonotypes[[index.control[[f]]]]=setdiff(setdiff(levels(factor(Clonal_Obs[[index.control[[f]]]]@meta.data$cdr3_na)),"unassigned"),combine.clonotypes)
      combine.clonotypes=union(combine.clonotypes,sample.clonotypes[[index.control[[f]]]])
      }
    }
  }
  
  #assign counts for ranking by UMO, then EBV, BKV - CMV
  print("Assigning counts for ranking...")
  NClonotypes=length(allclonotypes)
  Shared=rep(0,NClonotypes)
  Clonefreq=rep(0,NClonotypes)
  for(k in 1:NClonotypes){
    for(q in 1:conditions){
      for(d in 3:n.batch){
    if(length(intersect(allclonotypes[k],sample.clonotypes[[q*n.batch + d]]))>0){
      Clonefreq[k]=subset(Clonal_Obs[[q*n.batch + d]],cdr3_na==allclonotypes[k])@meta.data$countcln[1]
      if(length(intersect(allclonotypes[k],intersectclonotypes))>0){
        Shared[k]="yes"
      }
      else{
        Shared[k]="no"
      }
    }
        
      }
      }
  }
  
  Clonotypes_df = data.frame(cdr3=allclonotypes,shared=Shared,frequency=Clonefreq)
  Clonotypes_data = data.frame(cdr3=Clonotypes_df$cdr3[order(Clonotypes_df$frequency,decreasing=TRUE)],
                               shared=Clonotypes_df$shared[order(Clonotypes_df$frequency,decreasing=TRUE)],
                               frequency=Clonotypes_df$frequency[order(Clonotypes_df$frequency,decreasing=TRUE)])
  
  Clone=rep(0,length(allclonotypes))
  CD4_c=vector(mode = "list", length = conditions*n.batch)
  CD8_c=vector(mode = "list", length = conditions*n.batch)
  All_c=vector(mode = "list", length = conditions*n.batch)
  Total_c=vector(mode = "list", length = conditions*n.batch)
  TCR_seq=vector(mode = "list", length = conditions*n.batch)
  CD4df=vector(mode = "list", length = conditions*n.batch)
  CD8df=vector(mode = "list", length = conditions*n.batch)
  Alldf=vector(mode = "list", length = conditions*n.batch)
  for(q in 1:conditions){
    for(d in 3:n.batch){
      CD4_c[[q*n.batch + d]]=rep(0,length(allclonotypes))
      CD8_c[[q*n.batch + d]]=rep(0,length(allclonotypes))
      All_c[[q*n.batch + d]]=rep(0,length(allclonotypes))
      Total_c[[q*n.batch + d]]=rep(0,length(allclonotypes))
      TCR_seq[[q*n.batch + d]]=Clonal_Obs[[q*n.batch + d]]@meta.data$cdr3_na
      CD4df[[q*n.batch + d]]=data.frame(clonotype=subset(Clonal_Obs[[q*n.batch + d]],CD4cells=="1")@meta.data$cdr3_na,mark=subset(Clonal_Obs[[q*n.batch + d]],CD4cells=="1")@meta.data[[8+moipos]])
      CD8df[[q*n.batch + d]]=data.frame(clonotype=subset(Clonal_Obs[[q*n.batch + d]],CD8cells=="1")@meta.data$cdr3_na,mark=subset(Clonal_Obs[[q*n.batch + d]],CD8cells=="1")@meta.data[[8+moipos]])
      Alldf[[q*n.batch + d]]=data.frame(clonotype=Clonal_Obs[[q*n.batch + d]]@meta.data$cdr3_na,mark=Clonal_Obs[[q*n.batch + d]]@meta.data[[8+moipos]])
    }
  }
  
  print("Assigning proportional specificity...")
  for(k in 1:NClonotypes){
    Clone[k]=Clonotypes_data$cdr3[k]
    Ncellstotunstim=0
    Ncellsaboveunstim=0
    Ncellstotcmv_sct=0
    Ncellsabovecmv_sct=0
    Ncellstotebv_sct=0
    Ncellsaboveebv_sct=0
    Ncellstotbkv_sct=0
    Ncellsabovebkv_sct=0
    tempvec=subset(CD4df[[q*n.batch + d]],clonotype==Clonotypes_data$cdr3[k])$mark
    tempvectot=subset(Alldf[[q*n.batch + d]],clonotype==Clonotypes_data$cdr3[k])$mark
    if(length(intersect(Clonotypes_data$cdr3[k],TCR_seq[[q*n.batch + d]]))>0){
      Ncellstotunstim=length(tempvectot)
      Ncellsaboveunstim=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD4_c[[q*n.batch + d]][k]=0
      All_c[[q*n.batch + d]][k]=0
      if(Ncellstotunstim>0){
        
        CD4_c[[q*n.batch + d]][k]=(Ncellsaboveunstim/Ncellstotunstim)*100
        All_c[[q*n.batch + d]][k]=(Allaboveunstim/Ncellstotunstim)*100
        
      }
    }
    
    tempvec=subset(CD8df[[q*n.batch + d]],clonotype==Clonotypes_data$cdr3[k])$mark
    tempvectot=subset(Alldf[[q*n.batch + d]],clonotype==Clonotypes_data$cdr3[k])$mark
    if(length(intersect(Clonotypes_data$cdr3[k],TCR_seq[[q*n.batch + d]]))>0){
      Ncellstotunstim=length(tempvectot)
      Ncellsaboveunstim=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD8_c[[q*n.batch + d]][k]=0
      if(Ncellstotunstim>0){
        
        CD8_c[[q*n.batch + d]][k]=(Ncellsaboveunstim/Ncellstotunstim)*100
        
      }
    }
    
    Total_c[[q*n.batch + d]][k]=Ncellstotunstim
    Ncellstotunstim=0
  }
  
  Summarydf=data.frame(matrix(ncol = conditions*n.batch*4+2, nrow = length(Clone)))
  Summarydf[,1]=Clone
  for(q in 1:conditions){
    for(d in 1:n.batch){
  
      Summarydf[,4*((q*n.batch + d)-1)+1]=CD4_c[[q*n.batch + d]]
      Summarydf[,4*((q*n.batch + d)-1)+2]=CD8_c[[q*n.batch + d]]
      Summarydf[,4*((q*n.batch + d)-1)+3]=All_c[[q*n.batch + d]]
      Summarydf[,4*((q*n.batch + d)-1)+4]=Total_c[[q*n.batch + d]]
      
    }
  }
  
  Amino_Acid_seq=rep(0,length(Summarydf$Clonotype))
  AAs=vector(mode = "list", length = conditions*n.batch)
  NTs=vector(mode = "list", length = conditions*n.batch)
  
  for(q in 1:conditions){
    for(d in 1:n.batch){
      AAs=Clonal_Obs[[q*n.batch + d]]@meta.data$cdr3
      NTs=Clonal_Obs[[q*n.batch + d]]@meta.data$cdr3_na
    }}
  
  for(b in 1:conditions){
    for(d in 1:n.batch){
  for(q in 1:length(Amino_Acid_seq)){
    if(length(intersect(NTs[[b*n.batch + d]],Summarydf$Clonotype[q]))>0){
      idx=match(Summarydf$Clonotype[q],NTs[[b*n.batch + d]])
      Amino_Acid_seq[q]=AAs[[b*n.batch + d]][idx]
    }
    
    else{
      print("Error no matching TCR found")
    }
  }
    }}
 
  Summarydf[,4*(conditions*n.batch)+2]=Amino_Acid_seq
  names(Summarydf) <- names.spreadsheet
  write.csv(Summarydf,clonotype_path)
  
}
