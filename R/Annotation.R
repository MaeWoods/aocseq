#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' List cell types and the proportion of cells expressing a choice of marker genes
#' at normalized counts greater than the control
#'
#' This function imports a Seurat object list processed by the aocseq::CombineData
#' function and produces a spreadsheet of cell types raned by frequency.
#' The table displays the percentage of CD4 and CD8 cells expressing a choice of marker
#' genes at normalized counts greater than the control.
#'
#' @param cell.data A Seurat object pre-processed with aocseq::CombineData.
#' @param goi A marker of interest.
#' @param cell.type Clone identifier, e.g. cdr3_na, phenotype.
#' @param goi.threshold Quantile that determines high counts per cell.
#' @param TCR True or false if true, adds amino acid sequence of TCR if cell.type is set to cdr3_na.
#' @param index.control.input Index of the control sequencing sample for each condition.
#' @param conditions Number of experiments within ClonalData object
#' @param path String. Directory of the table output.
#' @param names.spreadsheet List. Control and sample names.
#' @param preset Bool. If set to 1, the control is part of the dataset. If 0, a threshold value for each condition must be set.
#' @param threshold.entry List. List of thresholds for each condition if preset is 0.
#' @return A data frame containing clonotypes and summary statistics of transcript expression.
#' @concept Annotation
#' @export
Annotatecell.types <- function(
    cell.data,
    goi,
    cell.type="cdr3_na",
    goi.threshold=.975,
    TCR=TRUE,
    index.control.input=1,
    conditions=1,
    path="",
    names.spreadsheet=c(-1),
    preset=1,
    threshold.entry=FALSE
){
  cell.typePos=0
  if(length(index.control.input)==1){
    index.control=rep(1,length(cell.data))
  }
  path=paste(path,paste(paste(goi,"/SummaryTable",sep=""),".csv",sep=""),sep="")
  print("path=")
  print(path)
  n.batch=length(cell.data)
  if(names.spreadsheet[1]==-1){
    names.spreadsheet=rep(" ",length(cell.data)*4)
    names.spreadsheet=append(paste("Clone ",cell.type,sep=""),names.spreadsheet)
    if(TCR==TRUE){
    names.spreadsheet=append(names.spreadsheet,"Clone (amino)")
    }
  }
  else{
    for(g in 1:length(cell.data)){
      if(g==1){
        names.spreadsheet=append(paste(names.spreadsheet[g],c(" CD4 (%)"," CD8 (%)"," (%)"," Cells"),sep=""),names.spreadsheet[(g+1):length(names.spreadsheet)])
      }
      else if(g<length(cell.data)){
        names.spreadsheet=append(append(names.spreadsheet[1:(4*(g-1))],paste(names.spreadsheet[4*(g-1)+1],c(" CD4 (%)"," CD8 (%)"," (%)"," Cells"),sep="")),
                                 names.spreadsheet[(4*(g-1)+g):length(names.spreadsheet)])
      }
      else{
        names.spreadsheet=append(names.spreadsheet[1:(4*(g-1))],paste(names.spreadsheet[4*(g-1)+1],c(" CD4 (%)"," CD8 (%)"," (%)"," Cells"),sep=""))
        
      }
    }
    names.spreadsheet=append(paste("Clone ",cell.type,sep=""),names.spreadsheet)
    if(TCR==TRUE){
    names.spreadsheet=append(names.spreadsheet,"Clone (amino)")
    }
  }
  goi_samp=vector(mode = "list", length = conditions*n.batch)
  cell.type_samp=vector(mode = "list", length = conditions*n.batch)
  cutoff=vector(mode = "list", length = conditions)
  #match goi
  for(q in 0:(conditions-1)){
    for(k in 1:n.batch){
      goi_samp[[q*n.batch + k]]=match(goi, row.names(cell.data[[q*n.batch + k]][["SCT"]]@data))
      cell.type_samp[[q*n.batch + k]]=match(cell.type, names(cell.data[[q*n.batch + k]]@meta.data))
      if((preset==1)&&(k==index.control[1])){
        cutoff[[(q+1)]]=quantile(cell.data[[index.control[(q+1)]]][["SCT"]]@data[goi_samp[[q*n.batch + k]],],goi.threshold)[[1]]
      }else if(preset==0){
        cutoff[[(q+1)]]=threshold.entry[q]
        
      }
    }
  }
  if(length(levels(factor(cell.type_samp)))>1){
    print("Error: cell.type meta data column not in same position for each sample")
  }
  cell.typePos=cell.type_samp[[1]]
  
  #Find intersection and union of all clonotypes
  allclonotypes = c()
  intersectclonotypes=c()
  if(length(cell.data)==1){
    allclonotypes = levels(factor(cell.data[[1]]@meta.data[[cell.typePos]]))
    intersectclonotypes=levels(factor(cell.data[[1]]@meta.data[[cell.typePos]]))
    
  }
  else{
    allclonotypes = union(cell.data[[1]]@meta.data[[cell.typePos]],cell.data[[2]]@meta.data[[cell.typePos]])
    intersectclonotypes=intersect(cell.data[[1]]@meta.data[[cell.typePos]],cell.data[[2]]@meta.data[[cell.typePos]])
    for(g in 1:length(cell.data)){
      if(g<3){
        
      }
      else{
        allclonotypes = union(allclonotypes,cell.data[[g]]@meta.data[[cell.typePos]])
        intersectclonotypes=intersect(intersectclonotypes,cell.data[[g]]@meta.data[[cell.typePos]])
        
      }
    }
    for(q in 0:(conditions-1)){
      for(k in 1:n.batch){
        allclonotypes = union(allclonotypes,cell.data[[q*n.batch + k]]@meta.data[[cell.typePos]])
        intersectclonotypes=intersect(intersectclonotypes,cell.data[[q*n.batch + k]]@meta.data[[cell.typePos]])
      }
    }
  }
  
  allclonotypes = setdiff(levels(factor(allclonotypes)),"unassigned")
  intersectclonotypes = setdiff(levels(factor(intersectclonotypes)),"unassigned")
  
  sample.clonotypes=vector(mode = "list", length = conditions*n.batch)
  
  if(preset==1){
    sample.clonotypes[[index.control[1]]]=setdiff(levels(factor(cell.data[[index.control[1]]]@meta.data[[cell.typePos]])),"unassigned")
    combine.clonotypes = sample.clonotypes[[index.control[1]]]
  }else{
    combine.clonotypes = setdiff(levels(factor(cell.data[[1]]@meta.data[[cell.typePos]])),"unassigned")
    combine.clonotypes = sample.clonotypes[[1]]
  }
  for(q in 0:(conditions-1)){
    for(k in 1:n.batch){
      if(((preset==1)&(index.control[1]!=k))|((preset==0)&(q>0))){
        sample.clonotypes[[q*n.batch + k]]=setdiff(setdiff(levels(factor(cell.data[[q*n.batch + k]]@meta.data[[cell.typePos]])),"unassigned"),combine.clonotypes)
        combine.clonotypes=union(combine.clonotypes,sample.clonotypes[[q*n.batch + k]])
      }
    }
    
  }
  
  #assign counts for ranking by UMO, then EBV, BKV - CMV
  print("Assigning counts for ranking...")
  NClonotypes=length(allclonotypes)
  Shared=rep(0,NClonotypes)
  Clonefreq=rep(0,NClonotypes)
  for(k in 1:NClonotypes){
    for(q in 0:(conditions-1)){
      for(d in 1:n.batch){
        if(length(intersect(allclonotypes[k],sample.clonotypes[[q*n.batch + d]]))>0){
          Clonefreq[k]=cell.data[[q*n.batch + d]]@meta.data$countcln[which(cell.data[[q*n.batch + d]]@meta.data[[cell.typePos]]==allclonotypes[k])][1]
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
  for(q in 0:(conditions-1)){
    for(d in 1:n.batch){
      CD4_c[[q*n.batch + d]]=rep(0,length(allclonotypes))
      CD8_c[[q*n.batch + d]]=rep(0,length(allclonotypes))
      All_c[[q*n.batch + d]]=rep(0,length(allclonotypes))
      Total_c[[q*n.batch + d]]=rep(0,length(allclonotypes))
      TCR_seq[[q*n.batch + d]]=cell.data[[q*n.batch + d]]@meta.data[[cell.typePos]]
      CD4df[[q*n.batch + d]]=data.frame(clonotype=subset(cell.data[[q*n.batch + d]],CD4cells=="1")@meta.data[[cell.typePos]],mark=subset(cell.data[[q*n.batch + d]],CD4cells=="1")@meta.data[[match(paste("Threshold_",goi,sep=""),names(cell.data[[q*n.batch + d]]@meta.data))]])
      CD8df[[q*n.batch + d]]=data.frame(clonotype=subset(cell.data[[q*n.batch + d]],CD8cells=="1")@meta.data[[cell.typePos]],mark=subset(cell.data[[q*n.batch + d]],CD8cells=="1")@meta.data[[match(paste("Threshold_",goi,sep=""),names(cell.data[[q*n.batch + d]]@meta.data))]])
      Alldf[[q*n.batch + d]]=data.frame(clonotype=cell.data[[q*n.batch + d]]@meta.data[[cell.typePos]],mark=cell.data[[q*n.batch + d]]@meta.data[[match(paste("Threshold_",goi,sep=""),names(cell.data[[q*n.batch + d]]@meta.data))]])
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
    for(q in 0:(conditions-1)){
      for(d in 1:n.batch){
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
    }
  }
  
  #Compute the median clone size (+1) across replicate samples
  AvgCloneSize=rep(0,NClonotypes)
  for(k in 1:NClonotypes){
    
    vec_Cs=0
    for(q in 0:(conditions-1)){
      for(d in 1:n.batch){
        vec_Cs=append(vec_Cs,Total_c[[q*n.batch + d]][k])
      }
    }
    vec_Cs=vec_Cs[-1]
    AvgCloneSize[k]=ceiling(mean(vec_Cs))
  }
  Summarydf=data.frame()
  if(TCR==TRUE){
  Summarydf=data.frame(matrix(ncol = ((conditions*n.batch)*4+2), nrow = length(Clone)))
  }
  else{
    Summarydf=data.frame(matrix(ncol = ((conditions*n.batch)*4+1), nrow = length(Clone)))
  }
  Summarydf[,1]=Clone
  for(q in 0:(conditions-1)){
    for(d in 1:n.batch){
      
      Summarydf[,1+(q*n.batch + d-1)*4+1]=CD4_c[[q*n.batch + d]]
      Summarydf[,1+(q*n.batch + d-1)*4+2]=CD8_c[[q*n.batch + d]]
      Summarydf[,1+(q*n.batch + d-1)*4+3]=All_c[[q*n.batch + d]]
      Summarydf[,1+(q*n.batch + d-1)*4+4]=Total_c[[q*n.batch + d]]
      
    }
  }
  names(Summarydf) <- names.spreadsheet
  Amino_Acid_seq=rep(0,length(Summarydf[,1]))
  AAs=vector(mode = "list", length = conditions*n.batch)
  NTs=vector(mode = "list", length = conditions*n.batch)
  
  
  for(q in 0:(conditions-1)){
    for(d in 1:n.batch){
      if(TCR==TRUE){
      AAs[[q*n.batch + d]]=cell.data[[q*n.batch + d]]@meta.data$cdr3
      }
      NTs[[q*n.batch + d]]=cell.data[[q*n.batch + d]]@meta.data[[cell.typePos]]
    }}
  
  if(TCR==TRUE){
  for(b in 0:(conditions-1)){
    for(d in 1:n.batch){
      for(q in 1:length(Amino_Acid_seq)){
        if(length(intersect(NTs[[b*n.batch + d]],Summarydf[,1][q]))>0){
          idx=match(Summarydf[,1][q],NTs[[b*n.batch + d]])
          Amino_Acid_seq[q]=AAs[[b*n.batch + d]][idx]
        }
        
        else{
          print(paste("No matching TCR found in sample",d))
        }
      }
    }}
  
  Summarydf[,4*(conditions*n.batch)+2]=Amino_Acid_seq
  }
  names(Summarydf) <- names.spreadsheet
  Summarydf$avg=AvgCloneSize
  IndexMedian=order(AvgCloneSize,decreasing=TRUE)
  Ncolumns=dim(Summarydf)[2]
  for(k in 1:Ncolumns){
    Summarydf[,k]=Summarydf[IndexMedian,k]
  }
  write.csv(Summarydf,path)
  return(Summarydf)
  
}
