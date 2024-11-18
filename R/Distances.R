#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Return the Rscore of single cells in a query dataset compared to reference
#'
#' This function will compute the response score (Rscore) for each single cell in samples that have been activated
#' and sequenced. For each cell in the reference, we construct the positive definite covariance matrix
#' of the reference data set signature.ref. Takes as input two Seurat objects, one containing
#' the query cells and another containing the reference. A gene list
#' used to construct the distance and output path for plotting. Other parameters are listed for debugging,
#' but can be left
#' as default values.
#' there are different functions for different genesets
#'
#' @param output.array A Seurat object containing cells that are to be assigned an Rscore.
#' @param normalized.array A matrix of query data normalized counts.
#' @param signature.ref Matrix of doubles. Reference data used to calculate the Rscore.
#' @param Glist Vector. List of genes in the gene signature.
#' @param scramble Bool. To asses the effect of the gene signature on the Rscore, select a random set of genes and compute their distance from the reference data.
#' @param verbose Print progress bars and output
#' @return A Seurat object containing metadata with the Rscore.
#' @concept annotation
#'
#' @export
ClassifyCells <- function(
    output.array,
    normalized.array,
    signature.ref,
    Glist,
    distance=0,
    scramble=FALSE,
    withSCT=FALSE,
    verbose = TRUE
){

  if(distance==0){
  Nspecific=dim(signature.ref)[2]
  GlistInd=match(Glist,row.names(signature.ref))
  if(scramble){
    GlistInd2=ceiling(runif(length(GlistInd),0,dim(normalized.array)[1]))
    print("scrambled genes")
  } else{
    GlistInd2=match(Glist,row.names(normalized.array))
    print("Genes matched...")
  }
  SignatureCells=signature.ref[GlistInd,]
  sigvec=rep(0,dim(SignatureCells)[1]*dim(SignatureCells)[2])
  for(s in 1:dim(SignatureCells)[1]){
    for(g in 1:dim(SignatureCells)[2]){
      sigvec[(s-1)*dim(SignatureCells)[2] + g] = SignatureCells[s,g]
    }
  }
  sigvec = as.double(sigvec)

  TestCells=normalized.array[GlistInd2,]
  Tvec=rep(0,dim(TestCells)[1]*dim(TestCells)[2])
  for(s in 1:dim(TestCells)[1]){
    for(g in 1:dim(TestCells)[2]){
      Tvec[(s-1)*dim(TestCells)[2] + g] = TestCells[s,g]
    }
  }
  Nsig=dim(SignatureCells)[2]
  Ngenes=length(Glist)
  Ncells=dim(TestCells)[2]
  SpecificityDistance=rep(0,dim(TestCells)[2])
  TestCells=Tvec
  dist=rep(0,Ncells)
  meansinp=rep(0,Ngenes)
  colsv=rep(0,Ngenes)
  covinpt=rep(0,Ngenes*Ngenes)
  sigtest=rep(0,Ngenes*(Ncells+1))

  SQ=GetMahalanobis(Nsig, Ncells, Ngenes, as.double(sigvec), as.double(TestCells),as.double(dist),as.double(meansinp),as.double(colsv),as.double(covinpt),as.double(sigtest))
  distarrayctrl=array( unlist( SQ ))
  if(withSCT){
    if(scramble){
      output.array=AddMetaData(output.array, distarrayctrl, col.name = 'MdistSCTscr')
    }
    else{
      output.array=AddMetaData(output.array, distarrayctrl, col.name = 'MdistSCT')
    }
  }
  else{
  if(scramble){
    output.array=AddMetaData(output.array, distarrayctrl, col.name = 'Mdistscr')
  }
  else{
  output.array=AddMetaData(output.array, distarrayctrl, col.name = 'Mdist')
  }
  }
 }
  else if(distance==2){
    Nspecific=dim(signature.ref)[2]
    GlistInd=match("IFNG",row.names(signature.ref))
    if(scramble){
      GlistInd2=ceiling(runif(length(GlistInd),0,dim(normalized.array)[1]))
      print("scrambled genes")
    } else{
      GlistInd2=match("IFNG",row.names(normalized.array))
      print("Genes matched...")
    }
    SignatureCells=signature.ref[GlistInd,]
    TestCells=normalized.array[GlistInd2,]
    lenretvec=length(unname(TestCells))
    lenref=length(unname(SignatureCells))
    distarrayctrl=rep(0,lenretvec)
    meansig=mean(unname(SignatureCells))
    stdsig=sd(unname(SignatureCells))
    for(j in 1:lenretvec){
      sumj=0
      sumj=((abs(TestCells[j]-meansig))/stdsig)
      distarrayctrl[j]=sumj
    }
    if(withSCT){
      if(scramble){
        output.array=AddMetaData(output.array, distarrayctrl, col.name = 'MdistSCTscr')
      }
      else{
        output.array=AddMetaData(output.array, distarrayctrl, col.name = 'MdistSCT')
      }
    }
    else{
      if(scramble){
    output.array=AddMetaData(output.array, distarrayctrl, col.name = 'Mdistscr')
      }
      else{
        output.array=AddMetaData(output.array, distarrayctrl, col.name = 'Mdist')
      }
    }
  }
  else if(distance==1){
    Nspecific=dim(signature.ref)[2]
    GlistInd=match(Glist,row.names(signature.ref))
    if(scramble){
      GlistInd2=ceiling(runif(length(GlistInd),0,dim(normalized.array)[1]))
      print("scrambled genes")
    } else{
      GlistInd2=match(Glist,row.names(normalized.array))
      print("Genes matched...")
    }
    SignatureCells=signature.ref[GlistInd,]
    TestCells=normalized.array[GlistInd2,]
    lenretvec=dim(TestCells)[2]
    lenref=dim(SignatureCells)[2]
    distarrayctrl=rep(0,lenretvec)
    for(j in 1:lenretvec){
      sumj=0
      for(g in 1:lenref){
        for(h in 1:length(GlistInd)){
          sumj=sumj+((abs(SignatureCells[h,g]-TestCells[h,j]))/lenref)
        }

      }
      distarrayctrl[j]=sumj
    }
    if(withSCT){
      if(scramble){
        output.array=AddMetaData(output.array, distarrayctrl, col.name = 'MdistSCTscr')
      }
      else{
        output.array=AddMetaData(output.array, distarrayctrl, col.name = 'MdistSCT')
      }
    }
    else{
      if(scramble){
        output.array=AddMetaData(output.array, distarrayctrl, col.name = 'Mdistscr')
      }
      else{
        output.array=AddMetaData(output.array, distarrayctrl, col.name = 'Mdist')
      }
    }
  }
  else if(distance==3){
    
    cell_types=levels(factor(output.array@meta.data$cdr3_na))[1:3]
    
    CMV_outliers<-CellTypeLoop(output.array,cell_types,signature.ref,Glist,
              num_trees=10,max_height=20,subsample_count=ncol(control_set)+1, cutoff=.75)
    
  }
  return(output.array)

}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Return a data frame of meta data associated with single cell data sets that
#' contains a distance from a reference
#'
#' This function will compute the distance for each single cell. For
#' each cell in the reference and then return it in a data frame along with all the additional
#' single cell metadata
#'
#' @param cell.data A Seurat object containing cells that are to be assigned a distance.
#' @param signature.ref Matrix of doubles. Reference data used to calculate the Rscore.
#' @param path.glist Directory. Directory of the gene list used for gene signature, stored in csv format.
#' @param distance Choice of model. 0-Mahalanobis 1-taxicab 2-z-score.
#' @param scramble Bool. To asses the effect of the gene signature on the Rscore, select a random set of genes and compute their distance from the reference data.
#' @param verbose Print error log.
#' @export
AddDistances <- function(
    cell.data,
    signature.ref,
    path.glist,
    distance=0,
    scramble=FALSE,
    verbose = TRUE
){

  Glist=setdiff(read.csv(path.Glist)$x,"")[1:25]
  G1l=match(Glist,row.names(cell.data@assays$RNA@counts))
  SPFlog1pPF=as.matrix(cell.data@assays$RNA@counts)
  PF=log(1+(t((t(SPFlog1pPF)/mean(colSums(SPFlog1pPF))))))
  PFlog1pPF=t(t(PF)/mean(colSums(PF)))
  fullPFlog1pPF=PFlog1pPF
  cell.arrayPFlog1pPF=fullPFlog1pPF[G1l,colnames(cell.data)]
  cell.arrayPFlog1pPF <- replace(cell.arrayPFlog1pPF, is.na(cell.arrayPFlog1pPF), 0)

  cell.data=ClassifyCells(cell.data,cell.arrayPFlog1pPF,signature.ref,Glist,distance)

  if(PRINT){
    write_xlsx(data.frame(cell.data@meta.data),file.path)
  }
  return(cell.data)


}

#' Return a data frame of meta data associatd with single cell data sets that
#' contains a distance from a reference
#'
#' This function will return the percentage of cells that are outliers when compared to a reference data set.
#'
#' @param test_data A Seurat object containing cells that are to be assigned a distance.
#' @param cell_types A list of the cell types (currently this is clonotypes).
#' @param control_set Normalized reference data.
#' @param genes Genes to pull from the query dataset.
#' @param num_trees Number of trees used for the isolation forest.
#' @param max_height Maximum height used for the isolation forest.
#' @param subsample_count Subsampling number for isolation forest.
#' @param cutoff Maximum cutoff.
#' @export
CellTypeLoop<-function(
    test_data, 
    cell_types,
    control_set, 
    genes, 
    num_trees, 
    max_height, 
    subsample_count=ncol(control_set)+1, 
    cutoff=.75){
  
  unique_clones<-levels(factor(cell_types))
  df<-data.frame(clones=unique_clones, outlier_fraction=vector(mode = 'numeric', length = length(unique_clones)))
  for (i in 1:length(unique_clones)) {
    print(paste0('processing ', unique_clones[i]))
    select_clone<-subset(test_data, cdr3_na==unique_clones[i])
    select_clone_mat<-select_clone@assays$SCT@data
    clone_result<-percent_outlier(select_clone_mat, control_set, genes, num_trees, max_height, subsample_count, cutoff)
    df[i,'outlier_fraction']<-clone_result$outlier_fraction
  }
  return(df)
}

myClonotypes=read.csv("/Users/maewoodsphd/mVSTManuscript/SupplementaryTables/final_mVSTIFNG/SummaryTable.csv")
test_data=ThreeHourStim[[2]]
cell_types=myClonotypes$Clone..nucleic.[10:13]
control_set=ThreeHourStim[[1]]@assays$SCT@data
genes=c("IFNG","GZMB")
num_trees=10
max_height=20
subsample_count=ncol(control_set)+1
cutoff=.75

#' This function determines what constitues an outlier and what does not in the isolation forest algorithm
#'
#' @param df A data frame containing isolation forest distances
#' @export
anomaly_score<-function(df){
  c<-2*(log(dim(df)[1]-1)+0.5772156649) - (2.0*(log(dim(df)[1]-1)/(log(dim(df)[1]*1.0))))
  df[,'anomaly_score']<-2^(-df[,'avg_height']/c)
  return(df)
}

#' Return a data frame of meta data associatd with single cell data sets that
#' contains a distance from a reference
#'
#' This function will return the percentage of cells that are outliers when compared to a reference data set.
#'
#' @param test_data A Seurat object containing cells that are to be assigned a distance.
#' @param cell_types A list of the cell types (currently this is clonotypes).
#' @param control_set Normalized reference data.
#' @param genes Genes to pull from the query dataset.
#' @param num_trees Number of trees used for the isolation forest.
#' @param max_height Maximum height used for the isolation forest.
#' @param subsample_count Subsampling number for isolation forest.
#' @param cutoff Maximum cutoff.
#' @export
percent_outlier<-function(
    test_set, 
    control_set, 
    genes, 
    num_trees, 
    max_height, 
    subsample_count=ncol(control_set)+1, 
    cutoff=.75){
  numcells<-ncol(test_set)
  test_set<-test_set[genes,]
  control_set<-control_set[genes,]
  num_outliers<-0
  anomaly_score_list<-c()
  height_list=0
  for (i in 1:numcells) {
    cell<-test_set[, i]
    if(length(genes)>1){
      working_set<-cbind(as.matrix(control_set), cell)
      
      test_df<-as.data.frame(working_set[1,])
      for (j in 2:length(genes)) {
        working_gene=genes[j]
        test_df[,working_gene]<-working_set[j,]
      }
    }
    else{
      working_set<-rbind(as.matrix(control_set), cell)
      test_df<-as.data.frame(working_set)
    }
    
    colnames(test_df)<-genes
    height_df<-Iso_forest(test_df, num_trees, max_height, subsample_count,kurtosis_param=TRUE)
    results_df<-anomaly_score(height_df)
    cell_results<-results_df[do.call(paste0, results_df[,1:length(genes), drop=FALSE]) == do.call(paste0, as.list(cell)), ]
    if(cell_results$anomaly_score>0.75){
      num_outliers=num_outliers+1
    }
    height_list<-append(height_list, cell_results$avg_height)
    anomaly_score_list<-append(anomaly_score_list, cell_results$anomaly_score)
  }
  height_list=height_list[-1]
  return(data.frame(
    clone_height=mean(height_list), 
    clone_AS=mean(anomaly_score_list), 
    outlier_fraction=num_outliers/numcells)
  )
}
