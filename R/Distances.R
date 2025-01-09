#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Cell classifier function that assigns a mean distance per cell type based on several metrics used for single cell RNA sequencing classification.
#'
#' This function will compute a classification for each single cell in samples that have been activated
#' and sequenced. Takes as input two Seurat objects, one containing
#' the query cells and another containing the reference. A gene list
#' used to construct the distance and output path for plotting and storing the mean distance per cell type. 
#' Other parameters are listed for debugging, but can be left as default values.
#'
#' @param output.array A Seurat object containing cells that are to be assigned an Rscore.
#' @param reference.data Matrix of doubles. Reference data used to calculate the Rscore.
#' @param gene.list Vector. List of genes in the gene signature.
#' @param cell.types Vector. List of cell types for classification.
#' @param distance Can be set to 0,1,2 or 3. Choice of 0 returns the Mahalanobis distance, 1 returns the taxicab distance, 2 returns the z-score and 3 returns the isolation forest distance. 
#' @param scramble Bool. To asses the effect of the gene signature on the Rscore, select a random set of genes and compute their distance from the reference data.
#' @param withSCT Print progress bars and output
#' @param ntrees Print progress bars and output
#' @param maxheight Print progress bars and output
#' @param PRINT True or false
#' @param file.path path to save metadata as spreadsheet if PRINT is set to TRUE
#' @param verbose Print progress bars and output
#' @return A Seurat object containing aocseq cell classification metadata.
#' @concept Single cell analysis
#'
#' @export
ClassifyCells <- function(
    output.array,
    reference.data,
    gene.list,
    cell.types=c(),
    distance=0,
    scramble=FALSE,
    withSCT=FALSE,
    ntrees=10,
    maxheight=20,
    PRINT=FALSE,
    file.path=".",
    verbose = TRUE
){

  if(distance==0){
    
    G1l=match(gene.list,row.names(output.array@assays$RNA@counts))
    SPFlog1pPF=as.matrix(output.array@assays$RNA@counts)
    PF=log(1+(t((t(SPFlog1pPF)/mean(colSums(SPFlog1pPF))))))
    PFlog1pPF=t(t(PF)/mean(colSums(PF)))
    fullPFlog1pPF=PFlog1pPF
    cell.arrayPFlog1pPF=fullPFlog1pPF[G1l,colnames(output.array)]
    normalized.array <- replace(cell.arrayPFlog1pPF, is.na(cell.arrayPFlog1pPF), 0)

  Nspecific=dim(reference.data)[2]
  gene.listInd=match(gene.list,row.names(reference.data))
  if(scramble){
    gene.listInd2=ceiling(runif(length(gene.listInd),0,dim(normalized.array)[1]))
    print("scrambled genes")
  } else{
    gene.listInd2=match(gene.list,row.names(normalized.array))
    print("Genes matched...")
  }
  SignatureCells=reference.data[gene.listInd,]
  sigvec=rep(0,dim(SignatureCells)[1]*dim(SignatureCells)[2])
  for(s in 1:dim(SignatureCells)[1]){
    for(g in 1:dim(SignatureCells)[2]){
      sigvec[(s-1)*dim(SignatureCells)[2] + g] = SignatureCells[s,g]
    }
  }
  sigvec = as.double(sigvec)

  TestCells=normalized.array[gene.listInd2,]
  Tref=rep(0,dim(TestCells)[1]*dim(TestCells)[2])
  for(s in 1:dim(TestCells)[1]){
    for(g in 1:dim(TestCells)[2]){
      Tref[(s-1)*dim(TestCells)[2] + g] = TestCells[s,g]
    }
  }
  Nsig=dim(SignatureCells)[2]
  Ngenes=length(gene.list)
  Ncells=dim(TestCells)[2]
  SpecificityDistance=rep(0,dim(TestCells)[2])
  TestCells=Tref
  dist=rep(0,Ncells)
  meansinp=rep(0,Ngenes)
  colsv=rep(0,Ngenes)
  covinpt=rep(0,Ngenes*Ngenes)
  sigtest=rep(0,Ngenes*(Ncells+1))

  SQ=GetMahalanobis(Nsig, Ncells, Ngenes, as.double(sigvec), as.double(TestCells),as.double(dist),as.double(meansinp),as.double(colsv),as.double(covinpt),as.double(covinpt),as.double(covinpt),as.double(covinpt),as.double(sigtest))
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
    
    G1l=match(gene.list,row.names(output.array@assays$RNA@counts))
    SPFlog1pPF=as.matrix(output.array@assays$RNA@counts)
    PF=log(1+(t((t(SPFlog1pPF)/mean(colSums(SPFlog1pPF))))))
    PFlog1pPF=t(t(PF)/mean(colSums(PF)))
    fullPFlog1pPF=PFlog1pPF
    cell.arrayPFlog1pPF=fullPFlog1pPF[G1l,colnames(output.array)]
    normalized.array <- replace(cell.arrayPFlog1pPF, is.na(cell.arrayPFlog1pPF), 0)
    
    Nspecific=dim(reference.data)[2]
    gene.listInd=match("IFNG",row.names(reference.data))
    if(scramble){
      gene.listInd2=ceiling(runif(length(gene.listInd),0,dim(normalized.array)[1]))
      print("scrambled genes")
    } else{
      gene.listInd2=match("IFNG",row.names(normalized.array))
      print("Genes matched...")
    }
    SignatureCells=reference.data[gene.listInd,]
    TestCells=normalized.array[gene.listInd2,]
    lenreTref=length(unname(TestCells))
    lenref=length(unname(SignatureCells))
    distarrayctrl=rep(0,lenreTref)
    meansig=mean(unname(SignatureCells))
    stdsig=sd(unname(SignatureCells))
    for(j in 1:lenreTref){
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
    
    G1l=match(gene.list,row.names(output.array@assays$RNA@counts))
    SPFlog1pPF=as.matrix(output.array@assays$RNA@counts)
    PF=log(1+(t((t(SPFlog1pPF)/mean(colSums(SPFlog1pPF))))))
    PFlog1pPF=t(t(PF)/mean(colSums(PF)))
    fullPFlog1pPF=PFlog1pPF
    cell.arrayPFlog1pPF=fullPFlog1pPF[G1l,colnames(output.array)]
    normalized.array <- replace(cell.arrayPFlog1pPF, is.na(cell.arrayPFlog1pPF), 0)
    
    Nspecific=dim(reference.data)[2]
    gene.listInd=match(gene.list,row.names(reference.data))
    if(scramble){
      gene.listInd2=ceiling(runif(length(gene.listInd),0,dim(normalized.array)[1]))
      print("scrambled genes")
    } else{
      gene.listInd2=match(gene.list,row.names(normalized.array))
      print("Genes matched...")
    }
    SignatureCells=reference.data[gene.listInd,]
    TestCells=normalized.array[gene.listInd2,]
    lenreTref=dim(TestCells)[2]
    lenref=dim(SignatureCells)[2]
    distarrayctrl=rep(0,lenreTref)
    for(j in 1:lenreTref){
      sumj=0
      for(g in 1:lenref){
        for(h in 1:length(gene.listInd)){
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
    
    cell.types=levels(factor(output.array@meta.data$celltype))
    
    cell_outliers<-CellTypeLoop(output.array,cell.types,reference.data,gene.list,
              ntrees,maxheight,subsample.count=ncol(reference.data)+1)
    
    output.array=AddMetaData(output.array, cell_outliers$outlier_fraction, col.name = 'Isoforest')
    
  }
  
  if(PRINT){
    write_xlsx(data.frame(output.array@meta.data),file.path)
  }
  
  return(output.array)

}

#' Return a data frame of meta data associated with single cell data sets that
#' contains a distance from a reference that is determined using the isolation forest algorithm
#'
#' This function will return the percentage of cells that are outliers when compared to a reference data set.
#'
#' @param cell.data A Seurat object containing cells that are to be assigned a distance.
#' @param cell.types A list of the cell types (currently this is clonotypes).
#' @param reference.data Normalized reference data.
#' @param genes Genes to pull from the query data set.
#' @param ntrees Number of trees used for the isolation forest.
#' @param maxheight Maximum height used for the isolation forest.
#' @param solver True or false If true an Rcpp implementation of the software is used
#' @param subsample.count Sub sampling number for isolation forest.
#' @return A data frame containing the mean normalized isolation forest score for each cell type.
#' @concept Single cell analysis
#'
#' @export
CellTypeLoop<-function(
    cell.data, 
    cell.types,
    reference.data, 
    genes, 
    ntrees, 
    maxheight, 
    solver=TRUE,
    subsample.count=ncol(reference.data)+1){
  
  if(solver){
  unique_types<-levels(factor(cell.types))
  df<-data.frame(clones=unique_types, outlier_fraction=vector(mode = 'numeric', length = length(unique_types)))
  for (i in 1:length(unique_types)) {
    
    ##Select only counts from signature genes
    cell.data=as.matrix(subset(cell.data,`celltype`==unique_types[i])@assays$SCT@counts[genes,])
    ##Set clonotype score will hold normalized heights for all cells in the clonotype
    numcells=dim(cell.data)[2]
    test_df<-as.data.frame(cbind(as.matrix(reference.data), cell.data[, 1]))
    TestCells=unname(as.matrix(test_df))
    Tref=rep(0,dim(TestCells)[1]*dim(TestCells)[2])
    Tsg=rep(0,dim(TestCells)[2])
    for(s in 1:dim(TestCells)[1]){
      for(g in 1:dim(TestCells)[2]){
        Tref[(s-1)*dim(TestCells)[2] + g] = test_df[s,g]
        Tsg[g] = test_df[1,g]
      }
    }
    
    allcellsF=unname(as.matrix(as.data.frame(as.matrix(cell.data))))
    Tcd=rep(0,dim(allcellsF)[1]*dim(allcellsF)[2])
    Tsg2=rep(0,numcells)
    for(s in 1:dim(allcellsF)[1]){
      for(g in 1:numcells){
        Tcd[(s-1)*numcells + g] = allcellsF[s,g]
        Tsg2[g] = 0
      }
    }
    
    TsgKurt=rep(0,n_genes)
    for(g in 1:n_genes){
      TsgKurt[g] = 0
    }
    
    dimN=dim(test_df)[2]
    maxkurtosis=0
    clone_result=isoForest(ntrees, n_genes,dimN,numcells, 
    as.double(Tref),as.double(Tsg),as.double(Tsg),as.double(Tsg), 
    as.double(Tsg),as.double(Tsg),as.double(TsgKurt),as.double(Tsg),
    as.double(Tsg),as.double(Tsg),as.double(Tsg),maxkurtosis,as.double(Tcd),as.double(Tsg2))
    
    SetClonotypeScores=clone_result
    IsodataFrame$isoF[f]=mean(clone_result)
  }
  df$outlier_fraction = IsodataFrame$isoF
  }
  else{
    for (i in 1:length(unique_types)) {
    print(paste0('processing ', unique_types[i]))
    select_clone<-subset(cell.data, cdr3_na==unique_types[i])
    select_clone_mat<-select_clone@assays$SCT@data
    clone_result<-percent_outlier(select_clone_mat, reference.data, genes, ntrees, maxheight, subsample.count)
    df[i,'outlier_fraction']<-clone_result$outlier_fraction
    }
  }
  return(df)
}

#' This function normalizes the data from the R implementation of the isolation forest
#'
#' @param df A data frame containing isolation forest distances
#' @return Normalized values
#' @concept Routine functions
#'
#' @export
NormalizationScore<-function(
    df){
  
  c<-2*(log(dim(df)[1]-1)+0.5772156649) - (2.0*(log(dim(df)[1]-1)/(log(dim(df)[1]*1.0))))
  df[,'normalization_score']<-2^(-df[,'avg_height']/c)
  return(df)
  
}

#' This function sets the cell.types as meta data.
#'
#' @param cell.data A Seurat Object.
#' @param cell.types Single cell metadata.
#' @return A Seurat Object.
#' @concept Annotation
#'
#' @export
SetCellType<-function(
    cell.data, 
    cell.types
    ){
  cell.data=AddMetaData(cell.data,cell.types,"celltype")
  return(cell.data)
}

#' Return a data frame of meta data associated with single cell data sets that
#' contains a distance from a reference
#'
#' This function will return the percentage of cells that are outliers when compared to a reference data set.
#'
#' @param cell.data A Seurat object containing cells that are to be assigned a distance.
#' @param reference.data Normalized reference data.
#' @param genes Genes to pull from the query dataset.
#' @param ntrees Number of trees used for the isolation forest.
#' @param maxheight Maximum height used for the isolation forest.
#' @param subsample.count Subsampling number for isolation forest.
#' @return A data frame with the percentage of cells that are outliers.
#' @concept Single cell analysis
#'
#' @export
percent_outlier<-function(
    cell.data, 
    reference.data, 
    genes, 
    ntrees, 
    maxheight,
    subsample.count=ncol(reference.data)+1){

    ################################################################################################
    ### This section of the code implements a version of the isolation forest written in Rscript ###
    ################################################################################################
  numcells<-ncol(cell.data)
  cell.data<-cell.data[genes,]
  reference.data<-reference.data[genes,]
  num_outliers<-0
  normalization_score_list<-c()
  height_list=0
  for (i in 1:numcells) {
    cell<-cell.data[, i]
    if(length(genes)>1){
      working_set<-cbind(as.matrix(reference.data), cell)
      
      test_df<-as.data.frame(working_set[1,])
      for (j in 2:length(genes)) {
        working_gene=genes[j]
        test_df[,working_gene]<-working_set[j,]
      }
    }
    else{
      working_set<-rbind(as.matrix(reference.data), cell)
      test_df<-as.data.frame(working_set)
    }
    
    colnames(test_df)<-genes
    height_df<-Iso_forest(test_df, ntrees, maxheight, subsample.count,kurtosis_param=TRUE)
    results_df<-NormalizationScore(height_df)
    cell_results<-results_df[do.call(paste0, results_df[,1:length(genes), drop=FALSE]) == do.call(paste0, as.list(cell)), ]
    if(cell_results$normalization_score>0.75){
      num_outliers=num_outliers+1
    }
    height_list<-append(height_list, cell_results$avg_height)
    normalization_score_list<-append(normalization_score_list, cell_results$normalization_score)
  }
  height_list=height_list[-1]
  return(data.frame(
    clone_height=mean(height_list), 
    clone_AS=mean(normalization_score_list), 
    outlier_fraction=num_outliers/numcells)
  )
  
}
