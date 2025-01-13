#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' These are functions to produce reference datasets to compute a response score of
#' activated or perturbed T cells.
#'
#' This function returns a matrix of normalized gene expression along the rows and cells
#' along the columns.
#'
#' @param cell.data List of Seurat objects. Seurat objects containing expression data processed.
#' by aocseq that will be used to score query T cell response.
#' @param gene.list List of genes to be included in the reference matrix.
#' @param cell.type.list Character list. List of cell types to include in the reference data.
#' @param celltype Cell type to access metadata of reference dataset.
#' @param n.inlist Number of genes with high expression in each cell included in the reference matrix.
#' @param cellcycle Adds Seurat cell cycle phase to reference data.
#' @param SCT Returns reference matrix in sctransform space.
#' @param verbose Print progress bars and output.
#' @return A reference matrix.
#' @concept Single cell analysis
#' @export
MakeReference <- function(
    cell.data,
    gene.list,
    cell.type.list,
    threshold=0.975,
    n.inlist=10,
    cellcycle=FALSE,
    SCT=FALSE,
    verbose=TRUE
){
  
  n.samples=length(cell.data)
  if(cellcycle){
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  for(k in 1:n.samples){
    cell.data[[k]] <- CellCycleScoring(cell.data[[k]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  }
  }
  Sig = vector(mode = "list", length = n.samples)
  SigSCT = vector(mode = "list", length = n.samples)
  cell.data.norm = vector(mode = "list", length = n.samples)

  for(g in 1:n.samples){
    cell.data[[g]]=SetCellType(cell.data[[g]],cell.data[[g]]@meta.data$cdr3_na,"celltype")
    cell.data.norm[[g]]=as.matrix(cell.data[[g]]@assays$RNA@counts)
    cell.data.norm[[g]]=log(1+(cell.data.norm[[g]]/mean(colSums(cell.data.norm[[g]]))))
    cell.data.norm[[g]]=cell.data.norm[[g]]/mean(colSums(cell.data.norm[[g]]))
  }
  for(g in 1:n.samples){
    numG=length(gene.list)
    cutoff=rep(0,length(gene.list))
    gene.listIndS=match(gene.list,row.names(cell.data.norm[[g]]))
    
    for(s in 1:numG){
      cutoff[s]=quantile(cell.data.norm[[g]][gene.listIndS[s],],threshold)
      if(cutoff[s]==0){
        print(paste(paste("Gene: ",row.names(cell.data.norm[[g]][gene.listIndS[s],]),sep="")," contains zero counts. Matrix construction failed. Choose a different gene with more counts",sep="")
        )
      }
    }
    
    Thresholds=matrix("unassigned",nrow=dim(cell.data[[g]])[2],ncol=length(gene.listIndS))
    for(s in 1:length(gene.listIndS)){
      vec1=cell.data.norm[[g]][gene.listIndS[s],]
      for(j in 1:dim(cell.data[[g]])[2]){
        if(vec1[j]>=cutoff[s]){
          Thresholds[j,s]="high"
        }
      }
      cell.data[[g]]=AddMetaData(cell.data[[g]], Thresholds[,s], col.name = paste("Signature_",gene.list[s],sep=""))
    }
    
    gene_meta=match(paste("Signature_",gene.list[1],sep=""), names(cell.data[[g]]@meta.data))
    MatrixOfValues=as.matrix(cell.data[[g]]@meta.data[gene_meta:(gene_meta+(numG-1))])
    NcellsSigmat=dim(MatrixOfValues)[[1]]
    SignatureCell=rep(0,NcellsSigmat)
    ReferenceCells=0
    for(h in 1:NcellsSigmat){
      if(length(subset(MatrixOfValues[h,],MatrixOfValues[h,]=="high"))>n.inlist){
        SignatureCell[h]=1
        ReferenceCells=1
      }
    }
    if(ReferenceCells==0){
      print(paste("No reference cells in the matrix. Matrix construction failed, choose a high value for n.inlist. Value is currently: ",n.inlist,sep=""))
    }
    
    cell.data[[g]]=AddMetaData(cell.data[[g]], SignatureCell, col.name = "SignatureCell")
    signature.ref <- subset(cell.data[[g]],(celltype %in% cell.type.list) & (SignatureCell==1))
    colnamesSig2=colnames(signature.ref)
    gene.listMat2=match(gene.list,row.names(cell.data.norm[[g]]))
    Sig[[g]]=cell.data.norm[[g]][gene.listMat2,colnamesSig2]
    if(SCT){
      SigSCT[[g]]=cell.data[[g]]@assays$RNA@data[gene.listMat2,colnamesSig2]
    }
  }
  
  if(SCT){
  ref.matrixSCT=SigSCT[[1:n.samples]]
    return(ref.matrixSCT)
  }
  else{
    ref.matrix=Sig[[1:n.samples]]
  return(ref.matrix)
  }
}

