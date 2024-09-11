#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' These are functions to produce reference datasets to compute a response score of
#' activated or perturbed T cells.
#'
#' This function returns a matrix of normalized gene expression along the rows and cells
#' along the columns.
#'
#' @param Perturbation List of Seurat objects. Seurat objects containing expression data processed.
#' by aocseq that will be used to score query T cell response.
#' @param path.glist path to the gene list of the gene signature.
#' @param TCR.list Character list. List of TCRs for clonotypes to include in the reference data.
#' @param sample.return Index of samples from the perturebation list that should be returned in the matrix.
#' @param n.inlist Number of genes from the gene signature that must have high expression for the cell to be included in the reference.
#' @param verbose Print progress bars and output.
#' @return A reference matrix.
#' @concept annotation.
#' @export
MakeReference <- function(
    Perturbation,
    path.glist,
    TCR.list,
    sample.return=3,
    n.inlist=10,
    verbose=TRUE
){

  n.samples=length(Perturbation)

  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  for(k in 1:n.samples){
  Perturbation[[k]] <- CellCycleScoring(Perturbation[[k]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  }
  Sig = vector(mode = "list", length = n.samples)
  Perturbation.norm = vector(mode = "list", length = n.samples)
#All tetramer sorted clonotypes

for(g in 1:n.samples){
  Perturbation.norm[[g]]=Perturbation[[g]]
  Perturbation.norm[[g]]=as.matrix(Perturbation[[g]]@assays$RNA@counts)
  Perturbation.norm[[g]]=log(1+(Perturbation.norm[[g]]/mean(colSums(Perturbation.norm[[g]]))))
  Perturbation.norm[[g]]=Perturbation.norm[[g]]/mean(colSums(Perturbation.norm[[g]]))
}
for(g in 1:n.samples){
Glist=setdiff(read.csv(path.glist)$x,"")[1:25]
cutoff=rep(0,length(Glist))
GlistIndS=match(Glist,row.names(Perturbation.norm[[g]]))

for(g in 1:25){
  cutoff[g]=quantile(Perturbation.norm[[g]][GlistIndS[g],],.915)
}

Thresholds=matrix("unassigned",nrow=dim(Perturbation[[g]])[2],ncol=length(GlistIndS))
for(s in 1:length(GlistIndS)){
  vec1=Perturbation1[GlistIndS[s],]
  for(j in 1:dim(Perturbation[[g]])[2]){
    if(vec1[j]>cutoff[s]){
      Thresholds[j,s]="high"
    }
  }
  Perturbation[[g]]=AddMetaData(Perturbation[[g]], Thresholds[,s], col.name = paste("Signature_",Glist[s],sep=""))
}

MatrixOfValues=as.matrix(Perturbation[[g]]@meta.data[19:43])
NcellsSigmat=dim(MatrixOfValues)[[1]]
SignatureCell=rep(0,NcellsSigmat)
for(h in 1:NcellsSigmat){
  if(length(subset(MatrixOfValues[h,],MatrixOfValues[h,]=="high"))>n.inlist){
    SignatureCell[h]=1
  }
}

Perturbation[[g]]=AddMetaData(Perturbation[[g]], SignatureCell, col.name = "SignatureCell")
signature.ref <- subset(Perturbation[[g]],(cdr3_na %in% TCR.list) & (SignatureCell==1)
                         & Signature_CCL4L2=="high"  & Signature_IFNG=="high")

colnamesSig2=colnames(signature.ref)
GlistMat2=match(Glist,row.names(Perturbation.norm[[g]]))
Sig[[g]]=Perturbation1[GlistMat2,colnamesSig2]
}

##Uses 1 and 3 for whole and RAD subsets for the TA - specific TCRs
Sigall_T=Sig[[sample.return]]

dim(Sigall_T)

return(Sigall_T)
}

