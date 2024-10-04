library(Seurat)
library(Matrix)
control<-CombineData('C:/Projects/aocseq/Examples/AOCSEQ_Example/SeqData/gex/Control',
                  c("IFNG","CRTAM","GZMB","CCR7","CCL1","SELL","CD28","CD27"),
                  'C:/Projects/aocseq/Examples/AOCSEQ_Example/SeqData/vdj/Control/all_contig_annotations.csv',
                  file.saved = 'control.rds',
                  sample.name = 'control_data',
                  nvariable_features = 500)

CMV<-CombineData('C:/Projects/aocseq/Examples/AOCSEQ_Example/SeqData/gex/CMV',
                 c("IFNG","CRTAM","GZMB","CCR7","CCL1","SELL","CD28","CD27"),
                 'C:/Projects/aocseq/Examples/AOCSEQ_Example/SeqData/vdj/CMV/all_contig_annotations.csv',
                 file.saved = 'CMV.rds',
                 sample.name = 'CMV_data',
                 nvariable_features = 500)

clonotype1<-subset(CMV[[1]], clonotype=='clonotype1')
clonotype1_mat<-clonotype1@assays$SCT@data
control_mat<-control[[1]]@assays$SCT@data

IndexToPointer <- function(j) {
  p <- vector(mode = 'integer', length = max(j) + 1)
  index <- seq.int(from = 2, to = length(x = p))
  for (i in seq_along(along.with = index)) {
    p[index[i]] <- sum(j <= i)
  }
  return(p)
}
PointerToIndex <- function(p) {
  dp <- diff(x = p)
  j <- rep.int(x = seq_along(along.with = dp), times = dp)
  return(j)
}
Transpose.dgCMatrix <- function(x, ...) {
  i.order <- order(slot(object = x, name = 'i'))
  return(sparseMatrix(
    i = PointerToIndex(p = slot(object = x, name = 'p'))[i.order],
    p = IndexToPointer(j = slot(object = x, name = 'i') + 1),
    x = slot(object = x, name = 'x')[i.order],
    dims = rev(x = dim(x = x)),
    dimnames = rev(x = dimnames(x = x)),
    giveCsparse = TRUE
  ))
}

cell<-clonotype1_mat[c('IFNG'),1]

rownames(cell)<-colnames(clonotype1_mat[,1])
test_set<-rbind(as.matrix(control_mat[c('IFNG'),]), cell)
colnames(test_set[1,])

test_df<-as.data.frame(test_set[1,])
test_df$GZMB<-test_set[2,]
colnames(test_df)<-c('IFNG', 'GZMB')
Iso_forest(test_df, 25, 20, 1000)

percent_outlier<-function(test_set, control_set, genes, num_trees, max_height, subsample_count=ncol(control_set)+1, cutoff=.75){
  numcells<-ncol(test_set)
  test_set<-test_set[genes,]
  control_set<-control_set[genes,]
  num_outliers<-0
  height_vec<-c()
  anomaly_score_list<-c()
  for (i in 1:numcells) {
    print(i)
    
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
    
    height_df<-Iso_forest_kurt(test_df, num_trees, max_height, subsample_count)
    
    results_df<-anomaly_score(height_df)
    cell_results<-results_df[do.call(paste0, results_df[,1:length(genes), drop=FALSE]) == do.call(paste0, as.list(cell)), ]
    
    if(cell_results$anomaly_score>0.75){
      num_outliers=num_outliers+1
    }
    height_list<-append(height_list, cell_results$avg_height)
    anomaly_score_list<-append(anomaly_score_list, cell_results$anomaly_score)
  }
  return(data.frame(clone_height=mean(height_list), clone_AS=mean(anomaly_score_list), outlier_fraction=num_outliers/numcells))
}
percent_outlier(clonotype1_mat, control_mat, genes=c('IFNG'), 5, 20, (ncol(control_mat)+1)/2)

SCPA_genes<-c('CD69',
              'TNF',
              'TNFSF9',
              'CCL4L2',
              'CCL3',
              'IFNG',
              'ZFP36L1',
              'CCL4',
              'BTG2',
              'EGR2',
              'NFKBID',
              'PHLDA1',
              'SLA',
              'XCL1',
              'GZMB',
              'MIR155HG',
              'FASLG',
              'PRNP',
              'IRF4',
              'SRGN')

clonotype_loop<-function(clono_data, control_set, genes, num_trees, max_height, subsample_count=ncol(control_set)+1, cutoff=.75){
  unique_clones<-unique(clono_data[[1]]$clonotype)
  df<-data.frame(clones=unique_clones, outlier_fraction=vector(mode = 'numeric', length = length(unique_clones)))
  for (i in 1:length(unique_clones)) {
    print(paste0('processing ', unique_clones[i]))
    select_clone<-subset(clono_data[[1]], clonotype==unique_clones[i])
    select_clone_mat<-select_clone@assays$SCT@data
    clone_result<-percent_outlier(select_clone_mat, control_set, genes, num_trees, max_height, subsample_count, cutoff)
    df[i,'outlier_fraction']<-clone_result
  }
  return(df)
}
CMV_outliers<-clonotype_loop(CMV, control_mat, SCPA_genes, 10, 20, ncol(control_mat)+1, cutoff)


