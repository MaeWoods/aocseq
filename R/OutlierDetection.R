IFNG<-FetchData(CellData[[1]], 'IFNG')
IFNG_GZMB<-FetchData(CellData[[1]], c('IFNG', 'GZMB'))
IFNG_df<-data.frame(data=IFNG$IFNG)

construct_tree<-function(data, current_height, max_height){
  if(current_height==0){
    output_df<<-unique(data)
    output_df$height<<-rep(-1, nrow(output_df))
  }
  sample_attribute<-data[sample(1:ncol(data), 1)]
  if(current_height == max_height|nrow(unique(data))==1){
    for (i in 1:nrow(unique(data))) {
      
      
      output_df$height[do.call(paste0, output_df[,1:ncol(data),drop=FALSE]) %in% do.call(paste0, unique(data)[i,1:ncol(data),drop=FALSE])]<<-current_height
      
       
      
      
    }
  }
  else{
    
    
    split_value<-runif(1, min=min(sample_attribute), max=max(sample_attribute))
    split_left<-data[sample_attribute[,1]<=split_value, ,drop=FALSE]
    split_right<-data[sample_attribute[,1]>split_value, , drop=FALSE]
    if(nrow(split_left)!=0){
      left<-construct_tree(split_left, current_height+1, max_height)
    }
    if(nrow(split_right)!=0){
      right<-construct_tree(split_right, current_height+1, max_height)
    }
    
  }
  if(all(output_df$height!=-1)){
    return(output_df)
  }
}

test<-construct_tree(IFNG, 0, 20)  


Iso_forest<-function(data, num_trees, max_height, subsample_count=nrow(data)){
  results_df<-unique(data)
  results_df$avg_height<-rep(-1, nrow(results_df))
  tree_counter<-0
  while(tree_counter<num_trees) {
    
    random_subsample<-data[sample(1:nrow(data), min(subsample_count, nrow(data))), ,drop=FALSE]
    tree<-construct_tree(random_subsample, 0, max_height)
    
    for(i in 1:nrow(tree)){
    matched_logical<-do.call(paste0, results_df[,1:ncol(data), drop=FALSE]) %in% do.call(paste0, tree[i,1:ncol(data), drop=FALSE])
    
    if(results_df$avg_height[matched_logical]==-1){
      
      results_df$avg_height[matched_logical]<-tree[i, 'height']
    }
    else{
      
      results_df$avg_height[matched_logical]<-mean(results_df$avg_height[matched_logical],tree[i, 'height'])
    }
    
    }
    tree_counter<-tree_counter+1
  
  }
  if(any(results_df$avg_height==-1)){
    print('Some data points were never sampled. Increase num_trees or subsample_count.')
  }
  return(results_df)
}
test<-Iso_forest(IFNG_GZMB, 100, 20, 1000)  

anomaly_score<-function(df){
  c<-2*(log(dim(df)[1]-1)+0.5772156649) - (2.0*(log(dim(df)[1]-1)/(log(dim(df)[1]*1.0))))
  df[,'anomaly_score']<-2^(-df[,'avg_height']/c)
  return(df)
}
marker.gene=c("IFNG","CRTAM","GZMB","CCR7","CCL1","SELL","CD28","CD27")

iForest_wrapper<-function(CombineData_output,
                          gene_list, 
                          num_trees, 
                          max_height, 
                          subsample_count=nrow(data)){
  gene_counts<-FetchData(CombineData_output,gene_list[1])
  for (i in 2:length(gene_list)) {
    gene_counts[,gene_list[i]]<-FetchData(CombineData_output, gene_list[i])
  }
  df<-Iso_forest(gene_counts,num_trees,max_height,subsample_count)
  output<-anomaly_score(df)
  return(output)
}
iForest_wrapper(CellData[[1]], marker.gene, 10, 20, 1000)
