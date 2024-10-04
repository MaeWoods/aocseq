library(moments)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# construct_tree
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Build a binary tree from a numerical data frame
#' 
#' This function will build a binary tree from a numerical data frame and return
#' a data frame containing each unique data point from the input data and the
#' height at which each data point becomes "isolated" (reaches an external node)
#' or the max_height.
#' 
#' @param data a numerical data frame
#' @param current_height a parameter which tracks the current height of a given construct_tree instance, allows for recursive calling of the function
#' @param max_height the maximum height the tree will grow to before stopping
#' 
#' @return a data frame containing the data for each unique point in the input data and the height at which that point was isolated
#' @export

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

construct_tree_kurt<-function(data, current_height, max_height){
  
  
  if(current_height==0){
    output_df<<-unique(data)
    output_df$ID<<-do.call(paste0, output_df[,1:ncol(data),drop=FALSE])
    output_df$height<<-rep(-1, nrow(output_df))
  }
  
  
  if(current_height == max_height|nrow(unique(data))==1){
    
    data$ID<-do.call(paste0, data)
    
    output_df$height[output_df$ID%in%data$ID]<<-current_height
    
      
      
      
    }
  
  else{
    data_pared<-data[unlist(lapply(data, function(x){length(unique(x))!=1}))]
    
    sample_attribute<-data_pared[kurtosis(data_pared)==max(kurtosis(data_pared), na.rm = TRUE)]
    
    split_value<-runif(1, min=min(sample_attribute), max=max(sample_attribute))
    split_left<-data[sample_attribute[,1]<=split_value, ,drop=FALSE]
    split_right<-data[sample_attribute[,1]>split_value, , drop=FALSE]
    
    if(nrow(split_left)!=0){
      left<-construct_tree_kurt(split_left, current_height+1, max_height)
    }
    if(nrow(split_right)!=0){
      right<-construct_tree_kurt(split_right, current_height+1, max_height)
    }
    
  }
  
  if(all(output_df$height!=-1)){
    return(output_df)
  }
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Iso_forest
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Build many trees on samples from data and average the heights
#' 
#' @param data a numerical data frame
#' @param num_trees the number of trees to build and average over
#' @param max_height the maximum height each tree will grow to before stopping
#' @param subsample_count the size of the random sample that a tree will be built on. Cannot be larger than the number of rows in the data
#' 
#' @return a data frame containing the data for each unique point in the input data and the average height at which that point was isolated
#' @export


Iso_forest<-function(data, num_trees, max_height, subsample_count=nrow(data)){
  results_df<-unique(data)
  results_df$avg_height<-rep(-1, nrow(results_df))
  tree_counter<-0
  while(tree_counter<num_trees) {
    print(tree_counter)
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

Iso_forest_kurt<-function(data, num_trees, max_height, subsample_count=nrow(data)){
  results_df<-unique(data)
  results_df$avg_height<-rep(-1, nrow(results_df))
  results_df$ID<-do.call(paste0, results_df[,1:ncol(data), drop=FALSE])
  tree_counter<-0
  while(tree_counter<num_trees) {
    #print(tree_counter)
    random_subsample<-data[sample(1:nrow(data), min(subsample_count, nrow(data))), ,drop=FALSE]
    tree<-construct_tree_kurt(random_subsample, 0, max_height)
    tree$ID<-do.call(paste0, tree[,1:ncol(data), drop=FALSE])
    inds <- results_df$ID%in%tree$ID
    results_df$avg_height[!is.na(inds)&results_df$avg_height==-1]<-na.omit(tree$height[inds])
    results_df$avg_height[!is.na(inds)&results_df$avg_height!=-1]<-rowMeans(cbind(results_df$avg_height[!is.na(inds)&results_df$avg_height!=-1],na.omit(tree$height[inds])))

    
    tree_counter<-tree_counter+1
    }
    
    
  
  if(any(results_df$avg_height==-1)){
    print('Some data points were never sampled. Increase num_trees or subsample_count.')
  }
  return(subset(results_df, select=-ID))
}

Iso_forest_kurt_par<-function(data, num_trees, max_height, subsample_count=nrow(data)){
  results_df<-unique(data)
  results_df$avg_height<-rep(-1, nrow(results_df))
  results_df$ID<-do.call(paste0, results_df[,1:ncol(data), drop=FALSE])
  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl)
  foreach(i=1:num_trees)%dopar% {
    #print(tree_counter)
    random_subsample<-data[sample(1:nrow(data), min(subsample_count, nrow(data))), ,drop=FALSE]
    tree<-construct_tree_kurt(random_subsample, 0, max_height)
    tree$ID<-do.call(paste0, tree[,1:ncol(data), drop=FALSE])
    inds <- results_df$ID%in%tree$ID
    results_df$avg_height[!is.na(inds)&results_df$avg_height==-1]<-na.omit(tree$height[inds])
    results_df$avg_height[!is.na(inds)&results_df$avg_height!=-1]<-rowMeans(cbind(results_df$avg_height[!is.na(inds)&results_df$avg_height!=-1],na.omit(tree$height[inds])))
    
    
    
  }
  
  parallel::stopCluster(cl)
  
  if(any(results_df$avg_height==-1)){
    print('Some data points were never sampled. Increase num_trees or subsample_count.')
  }
  return(subset(results_df, select=-ID))
}
 

anomaly_score<-function(df){
  c<-2*(log(dim(df)[1]-1)+0.5772156649) - (2.0*(log(dim(df)[1]-1)/(log(dim(df)[1]*1.0))))
  df[,'anomaly_score']<-2^(-df[,'avg_height']/c)
  return(df)
}


iForest_wrapper<-function(CombineData_output,
                          gene_list, 
                          num_trees, 
                          max_height, 
                          subsample_count=ncol(CombineData_output)){
  gene_counts<-FetchData(CombineData_output,gene_list[1])
  for (i in 2:length(gene_list)) {
    gene_counts[,gene_list[i]]<-FetchData(CombineData_output, gene_list[i])
  }
  df<-Iso_forest_kurt(gene_counts,num_trees,max_height,subsample_count)
  output<-anomaly_score(df)
  return(output)
}
Iso_forest_kurt_par(IFNG, 5,10)
iForest_wrapper(CMV[[1]], SCPA_genes, 5, 10, 2500)
IFNG<-FetchData(CMV[[1]], c('IFNG', 'GZMB'))
colnames(IFNG[kurtosis(IFNG)==max(kurtosis(IFNG))])
d<-kurtosis(IFNG)
d['GZMB']<-NaN
max(d, na.rm = TRUE)
max(c(1, NaN), na.rm = TRUE)
IFNG[unlist(lapply(IFNG, function(x){length(unique(x))!=1}))]

library(isotree)
isotree_model<-isolation.forest(IFNG, sample_size = 1000, output_dist = TRUE)
