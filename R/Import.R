##################################################
#' @import Seurat
#' @import matrixStats
#' @import readr
#'
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# EM algorithm
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#'
#' Calculate the posterior probabilities (soft labels) that each component
#' has to each data point. Iterative function for the expectation maximization (EM) algorithm.
#'
#' @title EStep
#' @param sd.vector Vector containing the standard deviations of each component.
#' @param sd.vector Vector containing the mean of each component.
#' @param alpha.vector Vector containing the mixing weights  of each component.
#' @return Named list containing the loglik and posterior.df.
#' @concept Routine functions
#' @export
EStep <- function(
    x, 
    mu.vector, 
    sd.vector, 
    alpha.vector) {
  
  comp1.prod <- dnorm(x, mu.vector[1], sd.vector[1]) * alpha.vector[1]
  comp2.prod <- dnorm(x, mu.vector[2], sd.vector[2]) * alpha.vector[2]
  sum.of.comps <- comp1.prod + comp2.prod
  comp1.post <- comp1.prod / sum.of.comps
  comp2.post <- comp2.prod / sum.of.comps

  sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
  sum.of.comps.ln.sum <- sum(na.omit(sum.of.comps.ln))

  return(list("loglik" = sum.of.comps.ln.sum,
       "posterior.df" = cbind(comp1.post, comp2.post)))
}

#' Maximization Step of the EM Algorithm
#'
#' Update the Component Parameters. Iterative function for the expectation maximization (EM) algorithm.
#'
#' @param x Input data.
#' @param posterior.df Posterior probability data.frame.
#' 
#' @return Named list containing the mean (mu), variance (var), and mixing weights (alpha) for each component.
#' @concept Routine functions
#' @export
MStep <- function(
    x, 
    posterior.df
){
  
  comp1.n <- sum(na.omit(posterior.df[, 1]))
  comp2.n <- sum(na.omit(posterior.df[, 2]))

  comp1.mu <- 1/comp1.n * sum(na.omit(posterior.df[, 1] * x))
  comp2.mu <- 1/comp2.n * sum(na.omit(posterior.df[, 2] * x))

  comp1.var <- sum(na.omit(posterior.df[, 1] * (x - comp1.mu)^2)) * 1/comp1.n
  comp2.var <- sum(na.omit(posterior.df[, 2] * (x - comp2.mu)^2)) * 1/comp2.n

  comp1.alpha <- comp1.n / length(x)
  comp2.alpha <- comp2.n / length(x)

  return(list("mu" = c(comp1.mu, comp2.mu),
       "var" = c(comp1.var, comp2.var),
       "alpha" = c(comp1.alpha, comp2.alpha)))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GMM demux
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Flexible implementation of the GMM demux algorithm
#'
#' Demultiplex arbitrary number of hashtags
#'
#' @param s.name This function creates a 10X genomics Seurat object. s.name sets the name of the project.
#' @param data Input data. 10x genomics Seurat object with hashtag antibodies.
#' @param hashtag.index List containing positions of hashtags to be de-multiplexed.
#' @param nameshashtags Sets the values of the sample for the new orig.ident to be added to the de-multiplexed Seurat object.
#' @param set.col.name Sets the column name to be added to the de-multiplexed Seurat object.
#' @param Seurat_Object True or false. Default is FALSE and assumes input data is a data frame or matrix.
#' 
#' @return A Seurat object that has been split by hastagged demultiplexing. 
#' @concept Statistical inference
#' @export
GMMDemux<-function(
    s.name,
    data, 
    hashtag.index, 
    nameshashtags, 
    set.col.name='Hashtags',
    Seurat_Object=FALSE
){


  
  denom_fn <- function(inputmat,nvar,pos) 
  {
    matrixStats::product((subset(inputmat,inputmat[,pos]>0)[,pos]))^(1/nvar)
  }
  
  Total_T_Cells=0
  pbmc.htos=0
  if(Seurat_Object==FALSE){
  Total_T_Cells=length(data$`Gene Expression`@p)
  pbmc.htos=data$`Antibody Capture`
  }
  else{
    Total_T_Cells=dim(data)[2]
    pbmc.htos=data[['HTO']]@counts
  }
  temp=t(as.matrix(pbmc.htos))+1
  tempqp=temp
  
  Denom_arr=c(denom_fn(temp,Total_T_Cells,1),
              denom_fn(temp,Total_T_Cells,2),
              denom_fn(temp,Total_T_Cells,3),
              denom_fn(temp,Total_T_Cells,4),
              denom_fn(temp,Total_T_Cells,5),
              denom_fn(temp,Total_T_Cells,6),
              denom_fn(temp,Total_T_Cells,7),
              denom_fn(temp,Total_T_Cells,8),
              denom_fn(temp,Total_T_Cells,9),
              denom_fn(temp,Total_T_Cells,10))
  
  clrfn <- function(j)
  {
    log(temp[,j]/Denom_arr[j])
  }
  
  tempqp=matrix(0, nrow = dim(temp)[1], ncol = dim(temp)[2]) 
  tempqp=lapply(list(1:10), clrfn)[[1]]
  tempqp[!is.finite(tempqp)] <- NA
  
  mu.cont<-c()
  var.cont<-c()
  alpha.cont<-c()
  #KMeans and model fitting loop
  for (i in 1:length(hashtag.index)) {

    tempqp[,hashtag.index[i]]<-log(temp[,hashtag.index[i]]/matrixStats::product((subset(temp,temp[,hashtag.index[i]]>0)[,hashtag.index[i]])^(1/Total_T_Cells)))
    wait <- tempqp[,hashtag.index[i]]

    wait.kmeans <- kmeans(na.omit(wait), 2)
    wait.kmeans.cluster <- wait.kmeans$cluster
    wait.df <- data.frame(x = na.omit(wait), num=rank(na.omit(wait),ties.method = "first"), cluster = wait.kmeans.cluster)
    wait.summary.df1=data.frame(cluster=c(1,2),
      mu=c(mean(subset(wait.df,cluster==1)$x),mean(subset(wait.df,cluster==2)$x)),
       variance = c(var(subset(wait.df,cluster==1)$x),var(subset(wait.df,cluster==2)$x)),
       std = c(sd(subset(wait.df,cluster==1)$x),sd(subset(wait.df,cluster==2)$x)),
       size=c(length(subset(wait.df,cluster==1)$x),length(subset(wait.df,cluster==2)$x)))
    wait.summary.df1$alpha=(wait.summary.df1$size)/sum(wait.summary.df1$size)

    print(wait.summary.df1)

    for (i in 1:50) {
      if (i == 1) {
        # Initialization
        e.step <- EStep(wait, wait.summary.df1[["mu"]], wait.summary.df1[["std"]],
                         wait.summary.df1[["alpha"]])
        m.step1 <- MStep(wait, e.step[["posterior.df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      }
      else {
        # Repeat E and M steps till convergence
        e.step <- EStep(wait, m.step1[["mu"]], sqrt(m.step1[["var"]]),
                         m.step1[["alpha"]])
        m.step1 <- MStep(wait, e.step[["posterior.df"]])
        loglik.vector <- c(loglik.vector, e.step[["loglik"]])

        loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
        if(loglik.diff < 1e-6) {
          break
        } else {
          cur.loglik <- e.step[["loglik"]]
        }
      }
    }
    mu.cont<-append(mu.cont, m.step1$mu)
    var.cont<-append(var.cont, m.step1$var)
    alpha.cont<-append(alpha.cont, m.step1$alpha)
  }

  Mhash<-matrix(0, nrow = dim(temp)[1], ncol = length(hashtag.index)+1)
  donorlabel=rep("unassigned",length(tempqp[,1]))
  samplelabels=rep("unassigned",length(tempqp[,1]))

  #hashtag probs loop
  for(j in 1:length(tempqp[,1])){
    high.cont<-c()
    low.cont<-c()
    for(k in hashtag.index){

      if(mu.cont[match(k,hashtag.index)*2]>mu.cont[(match(k,hashtag.index)*2)-1]){
        x=tempqp[j,k]
        low=alpha.cont[match(k,hashtag.index)*2-1]
        high=alpha.cont[match(k,hashtag.index)*2]
        xgivenPz_ihigh = dnorm(x, mean = mu.cont[match(k,hashtag.index)*2], sd = var.cont[match(k,hashtag.index)*2], log = FALSE)
        xgivenPz_ilow = dnorm(x, mean = mu.cont[match(k,hashtag.index)*2-1], sd = var.cont[match(k,hashtag.index)*2-1], log = FALSE)
        Px_i = xgivenPz_ilow*low + xgivenPz_ihigh*high
        K1Pz_ihighgiven_x = (xgivenPz_ihigh*high)/Px_i
        K1Pz_ilowgiven_x = (xgivenPz_ilow*low)/Px_i

      }
      else{
        x=tempqp[j,k]
        low=alpha.cont[match(k,hashtag.index)*2]
        high=alpha.cont[match(k,hashtag.index)*2-1]
        xgivenPz_ihigh = dnorm(x, mean = mu.cont[match(k,hashtag.index)*2-1], sd = var.cont[match(k,hashtag.index)*2-1], log = FALSE)
        xgivenPz_ilow = dnorm(x, mean = mu.cont[match(k,hashtag.index)*2], sd=var.cont[match(k,hashtag.index)*2], log = FALSE)
        Px_i = xgivenPz_ilow*low + xgivenPz_ihigh*high
        K1Pz_ihighgiven_x = (xgivenPz_ihigh*high)/Px_i
        K1Pz_ilowgiven_x = (xgivenPz_ilow*low)/Px_i

      }
      high.cont<-append(high.cont, K1Pz_ihighgiven_x)
      low.cont<-append(low.cont, K1Pz_ilowgiven_x)
    }

    probs<-vector(mode = 'numeric', length = length(high.cont)+1)
    #probs calculation
    for (i in 1:length(probs)) {
      if(i<length(probs)){
        probs[i]<-high.cont[i]*prod(low.cont[-i])
      }
      else{
        probs[i]<-1-sum(probs)
      }

    }

    #donor/hashtag assignment, maybe make this custom, e.g. not donor but a specified label
    for (i in 1:length(probs)) {
      if((max(probs)==probs[i])&(i<length(probs))){
        Mhash[j,i]=1
        donorlabel[j]=nameshashtags[i]
      }
      else if ((max(probs)==probs[i])&(i==length(probs))){
        Mhash[j,i]=1
        donorlabel[j]='none'
      }
    }
  }
  Full=matrix(NA, nrow = dim(temp)[1], ncol = dim(temp)[2])
  for(j in 1:length(temp[,1])){

    for(k in 1:length(temp[1,])){



      if(any(Mhash[j,]>0)){
        Full[j,k]=temp[j,k]
      }

      else{}

    }
  }

  colnames(Full) <- colnames(temp)
  row.names(Full) <- row.names(temp)
  Full1=na.omit(Full)
  donormat <- t(as.matrix(donorlabel))
  if(Seurat_Object==FALSE){
  colnames(donormat) <- row.names(temp)
  joint.bcs <- intersect(colnames(data$`Gene Expression`), colnames(t(Full1)))
  pbmc.htos <- t(Full1)[, joint.bcs]
  
  data=CreateSeuratObject(counts = data$`Gene Expression`,project = s.name)
  data[["HTO"]]=CreateAssayObject(counts = pbmc.htos)
  }
  else{
  data=AddMetaData(data, donorlabel, col.name = set.col.name)
  }

  for (i in 1:length(nameshashtags)) {
    for(j in 1:length(tempqp[,1])){
    if(donorlabel[j]==nameshashtags[i]){
      samplelabels[j]=nameshashtags[i]
    }
    }
  }
  data$orig.ident=samplelabels
  return(data)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CombineData
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Combine single cell data and annotate with cell labels based on functionality
#'
#' This function will read in 10X genomics outputs from cellranger and combine
#' multiple Seurat objects, assigning labels based on additional tables stored
#' in csv format. Takes as input gene expression from scRNAseq and VDJ enrichment.
#' A gene list used to estimate specificity can be chosen and is application specific.
#' Other parameters are listed for debugging, but can be left as default values.
#'
#' @param gex.path path to gene expression data in the format of cellranger feature barcode matrices.
#' @param marker.gene list of marker genes used for specificity analysis.
#' @param vdj.path path to VDJ expression data in the format of cellranger csv files.
#' @param threshold.cutoff list of marker genes used for specificity analysis.
#' @param file.saved directory for storing Seurat object RDS files.
#' @param index.control Index of the control sequencing sample.
#' @param sample.name Names that will be added to each sample orig.ident.
#' @param preset If set to 1, preset determines the threshold for cutoff if a control is included in the assay. If set to 0 a cutoff is set by the threshold.entry, which defines what value is high for all samples. If set to 2 a cutoff is set per sample.
#' @param threshold.entry Double. Sets the percentile above which gene expression is labelled as high.
#' @param demultiplex Boolean value used to indicate if the data was hashtagged and therefore requires to be split.
#' @param demultiplex.index Indexes to select columns with hashtag counts in the 10X genomics Seurat object assay. This depends on the 10X genomics chemistry and output from cellranger.
#' @param names.hashtag hashtag sample names.
#' @param n.ht.per.sample Number of hashtags used per sample (only set if hashtagging has been used).
#' @param tenX_conversion Parameter for 10X compatability.
#' @param nFeature_RNA_lower Threshold for the minimum number of unique gene identifiers detected in a single cell.
#' @param nFeature_RNA_upper Threshold for the maximum number of unique gene identifiers detected in a single cell.
#' @param nvariable_features Threshold for the maximum total number of unique gene identifiers detected in a single cell for dimension reduction.
#' @param percent.mt_upper Threshold for the maximum percentage of unique mitochondrial gene identifiers detected in a single cell.
#' @param verbose Print progress bars and output.
#' @param QC_plots Print progress bars and output.
#'
#' @return A Seurat object list containing metadata and VDJ annotations.
#' @concept Annotation & quality control
#' @export
CombineData <- function(
  gex.path,
  marker.gene,
  vdj.path=c(),
  threshold.cutoff=.975,
  file.saved="samples.rds",
  index.control=c(-1),
  sample.name=c(-1),
  preset=1,
  threshold.entry=0,
  demultiplex=FALSE,
  demultiplex.index=c(),
  nameshashtags=c(),
  n.ht.per.sample=1,
  tenX_conversion="true",
  nFeature_RNA_lower=100,
  nFeature_RNA_upper=10000,
  nvariable_features=3000,
  percent.mt_upper=5,
  data.input="raw",
  upperQ=.95,
  lowerQ=.05,
  verbose=TRUE,
  QC_plots=FALSE
){

  #Set undefined parameters
  if(index.control[1]==-1){
    if(demultiplex){
      index.control=rep(1,n.ht.per.sample)
    }
    else{
    index.control=rep(1,length(gex.path))
    }
  }
  if(sample.name[1]==-1){
    for(j in 1:length(gex.path)){
      lengthpath=nchar(gex.path[j])
      this.name=''
      k=1
        temp=substr(gex.path[j],(lengthpath-(k-1)),(lengthpath-(k-1)))
        while((temp!="/")&&(k<lengthpath)){
          temp=substr(gex.path[j],(lengthpath-(k-1)),(lengthpath-(k-1)))
          sample.name=as.character(1:length(gex.path))
          if(k==1){
            this.name=temp
          }
          else{
            this.name=paste(this.name,temp,sep="")
          }
          k=k+1
        }
      
      newtemp=''
      for(k in 1:nchar(temp)){
        if(k==1){
          newtemp=substr(temp,(nchar(temp)-(k-1)),(nchar(temp)-(k-1)))
        }else{
          newtemp=paste(newtemp,substr(temp,(nchar(temp)-(k-1)),(nchar(temp)-(k-1))),sep="")
        }
      }
      sample.name[j]=newtemp
    }
  }
  
    n.samples=length(gex.path)
    myvec1=list()
  
  Clonal_Obs= vector(mode = "list", length = n.samples)
  VDJ_Obs= vector(mode = "list", length = n.samples)
  tcrhash= vector(mode = "list", length = n.samples)
  cutoff= vector(mode = "list", length = n.samples)
  
  if(packageVersion("Seurat")<'5.0.0'){
    
    if(demultiplex){
      n.hashtag.samples=length(gex.path)
      n.samples=n.ht.per.sample*n.hashtag.samples
      hashtagdata= vector(mode = "list", length = n.hashtag.samples)
      for(q in 1:n.hashtag.samples){
        clonal.data <- Read10X(data.dir = gex.path[[q]])
        hashtagdata[[q]] <- GMM_demux(sample.name[((q-1)*n.ht.per.sample)+1],clonal.data, demultiplex.index[(((q-1)*n.ht.per.sample)+1):((q)*n.ht.per.sample)],nameshashtags[(((q-1)*n.ht.per.sample)+1):((q)*n.ht.per.sample)])
        
      }
      
      Clonal_Obs= vector(mode = "list", length = n.ht.per.sample*n.hashtag.samples)
      for(q in 1:n.hashtag.samples){
        for(k in 1:n.ht.per.sample){
          Idents(hashtagdata[[q]]) <- "orig.ident"
          Clonal_Obs[[(q-1)*n.ht.per.sample+k]] = subset(hashtagdata[[q]],orig.ident %in% nameshashtags[(q-1)*n.ht.per.sample+k])
          Clonal_Obs[[(q-1)*n.ht.per.sample+k]][["percent.mt"]] <- PercentageFeatureSet(Clonal_Obs[[(q-1)*n.ht.per.sample+k]], pattern = "^MT-")
          if(QC_plots){
            if(data.input=="raw"){
              Clonal_Obs[[k]]<- subset(Clonal_Obs[[k]], nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < percent.mt_upper)
            }
            nFeature_RNA_upperQ=quantile(Clonal_Obs[[(q-1)*n.ht.per.sample+k]]@meta.data$nFeature_RNA,upperQ)
            nFeature_RNA_lowerQ=quantile(Clonal_Obs[[(q-1)*n.ht.per.sample+k]]@meta.data$nFeature_RNA,lowerQ)
            nCount_RNA_upperQ=quantile(Clonal_Obs[[(q-1)*n.ht.per.sample+k]]@meta.data$nCount_RNA,upperQ)
            nCount_RNA_lowerQ=quantile(Clonal_Obs[[(q-1)*n.ht.per.sample+k]]@meta.data$nCount_RNA,lowerQ)
            percent.mt_upperQ=quantile(Clonal_Obs[[(q-1)*n.ht.per.sample+k]]@meta.data$percent.mt,upperQ)
            print(QCPlot(Clonal_Obs[[(q-1)*n.ht.per.sample+k]], nFeature_RNA_lower=nFeature_RNA_lowerQ,
                         nFeature_RNA_upper=nFeature_RNA_upperQ,
                         nCount_RNA_lower=nCount_RNA_lowerQ,
                         nCount_RNA_upper=nCount_RNA_upperQ,
                         percent.mt_upper=percent.mt_upperQ))
            }
          Clonal_Obs[[(q-1)*n.ht.per.sample+k]]<- subset(Clonal_Obs[[(q-1)*n.ht.per.sample+k]], nFeature_RNA > nFeature_RNA_lowerQ & nFeature_RNA < nFeature_RNA_upperQ & percent.mt < percent.mt_upperQ)
        }
      }
      rm(clonal.data)
    }else{
      for(k in 1:n.samples){
        print(paste(paste("Reading in gene expression for sample ",k,sep=""),".",sep=""))
        clonal.data <- Read10X(data.dir = gex.path[[k]])
        Clonal_Obs[[k]] <- CreateSeuratObject(counts = clonal.data,project = sample.name[k])
        Clonal_Obs[[k]][["percent.mt"]] <- PercentageFeatureSet(Clonal_Obs[[k]], pattern = "^MT-")
        if(QC_plots){
          if(data.input=="raw"){
          Clonal_Obs[[k]]<- subset(Clonal_Obs[[k]], nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < percent.mt_upper)
          }
          nFeature_RNA_upperQ=quantile(Clonal_Obs[[k]]@meta.data$nFeature_RNA,upperQ)
          nFeature_RNA_lowerQ=quantile(Clonal_Obs[[k]]@meta.data$nFeature_RNA,lowerQ)
          nCount_RNA_upperQ=quantile(Clonal_Obs[[k]]@meta.data$nCount_RNA,upperQ)
          nCount_RNA_lowerQ=quantile(Clonal_Obs[[k]]@meta.data$nCount_RNA,lowerQ)
          percent.mt_upperQ=quantile(Clonal_Obs[[k]]@meta.data$percent.mt,upperQ)
          print(QCPlot(Clonal_Obs[[k]], nFeature_RNA_lower=nFeature_RNA_lowerQ,
                       nFeature_RNA_upper=nFeature_RNA_upperQ,
                       nCount_RNA_lower=nCount_RNA_lowerQ,
                       nCount_RNA_upper=nCount_RNA_upperQ,
                       percent.mt_upper=percent.mt_upperQ))
          pdf(paste(gex.path[k],"/QC.pdf",sep=""),width=10,height=7)
          print(QCPlot(Clonal_Obs[[k]], nFeature_RNA_lower=nFeature_RNA_lowerQ,
                     nFeature_RNA_upper=nFeature_RNA_upperQ,
                     nCount_RNA_lower=nCount_RNA_lowerQ,
                     nCount_RNA_upper=nCount_RNA_upperQ,
                     percent.mt_upper=percent.mt_upperQ))
              dev.off()
          }
        Clonal_Obs[[k]]<- subset(Clonal_Obs[[k]], nFeature_RNA > nFeature_RNA_lowerQ & nFeature_RNA < nFeature_RNA_upperQ & percent.mt < percent.mt_upperQ)
        
      }
      rm(clonal.data)
    }
    
    for(k in 1:n.samples){
      RnaStoreUMO=Clonal_Obs[[k]]@assays$RNA@counts
      
      #VDJ list with barcodes
      if(length(vdj.path)>0){
        n.hashtag.samples=length(gex.path)
        print(paste(paste("Reading in VDJ for sample ",k,sep=""),".",sep=""))
        if(demultiplex){
          if(k==1){
            vdj.path.temp=vdj.path
            vdj.path=0
            for(q in 1:n.hashtag.samples){
              vdj.path=append(vdj.path,rep(vdj.path.temp[q],n.ht.per.sample))
            }
            vdj.path=vdj.path[-1]
          }
          VDJ_Obs[[k]] <- read.csv(vdj.path[[k]])
        }
        else{
          VDJ_Obs[[k]] <- read.csv(vdj.path[[k]])
        }
        
        tcrUMO_cell=subset(VDJ_Obs[[k]],productive==tenX_conversion & is_cell==tenX_conversion)
        TCRlistUMO=tcrUMO_cell$barcode
        joint.bcsUMO <- intersect(colnames(RnaStoreUMO), TCRlistUMO)
        #Remove barcodes from list that don't match VDJ
        tcrhash[[k]]=subset(VDJ_Obs[[k]],VDJ_Obs[[k]]$barcode %in% joint.bcsUMO)
        mvsts.UMO <- RnaStoreUMO[, joint.bcsUMO]
        Clonal_Obs[[k]]=Clonal_Obs[[k]][,joint.bcsUMO]
      }
      else{
        mvsts.UMO<-RnaStoreUMO
      }
      ####Genes of interest
      Gene_indUMO=match(marker.gene,row.names(mvsts.UMO))
      CD8_UMO=match("CD8A",row.names(mvsts.UMO))
      CD4_UMO=match("CD4",row.names(mvsts.UMO))
      CD8cells=rep(0,length(mvsts.UMO[Gene_indUMO[1],]))
      CD4cells=rep(0,length(mvsts.UMO[Gene_indUMO[1],]))
      
      vec1=mvsts.UMO[CD4_UMO,]
      vec2=mvsts.UMO[CD8_UMO,]
      
      for(g in 1:length(CD8cells)){
        if((vec2[g]==0)&(vec1[g]>0)){
          CD4cells[g]=1
        }
        else if((vec2[g]>0)&(vec1[g]==0)){
          CD8cells[g]=1
        }
        else{}
      }
      rm(vec1)
      rm(vec2)
      
      Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], CD8cells, col.name = 'CD8cells')
      Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], CD4cells, col.name = 'CD4cells')
      Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], rep(sample.name[k],length(CD4cells)), col.name = 'sampleref')
      
      rm(CD8cells)
      rm(CD4cells)
      #Check for errors - should you wish to plot
      #plot(log10(mvsts.UMO[Gene_indUMO[1],]),log10(mvsts.UMO[Gene_indUMO[2],]),xlab=paste(Gene_indUMO[3]," (log10 UMIs)",sep=""),ylab=paste(Gene_indUMO[2]," (log10 UMIs)",sep=""),col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))
      #plot(log10(mvsts.UMO[CD8_UMO,]),log10(mvsts.UMO[CD4_UMO,]),xlab="CD8",ylab="CD4",col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))
      print(paste(paste("Starting scTRansform for sample ",k,sep=""),"...",sep=""))
      Clonal_Obs[[k]] <- SCTransform(Clonal_Obs[[k]], vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n = nvariable_features)
      Gene_indUMO=match(marker.gene,row.names(Clonal_Obs[[k]][["SCT"]]@data))
      Gene_indUMORNA=match(marker.gene,row.names(Clonal_Obs[[k]][['RNA']]@counts))
      
      number.marker.genes=length(marker.gene)
      
      if(k %in% index.control){
        mask=which(is.na(Gene_indUMO)==TRUE)
      }
      
      #Set genes not present to be first gene (this value isn't used)
      Gene_indUMO[mask]=1
   
      ###################################################################################
      ### Determine if all genes are in sctransform data                              ###
      ###################################################################################
      if((k %in% index.control)&&(length(mask)>0)){
        for(s in 1:length(mask)){
          print(paste(paste("gene: ",marker.gene[mask[s]],sep=" ")," is not present in SCT",sep=""))
          print(paste(paste("Continuing with remaining gois and replacing cutoff for gene: ",marker.gene[mask[s]],sep=""), " with count data...",sep="")
          )
        }
      }
      ###################################################################################
      ### Create matrices d1 & d2 which are the SCT & RNA values of the marker genes  ###
      ###################################################################################
      d1=as.matrix(Clonal_Obs[[k]][["SCT"]]@data)[Gene_indUMO,]
      d2=as.matrix(Clonal_Obs[[k]][['RNA']]@counts)[Gene_indUMORNA,]
      alist <- 1:number.marker.genes
      
      ###################################################################################
      ### Start with default threshold setting 1, that assumes inclusion of a control ###
      ###################################################################################
      if((k %in% index.control)&&(preset==1)){
        if(length(alist)==1){
          cutoff[[k]]=quantile(unname(d1),threshold.cutoff)[[1]]
        }
        else if(length(alist)>1){
          myvec1=0
          if(length(mask)>0){
            for(j in 1:number.marker.genes){
              if(j %in% mask){
            myvec1=append(myvec1,quantile(unname(d2[j,]),threshold.cutoff)[[1]])
              }
              else{
                myvec1=append(myvec1,quantile(unname(d1[j,]),threshold.cutoff)[[1]])
              }
            }
            myvec1=myvec1[-1]
            cutoff[[k]]=myvec1
          }
          else{
            cutoff[[k]]=lapply(alist, function(alist) quantile(unname(d1),threshold.cutoff)[[1]])
            
          }}
      }else if(preset==0){
        ###################################################################################
        ### Preset = 0 sets a global value to define high expression of a gene          ###
        ###################################################################################
        cutoff[[k]]=rep(threshold.entry,number.marker.genes)
      }else if(preset==2){
        ###################################################################################
        ### Preset = 2 sets a threshold for each gene within each sample                ###
        ###################################################################################
        if(length(alist)==1){
          cutoff[[k]]=quantile(unname(d1),threshold.cutoff)[[1]]
        }
        else if(length(alist)>1){
          myvec1=0
          if(length(mask)>0){
            for(j in 1:number.marker.genes){
              if(j %in% mask){
                myvec1=append(myvec1,quantile(unname(d2[j,]),threshold.cutoff)[[1]])
              }
              else{
                myvec1=append(myvec1,quantile(unname(d1[j,]),threshold.cutoff)[[1]])
              }
            }
            myvec1=myvec1[-1]
            cutoff[[k]]=myvec1
          }
          else{
            cutoff[[k]]=lapply(alist, function(alist) quantile(unname(d1),threshold.cutoff)[[1]])
            
          }}
        }else{
        print(paste(paste("Sample: ",k,sep="")," not used to set the marker gene thresholds"))
      }
      
      ###################################################################################
      ### Create a threshold matrix that will be used for the marker gene metadata    ###
      ###################################################################################
      Thresholds=matrix("unassigned",nrow=dim(Clonal_Obs[[k]])[2],ncol=length(Gene_indUMO))
      
      if(length(mask)>0){
        for(s in 1:number.marker.genes){
          if(s %in% mask){
            vec1=as.matrix(Clonal_Obs[[k]][['RNA']]@counts)[Gene_indUMORNA[s],]
            for(j in 1:dim(Clonal_Obs[[k]])[2]){
              if(vec1[j]>cutoff[[index.control[k]]][s]){
                Thresholds[j,s]="high"
              }
            }
          }
          else{
          vec1=as.matrix(Clonal_Obs[[k]][["SCT"]]@data)[Gene_indUMO[s],]
          for(j in 1:dim(Clonal_Obs[[k]])[2]){
            if(vec1[j]>cutoff[[index.control[k]]][s]){
              
              Thresholds[j,s]="high"
            }
          }
          }
          Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], Thresholds[,s], col.name = paste("Threshold_",marker.gene[s],sep=""))
        }
      }
      else{
        for(s in 1:number.marker.genes){
          vec1=unname(as.matrix(Clonal_Obs[[k]][["SCT"]]@data)[Gene_indUMO[s],])
          for(j in 1:dim(Clonal_Obs[[k]])[2]){
            if(vec1[j]>cutoff[[index.control[k]]][s]){
              Thresholds[j,s]="high"
            }
          }
         Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], Thresholds[,s], col.name = paste("Threshold_",marker.gene[s],sep=""))
        }
      }
      
      rm(mvsts.UMO)
      rm(RnaStoreUMO)
      
    }
    
    ###############################################################
    ##-----------------------------------------------------------##
    ##         Annotate cells with clonotypes (TCRB)             ##
    ##-----------------------------------------------------------##
    ###############################################################
    if(length(vdj.path)>0){
      monoTRBNA=vector(mode = "list", length = n.samples)
      monoTRBAA=vector(mode = "list", length = n.samples)
      for(k in 1:n.samples){
        monoTRBNA[[k]]=rep("unassigned",dim(Clonal_Obs[[k]])[2])
        monoTRBAA[[k]]=rep("unassigned",dim(Clonal_Obs[[k]])[2])
        
        ##TRB loop for each single cell
        for(j in 1:dim(Clonal_Obs[[k]])[2]){
          
          set_umo=subset(tcrhash[[k]],barcode==colnames(Clonal_Obs[[k]])[j] & productive==tenX_conversion & chain=="TRB")
          TCRs=set_umo
          len=length(set_umo$reads)
          if(len==0){
            
          }
          else if(len==1){
            monoTRBNA[[k]][j]=set_umo$cdr3_nt
            monoTRBAA[[k]][j]=set_umo$cdr3
          }
          else{
            position= which.max(set_umo$reads)
            monoTRBNA[[k]][j]=set_umo$cdr3_nt[position]
            monoTRBAA[[k]][j]=set_umo$cdr3[position]
          }
          
        }
        
      }
      
      IntersectTCRs=monoTRBNA[[1]]
      if(n.samples>1){
        for(h in 2:n.samples){
          IntersectTCRs = intersect(IntersectTCRs,monoTRBNA[[h]])
        }
      }
      
      SizesC_Obj= vector(mode = "list", length = n.samples)
      IMTCRs= vector(mode = "list", length = n.samples)
      all.Cname= vector(mode = "list", length = n.samples)
      all.Size= vector(mode = "list", length = n.samples)
      for(q in 1:n.samples)
      {
        all.Cname[[q]]=rep("unassigned",dim(Clonal_Obs[[q]])[2])
        all.Size[[q]]=rep(0,dim(Clonal_Obs[[q]])[2])
      }
      for(h in 1:n.samples){
        Tailmono=setdiff(monoTRBNA[[h]],IntersectTCRs)
        IMTCRs[[h]]=c(IntersectTCRs,Tailmono)
        
        SizesC_Obj[[h]]= rep(0,length(IMTCRs))
        for(j in 1:length(IMTCRs[[h]])){
          SizesC_Obj[[h]][j]=length(subset(tcrhash[[h]],cdr3_nt==IMTCRs[[h]][j] & productive==tenX_conversion & chain=="TRB")$barcode)
        }
        
      }
      
      
      Nclono= vector(mode = "list", length = n.samples)
      for(q in 1:n.samples){
        OSM=order(SizesC_Obj[[q]],decreasing=TRUE)
        Nclono[[q]]=rep("blank",length(IMTCRs[[q]]))
        for(f in 1: length(IMTCRs[[q]])){
          if(f>length(IntersectTCRs)){
            Nclono[[q]][OSM[f]]=paste(paste(paste("clonotype",f,sep=""),"_",sep=""),q,sep="")
          }else{
            Nclono[[q]][OSM[f]]=paste("clonotype",f,sep="")
          }
          
        }
      }
      for(k in 1:n.samples){
        clonotypes.nucleic.acid=data.frame(TCRs=IMTCRs[[k]],clonoT=Nclono[[k]],freqT=SizesC_Obj[[k]])
        for(j in 1:dim(Clonal_Obs[[k]])[2]){
          
          if(monoTRBNA[[k]][j]=="unassigned"){
            
          }
          else{
            all.Size[[k]][j]=subset(clonotypes.nucleic.acid,TCRs==monoTRBNA[[k]][j])$freqT
            all.Cname[[k]][j]=subset(clonotypes.nucleic.acid,TCRs==monoTRBNA[[k]][j])$clonoT
          }
          
        }
        
        Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], all.Cname[[k]], col.name = 'clonotype')
        Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], all.Size[[k]], col.name = 'countcln')
      }
      
      #####################################################
      ## Now add nucleotide CDR3 and Amino Acid CDR3     ##
      #####################################################
      
      clono.nucleotide.seq= vector(mode = "list", length = n.samples)
      clono.aminoacid.seq= vector(mode = "list", length = n.samples)
      for(k in 1:n.samples){
        clono.nucleotide.seq[[k]]=rep("unassigned",dim(Clonal_Obs[[k]])[2])
        clono.aminoacid.seq[[k]]=rep("unassigned",dim(Clonal_Obs[[k]])[2])
        
        ##TRB loop
        for(j in 1:dim(Clonal_Obs[[k]])[2]){
          
          set_umo=subset(tcrhash[[k]],barcode==colnames(Clonal_Obs[[k]])[j] & productive==tenX_conversion & chain=="TRB")
          TCRs=set_umo
          len=length(set_umo$reads)
          if(len==0){
            
          }
          else if(len==1){
            mTRBNA=set_umo$cdr3_nt
            mTRBAA=set_umo$cdr3
            clono.nucleotide.seq[[k]][j]=mTRBNA
            clono.aminoacid.seq[[k]][j]=mTRBAA
          }
          else{
            position= which.max(set_umo$reads)
            mTRBNA=set_umo$cdr3_nt[position]
            mTRBAA=set_umo$cdr3[position]
            clono.nucleotide.seq[[k]][j]=mTRBNA
            clono.aminoacid.seq[[k]][j]=mTRBAA
          }
          
        }
        
        Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], clono.nucleotide.seq[[k]], col.name = 'cdr3_na')
        Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], clono.aminoacid.seq[[k]], col.name = 'cdr3')
      }
    }
    Dataset=Clonal_Obs
    saveRDS(Dataset,file.saved)
    return(Dataset)
  }
  else{

if(demultiplex){
  n.hashtag.samples=length(gex.path)
  n.samples=n.ht.per.sample*n.hashtag.samples
  hashtagdata= vector(mode = "list", length = n.hashtag.samples)
  for(q in 1:n.hashtag.samples){
    clonal.data <- Read10X(data.dir = gex.path[[q]])
hashtagdata[[q]] <- GMM_demux(sample.name[((q-1)*n.ht.per.sample)+1],clonal.data, demultiplex.index[(((q-1)*n.ht.per.sample)+1):((q)*n.ht.per.sample)],nameshashtags[(((q-1)*n.ht.per.sample)+1):((q)*n.ht.per.sample)])

}

Clonal_Obs= vector(mode = "list", length = n.ht.per.sample*n.hashtag.samples)
for(q in 1:n.hashtag.samples){
for(k in 1:n.ht.per.sample){
  Idents(hashtagdata[[q]]) <- "orig.ident"
  
Clonal_Obs[[(q-1)*n.ht.per.sample+k]] = subset(hashtagdata[[q]],orig.ident %in% nameshashtags[(q-1)*n.ht.per.sample+k])
Clonal_Obs[[(q-1)*n.ht.per.sample+k]][["percent.mt"]] <- PercentageFeatureSet(Clonal_Obs[[(q-1)*n.ht.per.sample+k]], pattern = "^MT-")
if(QC_plots){
  if(data.input=="raw"){
    Clonal_Obs[[k]]<- subset(Clonal_Obs[[k]], nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < percent.mt_upper)
  }
  nFeature_RNA_upperQ=quantile(Clonal_Obs[[(q-1)*n.ht.per.sample+k]]@meta.data$nFeature_RNA,upperQ)
  nFeature_RNA_lowerQ=quantile(Clonal_Obs[[(q-1)*n.ht.per.sample+k]]@meta.data$nFeature_RNA,lowerQ)
  nCount_RNA_upperQ=quantile(Clonal_Obs[[(q-1)*n.ht.per.sample+k]]@meta.data$nCount_RNA,upperQ)
  nCount_RNA_lowerQ=quantile(Clonal_Obs[[(q-1)*n.ht.per.sample+k]]@meta.data$nCount_RNA,lowerQ)
  percent.mt_upperQ=quantile(Clonal_Obs[[(q-1)*n.ht.per.sample+k]]@meta.data$percent.mt,upperQ)
  print(QCPlot(Clonal_Obs[[(q-1)*n.ht.per.sample+k]], nFeature_RNA_lower=nFeature_RNA_lowerQ,
               nFeature_RNA_upper=nFeature_RNA_upperQ,
               nCount_RNA_lower=nCount_RNA_lowerQ,
               nCount_RNA_upper=nCount_RNA_upperQ,
               percent.mt_upper=percent.mt_upperQ))
}
Clonal_Obs[[(q-1)*n.ht.per.sample+k]]<- subset(Clonal_Obs[[(q-1)*n.ht.per.sample+k]], subset = nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < percent.mt_upper)
}
}
rm(clonal.data)
}else{
  for(k in 1:n.samples){
    print(paste(paste("Reading in gene expression for sample ",k,sep=""),".",sep=""))
    clonal.data <- Read10X(data.dir = gex.path[[k]])
    clonal.data
  Clonal_Obs[[k]] <- CreateSeuratObject(counts = clonal.data,project = sample.name[k])
  Clonal_Obs[[k]][["percent.mt"]] <- PercentageFeatureSet(Clonal_Obs[[k]], pattern = "^MT-")
  if(QC_plots){
    if(data.input=="raw"){
      Clonal_Obs[[k]]<- subset(Clonal_Obs[[k]], nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < percent.mt_upper)
    }
    nFeature_RNA_upperQ=quantile(Clonal_Obs[[k]]@meta.data$nFeature_RNA,upperQ)
    nFeature_RNA_lowerQ=quantile(Clonal_Obs[[k]]@meta.data$nFeature_RNA,lowerQ)
    nCount_RNA_upperQ=quantile(Clonal_Obs[[k]]@meta.data$nCount_RNA,upperQ)
    nCount_RNA_lowerQ=quantile(Clonal_Obs[[k]]@meta.data$nCount_RNA,lowerQ)
    percent.mt_upperQ=quantile(Clonal_Obs[[k]]@meta.data$percent.mt,upperQ)
    print(QCPlot(Clonal_Obs[[k]], nFeature_RNA_lower=nFeature_RNA_lowerQ,
                 nFeature_RNA_upper=nFeature_RNA_upperQ,
                 nCount_RNA_lower=nCount_RNA_lowerQ,
                 nCount_RNA_upper=nCount_RNA_upperQ,
                 percent.mt_upper=percent.mt_upperQ))
    pdf(paste(gex.path[k],"/QC.pdf",sep=""),width=10,height=7)
    print(QCPlot(Clonal_Obs[[k]], nFeature_RNA_lower=nFeature_RNA_lowerQ,
                 nFeature_RNA_upper=nFeature_RNA_upperQ,
                 nCount_RNA_lower=nCount_RNA_lowerQ,
                 nCount_RNA_upper=nCount_RNA_upperQ,
                 percent.mt_upper=percent.mt_upperQ))
    dev.off()
  }
  Clonal_Obs[[k]]<- subset(Clonal_Obs[[k]], subset = nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < percent.mt_upper)

}
rm(clonal.data)
}
mask=0
for(k in 1:n.samples){
RnaStoreUMO=Clonal_Obs[[k]][['RNA']]$counts

#VDJ list with barcodes
if(length(vdj.path)>0){
  print(paste(paste("Reading in VDJ for sample ",k,sep=""),".",sep=""))
  n.hashtag.samples=length(gex.path)
  if(demultiplex){
  if(k==1){
  vdj.path.temp=vdj.path
  vdj.path=0
  for(q in 1:n.hashtag.samples){
    vdj.path=append(vdj.path,rep(vdj.path.temp[q],n.ht.per.sample))
  }
  vdj.path=vdj.path[-1]
  }
  VDJ_Obs[[k]] <- read.csv(vdj.path[[k]])
}
  else{
VDJ_Obs[[k]] <- read.csv(vdj.path[[k]])
}

tcrUMO_cell=subset(VDJ_Obs[[k]],productive==tenX_conversion & is_cell==tenX_conversion)
TCRlistUMO=tcrUMO_cell$barcode
joint.bcsUMO <- intersect(colnames(RnaStoreUMO), TCRlistUMO)
#Remove barcodes from list that don't match VDJ
tcrhash[[k]]=subset(VDJ_Obs[[k]],VDJ_Obs[[k]]$barcode %in% joint.bcsUMO)
mvsts.UMO <- RnaStoreUMO[, joint.bcsUMO]
Clonal_Obs[[k]]=Clonal_Obs[[k]][,joint.bcsUMO]
}
else{
  mvsts.UMO<-RnaStoreUMO
}
####Genes of interest
Gene_indUMO=match(marker.gene,row.names(mvsts.UMO))
CD8_UMO=match("CD8A",row.names(mvsts.UMO))
CD4_UMO=match("CD4",row.names(mvsts.UMO))
CD8cells=rep(0,length(mvsts.UMO[Gene_indUMO[1],]))
CD4cells=rep(0,length(mvsts.UMO[Gene_indUMO[1],]))

vec1=mvsts.UMO[CD4_UMO,]
vec2=mvsts.UMO[CD8_UMO,]

for(g in 1:length(CD8cells)){
  if((vec2[g]==0)&(vec1[g]>0)){
    CD4cells[g]=1
  }
  else if((vec2[g]>0)&(vec1[g]==0)){
    CD8cells[g]=1
  }
  else{}
}
rm(vec1)
rm(vec2)

Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], CD8cells, col.name = 'CD8cells')
Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], CD4cells, col.name = 'CD4cells')
Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], rep(sample.name[k],length(CD4cells)), col.name = 'sampleref')

rm(CD8cells)
rm(CD4cells)
#Check for errors - should you wish to plot
#plot(log10(mvsts.UMO[Gene_indUMO[1],]),log10(mvsts.UMO[Gene_indUMO[2],]),xlab=paste(Gene_indUMO[3]," (log10 UMIs)",sep=""),ylab=paste(Gene_indUMO[2]," (log10 UMIs)",sep=""),col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))
#plot(log10(mvsts.UMO[CD8_UMO,]),log10(mvsts.UMO[CD4_UMO,]),xlab="CD8",ylab="CD4",col=rgb(red=0, green = 0, blue = 0, alpha=0.5),pch=16,cex=1,xlim=c(0,4.5),ylim=c(0,4.5))
print(paste(paste("Starting scTRansform for sample ",k,sep=""),"...",sep=""))
Clonal_Obs[[k]] <- SCTransform(Clonal_Obs[[k]], vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n = nvariable_features)
Gene_indUMO=match(marker.gene,row.names(Clonal_Obs[[k]][["SCT"]]$data))
Gene_indUMORNA=match(marker.gene,row.names(Clonal_Obs[[k]][['RNA']]$counts))

number.marker.genes=length(marker.gene)

if(k %in% index.control){
  mask=which(is.na(Gene_indUMO)==TRUE)
}

#Set genes not present to be first gene (this value isn't used)
Gene_indUMO[mask]=1

###################################################################################
### Determine if all genes are in sctransform data                              ###
###################################################################################
if((k %in% index.control)&&(length(mask)>0)){
  for(s in 1:length(mask)){
    print(paste(paste("gene: ",marker.gene[mask[s]],sep=" ")," is not present in SCT",sep=""))
    print(paste(paste("Continuing with remaining gois and replacing cutoff for gene: ",marker.gene[mask[s]],sep=""), " with count data...",sep="")
    )
  }
}
###################################################################################
### Create matrices d1 & d2 which are the SCT & RNA values of the marker genes  ###
###################################################################################
d1=as.matrix(Clonal_Obs[[k]][["SCT"]]$data)[Gene_indUMO,]
d2=as.matrix(Clonal_Obs[[k]][['RNA']]$counts)[Gene_indUMORNA,]
alist <- 1:number.marker.genes

###################################################################################
### Start with default threshold setting 1, that assumes inclusion of a control ###
###################################################################################
if((k %in% index.control)&&(preset==1)){
  if(length(alist)==1){
    cutoff[[k]]=quantile(unname(d1),threshold.cutoff)[[1]]
  }
  else if(length(alist)>1){
    myvec1=0
    if(length(mask)>0){
      for(j in 1:number.marker.genes){
        if(j %in% mask){
          myvec1=append(myvec1,quantile(unname(d2[j,]),threshold.cutoff)[[1]])
        }
        else{
          myvec1=append(myvec1,quantile(unname(d1[j,]),threshold.cutoff)[[1]])
        }
      }
      myvec1=myvec1[-1]
      cutoff[[k]]=myvec1
    }
    else{
      cutoff[[k]]=lapply(alist, function(alist) quantile(unname(d1),threshold.cutoff)[[1]])
      
    }}
}else if(preset==0){
  ###################################################################################
  ### Preset = 0 sets a global value to define high expression of a gene          ###
  ###################################################################################
  cutoff[[k]]=rep(threshold.entry,number.marker.genes)
}else if(preset==2){
  ###################################################################################
  ### Preset = 2 sets a threshold for each gene within each sample                ###
  ###################################################################################
  if(length(alist)==1){
    cutoff[[k]]=quantile(unname(d1),threshold.cutoff)[[1]]
  }
  else if(length(alist)>1){
    myvec1=0
    if(length(mask)>0){
      for(j in 1:number.marker.genes){
        if(j %in% mask){
          myvec1=append(myvec1,quantile(unname(d2[j,]),threshold.cutoff)[[1]])
        }
        else{
          myvec1=append(myvec1,quantile(unname(d1[j,]),threshold.cutoff)[[1]])
        }
      }
      myvec1=myvec1[-1]
      cutoff[[k]]=myvec1
    }
    else{
      cutoff[[k]]=lapply(alist, function(alist) quantile(unname(d1),threshold.cutoff)[[1]])
      
    }}
}else{
  print(paste(paste("Sample: ",k,sep="")," not used to set the marker gene thresholds"))
}

###################################################################################
### Create a threshold matrix that will be used for the marker gene metadata    ###
###################################################################################
Thresholds=matrix("unassigned",nrow=dim(Clonal_Obs[[k]])[2],ncol=length(Gene_indUMO))

if(length(mask)>0){
  for(s in 1:number.marker.genes){
    if(s %in% mask){
      vec1=as.matrix(Clonal_Obs[[k]][['RNA']]$counts)[Gene_indUMORNA[s],]
      for(j in 1:dim(Clonal_Obs[[k]])[2]){
        if(vec1[j]>cutoff[[index.control[k]]][s]){
          Thresholds[j,s]="high"
        }
      }
    }
    else{
      vec1=as.matrix(Clonal_Obs[[k]][["SCT"]]$data)[Gene_indUMO[s],]
      for(j in 1:dim(Clonal_Obs[[k]])[2]){
        if(vec1[j]>cutoff[[index.control[k]]][s]){
          
          Thresholds[j,s]="high"
        }
      }
    }
    Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], Thresholds[,s], col.name = paste("Threshold_",marker.gene[s],sep=""))
  }
}
else{
  for(s in 1:number.marker.genes){
    vec1=unname(as.matrix(Clonal_Obs[[k]][["SCT"]]$data)[Gene_indUMO[s],])
    for(j in 1:dim(Clonal_Obs[[k]])[2]){
      if(vec1[j]>cutoff[[index.control[k]]][s]){
        Thresholds[j,s]="high"
      }
    }
    Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], Thresholds[,s], col.name = paste("Threshold_",marker.gene[s],sep=""))
  }
}

rm(mvsts.UMO)
rm(RnaStoreUMO)

}

###############################################################
##-----------------------------------------------------------##
##         Annotate cells with clonotypes (TCRB)             ##
##-----------------------------------------------------------##
###############################################################
if(length(vdj.path)>0){
monoTRBNA=vector(mode = "list", length = n.samples)
monoTRBAA=vector(mode = "list", length = n.samples)
for(k in 1:n.samples){
monoTRBNA[[k]]=rep("unassigned",dim(Clonal_Obs[[k]])[2])
monoTRBAA[[k]]=rep("unassigned",dim(Clonal_Obs[[k]])[2])

##TRB loop for each single cell
for(j in 1:dim(Clonal_Obs[[k]])[2]){

  set_umo=subset(tcrhash[[k]],barcode==colnames(Clonal_Obs[[k]])[j] & productive==tenX_conversion & chain=="TRB")
  TCRs=set_umo
  len=length(set_umo$reads)
  if(len==0){

  }
  else if(len==1){
    monoTRBNA[[k]][j]=set_umo$cdr3_nt
    monoTRBAA[[k]][j]=set_umo$cdr3
  }
  else{
    position= which.max(set_umo$reads)
    monoTRBNA[[k]][j]=set_umo$cdr3_nt[position]
    monoTRBAA[[k]][j]=set_umo$cdr3[position]
  }

}

}

IntersectTCRs=monoTRBNA[[1]]
if(n.samples>1){
for(h in 2:n.samples){
  IntersectTCRs = intersect(IntersectTCRs,monoTRBNA[[h]])
}
}

SizesC_Obj= vector(mode = "list", length = n.samples)
IMTCRs= vector(mode = "list", length = n.samples)
all.Cname= vector(mode = "list", length = n.samples)
all.Size= vector(mode = "list", length = n.samples)
for(q in 1:n.samples)
{
  all.Cname[[q]]=rep("unassigned",dim(Clonal_Obs[[q]])[2])
  all.Size[[q]]=rep(0,dim(Clonal_Obs[[q]])[2])
}
for(h in 1:n.samples){
Tailmono=setdiff(monoTRBNA[[h]],IntersectTCRs)
IMTCRs[[h]]=c(IntersectTCRs,Tailmono)

  SizesC_Obj[[h]]= rep(0,length(IMTCRs))
  for(j in 1:length(IMTCRs[[h]])){
    SizesC_Obj[[h]][j]=length(subset(tcrhash[[h]],cdr3_nt==IMTCRs[[h]][j] & productive==tenX_conversion & chain=="TRB")$barcode)
  }

}


Nclono= vector(mode = "list", length = n.samples)
for(q in 1:n.samples){
  OSM=order(SizesC_Obj[[q]],decreasing=TRUE)
  Nclono[[q]]=rep("blank",length(IMTCRs[[q]]))
for(f in 1: length(IMTCRs[[q]])){
  if(f>length(IntersectTCRs)){
    Nclono[[q]][OSM[f]]=paste(paste(paste("clonotype",f,sep=""),"_",sep=""),q,sep="")
  }else{
    Nclono[[q]][OSM[f]]=paste("clonotype",f,sep="")
  }

}
}
for(k in 1:n.samples){
clonotypes.nucleic.acid=data.frame(TCRs=IMTCRs[[k]],clonoT=Nclono[[k]],freqT=SizesC_Obj[[k]])
for(j in 1:dim(Clonal_Obs[[k]])[2]){

  if(monoTRBNA[[k]][j]=="unassigned"){

  }
  else{
    all.Size[[k]][j]=subset(clonotypes.nucleic.acid,TCRs==monoTRBNA[[k]][j])$freqT
    all.Cname[[k]][j]=subset(clonotypes.nucleic.acid,TCRs==monoTRBNA[[k]][j])$clonoT
  }

}

Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], all.Cname[[k]], col.name = 'clonotype')
Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], all.Size[[k]], col.name = 'countcln')
}

#####################################################
## Now add nucleotide CDR3 and Amino Acid CDR3     ##
#####################################################

clono.nucleotide.seq= vector(mode = "list", length = n.samples)
clono.aminoacid.seq= vector(mode = "list", length = n.samples)
for(k in 1:n.samples){
clono.nucleotide.seq[[k]]=rep("unassigned",dim(Clonal_Obs[[k]])[2])
clono.aminoacid.seq[[k]]=rep("unassigned",dim(Clonal_Obs[[k]])[2])

##TRB loop
for(j in 1:dim(Clonal_Obs[[k]])[2]){

  set_umo=subset(tcrhash[[k]],barcode==colnames(Clonal_Obs[[k]])[j] & productive==tenX_conversion & chain=="TRB")
  TCRs=set_umo
  len=length(set_umo$reads)
  if(len==0){

  }
  else if(len==1){
    mTRBNA=set_umo$cdr3_nt
    mTRBAA=set_umo$cdr3
    clono.nucleotide.seq[[k]][j]=mTRBNA
    clono.aminoacid.seq[[k]][j]=mTRBAA
  }
  else{
    position= which.max(set_umo$reads)
    mTRBNA=set_umo$cdr3_nt[position]
    mTRBAA=set_umo$cdr3[position]
    clono.nucleotide.seq[[k]][j]=mTRBNA
    clono.aminoacid.seq[[k]][j]=mTRBAA
  }

}

Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], clono.nucleotide.seq[[k]], col.name = 'cdr3_na')
Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], clono.aminoacid.seq[[k]], col.name = 'cdr3')
}
}
Dataset=Clonal_Obs
saveRDS(Dataset,file.saved)
return(Dataset)
}
}
