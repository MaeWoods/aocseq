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
#' has to each data point.
#'
#' @param sd.vector Vector containing the standard deviations of each component
#' @param sd.vector Vector containing the mean of each component
#' @param alpha.vector Vector containing the mixing weights  of each component
#' @return Named list containing the loglik and posterior.df
#' @export
e_step <- function(x, mu.vector, sd.vector, alpha.vector) {
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
#' Update the Component Parameters
#'
#' @param x Input data.
#' @param posterior.df Posterior probability data.frame.
#' @return Named list containing the mean (mu), variance (var), and mixing
#'   weights (alpha) for each component.
#' @export
m_step <- function(x, posterior.df) {
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
#' @param data Input data. 10x object with hashtags
#' @param hashtag_index list containing positions of hashtags to be demultiplexed
#' @param nameshashtags Subsets of the sample for new orig.ident
#' @param s.name Name of the sample
#' @return A Seurat object that is has
#' @export
GMM_demux<-function(s.name,data, hashtag_index, nameshashtags){

  #Total_T_Cells=length(data$`Gene Expression`@p)
  #pbmc.htos <- data$`Antibody Capture`
  Total_T_Cells=5000
  pbmc.htos=10

  temp=t(as.matrix(pbmc.htos))+1
  tempqp=temp
  mu.cont<-c()
  var.cont<-c()
  alpha.cont<-c()
  #KMeans and model fitting loop
  for (i in 1:length(hashtag_index)) {

    tempqp[,hashtag_index[i]]<-log(temp[,hashtag_index[i]]/matrixStats::product((subset(temp,temp[,hashtag_index[i]]>0)[,hashtag_index[i]])^(1/Total_T_Cells)))
    wait <- tempqp[,hashtag_index[i]]

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
        e.step <- e_step(wait, wait.summary.df1[["mu"]], wait.summary.df1[["std"]],
                         wait.summary.df1[["alpha"]])
        m.step1 <- m_step(wait, e.step[["posterior.df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      }
      else {
        # Repeat E and M steps till convergence
        e.step <- e_step(wait, m.step1[["mu"]], sqrt(m.step1[["var"]]),
                         m.step1[["alpha"]])
        m.step1 <- m_step(wait, e.step[["posterior.df"]])
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

  Mhash<-matrix(0, nrow = dim(temp)[1], ncol = length(hashtag_index)+1)
  donorlabel=rep("unassigned",length(tempqp[,1]))
  samplelabels=rep("unassigned",length(tempqp[,1]))

  #hashtag probs loop
  for(j in 1:length(tempqp[,1])){
    high.cont<-c()
    low.cont<-c()
    for(k in hashtag_index){

      if(mu.cont[match(k,hashtag_index)*2]>mu.cont[(match(k,hashtag_index)*2)-1]){
        x=tempqp[j,k]
        low=alpha.cont[match(k,hashtag_index)*2-1]
        high=alpha.cont[match(k,hashtag_index)*2]
        xgivenPz_ihigh = dnorm(x, mean = mu.cont[match(k,hashtag_index)*2], sd = var.cont[match(k,hashtag_index)*2], log = FALSE)
        xgivenPz_ilow = dnorm(x, mean = mu.cont[match(k,hashtag_index)*2-1], sd = var.cont[match(k,hashtag_index)*2-1], log = FALSE)
        Px_i = xgivenPz_ilow*low + xgivenPz_ihigh*high
        K1Pz_ihighgiven_x = (xgivenPz_ihigh*high)/Px_i
        K1Pz_ilowgiven_x = (xgivenPz_ilow*low)/Px_i

      }
      else{
        x=tempqp[j,k]
        low=alpha.cont[match(k,hashtag_index)*2]
        high=alpha.cont[match(k,hashtag_index)*2-1]
        xgivenPz_ihigh = dnorm(x, mean = mu.cont[match(k,hashtag_index)*2-1], sd = var.cont[match(k,hashtag_index)*2-1], log = FALSE)
        xgivenPz_ilow = dnorm(x, mean = mu.cont[match(k,hashtag_index)*2], sd=var.cont[match(k,hashtag_index)*2], log = FALSE)
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
      if(max(probs)==probs[i]&i<length(probs)){
        Mhash[j,i]=1
        donorlabel[j]=paste('donor_',as.character(i), sep = '')
      }
      else if (max(probs)==probs[i]&i==length(probs)){
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
  colnames(donormat) <- row.names(temp)
  joint.bcs <- intersect(colnames(data$`Gene Expression`), colnames(t(Full1)))
  pbmc.htos <- t(Full1)[, joint.bcs]
  data=CreateSeuratObject(counts = data$`Gene Expression`,project = s.name)
  data=AddMetaData(data, donorlabel, col.name = 'Hashtags')

  for (i in 1:length(nameshashtags)) {
    for(j in 1:length(tempqp[,1])){
    if(donorlabel[j]==paste('donor_',as.character(i), sep = '')){
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
#' @param gex.path path to gene expression data in the format of cellranger feature barcode matrices
#' @param vdj.path path to VDJ expression data in the format of cellranger csv files
#' @param marker.gene list of marker genes used for specificity analysis
#' @param threshold.cutoff list of marker genes used for specificity analysis
#' @param file.saved directory for storing Seurat object RDS files
#' @param index.control index of the control sequencing sample
#' @param preset If set to 1, preset determines the threshold for cutoff if a control is included in the assay
#' @param threshold.entry Double. Sets the percentile above which gene expression is labelled as high.
#' @param demultiplex Boolean value used to indicate if the data was hashtagged
#' @param hashtags Integer represents the number of hashtags used (if any)
#' @param verbose Print progress bars and output
#'
#' @return A Seurat object list containing metadata and VDJ annotations.
#' @concept annotation
#' @export
CombineData <- function(
  gex.path,
  marker.gene,
  vdj.path=c(),
  threshold.cutoff=.975,
  file.saved="samples.rds",
  index.control=1,
  n.samples=-1,
  c.index=-1,
  sample.name=-1,
  preset=1,
  threshold.entry=0,
  demultiplex=FALSE,
  n.hashtag.samples=1,
  n.samples.ht=1,
  tenX_conversion="true",
  hashtags=1,
  nFeature_RNA_lower=100,
  nFeature_RNA_upper=10000,
  nvariable_features=3000,
  percent.mt_upper=5,
  verbose=TRUE
){


  #Set undefined parameters
  if(c.index[1]==-1){
    c.index=rep(1,length(gex.path))
  }
  if(sample.name[1]==-1){
    sample.name=as.character(1:length(gex.path))
  }
  if(n.samples[1]==-1){
    n.samples=length(gex.path)
  }

Clonal_Obs= vector(mode = "list", length = n.samples)
VDJ_Obs= vector(mode = "list", length = n.samples)
tcrhash= vector(mode = "list", length = n.samples)
cutoff= vector(mode = "list", length = n.samples)

if(demultiplex){
  hashtagdata= vector(mode = "list", length = n.hashtag.samples)
  for(q in 1:n.hashtag.samples){
    clonal.data <- Read10X(data.dir = gex.path[[q]])
hashtagdata[[q]] <- GMM_demux(sample.name[((q-1)*n.samples.ht)+1],clonal.data, demultiplex.index[(((q-1)*n.samples.ht)+1):((q)*n.samples.ht)],nameshashtags[(((q-1)*n.samples.ht)+1):((q)*n.samples.ht)])

}

Clonal_Obs= vector(mode = "list", length = n.samples.ht*n.hashtag.samples)
for(q in 1:n.hashtag.samples){
for(k in 1:n.samples.ht){
  Idents(hashtagdata[[q]]) <- "orig.ident"
Clonal_Obs[[(q-1)*n.samples.ht+k]] = subset(hashtagdata[[q]],orig.ident %in% nameshashtags[(q-1)*n.samples.ht+k])
Clonal_Obs[[(q-1)*n.samples.ht+k]][["percent.mt"]] <- PercentageFeatureSet(Clonal_Obs[[(q-1)*n.samples.ht+k]], pattern = "^MT-")
Clonal_Obs[[(q-1)*n.samples.ht+k]]<- subset(Clonal_Obs[[(q-1)*n.samples.ht+k]], subset = nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < percent.mt_upper)
}
}
rm(clonal.data)
}else{
  for(k in 1:n.samples){
    print(paste(paste("Reading in gene expression for sample ",k,sep=""),".",sep=""))
    clonal.data <- Read10X(data.dir = gex.path[[k]])
  Clonal_Obs[[k]] <- CreateSeuratObject(counts = clonal.data,project = sample.name[k])
  Clonal_Obs[[k]][["percent.mt"]] <- PercentageFeatureSet(Clonal_Obs[[k]], pattern = "^MT-")
  Clonal_Obs[[k]]<- subset(Clonal_Obs[[k]], subset = nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < percent.mt_upper)

}
rm(clonal.data)
}

for(k in 1:n.samples){
RnaStoreUMO=Clonal_Obs[[k]]@assays$RNA@counts

#VDJ list with barcodes
print(paste(paste("Reading in VDJ for sample ",k,sep=""),".",sep=""))
if(length(vdj.path)>0){
if(demultiplex){
  if(k==1){
  vdj.path.temp=vdj.path
  vdj.path=0
  for(q in 1:n.hashtag.samples){
    vdj.path=append(vdj.path,rep(vdj.path.temp[q],n.samples.ht))
  }
  vdj.path=vdj.path[-1]
  }
  VDJ_Obs[[k]] <- read.csv(vdj.path[[k]])
}else{
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
d1=as.matrix(Clonal_Obs[[k]][["SCT"]]@data)[Gene_indUMO,]
alist <- 1:length(Gene_indUMO)
if((k %in% index.control)&&(preset==1)){
  if(length(Gene_indUMO)==1){
    cutoff[[k]]=lapply(alist, function(alist) quantile(unname(d1),threshold.cutoff)[[1]])
  }
  else{
cutoff[[k]]=lapply(alist, function(alist) quantile(d1[alist,],threshold.cutoff)[[1]])
}
}else if(preset==0){
  cutoff[[k]]=threshold.entry
}else if(preset==2){
  for(w in 1:n.samples){
    if(w==1){
    d1=as.matrix(Clonal_Obs[[k]][["SCT"]]@data)[Gene_indUMO,]
    }else{
      d1=cbind(d1,as.matrix(Clonal_Obs[[k]][["SCT"]]@data)[Gene_indUMO,])
    }
  }
  for(w in 1:n.samples){
  cutoff[[w]]=lapply(alist, function(alist) quantile(d1[alist,],threshold.cutoff)[[1]])
  }
  }else{
  print("cutoff")
  print(cutoff[[k]][1])
}

Thresholds=matrix("unassigned",nrow=dim(Clonal_Obs[[k]])[2],ncol=length(Gene_indUMO))
for(s in 1:length(Gene_indUMO)){
  vec1=Clonal_Obs[[k]][Gene_indUMO[s],][["SCT"]]@data
for(j in 1:dim(Clonal_Obs[[k]])[2]){
  if(vec1[j]>cutoff[[c.index[k]]][s]){
    Thresholds[j,s]="high"
  }
}
  Clonal_Obs[[k]]=AddMetaData(Clonal_Obs[[k]], Thresholds[,s], col.name = paste("Threshold_",marker.gene[s],sep=""))
}

rm(mvsts.UMO)
rm(RnaStoreUMO)

}

###############################################################
##-----------------------------------------------------------##
##         Annotate cells with clonotypes (TCRB)             ##
##-----------------------------------------------------------##
###############################################################
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
Dataset=Clonal_Obs
saveRDS(Dataset,file.saved)
return(Dataset)
}
