
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Functions to plot data stored in clonal objects
#' at normalized counts greater than the control
#' @import ggpubr
#' @import circlize
#'
#' This function will read in Seurat objects processed by aocseq::CombineData
#' and generate a ranked spreadsheet of cells grouped by the aocseq cell type meta data.
#' The table displays the percentage of CD4 and CD8 cells expressing a choice of marker
#' genes at normalized counts greater than the control.
#'
#' @param SummaryTable An aocseq summary table listing cell type frequencies.
#' @param names Names of samples that are being compared.
#' @param n.segments Number of segments to plot.
#' @param tune.gap Resolution of circos plot.
#' @param segment.names.spacer Parameter to control spacing of circos plot labels.
#' @param percentage.spacer.1 Parameter to control spacing of circos plot labels of sample 1.
#' @param percentage.spacer.2 Parameter to control spacing of circos plot labels of sample 2.
#' @param segment.col Colors of circos plot labels.
#' @param ribbon.hue Colors of circos plot ribbons.
#' @param save.dir save.dir to output file.
#' @return A circos plot of aocseq cell type frequencies.
#' @concept Visualization
#' @export
SegmentPlot <- function(
    SummaryTable,
    names=c("S1","S2"),
    segment.names=c("S1","S2"),
    n.segments=50,
    tune.gap=0.005,
    segment.names.spacer="16mm",
    percentage.spacer.1="14mm",
    percentage.spacer.2="13mm",
    segment.col=c("green","blue"),
    ribbon.hue="pink",
    save.dir="circos.pdf"
){
  s.one=match(names[1], names(SummaryTable))
  s.two=match(names[2], names(SummaryTable))
  SummaryTable=subset(SummaryTable,SummaryTable[,s.one]>0 |SummaryTable[,s.two]>0 )
  SummaryTable[,s.one]=SummaryTable[,s.one]/(sum(SummaryTable[,s.one]))
  SummaryTable[,s.two]=SummaryTable[,s.two]/(sum(SummaryTable[,s.two]))
  s1_indx=order(SummaryTable[,s.one],decreasing=TRUE)[1:n.segments]
  s2_indx=order(SummaryTable[,s.two],decreasing=TRUE)[1:n.segments]

  idex_s1=0
  idex_s2=0
  for(k in 1:n.segments){
    for(j in 1:n.segments){

      if(s1_indx[k]==s2_indx[j]){
        idex_s1=append(idex_s1,k)
        idex_s2=append(idex_s2,j)
      }

    }
  }
  idex_s1=idex_s1[-1]
  idex_s2=idex_s2[-1]

  pdf(save.dir,width=6.5,height=6.5)
  percentages=c(SummaryTable[,s.one][s1_indx],SummaryTable[,s.two][s2_indx])
  circos.clear()
  factors = 1:length(percentages)
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  circos.par(cell.padding = c(0, 0, 0, 0),gap.degree = tune.gap)
  circos.initialize(factors, xlim = cbind(rep(0, length(percentages)), percentages))
  circos.track(ylim = c(3, 4), track.height = 0.2, bg.border = NA)
  circos.track(ylim = c(0, 2), track.height = 0.15,
               bg.col = c(rep(segment.col[1],length(s1_indx)),rep(segment.col[2],length(s2_indx))), bg.border = NA)
  p_links1=idex_s1
  p_links2=idex_s2+length(s1_indx)
  for(i in 1:length(p_links1)){
    circos.link(p_links1[i], c(0,percentages[p_links1[i]]), p_links2[i], c(0,percentages[p_links2[i]]),
                col = rand_color(1, transparency = 0.5,hue="pink",luminosity =c("bright","light")), border = NA)

  }
  highlight.sector(1:1, track.index = 2, text.col=segment.col[1], text = paste(toString(signif((SummaryTable[,s.one][s1_indx[1]]),2)*100),"%"),
                   facing = "clockwise", niceFacing = TRUE, text.vjust = percentage.spacer.1, cex = 1,col = NA)
  highlight.sector(2:2, track.index = 2, text.col=segment.col[1], text = paste(toString(signif((SummaryTable[,s.one][s1_indx[2]]),2)*100),"%"),
                   facing = "clockwise", niceFacing = TRUE, text.vjust = percentage.spacer.1, cex = 1,col = NA)
  highlight.sector((length(s1_indx)+1):(length(s1_indx)+1), track.index = 2, text.col=segment.col[2], text = paste(toString(signif((SummaryTable[,s.two][s2_indx[1]]),2)*100),"%"),
                   facing = "clockwise", niceFacing = TRUE, text.vjust = percentage.spacer.2, cex = 1,col = NA)
  highlight.sector((length(s1_indx)+2):(length(s1_indx)+2), track.index = 2, text.col=segment.col[2], text = paste(toString(signif((SummaryTable[,s.two][s2_indx[2]]),2)*100),"%"),
                   facing = "clockwise", niceFacing = TRUE, text.vjust = percentage.spacer.2, cex = 1,col = NA)
  highlight.sector(1:length(s1_indx), track.index = 2, text.col=segment.col[1],text = segment.names[1],
                   facing = "bending.inside", niceFacing = TRUE, text.vjust = segment.names.spacer, cex = 1.3,col = NA)
  highlight.sector((length(s1_indx)):((length(s2_indx)+length(s1_indx))),text.col=segment.col[2], track.index = 2, text = segment.names[2],
                   facing = "bending.inside", niceFacing = TRUE, text.vjust = segment.names.spacer, cex = 1.3,col = NA)

  dev.off()
  
 
    percentages=c(SummaryTable[,s.one][s1_indx],SummaryTable[,s.two][s2_indx])
    circos.clear()
   
    factors = 1:length(percentages)
    par(mar = c(0.2, 0.2, 0.2, 0.2))
   
    circos.par(cell.padding = c(0, 0, 0, 0),gap.degree = tune.gap)
    circos.initialize(factors, xlim = cbind(rep(0, length(percentages)), percentages))
    circos.track(ylim = c(3, 4), track.height = 0.2, bg.border = NA)
    circos.track(ylim = c(0, 2), track.height = 0.15,
                 bg.col = c(rep(segment.col[1],length(s1_indx)),rep(segment.col[2],length(s2_indx))), bg.border = NA)
    p_links1=idex_s1
    p_links2=idex_s2+length(s1_indx)
    for(i in 1:length(p_links1)){
      circos.link(p_links1[i], c(0,percentages[p_links1[i]]), p_links2[i], c(0,percentages[p_links2[i]]),
                  col = rand_color(1, transparency = 0.5,hue="pink",luminosity =c("bright","light")), border = NA)
      
    }
    highlight.sector(1:1, track.index = 2, text.col=segment.col[1], text = paste(toString(signif((SummaryTable[,s.one][s1_indx[1]]),2)*100),"%"),
                     facing = "clockwise", niceFacing = TRUE, text.vjust = percentage.spacer.1, cex = 1,col = NA)
    highlight.sector(2:2, track.index = 2, text.col=segment.col[1], text = paste(toString(signif((SummaryTable[,s.one][s1_indx[2]]),2)*100),"%"),
                     facing = "clockwise", niceFacing = TRUE, text.vjust = percentage.spacer.1, cex = 1,col = NA)
    highlight.sector((length(s1_indx)+1):(length(s1_indx)+1), track.index = 2, text.col=segment.col[2], text = paste(toString(signif((SummaryTable[,s.two][s2_indx[1]]),2)*100),"%"),
                     facing = "clockwise", niceFacing = TRUE, text.vjust = percentage.spacer.2, cex = 1,col = NA)
    highlight.sector((length(s1_indx)+2):(length(s1_indx)+2), track.index = 2, text.col=segment.col[2], text = paste(toString(signif((SummaryTable[,s.two][s2_indx[2]]),2)*100),"%"),
                     facing = "clockwise", niceFacing = TRUE, text.vjust = percentage.spacer.2, cex = 1,col = NA)
    highlight.sector(1:length(s1_indx), track.index = 2, text.col=segment.col[1],text = segment.names[1],
                     facing = "bending.inside", niceFacing = TRUE, text.vjust = segment.names.spacer, cex = 1.3,col = NA)
    highlight.sector((length(s1_indx)):((length(s2_indx)+length(s1_indx))),text.col=segment.col[2], track.index = 2, text = segment.names[2],
                     facing = "bending.inside", niceFacing = TRUE, text.vjust = segment.names.spacer, cex = 1.3,col = NA)
    
  

}




#' Plot UMAP dimension reduction with aocseq meta data.
#'
#' This function will save dimension reduction UMAP plots generated by Seurat and overlay aocseq metadata.
#'
#' @param cell.data A Seurat object pre-processed with aocseq::CombineData.
#' @param save.dir save.dir to output file.
#' @param cell.types List of cell types for UMAP annotation.
#' @param ident.list Character array. Gene name, a marker of interest.
#' @param rseed Set seed for data recall.
#' @param ur.nfeatures Choose fewer dimension reduction parameters for rapid visualization of aocseq cell type metadata.
#' @param preload True or False. If integration anchors must be loaded, set this parameter to True, otherwise they will be computed.
#' @param preloadsave.dir save.dir to preload data, if preload is set to True.
#' @param verbose Print progress.
#'
#' @return return a Seurat UMAP plot with aocseq annotation.
#' slot
#' @concept Visualization
#' @export
UMAPReduce <- function(
    cell.data,
    save.dir=".",
    cell.types=c(),
    ident.list=c("orig.ident"),
    rseed=356,
    ur.nfeatures=500,
    preload=FALSE,
    preloadsave.dir="."){
  
  set.seed(rseed)
  # split the dataset into a list of 3 seurat objects
  
  # normalize and identify variable features for each dataset independently
  
  
  if(length(cell.data)>1){
    
    immune.list= vector(mode = "list", length = length(cell.data))
    for(k in 1:length(cell.data)){
      immune.list[[k]]=cell.data[[k]]
    }
    rm(cell.data)
    immune.list <- lapply(X = immune.list, FUN = function(x) {
      x <- Seurat::NormalizeData(x, do.scale=FALSE, nfeatures = ur.nfeatures)
      x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = ur.nfeatures)
    })
    features <- Seurat::SelectIntegrationFeatures(do.scale=FALSE, object.list = immune.list,nfeatures = ur.nfeatures,features.to.integrate=ur.nfeatures)
    if(preload){
      immune.anchors <- readRDS(preloadsave.dir)
      rm(immune.list)
    }
    else{
      immune.anchors <- Seurat::FindIntegrationAnchors(object.list = immune.list, anchor.features = features)
      
      rm(immune.list)
      saveRDS(immune.anchors,paste(save.dir,"/immune.anchors.rds",sep=""))
    }
    immune.integrated <- Seurat::IntegrateData(anchorset = immune.anchors,features.to.integrate=features)
    rm(immune.anchors)
    saveRDS(immune.integrated,paste(save.dir,"/immuneCombined.rds",sep=""))
    immune.integrated <- ScaleData(immune.integrated, verbose = FALSE)
    immune.integrated <- RunPCA(immune.integrated, npcs = 30, verbose = FALSE)
    immune.integrated <- RunUMAP(immune.integrated, reduction = "pca", dims = 1:30)
    #run this 2 together
    immune.integrated <- FindNeighbors(immune.integrated, reduction = "pca", dims = 1:30)
    immune.integrated <- FindClusters(immune.integrated, resolution = 0.5)
    
    for(k in 1:length(ident.list)){
      p=DimPlot(immune.integrated, reduction = "umap",split.by = ident.list[k])
      
      #immune.integrated=Idents(immune.integrated,indent.list[k])
      #DimPlot(immune.integrated, reduction = "umap",split.by = ident.list[k])
      pdf(paste(paste("UMAP",k,sep=""),".pdf",sep=""),width=10,height=5)
      print(DimPlot(immune.integrated, reduction = "umap",split.by = ident.list[k]))
      dev.off()
    }
    
    if(length(cell.types)>0){
      for(j in 1:length(cell.types)){
        pdf(paste(paste("UMAP_clonotype",j,sep=""),".pdf",sep=""),width=10,height=5)
        print(DimPlot(subset(immune.integrated,cdr3_na %in% cell.types[k]), reduction = "umap",split.by = indent.list[k]))
        dev.off()
      }
    }
    
  }
  else{
    x <- Seurat::NormalizeData(cell.data, do.scale=FALSE, nfeatures = ur.nfeatures)
    x <- Seurat::FindVariableFeatures(cell.data, selection.method = "vst", nfeatures = ur.nfeatures)
    
    immune.integrated <- ScaleData(x, verbose = FALSE)
    immune.integrated <- RunPCA(immune.integrated, npcs = 30, verbose = FALSE)
    immune.integrated <- RunUMAP(immune.integrated, reduction = "pca", dims = 1:30)
    #run this 2 together
    immune.integrated <- FindNeighbors(immune.integrated, reduction = "pca", dims = 1:30)
    immune.integrated <- FindClusters(immune.integrated, resolution = 0.5)
    
    for(k in 1:length(ident.list)){
      #immune.integrated=Idents(immune.integrated,indent.list[k])
      print(DimPlot(immune.integrated, reduction = "umap",split.by = ident.list[k]))
      pdf(paste(paste("UMAP",k,sep=""),".pdf",sep=""),width=10,height=5)
      print(DimPlot(immune.integrated, reduction = "umap",split.by = ident.list[k]))
      dev.off()
    }
    
  }
  
  
}



#' QC plot
#'
#' This function will import raw feature barcode matrix and plot the number of genes detected 
#' per cell, the number of reads detected per cell and the percentage of mitochondrial genes 
#' per cell. The results are displayed as QC thresholds to vizualize the cells that are removed from the data.
#'
#' @param cell.data A Seurat object pre-processed with aocseq::CombineData.
#' @param nFeature_RNA_lower Threshold for the minimum number of unique gene identifiers detected in a single cell.
#' @param nFeature_RNA_upper Threshold for the maximum number of unique gene identifiers detected in a single cell.
#' @param nCount_RNA_lower Threshold for the minimum number of total RNA reads detected in a single cell.
#' @param nCount_RNA_upper Threshold for the maximum number of total RNA reads detected in a single cell.
#' @param percent.mt_upper Threshold for the maximum percentage of unique mitochondrial gene identifiers detected in a single cell.
#' @param verbose Print progress
#'
#' @return return an assay containing predicted expression value in the data
#' slot
#' @concept Visualization
#' @export
QCPlot <- function(
    cell.data,
    nFeature_RNA_lower=100,
    nFeature_RNA_upper=10000,
    nCount_RNA_lower=100,
    nCount_RNA_upper=10000,
    percent.mt_upper=5
){
  
  ##Create a data frame
  combinedCounts=data.frame(nCount_RNA=cell.data@meta.data$nCount_RNA,
                            nFeature_RNA=cell.data@meta.data$nFeature_RNA,
                            percent.mt=cell.data@meta.data$percent.mt,
                            name=cell.data@meta.data$orig.ident,
                            ColSplitCount=ifelse(((cell.data@meta.data$nCount_RNA>nCount_RNA_upper)|
                                                  (cell.data@meta.data$nCount_RNA<nCount_RNA_lower)),"remove","use"),
                            ColSplitRNA=ifelse(((cell.data@meta.data$nFeature_RNA>nFeature_RNA_upper)|
                                                  (cell.data@meta.data$nFeature_RNA<nFeature_RNA_lower)),"remove","use"),
                            ColSplitMito=ifelse(((cell.data@meta.data$percent.mt>percent.mt_upper)),"remove","use"))
  
  positions=sample(length(cell.data@meta.data$nCount_RNA),5000)
  combinedCounts_low=data.frame(nCount_RNA=cell.data@meta.data$nCount_RNA[positions],
                            nFeature_RNA=cell.data@meta.data$nFeature_RNA[positions],
                            percent.mt=cell.data@meta.data$percent.mt[positions],
                            name=cell.data@meta.data$orig.ident[positions],
                            ColSplitCount=ifelse(((cell.data@meta.data$nCount_RNA[positions]>nCount_RNA_upper)|
                                                    (cell.data@meta.data$nCount_RNA[positions]<nCount_RNA_lower)),"remove","use"),
                            ColSplitRNA=ifelse(((cell.data@meta.data$nFeature_RNA[positions]>nFeature_RNA_upper)|
                                                  (cell.data@meta.data$nFeature_RNA[positions]<nFeature_RNA_lower)),"remove","use"),
                            ColSplitMito=ifelse(((cell.data@meta.data$percent.mt[positions]>percent.mt_upper)),"remove","use"))
  
  ##Plot the data frame
  p1=ggplot() + 
    geom_violin(data=combinedCounts,aes(x=factor(name), y=nFeature_RNA ), alpha = 0.6, fill="gray")+
    geom_point(data=combinedCounts_low,aes(x=factor(name), y=nFeature_RNA ,color=factor(ColSplitRNA)), alpha = 1 )+
    geom_hline(yintercept=nFeature_RNA_upper)+
    geom_hline(yintercept=nFeature_RNA_lower)+
    ylab("Count per cell") +
    xlab("Sample") +
    scale_color_manual("",values=c("red4","blue"))+
  ggtitle("Features") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 8),
          legend.position = "top",
          panel.background = element_blank())
  
  #print(p1)
  
  p2=ggplot() + 
    geom_violin(data=combinedCounts,aes(x=factor(name), y=nCount_RNA ),fill="gray", alpha = 0.6 )+
    geom_point(data=combinedCounts_low,aes(x=factor(name), y=nCount_RNA ,color=factor(ColSplitCount)), alpha = 1 )+
    geom_hline(yintercept=nCount_RNA_upper)+
    geom_hline(yintercept=nCount_RNA_lower)+
    ylab("Count per cell") +
    xlab("Sample") +
    scale_color_manual("",values=c("red4","blue"))+
  ggtitle("RNA count") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 8),
          legend.position = "top",
          panel.background = element_blank())
  
  #print(p2)
  
  p3=ggplot() + 
    geom_violin(data=combinedCounts,aes(x=factor(name), y=percent.mt ),fill="gray", alpha = 0.6 )+
    geom_point(data=combinedCounts_low,aes(x=factor(name), y=percent.mt ,color=factor(ColSplitMito)), alpha = 1 )+
    geom_hline(yintercept=percent.mt_upper)+
    ylab("Count per cell") +
    xlab("Sample") +
    scale_color_manual("",values=c("red4","blue"))+
  ggtitle("Percent mitochondria") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 8),
          legend.position = "top",
          panel.background = element_blank())
  
  ggpubr::ggarrange(p1,p2,p3,ncol=3)
  
}

#'
#' This function will take a set of Seurat objects and plot them as a heat map. Normalization options include no normalization - raw counts, log1p, sctransform, pflog1ppf
#'
#' @param cell.data A Seurat object pre-processed with aocseq::CombineData.
#' @param reference.data A reference matrix.
#' @param n.plot Number of cells in heatmap.
#' @param glist Character array. Gene name, a marker of interest.
#' @param save.dir Character array. Gene name, a marker of interest.
#' @param verbose Print progress.
#' @return Void.
#' @concept Visualization
#' 
#' @export
SaveHeatmap <- function(
    cell.data,
    reference.data,
    n.plot,
    glist,
    save.dir,
    verbose=TRUE
){
  
  cell.datainput=cell.data
  if(packageVersion("Seurat")<'5.0.0'){
  cell.data=as.matrix(cell.data@assays$RNA@counts)
  }
  else{
    cell.data=as.matrix(cell.data@assays$RNA$counts)
  }  
  PFcolsRoundone=mean(colSums(cell.data))
  cell.data=log((1+(cell.data)))
  PFcolsRoundtwo=mean(colSums(cell.data))
  cell.data=cell.data/PFcolsRoundtwo
  rm(PFcolsRoundone)
  rm(PFcolsRoundtwo)
  
  # colnamesSig2=colnames(signature2.ref)
  # glistMat2=match(glist,row.names(cell.data3))
  # Sig3_M=cell.data3[glistMat2,colnamesSig2]

  cutoff=rep(0,length(glist))
  glistIndSquer=match(glist,row.names(cell.data))
  
  glistIndS=match(glist,row.names(cell.data))
  
  for(g in 1:6){
    cutoff[g]=quantile(cell.data[glistIndS[g],],.915)
  }
  
  for(g in 7:7){
    cutoff[g]=quantile(cell.data[glistIndS[g],],.915)
  }
  
  for(g in 8:15){
    cutoff[g]=quantile(cell.data[glistIndS[g],],.915)
  }
  
  for(g in 16:20){
    cutoff[g]=quantile(cell.data[glistIndS[g],],.915)
  }
  
  
  Thresholds=matrix("unassigned",nrow=dim(cell.data)[2],ncol=length(glistIndS))
  for(s in 1:length(glistIndS)){
    vec1=cell.data[glistIndS[s],]
    for(j in 1:dim(cell.data)[2]){
      if(vec1[j]>cutoff[s]){
        Thresholds[j,s]="high"
      }
    }
    cell.datainput=AddMetaData(cell.datainput, Thresholds[,s], col.name = paste("Signature_",glist[s],sep=""))
  }
  
  MatrixOfValues=as.matrix(cell.datainput@meta.data[19:38])
  NcellsSigmat=dim(MatrixOfValues)[[1]]
  SignatureCell=rep(0,NcellsSigmat)
  for(h in 1:NcellsSigmat){
    if(length(subset(MatrixOfValues[h,],MatrixOfValues[h,]=="high"))>4){
      SignatureCell[h]=1
    }
  }
  
  
  cell.datainput=AddMetaData(cell.datainput, SignatureCell, col.name = "SignatureCell")
  signature5.ref <- subset(cell.datainput,(SignatureCell==1))
  
  colnamesSig5=colnames(signature5.ref)
  glistMat5=match(glist,row.names(cell.data))
  Sig_Unstim=cell.data[glistMat5,colnamesSig5]
  
  
  
  
  
  
  glist=setdiff(read.csv(save.dir.glist)$x,"")[1:20]
  cutoff=rep(0,length(glist))
  glistIndS=match(glist,row.names(Sig_Unstim))
  
  #reference.data
  mat11=Sig_Unstim[glistIndS,sample(1:dim(Sig_Unstim)[2],dim(reference.data)[2])]
  
  
  glist=glist[1:20]
  library(preprocessCore)
  col1Mat=cbind(reference.data[1:20,],mat11)
  
  col1Mat=normalize.quantiles(t(as.matrix(col1Mat)))
  colMat=as.matrix(col1Mat)
  
  
  
  #colMat=as.matrix(rbind(cbind(col1Mat,col2Mat),cbind(col3Mat,col4Mat)))
  NWgenes=glist
  HeatmapDataNW=data.frame(Name=glist,t(colMat))
  
  Ng_T=length(glist)
  Ng_T=1*length(glist)
  allcells=1
  allgenes=1:(length(glist))
  
  lengthA=n.plot
  lengthB=n.plot
  
  
  #Copy these subsets to a data frame so you can write as a matrix
  xv=rep(1:Ng_T,lengthA)
  vy=0
  yyv=0
  colsep=0
  for(j in 1:lengthA){
    vy=append(vy,HeatmapDataNW[,1+j])
    yyv=append(yyv,rep(j,Ng_T))
    colsep=append(colsep,1)
  }
  vy=vy[-1]
  yyv=yyv[-1]
  colsep=colsep[-1]
  #CD8RA
  xv=append(xv,rep(1:Ng_T,lengthB))
  for(j in 1:lengthB){
    vy=append(vy,HeatmapDataNW[,1+lengthA+j])
    yyv=append(yyv,rep(j+lengthA+5,Ng_T))
    colsep=append(colsep,2)
  }
  
  
  
  HeatDF=data.frame(x=xv,y=yyv,cc=vy,sep=colsep)
  
  pdf(save.dir,width=9,height=4)  
  print(ggplot(HeatDF,aes(x=y,y=101-x,fill=cc^(1/4),colour=cc^(1/4)))+
    geom_tile(size=0.1)+
    labs(fill="log(nCount)")+
    scale_colour_gradientn(colours=c("black", "cornflowerblue","blue4"))+
    scale_fill_gradientn(colours=c("black","midnightblue","navy","darkblue","blue4", "cornflowerblue","mediumspringgreen","white"))+
    # annotate("text",label=HeatmapDataNW[,1],x=rep(0.2,20),y=0:24,size=30,fontface="bold",color="black")+
    #  annotate("text",label=c("NTR","CLL1","CD70","BiCAR"),x=c(1,2,3,4),y=rep(-0.8,4),size=30,fontface="bold",color="black")+
    ggtitle("PFlog1pPF")+
    ylab(" ")+
    xlab("Construct")+
    theme(axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=9),
          legend.text = element_text(size=5),
          legend.key.width = unit(0.5, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.title = element_text(size=9),
          plot.title = element_text(size=9),
          panel.background = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()))
  dev.off()
}
