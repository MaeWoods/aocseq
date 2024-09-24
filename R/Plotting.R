#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Functions to plot data stored in clonal objects
#' at normalized counts greater than the control
#' @import ggpubr
#'
#' This function will read in Seurat objects processed by aocseq::CombineData
#' and generate a ranked spreadsheet of clonotypes grouped by CDR3 beta sequence.
#' The table displays the percentage of CD4 and CD8 cells expressing a choice of marker
#' genes at normalized counts greater than the control.
#'
#' @param SummaryTable A Seurat object pre-processed with aocseq::CombineData.
#' @param names A marker of interest.
#' @param n.segments Threshold of the goi normalized control counts per cell.
#' @param tune.gap List. Index of the control sequencing sample for each condition.
#' @param segment.names.spacer Number of experiments within ClonalData object
#' @param percentage.spacer.1 String. Directory of the table output.
#' @param percentage.spacer.2 List. Control and sample names.
#' @param segment.col Bool. If set to 1, the control is part of the dataset. If 0, a threshold
#' value for each condition must be set.
#' @param ribbon.hue List. List of thresholds for each condition if preset is 0.
#' @param path Path to outpt file.
#' @return A data frame containing clonotypes and summary statistics of transcript expression.
#' @concept Annotation
#' @export
segmentPlot <- function(
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
    path="circos.pdf"
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

  pdf(path,width=6.5,height=6.5)
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

}




#' Export differential gene expression from several different samples and phenotypes
#'
#' This function will take a set of annotation spreadsheets and export differential
#' expression between different cell types labelled by their marker genes
#' Currently the phenotypic separation is CD4 and CD8, but the reader is encouraged
#' to generate as many subsets as needed for downstream analysis
#'
#' @param Clonal_Obs A Seurat object pre-processed with aocseq::CombineData.
#' @param clonotype.path Character array. Directory of an aocseq clonotype annotation table.
#' @param save.dir Directory for storing differentially expressed genes.
#' @param goi Character array. Gene name, a marker of interest.
#' @param verbose Print progress
#'
#' @return return an assay containing predicted expression value in the data
#' slot
#' @concept integration
#' @export
UMAPReduce <- function(
    Clonal_Obs,
    save.dir=".",
    clonotypes=c(),
    ident.list=c("orig.ident"),
    rseed=356,
    ur.nfeatures=500,
    preload=FALSE,
    preloadpath="."){
  
  
  set.seed(rseed)
  # split the dataset into a list of 3 seurat objects
  
  # normalize and identify variable features for each dataset independently
  
  
  if(length(Clonal_Obs)>1){
    
    immune.list= vector(mode = "list", length = length(Clonal_Obs))
    for(k in 1:length(Clonal_Obs)){
      immune.list[[k]]=Clonal_Obs[[k]]
    }
    rm(Clonal_Obs)
    immune.list <- lapply(X = immune.list, FUN = function(x) {
      x <- Seurat::NormalizeData(x, do.scale=FALSE, nfeatures = ur.nfeatures)
      x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = ur.nfeatures)
    })
    features <- Seurat::SelectIntegrationFeatures(do.scale=FALSE, object.list = immune.list,nfeatures = ur.nfeatures,features.to.integrate=ur.nfeatures)
    if(preload){
      immune.anchors <- readRDS(preloadpath)
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
    
    if(length(clonotypes)>0){
      for(j in 1:length(clonotypes)){
        pdf(paste(paste("UMAP_clonotype",j,sep=""),".pdf",sep=""),width=10,height=5)
        print(DimPlot(subset(immune.integrated,cdr3_na %in% clonotypes[k]), reduction = "umap",split.by = indent.list[k]))
        dev.off()
      }
    }
    
  }
  else{
    x <- Seurat::NormalizeData(Clonal_Obs, do.scale=FALSE, nfeatures = ur.nfeatures)
    x <- Seurat::FindVariableFeatures(Clonal_Obs, selection.method = "vst", nfeatures = ur.nfeatures)
    
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
#' This function will take a set of annotation spreadsheets and export differential
#' expression between different cell types labelled by their marker genes
#' Currently the phenotypic separation is CD4 and CD8, but the reader is encouraged
#' to generate as many subsets as needed for downstream analysis
#'
#' @param Clonal_Obs A Seurat object pre-processed with aocseq::CombineData.
#' @param clonotype.path Character array. Directory of an aocseq clonotype annotation table.
#' @param save.dir Directory for storing differentially expressed genes.
#' @param goi Character array. Gene name, a marker of interest.
#' @param verbose Print progress
#'
#' @return return an assay containing predicted expression value in the data
#' slot
#' @concept integration
#' @export
QCPlot <- function(
    Clonal_Obs,
    nFeature_RNA_lower=100,
    nFeature_RNA_upper=10000,
    nCount_RNA_lower=100,
    nCount_RNA_upper=10000,
    percent.mt_upper=5
){
  
  ##Create a data frame
  combinedCounts=data.frame(nCount_RNA=Clonal_Obs@meta.data$nCount_RNA,
                            nFeature_RNA=Clonal_Obs@meta.data$nFeature_RNA,
                            percent.mt=Clonal_Obs@meta.data$percent.mt,
                            name=Clonal_Obs@meta.data$orig.ident,
                            ColSplitCount=ifelse(((Clonal_Obs@meta.data$nCount_RNA>nCount_RNA_upper)|
                                                  (Clonal_Obs@meta.data$nCount_RNA<nCount_RNA_lower)),"remove","use"),
                            ColSplitRNA=ifelse(((Clonal_Obs@meta.data$nFeature_RNA>nFeature_RNA_upper)|
                                                  (Clonal_Obs@meta.data$nFeature_RNA<nFeature_RNA_lower)),"remove","use"),
                            ColSplitMito=ifelse(((Clonal_Obs@meta.data$percent.mt>percent.mt_upper)),"remove","use"))
  
  ##Plot the data frame
  par(mfrow = c(3, 1))
  p1=ggplot() + 
    geom_violin(data=combinedCounts,aes(x=factor(name), y=nFeature_RNA ), alpha = 0.6, fill="gray")+
    geom_point(data=combinedCounts,aes(x=factor(name), y=nFeature_RNA ,color=factor(ColSplitRNA)), alpha = 1 )+
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
    geom_violin(data=combinedCounts,aes(x=factor(name), y=nCount_RNA ),fill="gray", alpha = 0.3 )+
    geom_point(data=combinedCounts,aes(x=factor(name), y=nCount_RNA ,color=factor(ColSplitCount)), alpha = 1 )+
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
    geom_violin(data=combinedCounts,aes(x=factor(name), y=percent.mt ), alpha = 0.6 )+
    geom_point(data=combinedCounts,aes(x=factor(name), y=percent.mt ,color=factor(ColSplitMito)), alpha = 1 )+
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
  
  print(ggpubr::ggarrange(p1,p2,p3,ncol=3), common.legend = TRUE, legend = "top")
  
}


SaveHeatmap <- function(
    Resting.cells,
    Act.ref,
    CellsInRef,
    path.glist,
    heatmap.path,
    verbose=TRUE
){
  
  Resting.cellsinput=Resting.cells
  
  Resting.cells=as.matrix(Resting.cells@assays$RNA@counts)
  PFcolsRoundone=mean(colSums(Resting.cells))
  Resting.cells=log((1+(Resting.cells)))
  PFcolsRoundtwo=mean(colSums(Resting.cells))
  Resting.cells=Resting.cells/PFcolsRoundtwo
  rm(PFcolsRoundone)
  rm(PFcolsRoundtwo)
  
  # colnamesSig2=colnames(signature2.ref)
  # GlistMat2=match(Glist,row.names(Resting.cells3))
  # Sig3_M=Resting.cells3[GlistMat2,colnamesSig2]
  
  
  Glist=setdiff(read.csv(path.glist)$x,"")[1:20]
  cutoff=rep(0,length(Glist))
  GlistIndSquer=match(Glist,row.names(Resting.cells))
  
  GlistIndS=match(Glist,row.names(Resting.cells))
  
  for(g in 1:6){
    cutoff[g]=quantile(Resting.cells[GlistIndS[g],],.915)
  }
  
  for(g in 7:7){
    cutoff[g]=quantile(Resting.cells[GlistIndS[g],],.915)
  }
  
  for(g in 8:15){
    cutoff[g]=quantile(Resting.cells[GlistIndS[g],],.915)
  }
  
  for(g in 16:20){
    cutoff[g]=quantile(Resting.cells[GlistIndS[g],],.915)
  }
  
  
  Thresholds=matrix("unassigned",nrow=dim(Resting.cells)[2],ncol=length(GlistIndS))
  for(s in 1:length(GlistIndS)){
    vec1=Resting.cells[GlistIndS[s],]
    for(j in 1:dim(Resting.cells)[2]){
      if(vec1[j]>cutoff[s]){
        Thresholds[j,s]="high"
      }
    }
    Resting.cellsinput=AddMetaData(Resting.cellsinput, Thresholds[,s], col.name = paste("Signature_",Glist[s],sep=""))
  }
  
  MatrixOfValues=as.matrix(Resting.cellsinput@meta.data[19:38])
  NcellsSigmat=dim(MatrixOfValues)[[1]]
  SignatureCell=rep(0,NcellsSigmat)
  for(h in 1:NcellsSigmat){
    if(length(subset(MatrixOfValues[h,],MatrixOfValues[h,]=="high"))>4){
      SignatureCell[h]=1
    }
  }
  
  
  Resting.cellsinput=AddMetaData(Resting.cellsinput, SignatureCell, col.name = "SignatureCell")
  signature5.ref <- subset(Resting.cellsinput,(SignatureCell==1))
  
  colnamesSig5=colnames(signature5.ref)
  GlistMat5=match(Glist,row.names(Resting.cells))
  Sig_Unstim=Resting.cells[GlistMat5,colnamesSig5]
  
  
  
  
  
  
  Glist=setdiff(read.csv(path.glist)$x,"")[1:20]
  cutoff=rep(0,length(Glist))
  GlistIndS=match(Glist,row.names(Sig_Unstim))
  
  #Act.ref
  mat11=Sig_Unstim[GlistIndS,sample(1:dim(Sig_Unstim)[2],dim(Act.ref)[2])]
  
  
  Glist=Glist[1:20]
  library(preprocessCore)
  col1Mat=cbind(Act.ref[1:20,],mat11)
  
  col1Mat=normalize.quantiles(t(as.matrix(col1Mat)))
  colMat=as.matrix(col1Mat)
  
  
  
  #colMat=as.matrix(rbind(cbind(col1Mat,col2Mat),cbind(col3Mat,col4Mat)))
  NWgenes=Glist
  HeatmapDataNW=data.frame(Name=Glist,t(colMat))
  
  Ng_T=length(Glist)
  Ng_T=1*length(Glist)
  allcells=1
  allgenes=1:(length(Glist))
  
  lengthA=CellsInRef
  lengthB=CellsInRef
  
  
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
  
  pdf(heatmap.path,width=9,height=4)  
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
