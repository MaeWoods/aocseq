Rcpp::sourceCpp("GetMahalanobis.cpp")

CD8MaxClonotypeUpper=quantile(subset(Summarydf,(Summarydf[,3]/100)*Summarydf[,5]>0 & Summarydf[,5]>0)[,3],1)[[1]]
CD8MaxClonotype=quantile(subset(Summarydf,(Summarydf[,3]/100)*Summarydf[,5]>0 & Summarydf[,5]>0)[,3],lower_lim)[[1]]
TCRActCD8=list(subset(Summarydf,(Summarydf[,3]>=CD8MaxClonotype) & (Summarydf[,3]<=CD8MaxClonotypeUpper))$cdr3_na)

#SignatureCellsI = subset(Clonal_T,(cdr3_na %in% TCRAct[[1]]) & (Threshold_IFNG=="EBV_high") & (CD8cells==1))


SixHourStim=readRDS("../RDS/ThreeHourStim.rds")
Clonal_T=SixHourStim[[3]]
Glist=read_csv("../SupplementaryTables/SummaryGenes/Genes_upper_3hour.csv")$x

ClonotypeTable_path = "../SupplementaryTables/SummaryIFNG_2_5%.csv"
Summarydf2=read.csv(ClonotypeTable_path)
Summarydf <- Summarydf2[,c(2,9,10,11,17)]
names(Summarydf) <- c("cdr3_na","CD43hrCMV","CD83hrCMV",
                      "3hrCMV",
                      "N3hrCMV")


upper_lim=0.99999999
lower_lim=0.8
bystander_lim=0.4

TraceSignature <- function(Clonal_T,Summarydf,path_Glist,upsetnames,upset_path,combination_path,SignatureCellsI){
  
Glist=read_csv(path_Glist)$x

##CD8 high IFNG for plotting
CD8MaxClonotypeUpper=quantile(subset(Summarydf,(Summarydf[,3]/100)*Summarydf[,5]>0 & Summarydf[,5]>0)[,3],1)[[1]]
CD8MaxClonotype=quantile(subset(Summarydf,(Summarydf[,3]/100)*Summarydf[,5]>0 & Summarydf[,5]>0)[,3],lower_lim)[[1]]
TCRActCD8=list(subset(Summarydf,(Summarydf[,3]>=CD8MaxClonotype) & (Summarydf[,3]<=CD8MaxClonotypeUpper))$cdr3_na)

Nspecific=dim(SignatureCellsI)[2]
GlistInd=match(Glist,row.names(SignatureCellsI))
if(scramble==TRUE){
  GlistInd2=ceiling(runif(length(GlistInd),0,dim(Clonal_T)[1]))
}
else{
GlistInd2=match(Glist,row.names(Clonal_T))
}
SignatureCells=SignatureCellsI[["SCT"]]@data[GlistInd,]
sigvec=rep(0,dim(SignatureCells)[1]*dim(SignatureCells)[2])
for(s in 1:dim(SignatureCells)[1]){
  for(g in 1:dim(SignatureCells)[2]){
    sigvec[(s-1)*dim(SignatureCells)[2] + g] = SignatureCells[s,g]
  }
}
sigvec = as.double(sigvec)

TestCells=Clonal_T[["SCT"]]@data[GlistInd2,]
Tvec=rep(0,dim(TestCells)[1]*dim(TestCells)[2])
for(s in 1:dim(TestCells)[1]){
  for(g in 1:dim(TestCells)[2]){
    Tvec[(s-1)*dim(TestCells)[2] + g] = TestCells[s,g]
  }
}
Nsig=dim(SignatureCells)[2]
Ngenes=length(Glist)
Ncells=dim(TestCells)[2]
SpecificityDistance=rep(0,dim(TestCells)[2]) 
TestCells=Tvec
dist=rep(0,Ncells)
meansinp=rep(0,Ngenes)
colsv=rep(0,Ngenes)
covinpt=rep(0,Ngenes*Ngenes)
sigtest=rep(0,Ngenes*(Ncells+1))

SQ=GetMahalanobis(Nsig, Ncells, Ngenes, as.double(sigvec), as.double(TestCells),as.double(dist),as.double(meansinp),as.double(colsv),as.double(covinpt),as.double(sigtest))


distarray=array( unlist( SQ ))

Clonal_T=AddMetaData(Clonal_T, distarray, col.name = 'Mdist')

CellIDsAgSpec=colnames(subset(Clonal_T,(cdr3_na %in% TCRAct[[1]]) & (Threshold_IFNG=="BKV_high") & (CD8cells==1)))

plotdfEBV=data.frame(CellBCs=colnames(Clonal_T),MahalanobisDistance=Clonal_T@meta.data$Mdist,SpecificityByMarker=rep(0,length(colnames(Clonal_T))))

for(j in 1:length(colnames(Clonal_T))){
  if(length(intersect(plotdfEBV$CellBCs[j],CellIDsAgSpec)>0)){
    plotdfEBV$SpecificityByMarker[j]=1
  }
}


ggplot(plotdfEBV, aes(y = MahalanobisDistance, x = SpecificityByMarker, fill = SpecificityByMarker, color = SpecificityByMarker, group = SpecificityByMarker)) +
  geom_boxplot(alpha = 0.1) +
  geom_point()

UpperQuantileEBV=quantile(subset(plotdfEBV,SpecificityByMarker==1)$MahalanobisDistance,.5)[[1]]
BarcodesForMdistCells=subset(plotdfEBV,MahalanobisDistance<UpperQuantileEBV)$CellBCs

TCRsdataframe=data.frame(TCR1=Clonal_T@meta.data$cdr3_na,BCsTCR=colnames(Clonal_T),spec=plotdfEBV$SpecificityByMarker,
                       MD=plotdfEBV$MahalanobisDistance  )

AllAgSpecTCRs=subset(TCRsdataframe,BCsTCR %in% BarcodesForMdistCells)$TCR1
AllSpec=subset(TCRsdataframe,BCsTCR %in% BarcodesForMdistCells)$spec

AllAgSpecTCRsZero=subset(TCRsdataframe,(BCsTCR %in% BarcodesForMdistCells) & (spec==0))$TCR1


FinalDF=data.frame(TCRall=Summarydf$cdr3_na,MdistCat=rep("No",length(Summarydf$cdr3_na)))

for(j in 1:length(FinalDF$TCRall)){
  if(length(intersect(FinalDF$TCRall[j],AllAgSpecTCRs))>0){
    FinalDF$MdistCat[j]="Yes"
  }
}
Summarydf$CategoryEBV=FinalDF$MdistCat
write.csv(Summarydf,"../SupplementaryTables/CMVCellsPredictedIFNG_2.csv")
}
