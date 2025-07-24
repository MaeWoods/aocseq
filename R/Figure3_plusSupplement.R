require(ggplot2)
require(pscl)
require(MASS)
require(boot)

library(ggplot2)
library(pscl)
library(MASS)
library(boot)
library(SiMRiv)
library(Seurat)
library(subplex)

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#######################
###     Figure3A    ###
#######################
Dataset1=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/ThreeHourStim.rds")[[1]]
IFNG_pos1=match("IFNG",row.names(Dataset1@assays$RNA@counts))
Dataset2=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/SixHourStim.rds")[[1]]
IFNG_pos2=match("IFNG",row.names(Dataset2@assays$RNA@counts))
inf_data_control = data.frame(ifng=c(Dataset1@assays$RNA@counts[IFNG_pos1,],
                                 Dataset2@assays$RNA@counts[IFNG_pos2,]))

rm(Dataset1)
rm(Dataset2)

Dataset1=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/ThreeHourStim.rds")[[2]]
IFNG_pos1=match("IFNG",row.names(Dataset1@assays$RNA@counts))
Dataset2=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/ThreeHourStim.rds")[[3]]
IFNG_pos2=match("IFNG",row.names(Dataset2@assays$RNA@counts))
inf_data_stim_3hour = data.frame(ifng=c(Dataset1@assays$RNA@counts[IFNG_pos1,],
                                     Dataset2@assays$RNA@counts[IFNG_pos2,]))

rm(Dataset1)
rm(Dataset2)

Dataset1=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/SixHourStim.rds")[[2]]
IFNG_pos1=match("IFNG",row.names(Dataset1@assays$RNA@counts))
Dataset2=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/SixHourStim.rds")[[3]]
IFNG_pos2=match("IFNG",row.names(Dataset2@assays$RNA@counts))
inf_data_stim_6hour = data.frame(ifng=c(Dataset1@assays$RNA@counts[IFNG_pos1,],
                                        Dataset2@assays$RNA@counts[IFNG_pos2,]))

inf_data_stim = data.frame(ifng=c(inf_data_stim_6hour$ifng,inf_data_stim_3hour$ifng))

Ndat = length(inf_data_control$ifng)
ifng_data_mix=inf_data_control$ifng
fr <- function(x) {   ## MLE of negative binomial distribution
  x1 <- x[1]
  x2 <- x[2]
  -1*sum(log(dnbinom(ifng_data_mix, size=x1, prob=x2, log = FALSE)))  
}
fr_zi <- function(x) {   ## MLE of zero inflated negative binomial distribution
  #parameters
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  -(sum(log(extraDistr::dzinb(ifng_data_mix,x1,x2,x3, log = FALSE))))
}
if(sum(ifng_data_mix)>0){
  outputT3=optim(par=c(10,0.06), fr,method = "BFGS")
  outputT2_zi=optim(par=c(10,0.06,0.5), fr_zi,method = "Nelder-Mead")
  x1o=outputT2_zi$par[1]
  x2o=outputT2_zi$par[2]
  x3o=outputT2_zi$par[3]
  s1=2*3-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
  s1_bic=3*log(length(ifng_data_mix))-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
  x1o=outputT3$par[1]
  x2o=outputT3$par[2]
  s2=2*2-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
  s2_bic=2*log(length(ifng_data_mix))-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
  AIC_score_control=s1-s2
  BIC_score_control=s1_bic-s2_bic
  if((s1-s2)<0){
    print(j)
    print(s1-s2)
  }
  if((s1_bic-s2_bic)<0){
    print("bic")
    print(s1_bic-s2_bic)
  }
}


Ndat = length(inf_data_stim$ifng)
ifng_data_mix=inf_data_stim$ifng
fr <- function(x) {   ## MLE of negative binomial distribution
  x1 <- x[1]
  x2 <- x[2]
  -1*sum(log(dnbinom(ifng_data_mix, size=x1, prob=x2, log = FALSE)))  
}
fr_zi <- function(x) {   ## MLE of zero inflated negative binomial distribution
  #parameters
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  
  -(sum(log(extraDistr::dzinb(ifng_data_mix,x1,x2,x3, log = FALSE))))
  
}
if(sum(ifng_data_mix)>0){
  outputT3=optim(par=c(10,0.06), fr,method = "BFGS")
  outputT2_zi=optim(par=c(10,0.06,0.5), fr_zi,method = "Nelder-Mead")
  x1o=outputT2_zi$par[1]
  x2o=outputT2_zi$par[2]
  x3o=outputT2_zi$par[3]
  s1=2*3-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
  s1_bic=3*log(length(ifng_data_mix))-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
  x1o=outputT3$par[1]
  x2o=outputT3$par[2]
  s2=2*2-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
  s2_bic=2*log(length(ifng_data_mix))-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
  AIC_score_stim=s1-s2
  BIC_score_stim=s1_bic-s2_bic
  if((s1-s2)<0){
    print(j)
    print(s1-s2)
  }
  if((s1_bic-s2_bic)<0){
    print("bic")
    print(s1_bic-s2_bic)
  }
}


Dataset1=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/ThreeHourStim.rds")[[1]]
uniqueClonotypes=levels(factor(Dataset1@meta.data$cdr3_na))
CC_seq=Dataset1@meta.data$cdr3_na
CC_size=Dataset1@meta.data$countcln
Lengths=0
sampls=0
df_c=data.frame(seq=CC_seq,si=CC_size)
for(k in 1:length(uniqueClonotypes)){
    Lengths=append(Lengths,length(subset(df_c,seq==uniqueClonotypes[k])$seq))
    sampls=append(sampls,1)
}
rm(Dataset1)

Dataset1=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/ThreeHourStim.rds")[[2]]
uniqueClonotypes=levels(factor(Dataset1@meta.data$cdr3_na))
CC_seq=Dataset1@meta.data$cdr3_na
CC_size=Dataset1@meta.data$countcln
df_c=data.frame(seq=CC_seq,si=CC_size)
for(k in 1:length(uniqueClonotypes)){
  Lengths=append(Lengths,length(subset(df_c,seq==uniqueClonotypes[k])$seq))
  sampls=append(sampls,2)
}
rm(Dataset1)

Dataset1=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/ThreeHourStim.rds")[[3]]
uniqueClonotypes=levels(factor(Dataset1@meta.data$cdr3_na))
CC_seq=Dataset1@meta.data$cdr3_na
CC_size=Dataset1@meta.data$countcln
df_c=data.frame(seq=CC_seq,si=CC_size)
for(k in 1:length(uniqueClonotypes)){
  Lengths=append(Lengths,length(subset(df_c,seq==uniqueClonotypes[k])$seq))
  sampls=append(sampls,3)
}
rm(Dataset1)

Lengths=Lengths[-1]
sampls=sampls[-1]
LC=data.frame(len=Lengths,stim=sampls)
options(scipen = 999)
C=ggplot() + # variant becomes a classifying factor
  geom_histogram(data=LC,aes(x=log10(len+1), y=..count..,fill= factor(stim) ,color=factor(stim)),size=0.5,  bins=12,alpha = 0.6)+
  ylab("Number of cells") +
  xlab("log clonotype size") +
    #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
  scale_fill_manual(values =c("lightgray","blue","orange","gray","blue3","tan1"))+
  scale_color_manual(values =c("lightgray","blue","orange","gray","blue3","tan1"))+
  #xlim(-0.5,8)+
  # ylim(0,1000)+
  ggtitle("Clonotype size: mVSTs") +
  theme(axis.line = element_line(colour = "black"),
        #  legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

pdf("/Users/maewoodsphd/mVSTManuscript/Figure3c.pdf",width=5,height=8)
C
dev.off()


###############################################################################
###     Plot negative binomial distributions but with AIC and BOC scores    ###
###############################################################################

Ndat = length(inf_data_control$ifng)
ifng_data_mix=inf_data_control$ifng
outputT3=optim(par=c(10,0.06), fr,method = "BFGS")
outputT2_zi=optim(par=c(10,0.06,0.5), fr_zi,method = "Nelder-Mead")
x1o=outputT2_zi$par[1]
x2o=outputT2_zi$par[2]
x3o=outputT2_zi$par[3]
s1=2*3-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
s1_bic=3*log(length(ifng_data_mix))-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
x1o=outputT3$par[1]
x2o=outputT3$par[2]
s2=2*2-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
s2_bic=2*log(length(ifng_data_mix))-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
AIC_score_control=s1-s2
BIC_score_control=s1_bic-s2_bic


x1o=outputT3$par[1]
x2o=outputT3$par[2]
xlnN1 = 0
retval = 0

#rectangle rule integration
for(j in 1:length(inf_data_control$ifng)){
  print(j)
  samplexv = ceiling(runif(1,-1,max(inf_data_control$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
  newx=dnbinom(samplexv, size=x1o, prob=x2o, log = FALSE)
  # print(newx)
  rv=runif(1,0,1)
  while(rv>newx){
    samplexv = ceiling(runif(1,-1,max(inf_data_control$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
    newx=dnbinom(samplexv, size=x1o, prob=x2o, log = FALSE)
    rv=runif(1,0,1)
    #print(newx)
  }
  if(j==1){
    xlnN12=samplexv
  }
  else{
    xlnN12=append(xlnN12,samplexv)
  }
}
xNB <- data.frame(xln1NB=xlnN12)
grob <- grobTree(textGrob(paste("M1-M2 AIC:",sprintf(AIC_score_control, fmt = '%#.2f')), x=0.6,  y=0.95, hjust=0,
                          gp=gpar(col="red", fontsize=10, fontface="italic")))
grob_bic <- grobTree(textGrob(paste("M1-M2 BIC:",round(BIC_score_control, digits=2)), x=0.6,  y=0.85, hjust=0,
                              gp=gpar(col="red", fontsize=10, fontface="italic")))
A11=ggplot() + # variant becomes a classifying factor
  geom_histogram(data=inf_data_control,aes(x=log(ifng+1), y=..count.. ), bins=40, alpha = 0.4 ,fill= "purple4",color="black")+
  geom_histogram(data=xNB,aes(x=log(xNB$xln1NB+1), y=..count.. ), bins=40, alpha = 0.6 , fill= "gray",color="black")+ # add a scatterplot; constant size, shape/fill depends on variant # add a scatterplot; constant size, shape/fill depends on variant
  ylab("Number of cells") +
  geom_vline(xintercept=1.098612)+
  xlab("log1p Interferon gamma counts") +
  annotation_custom(grob)+
  annotation_custom(grob_bic)+
  scale_y_log10()+
  xlim(-0.5,8)+
 # ylim(0,1000)+
  ggtitle("Negative binomial model: Control mVSTs") +
  theme(axis.line = element_line(colour = "black"),
        #  legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

Ndat = length(inf_data_stim$ifng)
ifng_data_mix=inf_data_stim$ifng
outputT3=optim(par=c(10,0.06), fr,method = "BFGS")
outputT2_zi=optim(par=c(10,0.06,0.5), fr_zi,method = "Nelder-Mead")
x1o=outputT2_zi$par[1]
x2o=outputT2_zi$par[2]
x3o=outputT2_zi$par[3]
s1=2*3-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
s1_bic=3*log(length(ifng_data_mix))-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
x1o=outputT3$par[1]
x2o=outputT3$par[2]
s2=2*2-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
s2_bic=2*log(length(ifng_data_mix))-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
AIC_score_stim=s1-s2
BIC_score_stim=s1_bic-s2_bic


x1o=outputT3$par[1]
x2o=outputT3$par[2]
xlnN12 = 0
retval = 0
#rectangle rule integration
for(j in 1:length(inf_data_stim$ifng)){
  print(j)
  samplexv = ceiling(runif(1,-1,max(inf_data_stim$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
  newx=dnbinom(samplexv, size=x1o, prob=x2o, log = FALSE)
  # print(newx)
  rv=runif(1,0,1)
  while(rv>newx){
    samplexv = ceiling(runif(1,-1,max(inf_data_stim$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
    newx=dnbinom(samplexv, size=x1o, prob=x2o, log = FALSE)
    rv=runif(1,0,1)
    #print(newx)
  }
  if(j==1){
    xlnN12=samplexv
  }
  else{
    xlnN12=append(xlnN12,samplexv)
  }
}


xNB2 <- data.frame(xln1NB2=xlnN12)
grob <- grobTree(textGrob(paste("M1-M2 AIC:",sprintf(AIC_score_stim, fmt = '%#.2f')), x=0.6,  y=0.95, hjust=0,
                          gp=gpar(col="red", fontsize=10, fontface="italic")))
grob_bic <- grobTree(textGrob(paste("M1-M2 BIC:",round(BIC_score_stim, digits=2)), x=0.6,  y=0.85, hjust=0,
                          gp=gpar(col="red", fontsize=10, fontface="italic")))
A12=ggplot() + # variant becomes a classifying factor
  geom_histogram(data=inf_data_stim,aes(x=log(ifng +1), y=..count.. ), bins=30, alpha = 0.4 ,fill= "purple4",color="black")+
  geom_histogram(data=xNB2,aes(x=log(xln1NB2+1), y=..count.. ), bins=30, alpha = 0.6 , fill= "gray",color="black")+ # add a scatterplot; constant size, shape/fill depends on variant # add a scatterplot; constant size, shape/fill depends on variant
  ylab("Number of cells") +
  xlab("log1p Interferon gamma counts") +
  #ylim(0,1000)+
  xlim(-0.5,8)+
  ggtitle("Negative binomial model: Stimulated mVSTs") +
  annotation_custom(grob)+
  annotation_custom(grob_bic)+
  scale_y_log10()+
  theme(axis.line = element_line(colour = "black"),
        #  legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

pdf("Figure3a.pdf",width=8,height=4)
grid.arrange(A11,A12, ncol = 2, nrow = 1, heights=c(0.5), widths=c(0.5,0.5))
dev.off()

##################################################
###     Corresponding supplementary table    #####
##################################################

stims_n=c("CAD-trio","CAD-tet")
orig.i_n=c("CAD-trio","CAD-tet")
freq_n=c("0.1 min.pct","0.25 min.pct","0.5 min.pct","0.75 min.pct","0.9 min.pct","0.99 min.pct")
DF_names=rep(0,10)
for(k in 1:1){
  for(j in 1:5){
    DF_names[(k-1)*5+j]=paste(stims_n[k],freq_n[j],sep="")
  }}
min_cellfraction=c(0.25,0.5,0.75,0.9,0.99)
DEsordered=vector(mode = "list", length = 14)
#clustervector=c(0,1,2,3,4,5,6,8,9,10,11,12,13,14,15,17)
#clustervector=c(0,1,2,3,4,5,6,7,8,9,12,13,14,16,17)
clustervector=c(17)
#clustervector=c(18)
for(cnum in clustervector){
  for(k in 1:1){
    for(j in 1:length(min_cellfraction)){
      DEGmarkers <- FindMarkers(subset(spatial.merge,orig.ident==orig.i_n[k]), group.by="seurat_clusters",ident.1 = c(cnum), logfc.threshold=0,
                                min.cells.feature = 1, min.cells.group = 1, test.use="bimod",min.pct = min_cellfraction[j])
      DEsordered[[(k-1)*5+j]]=data.frame(gene=row.names(DEGmarkers)[order(DEGmarkers$avg_log2FC,decreasing=TRUE)],
                                         `fold change`=DEGmarkers$avg_log2FC[order(DEGmarkers$avg_log2FC,decreasing=TRUE)],
                                         `p_value`=DEGmarkers$p_val_adj[order(DEGmarkers$avg_log2FC,decreasing=TRUE)])
      DEsordered[[(k-1)*5+j]]=subset(DEsordered[[(k-1)*5+j]], `p_value`<0.05)
      colnames(DEsordered[[(k-1)*5+j]]) <- gsub("\\.", " ", colnames(DEsordered[[(k-1)*5+j]]))
      colnames(DEsordered[[(k-1)*5+j]]) <- gsub("_", "-", colnames(DEsordered[[(k-1)*5+j]]))
    }
    print(k)
  }
  names(DEsordered) <- DF_names
  write_xlsx(DEsordered[1:5], path = paste(paste("CADTrioDEGs/Cluster",cnum,sep=""),".xlsx",sep=""))
}

##################################################
###     Corresponding supplementary figures    ###
##################################################

Dataset1=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/Tetramer.rds")[[1]]
Dataset2=readRDS("/Users/maewoodsphd/mVSTManuscript/RDS/Tetramer.rds")[[3]]
IFNG_pos1=match("IFNG",row.names(Dataset1@assays$RNA@counts))
IFNG_pos2=match("IFNG",row.names(Dataset2@assays$RNA@counts))
inf_data_control = data.frame(ifng=c(Dataset1@assays$RNA@counts[IFNG_pos1,]))
inf_data_stim = data.frame(ifng=c(Dataset2@assays$RNA@counts[IFNG_pos2,]))

Ndat = length(inf_data_control$ifng)
ifng_data_mix=inf_data_control$ifng
outputT3=optim(par=c(10,0.06), fr,method = "BFGS")
outputT2_zi=optim(par=c(10,0.06,0.5), fr_zi,method = "Nelder-Mead")
x1o=outputT2_zi$par[1]
x2o=outputT2_zi$par[2]
x3o=outputT2_zi$par[3]
s1=2*3-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
s1_bic=3*log(length(ifng_data_mix))-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
x1o=outputT3$par[1]
x2o=outputT3$par[2]
s2=2*2-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
s2_bic=2*log(length(ifng_data_mix))-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
AIC_score_control=s1-s2
BIC_score_control=s1_bic-s2_bic


x1o=outputT3$par[1]
x2o=outputT3$par[2]
xlnN1 = 0
retval = 0

#rectangle rule integration
for(j in 1:length(inf_data_control$ifng)){
  print(j)
  samplexv = ceiling(runif(1,-1,max(inf_data_control$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
  newx=dnbinom(samplexv, size=x1o, prob=x2o, log = FALSE)
  # print(newx)
  rv=runif(1,0,1)
  while(rv>newx){
    samplexv = ceiling(runif(1,-1,max(inf_data_control$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
    newx=dnbinom(samplexv, size=x1o, prob=x2o, log = FALSE)
    rv=runif(1,0,1)
    #print(newx)
  }
  if(j==1){
    xlnN12=samplexv
  }
  else{
    xlnN12=append(xlnN12,samplexv)
  }
}
xNB <- data.frame(xln1NB=xlnN12)
grob <- grobTree(textGrob(paste("M1-M2 AIC:",sprintf(AIC_score_control, fmt = '%#.2f')), x=0.6,  y=0.95, hjust=0,
                          gp=gpar(col="red", fontsize=10, fontface="italic")))
grob_bic <- grobTree(textGrob(paste("M1-M2 BIC:",round(BIC_score_control, digits=2)), x=0.6,  y=0.85, hjust=0,
                              gp=gpar(col="red", fontsize=10, fontface="italic")))
A13=ggplot() + # variant becomes a classifying factor
  geom_histogram(data=inf_data_control,aes(x=log(ifng+1), y=..count.. ), bins=40, alpha = 0.4 ,fill= "purple4",color="black")+
  geom_histogram(data=xNB,aes(x=log(xNB$xln1NB+1), y=..count.. ), bins=40, alpha = 0.6 , fill= "gray",color="black")+ # add a scatterplot; constant size, shape/fill depends on variant # add a scatterplot; constant size, shape/fill depends on variant
  ylab("Number of cells") +
  xlab("log1p Interferon gamma counts") +
  annotation_custom(grob)+
  annotation_custom(grob_bic)+
  scale_y_log10()+
  xlim(-0.5,8)+
  # ylim(0,1000)+
  ggtitle("M1: Stimulated EBVSTs") +
  theme(axis.line = element_line(colour = "black"),
        #  legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

Ndat = length(inf_data_stim$ifng)
ifng_data_mix=inf_data_stim$ifng
outputT3=optim(par=c(10,0.06), fr,method = "BFGS")
outputT2_zi=optim(par=c(10,0.06,0.5), fr_zi,method = "Nelder-Mead")
x1o=outputT2_zi$par[1]
x2o=outputT2_zi$par[2]
x3o=outputT2_zi$par[3]
s1=2*3-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
s1_bic=3*log(length(ifng_data_mix))-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
x1o=outputT3$par[1]
x2o=outputT3$par[2]
s2=2*2-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
s2_bic=2*log(length(ifng_data_mix))-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
AIC_score_stim=s1-s2
BIC_score_stim=s1_bic-s2_bic


x1o=outputT2_zi$par[1]
x2o=outputT2_zi$par[2]
x3o=outputT2_zi$par[3]
xlnN12 = 0
retval = 0
#rectangle rule integration
for(j in 1:length(inf_data_stim$ifng)){
  print(j)
  samplexv = ceiling(runif(1,-1,max(inf_data_stim$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
  newx=extraDistr::dzinb(samplexv, size=x1o, prob=x2o,pi=x3o, log = FALSE)
  # print(newx)
  rv=runif(1,0,1)
  while(rv>newx){
    samplexv = ceiling(runif(1,-1,max(inf_data_stim$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
    newx=extraDistr::dzinb(samplexv, size=x1o, prob=x2o,pi=x3o, log = FALSE)
    rv=runif(1,0,1)
    #print(newx)
  }
  if(j==1){
    xlnN12=samplexv
  }
  else{
    xlnN12=append(xlnN12,samplexv)
  }
}


xNB2 <- data.frame(xln1NB2=xlnN12)
grob <- grobTree(textGrob(paste("M1-M2 AIC:",sprintf(AIC_score_stim, fmt = '%#.2f')), x=0.2,  y=0.95, hjust=0,
                          gp=gpar(col="red", fontsize=10, fontface="italic")))
grob_bic <- grobTree(textGrob(paste("M1-M2 BIC:",round(BIC_score_stim, digits=2)), x=0.6,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=10, fontface="italic")))
A14=ggplot() + # variant becomes a classifying factor
  geom_histogram(data=inf_data_stim,aes(x=log(ifng +1), y=..count.. ), bins=30, alpha = 0.4 ,fill= "purple4",color="black")+
  geom_histogram(data=xNB2,aes(x=log(xln1NB2+1), y=..count.. ), bins=30, alpha = 0.6 , fill= "gray",color="black")+ # add a scatterplot; constant size, shape/fill depends on variant # add a scatterplot; constant size, shape/fill depends on variant
  ylab("Number of cells") +
  xlab("log1p Interferon gamma counts") +
  #ylim(0,1000)+
  xlim(-0.5,8)+
  ggtitle("M2: Stimulated CD45RA depleted EBVSTs") +
  annotation_custom(grob)+
  annotation_custom(grob_bic)+
  scale_y_log10()+
  theme(axis.line = element_line(colour = "black"),
        #  legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

pdf("SupplementForFigure3a.pdf",width=10,height=4)
grid.arrange(A13,A14, ncol = 2, nrow = 1, heights=c(0.5), widths=c(0.5,0.5))
dev.off()



CloneList=TCRsspecific

for(j in 1:length(CloneList)){
  if(length(intersect(CloneList[j],Dataset1@meta.data$cdr3_na))>0){
  if(length(c(subset(Dataset1,cdr3_na %in% CloneList[j])@assays$RNA@counts[IFNG_pos1,]))>15){
inf_data = data.frame(ifng=c(subset(Dataset3)@assays$RNA@counts[IFNG_pos3,]))
Ndat = length(inf_data$ifng)
ifng_data_mix=inf_data$ifng
fr <- function(x) {   ## MLE of negative binomial distribution
  x1 <- x[1]
  x2 <- x[2]
  -1*sum(log(dnbinom(ifng_data_mix, size=x1, prob=x2, log = FALSE)))  
}
fr_zi <- function(x) {   ## MLE of zero inflated negative binomial distribution
  #parameters
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  
  -(sum(log(extraDistr::dzinb(ifng_data_mix,x1,x2,x3, log = FALSE))))
  
}
if(sum(ifng_data_mix)>0){
outputT3=optim(par=c(10,0.06), fr,method = "BFGS")
outputT2_zi=optim(par=c(10,0.06,0.5), fr_zi,method = "Nelder-Mead")
x1o=outputT2_zi$par[1]
x2o=outputT2_zi$par[2]
x3o=outputT2_zi$par[3]
s1=2*3-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
s1_bic=3*log(length(ifng_data_mix))-2*(sum(log(extraDistr::dzinb(ifng_data_mix,x1o,x2o,x3o, log = FALSE))))
x1o=outputT3$par[1]
x2o=outputT3$par[2]
s2=2*2-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
s2_bic=2*log(length(ifng_data_mix))-2*(sum(log(dnbinom(ifng_data_mix,x1o,x2o, log = FALSE))))
if((s1-s2)<0){
print(j)
print(s1-s2)
}
if((s1_bic-s2_bic)<0){
  print("bic")
  print(s1_bic-s2_bic)
}
}
  }
  }
}


inf_data = data.frame(ifng=c(subset(Dataset3,cdr3_na %in% CloneList[66])@assays$RNA@counts[IFNG_pos3,]))
Ndat = length(inf_data$ifng)
ifng_data_mix=inf_data$ifng

fr <- function(x) {   ## MLE of negative binomial distribution
  x1 <- x[1]
  x2 <- x[2]
  -1*sum(log(dnbinom(inf_data_all_3_6$ifng, size=x1, prob=x2, log = FALSE)))  
}
outputT1=optim(par=c(1,0.3), fr,method = "L-BFGS-B")
outputT2=optim(par=c(2,0.1), fr,method = "Nelder-Mead")
outputT3=optim(par=c(10,0.06), fr,method = "BFGS")
outputT4=optim(par=c(1,0.3), fr,method = "CG",control=list(maxit=10000))
outputT5=optim(par=c(1,0.3), fr,method = "SANN",control=list(maxit=50000))

fr_zi <- function(x) {   ## MLE of zero inflated negative binomial distribution
  #parameters
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  
  -(sum(log(extraDistr::dzinb(inf_data_all_3_6$ifng,x1,x2,x3, log = FALSE))))
  
}
outputT1_zi=optim(par=c(1,0.1,0.1), fr_zi,method = "L-BFGS-B")
outputT2_zi=optim(par=c(50,0.1,0.9), fr_zi,method = "Nelder-Mead")
outputT3_zi=optim(par=c(2,0.2,0.1), fr_zi,method = "BFGS")
outputT4_zi=optim(par=c(20,0.2,0.1), fr_zi,method = "CG",control=list(maxit=5000))
outputT5_zi=optim(par=c(20,0.2,0.1), fr_zi,method = "SANN",control=list(maxit=5000))

#############Plot results
x1o=outputT2_zi$par[1]
x2o=outputT2_zi$par[2]
x3o=outputT2_zi$par[3]
xlnN1 = 0
retval = 0
#rectangle rule integration

#for(j in 1:length(log10(subset(ins_df,InsTrans==0)$InsSIns))){
for(j in 1:length(inf_data$ifng)){
  print(j)
  samplexv = ceiling(runif(1,-1,max(inf_data$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
  newx=(extraDistr::dzinb(samplexv, size=x1o, prob=x2o,pi=x3o, log = FALSE))
  # print(newx)
  rv=runif(1,0,1)
  while(rv>newx){
    samplexv = ceiling(runif(1,-1,max(inf_data$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
    newx=(extraDistr::dzinb(samplexv, size=x1o, prob=x2o,pi=x3o, log = FALSE))
    rv=runif(1,0,1)
    # print(newx)
  }
  if(j==1){
    xlnN1=samplexv
  }
  else{
    xlnN1=append(xlnN1,samplexv)
  }
}
xNB <- data.frame(xln1NB=xlnN1)
A11=ggplot() + # variant becomes a classifying factor
  geom_histogram(data=inf_data,aes(x=log(inf_data$ifng+1), y=..count.. ), bins=40, alpha = 0.6 ,fill= "blue")+
  geom_histogram(data=xNB,aes(x=log(xNB$xln1NB+1), y=..count.. ), bins=40, alpha = 0.4 , fill= "gray",color="black")+ # add a scatterplot; constant size, shape/fill depends on variant # add a scatterplot; constant size, shape/fill depends on variant
  ylab("Number of cells") +
  xlab("log1p Interferon gamma counts") +
  ylim(0,1000)+
  ggtitle("Zero inflated negative binomial model") +
  theme(axis.line = element_line(colour = "black"),
        #  legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



x1o=outputT3$par[1]
x2o=outputT3$par[2]
xlnN12 = 0
retval = 0
#rectangle rule integration
for(j in 1:length(inf_data$ifng)){
  print(j)
  samplexv = ceiling(runif(1,-1,max(inf_data$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
  newx=dnbinom(samplexv, size=x1o, prob=x2o, log = FALSE)
  # print(newx)
  rv=runif(1,0,1)
  while(rv>newx){
    samplexv = ceiling(runif(1,-1,max(inf_data$ifng)))#inf_data$ifng[ceiling(runif(1,0,length(inf_data$ifng)))]
    newx=dnbinom(samplexv, size=x1o, prob=x2o, log = FALSE)
    rv=runif(1,0,1)
     #print(newx)
  }
  if(j==1){
    xlnN12=samplexv
  }
  else{
    xlnN12=append(xlnN12,samplexv)
  }
}


xNB2 <- data.frame(xln1NB2=xlnN12)
A12=ggplot() + # variant becomes a classifying factor
  geom_histogram(data=inf_data,aes(x=log(ifng +1), y=..count.. ), bins=40, alpha = 0.6 ,fill= "blue")+
  geom_histogram(data=xNB2,aes(x=log(xln1NB2+1), y=..count.. ), bins=40, alpha = 0.4 , fill= "gray",color="black")+ # add a scatterplot; constant size, shape/fill depends on variant # add a scatterplot; constant size, shape/fill depends on variant
  ylab("Number of cells") +
  xlab("log1p Interferon gamma counts") +
  ylim(0,1000)+
  ggtitle("Negative binomial model") +
  theme(axis.line = element_line(colour = "black"),
        #  legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

A13=ggplot() + # variant becomes a classifying factor
  geom_histogram(data=inf_data_all_3_6,aes(x=log(inf_data_all_3_6$ifng+1), y=..count.. ), bins=40, alpha = 0.6 ,fill= "blue")+
  #geom_histogram(data=xNB2,aes(x=log(xNB2$xln1NB2+1), y=..count.. ), bins=10, alpha = 0.4 , fill= "gray",color="black")+ # add a scatterplot; constant size, shape/fill depends on variant # add a scatterplot; constant size, shape/fill depends on variant
  ylab("Number of cells") +
  xlab("log1p Interferon gamma counts") +
  ggtitle("Negative binomial model") +
  theme(axis.line = element_line(colour = "black"),
        #  legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

pdf("IFNG_all.pdf",width=8,height=4)
grid.arrange(A13,A13, ncol = 2, nrow = 1, heights=c(0.5), widths=c(0.5,0.5))
dev.off()


pdf("Tetramer_RAD_allCells.pdf",width=8,height=4)
grid.arrange(A11,A12, ncol = 2, nrow = 1, heights=c(0.5), widths=c(0.5,0.5))
dev.off()



#############Make UPset plots
SummarydfmVST=read.csv("/Users/maewoodsphd/LegacymVSTmanuscriptCode/SupplementaryTables/SummaryIFNG_2_5%.csv")

ThreehrTimeUS=levels(factor(subset(SummarydfmVST,
                                   SummarydfmVST$X3hr.Total.cells.US>0)$Clone))
ThreehrTimeCMV=levels(factor(subset(SummarydfmVST,SummarydfmVST$X3hr.Total.cells.CMV>0)$Clone))
ThreehrTimeEBV=levels(factor(subset(SummarydfmVST,SummarydfmVST$X3hr.Total.cells.EBV>0)$Clone))

listInput <- list(`Control` = ThreehrTimeUS, 
                  `CMV` = ThreehrTimeCMV,
                  `EBV` = ThreehrTimeEBV)
pdf("../../Figures/Manuscript/PDFs/Figure3b3hour.pdf",width=5,height=4)
upset(fromList(listInput), order.by = "freq")
dev.off()
