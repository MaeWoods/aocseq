saveRDS(ThreeHourStim,file="../RDS/ThreeHourStim.rds")
saveRDS(SixHourStim,file="../RDS/SixHourStim.rds")

ThreeHourStim=readRDS("../RDS/ThreeHourStim.rds")
SixHourStim=readRDS("../RDS/SixHourStim.rds")

markerGene=c("IFNG","CRTAM","CD69","TNF","CD70","TNFRSF9")

moi="TNFRSF9"
clonotype_path="../SupplementaryTables/SummaryTNFRSF9_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=6
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)


moi="CD70"
clonotype_path="../SupplementaryTables/SummaryCD70_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=5
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)

moi="TNF"
clonotype_path="../SupplementaryTables/SummaryTNF_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=4
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)


moi="CD69"
clonotype_path="../SupplementaryTables/SummaryCD69_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=3
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)


moi="CRTAM"
clonotype_path="../SupplementaryTables/SummaryCRTAM_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=2
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)


moi="IFNG"
clonotype_path="../SupplementaryTables/SummaryIFNG_2_5%.csv"
cell_path="../SupplementaryTables/Cells.csv"
moipos=1
AnnotateClonotypes(ThreeHourStim,SixHourStim,moi,clonotype_path,cell_path,moipos)


AnnotateClonotypes <- function(ThreeHour,SixHour,moi,clonotype_path,cell_path,moipos){
  
  #cell percentages t1
  unstim_t1 <- ThreeHour[[4]]
  moiunstim_t1=match(moi, row.names(unstim_t1[["SCT"]]@data))
  cmv_sct_t1 <- ThreeHour[[1]]
  moicmv_t1=match(moi, row.names(cmv_sct_t1[["SCT"]]@data))
  ebv_sct_t1 <- ThreeHour[[2]]
  moiebv_t1=match(moi, row.names(ebv_sct_t1[["SCT"]]@data))
  bkv_sct_t1 <- ThreeHour[[3]]
  moibkv_t1=match(moi, row.names(bkv_sct_t1[["SCT"]]@data))
  
  #cell percentages t2
  unstim_t2 <- SixHour[[4]]
  moiunstim_t2=match(moi, row.names(unstim_t2[["SCT"]]@data))
  cmv_sct_t2 <- SixHour[[1]]
  moicmv_t2=match(moi, row.names(cmv_sct_t2[["SCT"]]@data))
  ebv_sct_t2 <- SixHour[[2]]
  moiebv_t2=match(moi, row.names(ebv_sct_t2[["SCT"]]@data))
  bkv_sct_t2 <- SixHour[[3]]
  moibkv_t2=match(moi, row.names(bkv_sct_t2[["SCT"]]@data))
  
  cutoff_t1=quantile(unstim_t1[["SCT"]]@data[moiunstim_t1,],.975)[[1]]
  cutoff_t2=quantile(unstim_t2[["SCT"]]@data[moiunstim_t2,],.975)[[1]]
  
  allclonotypes = setdiff(levels(factor(union(cmv_sct_t2@meta.data$cdr3_na,union(ebv_sct_t2@meta.data$cdr3_na,union(unstim_t2@meta.data$cdr3_na,union(bkv_sct_t2@meta.data$cdr3_na,union(unstim_t1@meta.data$cdr3_na,union(cmv_sct_t1@meta.data$cdr3_na,union(bkv_sct_t1@meta.data$cdr3_na,ebv_sct_t1@meta.data$cdr3_na))))))))),"unassigned")
  
  intersectclonotypes = setdiff(levels(factor(intersect(intersect(intersect(intersect(bkv_sct_t1@meta.data$cdr3_na,
                                                                                      ebv_sct_t1@meta.data$cdr3_na),
                                                                            cmv_sct_t1@meta.data$cdr3_na),
                                                                  unstim_t1@meta.data$cdr3_na),intersect(intersect(intersect(bkv_sct_t2@meta.data$cdr3_na,
                                                                                                                             ebv_sct_t2@meta.data$cdr3_na),
                                                                                                                   cmv_sct_t2@meta.data$cdr3_na),
                                                                                                         unstim_t2@meta.data$cdr3_na)))),"unassigned")
  
  UMOclonotypes_t1=setdiff(levels(factor(unstim_t1@meta.data$cdr3_na)),"unassigned")
  EBVclonotypes_t1=setdiff(setdiff(levels(factor(ebv_sct_t1@meta.data$cdr3_na)),"unassigned"),UMOclonotypes_t1)
  BKVclonotypes_t1=setdiff(setdiff(levels(factor(bkv_sct_t1@meta.data$cdr3_na)),"unassigned"),union(EBVclonotypes_t1,UMOclonotypes_t1))
  CMVclonotypes_t1=setdiff(setdiff(levels(factor(cmv_sct_t1@meta.data$cdr3_na)),"unassigned"),union(BKVclonotypes_t1,union(EBVclonotypes_t1,UMOclonotypes_t1)))
  
  UMOclonotypes_t2=setdiff(setdiff(levels(factor(unstim_t2@meta.data$cdr3_na)),"unassigned"),union(CMVclonotypes_t1,
                                                                                                   union(BKVclonotypes_t1,union(EBVclonotypes_t1,UMOclonotypes_t1))))
  EBVclonotypes_t2=setdiff(setdiff(levels(factor(ebv_sct_t2@meta.data$cdr3_na)),"unassigned"),union(UMOclonotypes_t2,
                                                                                                    union(CMVclonotypes_t1,
                                                                                                          union(BKVclonotypes_t1,union(EBVclonotypes_t1,UMOclonotypes_t1)))))
  BKVclonotypes_t2=setdiff(setdiff(levels(factor(bkv_sct_t2@meta.data$cdr3_na)),"unassigned"),union(EBVclonotypes_t2,
                                                                                                    union(UMOclonotypes_t2,
                                                                                                          union(CMVclonotypes_t1,
                                                                                                                union(BKVclonotypes_t1,union(EBVclonotypes_t1,UMOclonotypes_t1))))))
  CMVclonotypes_t2=setdiff(setdiff(levels(factor(cmv_sct_t2@meta.data$cdr3_na)),"unassigned"),union(BKVclonotypes_t2,
                                                                                                    union(EBVclonotypes_t2,union(UMOclonotypes_t2,union(CMVclonotypes_t1,union(BKVclonotypes_t1,union(EBVclonotypes_t1,UMOclonotypes_t1)))))))
  
  
  #assign counts for ranking by UMO, then EBV, BKV - CMV
  print("Assigning counts for ranking...")
  NClonotypes=length(allclonotypes)
  Shared=rep(0,NClonotypes)
  Clonefreq=rep(0,NClonotypes)
  for(k in 1:NClonotypes){
    if(length(intersect(allclonotypes[k],UMOclonotypes_t1))>0){
      Clonefreq[k]=subset(unstim_t1,cdr3_na==allclonotypes[k])@meta.data$countcln[1]
      if(length(intersect(allclonotypes[k],intersectclonotypes))>0){
        Shared[k]="yes"
      }
      else{
        Shared[k]="no"
      }
    }
    else if(length(intersect(allclonotypes[k],EBVclonotypes_t1))>0){
      Clonefreq[k]=subset(ebv_sct_t1,cdr3_na==allclonotypes[k])@meta.data$countcln[1]
      if(length(intersect(allclonotypes[k],intersectclonotypes))>0){
        Shared[k]="yes"
      }
      else{
        Shared[k]="no"
      }
    }
    else if(length(intersect(allclonotypes[k],BKVclonotypes_t1))>0){
      Clonefreq[k]=subset(bkv_sct_t1,cdr3_na==allclonotypes[k])@meta.data$countcln[1]
      if(length(intersect(allclonotypes[k],intersectclonotypes))>0){
        Shared[k]="yes"
      }
      else{
        Shared[k]="no"
      }
    }
    else if(length(intersect(allclonotypes[k],CMVclonotypes_t1))>0){
      Clonefreq[k]=subset(cmv_sct_t1,cdr3_na==allclonotypes[k])@meta.data$countcln[1]
      if(length(intersect(allclonotypes[k],intersectclonotypes))>0){
        Shared[k]="yes"
      }
      else{
        Shared[k]="no"
      }
    }
    else if(length(intersect(allclonotypes[k],UMOclonotypes_t2))>0){
      Clonefreq[k]=subset(unstim_t2,cdr3_na==allclonotypes[k])@meta.data$countcln[1]
      if(length(intersect(allclonotypes[k],intersectclonotypes))>0){
        Shared[k]="yes"
      }
      else{
        Shared[k]="no"
      }
    }
    else if(length(intersect(allclonotypes[k],EBVclonotypes_t2))>0){
      Clonefreq[k]=subset(ebv_sct_t2,cdr3_na==allclonotypes[k])@meta.data$countcln[1]
      if(length(intersect(allclonotypes[k],intersectclonotypes))>0){
        Shared[k]="yes"
      }
      else{
        Shared[k]="no"
      }
    }
    else if(length(intersect(allclonotypes[k],BKVclonotypes_t2))>0){
      Clonefreq[k]=subset(bkv_sct_t2,cdr3_na==allclonotypes[k])@meta.data$countcln[1]
      if(length(intersect(allclonotypes[k],intersectclonotypes))>0){
        Shared[k]="yes"
      }
      else{
        Shared[k]="no"
      }
    }
    else if(length(intersect(allclonotypes[k],CMVclonotypes_t2))>0){
      Clonefreq[k]=subset(cmv_sct_t2,cdr3_na==allclonotypes[k])@meta.data$countcln[1]
      if(length(intersect(allclonotypes[k],intersectclonotypes))>0){
        Shared[k]="yes"
      }
      else{
        Shared[k]="no"
      }
    }
  }
  
  Clonotypes_df = data.frame(cdr3=allclonotypes,shared=Shared,frequency=Clonefreq)
  Clonotypes_data = data.frame(cdr3=Clonotypes_df$cdr3[order(Clonotypes_df$frequency,decreasing=TRUE)],
                               shared=Clonotypes_df$shared[order(Clonotypes_df$frequency,decreasing=TRUE)],
                               frequency=Clonotypes_df$frequency[order(Clonotypes_df$frequency,decreasing=TRUE)])
  
  Clone=rep(0,length(allclonotypes))
  CD4CMV_t1=rep(0,length(allclonotypes))
  CD4EBV_t1=rep(0,length(allclonotypes))
  CD4BKV_t1=rep(0,length(allclonotypes))
  CD4US_t1=rep(0,length(allclonotypes))
  CD8CMV_t1=rep(0,length(allclonotypes))
  CD8EBV_t1=rep(0,length(allclonotypes))
  CD8BKV_t1=rep(0,length(allclonotypes))
  CD8US_t1=rep(0,length(allclonotypes))
  AllCMV_t1=rep(0,length(allclonotypes))
  AllEBV_t1=rep(0,length(allclonotypes))
  AllBKV_t1=rep(0,length(allclonotypes))
  AllUS_t1=rep(0,length(allclonotypes))
  TotalBKV_t1=rep(0,length(allclonotypes))
  TotalCMV_t1=rep(0,length(allclonotypes))
  TotalEBV_t1=rep(0,length(allclonotypes))
  Totalunstim_t1=rep(0,length(allclonotypes))
  
  CD4CMV_t2=rep(0,length(allclonotypes))
  CD4EBV_t2=rep(0,length(allclonotypes))
  CD4BKV_t2=rep(0,length(allclonotypes))
  CD4US_t2=rep(0,length(allclonotypes))
  CD8CMV_t2=rep(0,length(allclonotypes))
  CD8EBV_t2=rep(0,length(allclonotypes))
  CD8BKV_t2=rep(0,length(allclonotypes))
  CD8US_t2=rep(0,length(allclonotypes))
  AllCMV_t2=rep(0,length(allclonotypes))
  AllEBV_t2=rep(0,length(allclonotypes))
  AllBKV_t2=rep(0,length(allclonotypes))
  AllUS_t2=rep(0,length(allclonotypes))
  TotalBKV_t2=rep(0,length(allclonotypes))
  TotalCMV_t2=rep(0,length(allclonotypes))
  TotalEBV_t2=rep(0,length(allclonotypes))
  Totalunstim_t2=rep(0,length(allclonotypes))
  
  
  print("Assigning proportional specificity...")
  TCRsunstimT1=unstim_t1@meta.data$cdr3_na
  TCRcmvT1=cmv_sct_t1@meta.data$cdr3_na
  TCRebvT1=ebv_sct_t1@meta.data$cdr3_na
  TCRbkvT1=bkv_sct_t1@meta.data$cdr3_na
  TCRsunstimT2=unstim_t2@meta.data$cdr3_na
  TCRcmvT2=cmv_sct_t2@meta.data$cdr3_na
  TCRebvT2=ebv_sct_t2@meta.data$cdr3_na
  TCRbkvT2=bkv_sct_t2@meta.data$cdr3_na
  CD4dfunstim_t1=data.frame(clonotype=subset(unstim_t1,CD4cells=="1")@meta.data$cdr3_na,mark=subset(unstim_t1,CD4cells=="1")@meta.data[[8+moipos]])
  CD4dfcmv_t1=data.frame(clonotype=subset(cmv_sct_t1,CD4cells=="1")@meta.data$cdr3_na,mark=subset(cmv_sct_t1,CD4cells=="1")@meta.data[[8+moipos]])
  CD4dfebv_t1=data.frame(clonotype=subset(ebv_sct_t1,CD4cells=="1")@meta.data$cdr3_na,mark=subset(ebv_sct_t1,CD4cells=="1")@meta.data[[8+moipos]])
  CD4dfbkv_t1=data.frame(clonotype=subset(bkv_sct_t1,CD4cells=="1")@meta.data$cdr3_na,mark=subset(bkv_sct_t1,CD4cells=="1")@meta.data[[8+moipos]])
  
  CD8dfunstim_t1=data.frame(clonotype=subset(unstim_t1,CD8cells=="1")@meta.data$cdr3_na,mark=subset(unstim_t1,CD8cells=="1")@meta.data[[8+moipos]])
  CD8dfcmv_t1=data.frame(clonotype=subset(cmv_sct_t1,CD8cells=="1")@meta.data$cdr3_na,mark=subset(cmv_sct_t1,CD8cells=="1")@meta.data[[8+moipos]])
  CD8dfebv_t1=data.frame(clonotype=subset(ebv_sct_t1,CD8cells=="1")@meta.data$cdr3_na,mark=subset(ebv_sct_t1,CD8cells=="1")@meta.data[[8+moipos]])
  CD8dfbkv_t1=data.frame(clonotype=subset(bkv_sct_t1,CD8cells=="1")@meta.data$cdr3_na,mark=subset(bkv_sct_t1,CD8cells=="1")@meta.data[[8+moipos]])
  
  dfunstim_t1=data.frame(clonotype=unstim_t1@meta.data$cdr3_na,mark=unstim_t1@meta.data[[8+moipos]])
  dfcmv_t1=data.frame(clonotype=cmv_sct_t1@meta.data$cdr3_na,mark=cmv_sct_t1@meta.data[[8+moipos]])
  dfebv_t1=data.frame(clonotype=ebv_sct_t1@meta.data$cdr3_na,mark=ebv_sct_t1@meta.data[[8+moipos]])
  dfbkv_t1=data.frame(clonotype=bkv_sct_t1@meta.data$cdr3_na,mark=bkv_sct_t1@meta.data[[8+moipos]])
  
  
  CD4dfunstim_t2=data.frame(clonotype=subset(unstim_t2,CD4cells=="1")@meta.data$cdr3_na,mark=subset(unstim_t2,CD4cells=="1")@meta.data[[8+moipos]])
  CD4dfcmv_t2=data.frame(clonotype=subset(cmv_sct_t2,CD4cells=="1")@meta.data$cdr3_na,mark=subset(cmv_sct_t2,CD4cells=="1")@meta.data[[8+moipos]])
  CD4dfebv_t2=data.frame(clonotype=subset(ebv_sct_t2,CD4cells=="1")@meta.data$cdr3_na,mark=subset(ebv_sct_t2,CD4cells=="1")@meta.data[[8+moipos]])
  CD4dfbkv_t2=data.frame(clonotype=subset(bkv_sct_t2,CD4cells=="1")@meta.data$cdr3_na,mark=subset(bkv_sct_t2,CD4cells=="1")@meta.data[[8+moipos]])
  
  CD8dfunstim_t2=data.frame(clonotype=subset(unstim_t2,CD8cells=="1")@meta.data$cdr3_na,mark=subset(unstim_t2,CD8cells=="1")@meta.data[[8+moipos]])
  CD8dfcmv_t2=data.frame(clonotype=subset(cmv_sct_t2,CD8cells=="1")@meta.data$cdr3_na,mark=subset(cmv_sct_t2,CD8cells=="1")@meta.data[[8+moipos]])
  CD8dfebv_t2=data.frame(clonotype=subset(ebv_sct_t2,CD8cells=="1")@meta.data$cdr3_na,mark=subset(ebv_sct_t2,CD8cells=="1")@meta.data[[8+moipos]])
  CD8dfbkv_t2=data.frame(clonotype=subset(bkv_sct_t2,CD8cells=="1")@meta.data$cdr3_na,mark=subset(bkv_sct_t2,CD8cells=="1")@meta.data[[8+moipos]])
  
  dfunstim_t2=data.frame(clonotype=unstim_t2@meta.data$cdr3_na,mark=unstim_t2@meta.data[[8+moipos]])
  dfcmv_t2=data.frame(clonotype=cmv_sct_t2@meta.data$cdr3_na,mark=cmv_sct_t2@meta.data[[8+moipos]])
  dfebv_t2=data.frame(clonotype=ebv_sct_t2@meta.data$cdr3_na,mark=ebv_sct_t2@meta.data[[8+moipos]])
  dfbkv_t2=data.frame(clonotype=bkv_sct_t2@meta.data$cdr3_na,mark=bkv_sct_t2@meta.data[[8+moipos]])
  
  
  for(k in 1:NClonotypes){
    Clone[k]=Clonotypes_data$cdr3[k]
    Ncellstotunstim=0
    Ncellsaboveunstim=0
    Ncellstotcmv_sct=0
    Ncellsabovecmv_sct=0
    Ncellstotebv_sct=0
    Ncellsaboveebv_sct=0
    Ncellstotbkv_sct=0
    Ncellsabovebkv_sct=0
    tempvec=subset(CD4dfunstim_t1,clonotype==Clonotypes_data$cdr3[k])$mark
    tempvectot=subset(dfunstim_t1,clonotype==Clonotypes_data$cdr3[k])$mark
    if(length(intersect(Clonotypes_data$cdr3[k],TCRsunstimT1))>0){
      Ncellstotunstim=length(tempvectot)
      Ncellsaboveunstim=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD4US_t1[k]=0
      AllUS_t1[k]=0
      if(Ncellstotunstim>0){
        
          CD4US_t1[k]=(Ncellsaboveunstim/Ncellstotunstim)*100
          AllUS_t1[k]=(Allaboveunstim/Ncellstotunstim)*100
        
      }
    }
    
    tempvec=subset(CD8dfunstim_t1,clonotype==Clonotypes_data$cdr3[k])$mark
    tempvectot=subset(dfunstim_t1,clonotype==Clonotypes_data$cdr3[k])$mark
    if(length(intersect(Clonotypes_data$cdr3[k],TCRsunstimT1))>0){
      Ncellstotunstim=length(tempvectot)
      Ncellsaboveunstim=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD8US_t1[k]=0
      if(Ncellstotunstim>0){
        
          CD8US_t1[k]=(Ncellsaboveunstim/Ncellstotunstim)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRcmvT1))>0){
      tempvec=subset(CD4dfcmv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfcmv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotcmv_sct=length(tempvectot)
      Ncellsabovecmv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD4CMV_t1[k]=0
      AllCMV_t1[k]=0
      if(Ncellstotcmv_sct>0){
        
          CD4CMV_t1[k]=(Ncellsabovecmv_sct/Ncellstotcmv_sct)*100
          AllCMV_t1[k]=(Allaboveunstim/Ncellstotcmv_sct)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRcmvT1))>0){
      tempvec=subset(CD8dfcmv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfcmv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotcmv_sct=length(tempvectot)
      Ncellsabovecmv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      CD8CMV_t1[k]=0
      if(Ncellstotcmv_sct>0){
        
          CD8CMV_t1[k]=(Ncellsabovecmv_sct/Ncellstotcmv_sct)*100
        
      }
    }
    
    
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRebvT1))>0){
      tempvec=subset(CD4dfebv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfebv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotebv_sct=length(tempvectot)
      Ncellsaboveebv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD4EBV_t1[k]=0
      AllEBV_t1[k]=0
      if(Ncellstotebv_sct>0){
       
          CD4EBV_t1[k]=(Ncellsaboveebv_sct/Ncellstotebv_sct)*100
          AllEBV_t1[k]=(Allaboveunstim/Ncellstotebv_sct)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRebvT1))>0){
      tempvec=subset(CD8dfebv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfebv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotebv_sct=length(tempvectot)
      Ncellsaboveebv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      CD8EBV_t1[k]=0
      if(Ncellstotebv_sct>0){
       
          CD8EBV_t1[k]=(Ncellsaboveebv_sct/Ncellstotebv_sct)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRbkvT1))>0){
      tempvec=subset(CD4dfbkv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfbkv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotbkv_sct=length(tempvectot)
      Ncellsabovebkv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD4BKV_t1[k]=0
      AllBKV_t1[k]=0
      if(Ncellstotbkv_sct>0){
          CD4BKV_t1[k]=(Ncellsabovebkv_sct/Ncellstotbkv_sct)*100
          AllBKV_t1[k]=(Allaboveunstim/Ncellstotbkv_sct)*100
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRbkvT1))>0){
      tempvec=subset(CD8dfbkv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfbkv_t1,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotbkv_sct=length(tempvectot)
      Ncellsabovebkv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      CD8BKV_t1[k]=0
      if(Ncellstotbkv_sct>0){
        
          CD8BKV_t1[k]=(Ncellsabovebkv_sct/Ncellstotbkv_sct)*100
        
      }
    }
    
    TotalBKV_t1[k] = Ncellstotbkv_sct
    TotalCMV_t1[k]=Ncellstotcmv_sct
    TotalEBV_t1[k]=Ncellstotebv_sct
    Totalunstim_t1[k]=Ncellstotunstim
    
    
    Ncellstotunstim=0
    Ncellsaboveunstim=0
    Ncellstotcmv_sct=0
    Ncellsabovecmv_sct=0
    Ncellstotebv_sct=0
    Ncellsaboveebv_sct=0
    Ncellstotbkv_sct=0
    Ncellsabovebkv_sct=0
    if(length(intersect(Clonotypes_data$cdr3[k],TCRsunstimT2))>0){
      tempvec=subset(CD4dfunstim_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfunstim_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotunstim=length(tempvectot)
      Ncellsaboveunstim=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD4US_t2[k]=0
      AllUS_t2[k]=0
      if(Ncellstotunstim>0){
        
          CD4US_t2[k]=(Ncellsaboveunstim/Ncellstotunstim)*100
          AllUS_t2[k]=(Allaboveunstim/Ncellstotunstim)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRsunstimT2))>0){
      tempvec=subset(CD8dfunstim_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfunstim_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotunstim=length(tempvectot)
      Ncellsaboveunstim=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      CD8US_t2[k]=0
      if(Ncellstotunstim>0){
        
          CD8US_t2[k]=(Ncellsaboveunstim/Ncellstotunstim)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRcmvT2))>0){
      tempvec=subset(CD4dfcmv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfcmv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotcmv_sct=length(tempvectot)
      Ncellsabovecmv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD4CMV_t2[k]=0
      AllCMV_t2[k]=0
      if(Ncellstotcmv_sct>0){
        
          CD4CMV_t2[k]=(Ncellsabovecmv_sct/Ncellstotcmv_sct)*100
          AllCMV_t2[k]=(Allaboveunstim/Ncellstotcmv_sct)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRcmvT2))>0){
      tempvec=subset(CD8dfcmv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfcmv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotcmv_sct=length(tempvectot)
      Ncellsabovecmv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      CD8CMV_t2[k]=0
      if(Ncellstotcmv_sct>0){
        
          CD8CMV_t2[k]=(Ncellsabovecmv_sct/Ncellstotcmv_sct)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRebvT2))>0){
      tempvec=subset(CD4dfebv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfebv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotebv_sct=length(tempvectot)
      Ncellsaboveebv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD4EBV_t2[k]=0
      AllEBV_t2[k]=0
      if(Ncellstotebv_sct>0){
        
          CD4EBV_t2[k]=(Ncellsaboveebv_sct/Ncellstotebv_sct)*100
          AllEBV_t2[k]=(Allaboveunstim/Ncellstotebv_sct)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRebvT2))>0){
      tempvec=subset(CD8dfebv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfebv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotebv_sct=length(tempvectot)
      Ncellsaboveebv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      CD8EBV_t2[k]=0
      if(Ncellstotebv_sct>0){
       
          CD8EBV_t2[k]=(Ncellsaboveebv_sct/Ncellstotebv_sct)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRbkvT2))>0){
      tempvec=subset(CD4dfbkv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfbkv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotbkv_sct=length(tempvectot)
      Ncellsabovebkv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      Allaboveunstim=length(tempvectot)-length(subset(tempvectot,tempvectot=="unassigned"))
      CD4BKV_t2[k]=0
      AllBKV_t2[k]=0
      if(Ncellstotbkv_sct>0){
        
          CD4BKV_t2[k]=(Ncellsabovebkv_sct/Ncellstotbkv_sct)*100
          AllBKV_t2[k]=(Allaboveunstim/Ncellstotbkv_sct)*100
        
      }
    }
    
    if(length(intersect(Clonotypes_data$cdr3[k],TCRbkvT2))>0){
      tempvec=subset(CD8dfbkv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      tempvectot=subset(dfbkv_t2,clonotype==Clonotypes_data$cdr3[k])$mark
      Ncellstotbkv_sct=length(tempvectot)
      Ncellsabovebkv_sct=length(tempvec)-length(subset(tempvec,tempvec=="unassigned"))
      CD8BKV_t2[k]=0
      if(Ncellstotbkv_sct>0){
        
          CD8BKV_t2[k]=(Ncellsabovebkv_sct/Ncellstotbkv_sct)*100
        
      }
    }
    
    TotalBKV_t2[k] = Ncellstotbkv_sct
    TotalCMV_t2[k]=Ncellstotcmv_sct
    TotalEBV_t2[k]=Ncellstotebv_sct
    Totalunstim_t2[k]=Ncellstotunstim
    
  }
  
  Summarydf=data.frame(Clonotype=Clone,CD4CMV3hr=CD4CMV_t1,CD8CMV3hr=CD8CMV_t1,AllCMV3hr=AllCMV_t1,CD4EBV3hr=CD4EBV_t1,CD8EBV3hr=CD8EBV_t1,AllEBV3hr=AllEBV_t1,CD4BKV3hr=CD4BKV_t1,
                       CD8BKV3hr=CD8BKV_t1,AllBKV3hr=AllBKV_t1,CD4US3hr=CD4US_t1,CD8US3hr=CD8US_t1,AllUS3hr=AllUS_t1,TotalcellsCMV3hr=TotalCMV_t1,
                       TotalcellsEBV3hr=TotalEBV_t1,TotalcellsBKV3hr=TotalBKV_t1,Totalcellsunstim3hr=Totalunstim_t1,
                       CD4CMV6hr=CD4CMV_t2,CD8CMV6hr=CD8CMV_t2,AllCMV6hr=AllCMV_t2,CD4EBV6hr=CD4EBV_t2,CD8EBV6hr=CD8EBV_t2,AllEBV6hr=AllEBV_t2,
                       CD4BKV6hr=CD4BKV_t2,CD8BKV6hr=CD8BKV_t2,AllBKV6hr=AllBKV_t2,CD4US6hr=CD4US_t2,CD8US6hr=CD8US_t2,AllUS6hr=AllUS_t2,
                       TotalcellsCMV6hr=TotalCMV_t2,
                       TotalcellsEBV6hr=TotalEBV_t2,TotalcellsBKV6hr=TotalBKV_t2,Totalcellsunstim6hr=Totalunstim_t2)
  
  Amino_Acid_seq=rep(0,length(Summarydf$Clonotype))
  
  AAsEBVT1 = ebv_sct_t1@meta.data$cdr3
  AAsEBVT2 = ebv_sct_t2@meta.data$cdr3
  AAsCMVT1 = cmv_sct_t1@meta.data$cdr3
  AAsCMVT2 = cmv_sct_t2@meta.data$cdr3
  AAsBKVT1 = bkv_sct_t1@meta.data$cdr3
  AAsBKVT2 = bkv_sct_t2@meta.data$cdr3
  AAsUnstimT1 = unstim_t1@meta.data$cdr3
  AAsUnstimT2 = unstim_t2@meta.data$cdr3
  
  NTsEBVT1 = ebv_sct_t1@meta.data$cdr3_na
  NTsEBVT2 = ebv_sct_t2@meta.data$cdr3_na
  NTsCMVT1 = cmv_sct_t1@meta.data$cdr3_na
  NTsCMVT2 = cmv_sct_t2@meta.data$cdr3_na
  NTsBKVT1 = bkv_sct_t1@meta.data$cdr3_na
  NTsBKVT2 = bkv_sct_t2@meta.data$cdr3_na
  NTsUnstimT1 = unstim_t1@meta.data$cdr3_na
  NTsUnstimT2 = unstim_t2@meta.data$cdr3_na
  
  for(q in 1:length(Amino_Acid_seq)){
    if(length(intersect(NTsEBVT1,Summarydf$Clonotype[q]))>0){
      idx=match(Summarydf$Clonotype[q],NTsEBVT1)
      Amino_Acid_seq[q]=AAsEBVT1[idx]
    }
    else if(length(intersect(NTsEBVT2,Summarydf$Clonotype[q]))>0){
      idx=match(Summarydf$Clonotype[q],NTsEBVT2)
      Amino_Acid_seq[q]=AAsEBVT2[idx]
    }
    else if(length(intersect(NTsCMVT1,Summarydf$Clonotype[q]))>0){
      idx=match(Summarydf$Clonotype[q],NTsCMVT1)
      Amino_Acid_seq[q]=AAsCMVT1[idx]
    }
    else if(length(intersect(NTsCMVT2,Summarydf$Clonotype[q]))>0){
      idx=match(Summarydf$Clonotype[q],NTsCMVT2)
      Amino_Acid_seq[q]=AAsCMVT2[idx]
    }
    else if(length(intersect(NTsBKVT1,Summarydf$Clonotype[q]))>0){
      idx=match(Summarydf$Clonotype[q],NTsBKVT1)
      Amino_Acid_seq[q]=AAsBKVT1[idx]
    }
    else if(length(intersect(NTsBKVT2,Summarydf$Clonotype[q]))>0){
      idx=match(Summarydf$Clonotype[q],NTsBKVT2)
      Amino_Acid_seq[q]=AAsBKVT2[idx]
    }
    else if(length(intersect(NTsUnstimT1,Summarydf$Clonotype[q]))>0){
      idx=match(Summarydf$Clonotype[q],NTsUnstimT1)
      Amino_Acid_seq[q]=AAsUnstimT1[idx]
    }
    else if(length(intersect(NTsUnstimT2,Summarydf$Clonotype[q]))>0){
      idx=match(Summarydf$Clonotype[q],NTsUnstimT2)
      Amino_Acid_seq[q]=AAsUnstimT2[idx]
    }
    else{
      print("Error no matching TCR found")
    }
  }
  
  names(Summarydf) <- c("Clone","3hr CD4 high CMV (%)","3hr CD8 high CMV (%)","3hr cells high CMV (%)"
                        ,"3hr CD4 high EBV (%)","3hr CD8 high EBV (%)","3hr cells high EBV (%)","3hr CD4 high BKV (%)",
                        "3hr CD8 high BKV (%)","3hr cells high BKV (%)", "3hr CD4 high US (%)",
                        "3hr CD8 high US (%)","3hr cells high US (%)","3hr Total cells CMV","3hr Total cells EBV",
                        "3hr Total cells BKV","3hr Total cells US",
                        "6hr CD4 high CMV (%)","6hr CD8 high CMV (%)","6hr cells high CMV (%)",
                        "6hr CD4 high EBV (%)","6hr CD8 high EBV (%)","6hr cells high EBV (%)","6hr CD4 high BKV (%)",
                        "6hr CD8 high BKV (%)","6hr cells high BKV (%)","6hr CD4 high US (%)",
                        "6hr CD8 high US (%)","6hr cells high US (%)","6hr Total cells CMV","6hr Total cells EBV",
                        "6hr Total cells BKV","6hr Total cells US")
  
  Summarydf$Amino_acid=Amino_Acid_seq
  write.csv(Summarydf,clonotype_path)
  
  
  
}
