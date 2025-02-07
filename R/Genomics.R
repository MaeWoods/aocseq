#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Code to pull the full length TCR sequence from 10X genomics all_contig_annotations.json file
#'
#' This function will import a set of single cell barcodes, extract the full length TCR transcript details
#' and output the details ranked by the cell with the greatest functional score.
#'
#' @param barcodes Barcodes used to identify cells in full contig annotations JSON file.
#' @param contig.annotations contig_annotations.csv file.
#' @param json.path Directory where JSON file is saved.
#' @param save.dir Directory where full length TCR transcripts will be saved.
#' @param score Orders the cells be classification score.
#' @param ranked.sheet False. If true implements and saves the cells by their characteristic rank.
#' @param verbose Print progress bars and output
#' @concept Genomics
#'
#' @export
GetTCRs <- function(
    barcodes,
    contig.annotations,
    json.path=".",
    save.dir=".",
    score="",
    ranked.sheet=FALSE,
    verbose = TRUE
){

  result <- fromJSON(file = json.path)
  for(s in 1:length(barcodes)){
    CHR1 <- lapply(result, function(x) { x$barcode == barcodes[s] })

    len.cells=dim(contig.annotations)[1]
    indicies=0
    for(k in 1:len.cells){
      if(length(CHR1[[k]])==0){

      }
      else{

        if(CHR1[k]==TRUE){
          indicies=append(indicies,k)
        }
      }
    }
    indicies=indicies[-1]
  length.fragments=length(indicies)

  frag.number = vector(mode = "list", length = length.fragments)
    for(k in 1:length.fragments){
      frag.number[[k]]=result[[indicies[k]]]
      capture.output(frag.number[[k]], file = paste(save.dir,paste(paste(paste("/Cell",s,sep=""),"_Chain",sep=""),k,sep=""),".txt",sep=""), append = TRUE)
    }

  }

}

#' This function will import a set of single cell barcodes, extract the full length TCR transcript details
#' and output the details ranked by the cell with the greatest functional score.
#'
#' @param cell.data single cell RNA sequencing data processed by aocseq.
#' @param gene gene list for UMI ranked cell types.
#' @param meta.data select a metadata column to rank cell types.
#' @param data.name name of metadata used for ranking.
#' @param threshold name of metadata used for ranking.
#' @param save.dir Path to the location where the output tables of ranked TCRs will be stored.
#' @param verbose Print progress bars and output.
#' @concept Genomics
#'
#' @export
RankTCRs <- function(
    cell.data,
    gene="IFNG",
    meta.data=FALSE,
    data.name="Mdist",
    threshold=c(-1,-1),
    save.dir=".",
    verbose = TRUE
){
  gene_id=0
  RnaStoreUMO=matrix(data=NA,nrow=dim(cell.data[['RNA']])[1],ncol=dim(cell.data[['RNA']])[2])
  if(packageVersion("Seurat")<'5.0.0'){
    RnaStoreUMO=cell.data[['RNA']]$counts
    gene_id=match(gene,row.names(cell.data[['RNA']]$counts))
  }
  else{
    RnaStoreUMO=cell.data[['RNA']]@counts
    gene_id=match(gene,row.names(cell.data[['RNA']]@counts))
  }
    clones=levels(factor(cell.data@meta.data$cdr3_na))
    clones.as.cells=cell.data@meta.data$cdr3_na
    size.clones=rep(0,length(clones))
    for(int k in 1:length(clones)){
      size.clones[k]=length(which(clones.as.cells==clones[k]))
    }
    out_df_pre=data.frame(barcodes=rep(0,length(clones)),
                          cdr3=clones,
                          clone.size=size.clones,
                          max_ele=rep(0,length(clones)),
                          mean_ele=rep(0,length(clones)))
    cell.meta=cell.data@meta.data
    subject.meta.array=rep(0,length(cell.data@meta.data[[1,]))
    if(meta.data==TRUE){
      pos=match(data.name,names(cell.meta))
      subject.meta.array=cell.meta[pos,]
    }
    else{
      subject.meta.array=RnaStoreUMO[gene_id,]
    }
    df_t=data.frame(cac=clones.as.cells,val=subject.meta.array)
    for(j in 1:length(clones)){
      out_df_pre$max_ele[j]=max(subset(cell.meta,cac==clones[j])$val)
      out_df_pre$mean_ele[j]=mean(subset(cell.meta,cac==clones[j])$val)
    }
    orderbyfactor=order(out_df_pre$max_ele,decreasing = TRUE)
    out_df=data.frame(barcodes=out_df_pre$barcodes[orderbyfactor],
                      clonotypes=out_df_pre$clones[orderbyfactor],
                      cdr3=out_df_pre$cdr3[orderbyfactor],
                      max_ele=out_df_pre$max_ele[orderbyfactor],
                      mean_ele=out_df_pre$mean_ele[orderbyfactor])
    write.csv(out_df,paste(save.dir,"/tcrs_max.csv",sep=""))
    orderbyfactor=order(out_df_pre$mean_ele,decreasing = TRUE)
    out_df=data.frame(barcodes=out_df_pre$barcodes[orderbyfactor],
                      clonotypes=out_df_pre$clones[orderbyfactor],
                      cdr3=out_df_pre$cdr3[orderbyfactor],
                      max_ele=out_df_pre$max_ele[orderbyfactor],
                      mean_ele=out_df_pre$mean_ele[orderbyfactor])
    write.csv(out_df,paste(save.dir,"/tcrs_avg.csv",sep=""))
}

#' This function will import a set of single cell barcodes, extract the full length TCR transcript details
#' and output the details ranked by the cell with the greatest functional score.
#'
#' @param n.cells Number of cells processed.
#' @param check_count Path to output files from alignr::GetTCRs
#' @param path Path to output files from alignr::GetTCRs
#' @param save.dir Directory where full length receptor sequences will be saved.
#' @param TRB.for Forward primer sequence for TCR beta chain
#' @param TRB.rev Reverse primer sequence for TCR beta chain
#' @param TRA.for Forward primer sequence for TCR alpha chain
#' @param TRA.rev Reverse primer sequence for TCR alpha chain
#' @param verbose Print progress bars and output
#' @concept Genomics
#'
#' @export
MakeReceptor <- function(
    n.cells,
    check_count,
    path=".",
    save.dir=".",
    TRB.for="A",
    TRB.rev="A",
    TRA.for="A",
    TRA.rev="A"
){

  TRAC="XIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSN
SAVAWSNKSDFACANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGF
RILLLKVAGFNLLMTLRLWSS"

  TRBC2="XDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVSTDP
QPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQI
VSAEAWGRADCGFTSESYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDSRG"

  TRBC1="XDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVSTDP
QPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQI
VSAEAWGRADCGFTSVSYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDF"

  curl_download("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
                paste(path,"HUMAN_9606_idmapping.dat.gz",sep=""),mode="w",quiet=FALSE)

  for(j in 1:n.cells){

    flankTRAC=substring(TRAC,2,12)
    flankTRBC1=substring(TRBC1,2,12)
    flankTRBC2=substring(TRBC2,2,12)
    FlankName="TBC"
    n.chains=list.files(path, pattern=paste(paste("Cell",j,sep=""),"_",sep=""))

    for(k in n.chains){
      JSON_output=read.csv(paste(path,k,sep=""),header=FALSE)
      names(JSON_output) <- "X"
      if(length(grep("TRA",t(JSON_output$X)))>0){
        ###This means the chain is an alpha chain
        Flank=flankTRAC
        FlankName="TCRAlpha"
      }else{
        if(length(grep("TRBC1",t(JSON_output$X)))>0){
          Flank=flankTRBC1
          FlankName="TCRBetaC1"
        }
        else if(length(grep("TRBC2",t(JSON_output$X)))>0){
          Flank=flankTRBC2
          FlankName="TCRBetaC2"
        }
      }
        ###This means the chain is an alpha chain
        scheck=length(grep("TRUE",t(JSON_output$X)[grep("high_confidence",t(JSON_output$X))+1]))+
          length(grep("TRUE",t(JSON_output$X)[grep("is_cell",t(JSON_output$X))+1]))+
          length(grep("TRUE",t(JSON_output$X)[grep("productive",t(JSON_output$X))+1]))+
          length(grep("TRUE",t(JSON_output$X)[grep("filtered",t(JSON_output$X))+1]))+
          length(grep("TRUE",t(JSON_output$X)[grep("is_gex_cell",t(JSON_output$X))+1]))+
          length(grep("TRUE",t(JSON_output$X)[grep("is_asm_cell",t(JSON_output$X))+1]))+
          length(grep("TRUE",t(JSON_output$X)[grep("full_length",t(JSON_output$X))+1]))

        if(scheck==check_count){
          ##Only proceed if a high confidence cell and productive receptor
          gene.name=strsplit(t(JSON_output$X)[grep("gene_name",t(JSON_output$X))+1][1],split = " ")[[1]][2]
          if(grepl("/",gene.name)){
            gene.name=paste(strsplit(gene.name,split="/")[[1]][1],strsplit(gene.name,split="/")[[1]][2],sep="")
          }
          locusID = match(gene.name,read.table(gzfile("HUMAN_9606_idmapping.dat.gz"), fill = TRUE)$V3)
          PiD = read.table(gzfile("HUMAN_9606_idmapping.dat.gz"), fill = TRUE)$V1[locusID]
          curl_download(paste(paste("https://rest.uniprot.org/uniprotkb/",PiD,sep=""),
                              ".fasta",sep=""),paste(gene.name,".txt",sep=""),mode="w",quiet=FALSE)
          genefasta=read.csv(paste(gene.name,".txt",sep=""))
          names(genefasta) <- "X"
          genesubstring=substring(genefasta$X[1],1,40)
          aasequence=strsplit(t(JSON_output$X)[grep("aa_sequence",t(JSON_output$X))+1],split = " ")[[1]][2]

          ###If no match, select fewer amino acids
          s=1
          while(!grepl(genesubstring,aasequence)){
            genesubstring=substr(genesubstring,1,nchar(genesubstring)-s)
          }

          if(grepl(genesubstring,aasequence)){

            aasequencestep=substr(aasequence,2,nchar(aasequence))
            s=1
            while(grepl(genesubstring,substr(aasequencestep,s,nchar(aasequencestep)))){
              s=s+1
              aasequencestep=substr(aasequencestep,s,nchar(aasequencestep))
            }
            ###Remove excess from receptor head
            if(s>1){
              aasequence=substr(aasequence,s+1,nchar(aasequence))
            }

            #####Now focus on the r.h.s of the sequence
            s=0
            while(!grepl(substr(Flank,1,nchar(Flank)-s),aasequence)){
              s=s+1
            }
            if(s>0){
              Flank=substr(Flank,1,nchar(Flank)-s)
            }

            aasequencestep=aasequence
            s=1
            while(grepl(Flank,substr(aasequencestep,1,nchar(aasequencestep)-s))){
              s=s+1
              aasequencestep=substr(aasequencestep,1,nchar(aasequencestep)-s)
            }
            ###Remove excess from receptor tail
            if(s>1){
              aasequence=substr(aasequencestep,1,nchar(aasequencestep)-(s+nchar(Flank)))
            }


          }

          ###Save chain to same folder
          write.csv(aasequence,paste(save.dir,paste(FlankName,k),sep=""))
        }
      }

    }

}
