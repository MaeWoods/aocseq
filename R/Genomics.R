#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Code to pull the full length TCR sequence from 10X genomics all_contig_annotations.json file
#'
#' This function will import a set of single cell barcodes, extract the full length TCR transcript details
#' and output the details ranked by the cell with the greatest functional score.
#'
#' @param barcodes Barcodes used to identify cells in full contig annotations JSON file.
#' @param contig.annoations contig_annotations.csv file.
#' @param json.path Directory where JSON file is saved.
#' @param save.dir Directory where full length TCR transcripts will be saved.
#' @param score Orders the cells be classification score. 
#' @param verbose Print progress bars and output
#' @return Void
#' @concept Genomics
#'
#' @export
GetTCR <- function(
    barcodes,
    contig.annoations,
    json.path=".",
    save.dir=".",
    score="",
    verbose = TRUE
){
  
  result <- fromJSON(file = json.path)
  for(s in 1:length(barcodes)){
    CHR1 <- lapply(result, function(x) { x$barcode == barcodes[s] })
    
    len.cells=dim(contig.annoations)[1]
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
      capture.output(frag.number[[k]], file = paste(save.dir,paste(paste(paste("Cell",s,sep=""),"_Chain",sep=""),k,sep=""),".txt",sep=""), append = TRUE)
    }
  
  }

}



