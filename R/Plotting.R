#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Functions to plot data stored in clonal objects
#' at normalized counts greater than the control
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
#' @param path
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
