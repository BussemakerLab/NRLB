#' Plot affinities of a sequence
#' 
#' @import ggplot2
#' 
#' @param genomicSequence
#' @param fits
#' @param index
#' @param mode
#' @param rc
#' @param nPeaks
#' @param annotate
#' @param rescale
#' @param genomicSequence2
#' @return
#' 
#' @examples
#'
#' @export
#' 
nrlb.plot.score.genome = function(genomicSequence, fits, index, mode=NULL, rc=FALSE, nPeaks=NULL, annotate=FALSE, rescale=NULL, 
                             genomicSequence2=NULL) {
  seq.string = deparse(substitute(genomicSequence))
  fit.string = deparse(substitute(fits))
  score = score.genome(genomicSequence, fits, index, mode, rc)
  if (fits[[1]]$k[index]=="Multi") {
    k = length(fits[[2]][[index]][[mode]]$NB)/4
  } else {
    k = as.numeric(fits[[1]]$k[index])
  }
  if (is.null(rescale)) {
    rescale = max.seq(fits, index, mode)$MaxAffinity
  }
  score = score/rescale
  max.score = max(score)
  #compare two tracks
  if (!is.null(genomicSequence2)) {
    if (length(genomicSequence2)!=length(genomicSequence)) {
      stop("Input genomic sequences are not of the same length.")
    }
    score2 = score.genome(genomicSequence2, fits, index, mode, rc)/rescale
    max.score = max(score, score2)
  }
  #create proper scale
  multiplier = 1
  while(max.score*multiplier<10) {
    multiplier = multiplier*10
  }
  #round up to nearest integer
  upper.bound = ceiling(max.score*multiplier)
  #find optimal step size
  divs = signif(seq(-upper.bound, upper.bound, length.out = 11)/multiplier, digits=2)
  #find top sites and rank them
  if (!is.null(nPeaks)) {
    idx = c(score[1,], score[2,])
    idx = cbind(idx, c(1:ncol(score), -(1:ncol(score))))
    idx = idx[order(-idx[,1]),]
    if (nPeaks>nrow(idx)) {
      nPeaks = nrow(idx)
    }
    Sequence = character(nPeaks)
    for (i in 1:nPeaks) {
      Sequence[i] = as.character(genomicSequence[abs(idx[i,2]):(abs(idx[i,2])+k-1)])
    }
    peaks = data.frame(Affinity=idx[1:nPeaks,1], Position=idx[1:nPeaks,2], Sequence)
    cat(paste0(fit.string," scores in ",seq.string,"\n"))
    print(peaks)
  }
  
  df = data.frame(Group=as.character(rep(1:2,each=ncol(score))), Position=c(1:ncol(score),1:ncol(score)), 
                  Affinity=c(score[1,], -score[2,]))
  if(!is.null(genomicSequence2)) {
    df = rbind(df, data.frame(Group=as.character(rep(3:4,each=ncol(score2))), Position=c(1:ncol(score2),1:ncol(score2)), 
                              Affinity=c(score2[1,], -score2[2,])))
    seq.string = paste0(seq.string, " and ", deparse(substitute(genomicSequence2)))
  }
  
  p = ggplot(df, aes(x=Position, y=Affinity, colour=Group)) +
    theme_bw() +
    geom_line() +
    ylab("Relative Affinity") +
    coord_fixed(ylim=c(-max.score,max.score)) +
    scale_y_continuous(breaks=divs, labels=format(abs(divs), nsmall=2))+
    theme(aspect.ratio=1, text=element_text(size=17, family="Helvetica")) +
    theme(legend.title=element_blank(), legend.position=c(0,0), legend.justification=c(0,0)) +
    theme(axis.title.x = element_text(vjust=0), axis.title.y = element_text(vjust=0)) +
    labs(title = paste0(fit.string," scores in ",seq.string))
  if (is.null(genomicSequence2)) {
    p = p+scale_color_manual(values=c("#000000", "#FF0000"), labels=c("Forward", "Reverse"))
  } else {
    p = p+scale_color_manual(values=c("#000000", "#FF0000", "blue", "green"), 
                             labels=c(paste0(deparse(substitute(genomicSequence)), " Forward"), 
                                      paste0(deparse(substitute(genomicSequence)), " Reverse"),
                                      paste0(deparse(substitute(genomicSequence2)), " Forward"), 
                                      paste0(deparse(substitute(genomicSequence2)), " Reverse"))) 
  }
  if (annotate && !is.null(nPeaks)) {
    xpos = abs(peaks$Position)+ncol(score)*.015
    ypos = sign(peaks$Position)*peaks$Affinity+max.score*.015
    p = p + annotate("text", x=xpos, y=ypos, label=as.character(1:nPeaks))
  }
  return(p)
}
