#' R0 QQ Plot
#' 
#' @param fileName TO_BE_ADDED
#' @param lowerBound TO_BE_ADDED
#' @param upperBound TO_BE_ADDED
#' @param divSize TO_BE_ADDED
#' @param completeCoverage TO_BE_ADDED
#' @param compare TO_BE_ADDED
#' @param ylim TO_BE_ADDED
#' @param xlim TO_BE_ADDED
#' @return TO_BE_ADDED
#' 
#' @examples
#'
#' @export
#' 
plot.R0.qq = function(fileName, lowerBound, upperBound, divSize, completeCoverage=TRUE, compare=NULL, ylim=NULL, xlim=NULL) {
  counts = read.table(fileName, header=FALSE)
  totCount = colSums(counts)
  nComparisons = nrow(compare)
  if (is.null(nComparisons)) {
    nComparisons = 1
    compare = t(as.matrix(compare))
  }
  if (completeCoverage) {
    totCount = sum(totCount*c(0:(ncol(counts)-1)))
  } else {
    totCount = sum(totCount*c(1:ncol(counts)))
  }
  x = seq(lowerBound, upperBound, by=divSize)
  #truncate zero rows
  all.zeros = apply(counts, 1, function(row) all(row==0))
  #find first and last non-zero elements
  for (i in 1:length(all.zeros)) {
    if (!all.zeros[i]){
      low = i
      break
    }
  }
  for (i in length(all.zeros):1) {
    if (!all.zeros[i]){
      high = i
      break
    }
  }
  x = x[low:high]
  #plot all combinations
  if (is.null(compare)) {
    compare = t(combn(c(1:ncol(counts)), m=2))
  } else if (completeCoverage) {
    compare = compare+1
  }
  #build comparisons
  ratios = NULL
  for (i in 1:nComparisons) {
    ratios = cbind(ratios, log10(counts[,compare[i,2]]/counts[,compare[i,1]]))
  }
  ratios = ratios[low:high,]
  #plot ratios
  colors = 1:nComparisons
  if (is.null(ylim)) {
    if (is.null(xlim)) {
      #ylim and xlim null
      matplot(x, ratios, type="l", lty=1, col=colors, xlab="log(w)", ylab="Log Ratio", main=fileName)
    } else {
      #ylim null and xlim not null
      matplot(x, ratios, type="l", lty=1, col=colors, xlab="log(w)", ylab="Log Ratio", main=fileName, xlim=xlim)
    }
  } else if (is.null(xlim)) {
    #ylim not null but xlim is null
    matplot(x, ratios, type="l", lty=1, col=colors, xlab="log(w)", ylab="Log Ratio", main=fileName, ylim=ylim)
  } else {
    #ylim and xlim not null
    matplot(x, ratios, type="l", lty=1, col=colors, xlab="log(w)", ylab="Log Ratio", main=fileName, xlim=xlim, ylim=ylim)
  }
  if (completeCoverage) {
    compare = compare-1
  }
  for (i in 1:nComparisons) {
    lines(x, .poisson.ratio(10^x*totCount, compare[i,1], compare[i,2]), col="grey", lty=2)
  }
  #Legend
  legend.text = NULL
  for (i in 1:nComparisons) {
    legend.text = c(legend.text, paste0(compare[i,], collapse="-"))
  }
  legend("bottomright", legend=legend.text, lty=1, col=colors, bty="n", ncol=ceiling(nrow(compare)/3), title="Rounds", cex=.75)
  locs = par('usr')
  text(x=.025*(locs[2]-locs[1])+locs[1], y=(.975*(locs[4]-locs[3])+locs[3]), labels=paste0("Total Reads: ", totCount), pos=4, cex=.75)
}

.poisson.ratio = function(x, countA, countB) {
  log10(x^(countB-countA)*factorial(countA)/factorial(countB))
}
