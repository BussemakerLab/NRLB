#' Plot likelihoods
#' 
#' @param fits
#' @param l
#' @param test
#' @return
#' @note does not work with multi-mode fits (notion of window is unknown)
#' 
#' @examples
#'
#' @export
#' 
plot.likelihoods = function(fits, l, test=FALSE) {
  old = list(fits$Information)
  for (colIdx in c(5, 7:12)) {
    new = NULL
    for (currList in 1:length(old)) {
      new = c(new, by(old[[currList]], old[[currList]][,colIdx], function(x) x))
    }
    old = new
  }
  tot.divs = length(old)
  likelihoods = vector(mode="list", length=tot.divs)
  fit.types = NULL
  for (i in 1:tot.divs) {
    topRow = old[[i]][1,]
    dat.string = c(topRow$k,sapply(7:12, function (x) if (as.logical(topRow[x])) {names(topRow)[x]} else {NA}))
    dat.string = dat.string[!is.na(dat.string)]
    fit.types = c(fit.types, paste(dat.string,collapse="-"))
    output = matrix(nrow=nrow(old[[i]]), ncol=3)
    output[,1] = l+2*old[[i]]$Flank-as.numeric(old[[i]]$k)+1
    output[,2] = old[[i]]$Shift
    if (test) {
      output[,3] = old[[i]]$TestLPerRead
    } else {
      output[,3] = old[[i]]$TrainLPerRead
    }
    likelihoods[[i]] = output
  }
  window.range = range(l-as.numeric(fits$Information$k)+2*fits$Information$Flank+1)
  if (test) {
    lik.range = range(fits$Information$TestLPerRead)
    plot(x = 0, y=0, type="n", ylim=lik.range, xlim=window.range, xlab="Windows", ylab="Testing -Log Likelihood")
  } else {
    lik.range = range(fits$Information$TrainLPerRead)
    plot(x = 0, y=0, type="n", ylim=lik.range, xlim=window.range, xlab="Windows", ylab="Training -Log Likelihood")
  }
  col.type = NULL
  lty.type = NULL
  for (i in 1:tot.divs) {
    col.type = 1 #c(col.type, (i %% 3)+1)
    lty.type = c(lty.type, i) #c(lty.type, (1+((i-1) %/% 3)))
    curr.break = by(likelihoods[[i]], likelihoods[[i]][,2], function(x) x)
    for (b in 1:length(curr.break)) {
      lines(x=curr.break[[b]][,1], y=curr.break[[b]][,3], col="grey", lty=4)
      points(x=curr.break[[b]][,1], y=curr.break[[b]][,3], col=1, pch=i, cex=.75)
    }
  }
  legend("topright", legend=fit.types, col=col.type, pch=lty.type, bty="n", ncol = 2)
}
