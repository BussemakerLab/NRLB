#' Score genome with a fit
#' 
#' @import Biostrings
#' 
#' @param genomicSequence TO_BE_ADDED
#' @param fits TO_BE_ADDED
#' @param index TO_BE_ADDED
#' @param mode TO_BE_ADDED
#' @param rc TO_BE_ADDED
#' @return TO_BE_ADDED
#' 
#' @examples
#'
#' @export
#' 
score.genome = function(genomicSequence, fits, index, mode=NULL, rc=FALSE, model=NULL) {
  if(!is.null(model)) {fits=model$fits; index=model$index; mode=model$mode}
  #Create sequences for rapid scoring
  fSeq = abs(toComplex(genomicSequence, c(A=1, C=2, G=3, T=4, N=0)))
  rSeq = abs(toComplex(reverseComplement(genomicSequence), c(A=1, C=2, G=3, T=4, N=0)))
  charFSeq = as.character(genomicSequence)
  charRSeq = as.character(reverseComplement(genomicSequence))
  l = length(genomicSequence)
  fit = fits[[2]][[index]]
  #Check to see if input is a multi-mode model
  if (fits[[1]]$k[index]=="Multi") {
    if (is.null(mode)) {
      stop("Multi-Mode Fit Detected: Mode Index Required")
    } else {
      fit = fit[[mode]]
      k = length(fit$NB)/4
    }
  } else {
    k = as.numeric(fits[[1]]$k[index])
  }
  adjK = k-1
  adjL = l+1
  nuc = fit$NB
  dim(nuc) = c(4,k)
  if (is.null(fit$DB)) {
    isDinuc = FALSE
    dinuc = NULL
  } else {
    isDinuc = TRUE
    dinuc = fit$DB
    dim(dinuc) = c(16, k-1)
    rownames(dinuc) = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
  }
  score = sapply(1:(l-adjK), FUN=function(x) .fastScore(substr(charFSeq, x, x+adjK), substr(charRSeq, adjL-x-adjK, adjL-x), 
                                                        fSeq[x:(x+adjK)], rSeq[(adjL-x-adjK):(adjL-x)], k, isDinuc, nuc, dinuc))
  if(rc) {
    score = score[c(2, 1),]
  }
  return(score)
}

#fast sequence scorer, optimized to work with score.genome
.fastScore = function(charFSeq, charRSeq, fSeq, rSeq, k, isDinuc, nuc, dinuc) {
  if (all(fSeq!=0)) {
    fTotal = 0
    rTotal = 0
    if (isDinuc) {
      for (j in 1:k) {
        fTotal = fTotal+nuc[fSeq[j],j]
        rTotal = rTotal+nuc[rSeq[j],j]
        if (j<k) {
          fTotal = fTotal+as.numeric(dinuc[substr(charFSeq, j, j+1),j])
          rTotal = rTotal+as.numeric(dinuc[substr(charRSeq, j, j+1),j])
        }
      }
    } else {
      for (j in 1:k) {
        fTotal = fTotal+nuc[fSeq[j],j]
        rTotal = rTotal+nuc[rSeq[j],j]
      }
    }
    output = c(exp(fTotal), exp(rTotal))
    return(output)
  } else {
    return(c(NA, NA))
  }
}
