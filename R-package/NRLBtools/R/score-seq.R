#' Score DNA sequence using a model
#' 
#' @param sequence TO_BE_ADDED
#' @param rSequence TO_BE_ADDED
#' @param k TO_BE_ADDED
#' @param isDinuc TO_BE_ADDED
#' @param nuc TO_BE_ADDED
#' @param dinuc TO_BE_ADDED
#' @return TO_BE_ADDED
#' @note requires explicit nuc, dinuc formatting (mode independent)
#' 
#' @examples
#'
#' @export
#' 
score.seq = function(sequence, rSequence, k, isDinuc, nuc, dinuc) {
  if (isDinuc) {
    fSeq = abs(toComplex(sequence, c(A=1, C=2, G=3, T=4)))
    rSeq = abs(toComplex(rSequence, c(A=1, C=2, G=3, T=4)))
    fTotal = 0
    rTotal = 0
    for (j in 1:k) {
      fTotal = fTotal+nuc[fSeq[j],j]
      rTotal = rTotal+nuc[rSeq[j],j]
      if (j<k) {
        fTotal = fTotal+as.numeric(dinuc[as.character(sequence[j:(j+1)]),j])
        rTotal = rTotal+as.numeric(dinuc[as.character(rSequence[j:(j+1)]),j])
      }
    }
  } else {
    fSeq = abs(toComplex(sequence, c(A=1, C=2, G=3, T=4)))
    rSeq = abs(toComplex(rSequence, c(A=1, C=2, G=3, T=4)))
    fTotal = 0
    rTotal = 0
    for (j in 1:k) {
      fTotal = fTotal+nuc[fSeq[j],j]
      rTotal = rTotal+nuc[rSeq[j],j]
    }
  }
  output = c(exp(fTotal), exp(rTotal))
  return(output)
}