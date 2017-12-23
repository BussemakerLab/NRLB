#' Find the highest affinity sequence for a model
#' 
#' @param fits TO_BE_ADDED
#' @param index TO_BE_ADDED
#' @param mode TO_BE_ADDED
#' @return TO_BE_ADDED
#' 
#' @note uses dynamic programming
#' 
#' @examples
#'
#' @export
#' 
optimal.site <- function(fits, index, mode = NULL) {
  #Get betas and transform them into a matrix
  fit = fits[[2]][[index]]
  k = fits[[1]]$k[index]
  isMulti = (k=="Multi")
  if (isMulti) {
    output = NULL
    #Loop over all modes
    if (is.null(mode)) {
      modes = 1:length(fit)
    } else {
      modes = mode
    }
    for (currMode in modes) {
      if (names(fit)[currMode]=="NSBinding") {
        next
      }
      nuc = fit[[currMode]]$NB
      k = length(nuc)/4
      dim(nuc) = c(4, k)
      if (is.null(fit[[currMode]]$DB)) {
        isDinuc = FALSE
        dinuc = NULL
      } else {
        isDinuc = TRUE
        dinuc = fit[[currMode]]$DB
        dim(dinuc) = c(16, k-1)
        rownames(dinuc) = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
      }
      output = rbind(output, .maxSeqHelper(isDinuc, k, nuc, dinuc))
    }
    return(output)
  } else {
    nuc = fit$NB
    k = as.numeric(k)
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
    return(.maxSeqHelper(isDinuc, k, nuc, dinuc))
  }
}

.maxSeqHelper = function(isDinuc, k, nuc, dinuc) {
  #Initialize loop (position 1)
  max.list = nuc[,1] #(A C G T)
  char.list= c("A", "C", "G", "T")
  temp.list= char.list
  curr.list = matrix(data=0, 4, 4)
  #Loop over all positions
  for (currPos in 2:k) {
    #Loop over all previous bases
    for (prevBase in 1:4) {
      curr.list[,prevBase] = max.list[prevBase] + nuc[, currPos]
      if (isDinuc) {
        curr.list[,prevBase] = curr.list[,prevBase] + dinuc[((prevBase-1)*4+1):(prevBase*4), (currPos-1)]
      }
    }
    for (currBase in 1:4) {
      max.list[currBase] = max(curr.list[currBase,])
      temp.list[currBase] = paste0(char.list[which.max(curr.list[currBase,])], c("A", "C", "G", "T")[currBase])
    }
    char.list = temp.list
  }
  return(data.frame(score=exp(max(max.list)), seq=char.list[which.max(max.list)], stringsAsFactors = FALSE))
}
