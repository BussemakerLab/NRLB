#' Quick Logo
#' 
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
ql = function(fits, index, mode=NULL, rc=FALSE) {
  if (fits[[1]]$k[index]=="Multi") {
    if (is.null(mode)) {
      mode = which.max(aff(m, index))
    }
  }
  logo(fits, index = index, mode = mode, rc = rc, fit.name=deparse(substitute(fits)))
}
