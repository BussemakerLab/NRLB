#' Quick Logo
#' 
#' @param fits
#' @param index
#' @param mode
#' @param rc
#' @return
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
