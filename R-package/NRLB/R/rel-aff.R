#' Calculate relative affinities for a fit
#' 
#' @param fits
#' @param index
#' @return
#' 
#' @examples
#'
#' @export
#' 
aff = function(fits, index) {
  modes = max.seq(fits, index)$MaxAffinity
  nsb   = exp(fits[[1]]$NSBind[index])
  tot   = c(modes, nsb)
  tot   = tot/max(tot)
  return(log(tot))
}
