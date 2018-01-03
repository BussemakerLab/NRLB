#' Footprint size associated with model
#' 
#' @param model NRLB model
#'
#' @return Footprint size in bp
#' 
#' @examples
#' 
#' m = NRLBtools::hox.models()$ExdScr
#' NRLBtools::footprint.size(m)
#'
#' @export
#' 
footprint.size = function(model) {
  nchar(NRLBtools::optimal.site(model$fits, model$index, model$mode)$seq)
}