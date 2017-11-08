#' Order fit results
#' 
#' @param fits
#' @param n
#' @return
#' 
#' @examples
#'
#' @export
#' 
o = function(fits, n=5) {
  info = fits[[1]]
  info = info[,-c(2, 3, 4, 17, 19, 20, 21)]
  info$k[info$k=="Multi"] = "-1"
  info = info[order(info$TrainLPerRead,-as.numeric(info$k), -info$Di, -as.numeric(as.POSIXct(info$Time))),]
  info$k[info$k=="-1"] = "Multi"
  return(info[1:n,])
}
