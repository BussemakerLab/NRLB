#' Use LogoGenerator tool from REDUCE Suite
#' 
#' @param output Directory in which output files are stored. The default is \code{tempdir()}.
#' @return An optional error code
#' @export
LogoGenerator = function(output, file, logo, title, ymin, ymax) {
  LogoGenerator_Cpp(output, file, logo, title, ymin, ymax);
}