#' Collection of NRLB models for Hox complexes
#'
#' @return Nested list containing model data and metadata
#'
#' @examples
#' 
#' models = hox.models()
#' names(models)
#' 
#' @export
#' 
hox.models = function() {
  df = read.table(system.file("extdata", "hox-models.tsv", package = "NRLBtools"), 
                  stringsAsFactors = FALSE)
  df$file = system.file("extdata/fitsets", df$file, package = "NRLBtools")
  df$fits = lapply(df$file, file.parser)
  nested.list = apply(df, 1, as.list)
  
  for(i in which(is.na(sapply(nested.list, function(x) x$mode)))) nested.list[[i]]$mode = NULL # turn NA into NULL

  return(nested.list)
}