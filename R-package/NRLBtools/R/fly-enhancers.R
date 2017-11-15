#' Define the set fly enhancers discussed in NRLB paper
#' 
#' @return Nested list containing DNA sequence and other information
#'
#' @examples
#' 
#' 
#' @export
#' 
fly.enhancers = function() {
  df = read.table(system.file("extdata", "fly-enhancers.tsv", package = "NRLBtools"),
                  stringsAsFactors = FALSE)
  genome = BSgenome.Dmelanogaster.UCSC.dm3::Dmelanogaster
  df$seq = as.list(
    IRanges::Views(genome,
                   GenomicRanges::GRanges(df$chr,
                                          IRanges::IRanges(df$start, df$end),
                                          strand = "+")
                   )
    )
  nested.list = apply(df, 1, as.list)
  return(nested.list)
}
