

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


#' Fetch pre-computed genomic score profiles
#' 
#' @param genome Complete genome sequence in BSgenome format
#' @param url.prefix Prefix of URL containing precomputed profiles
#' @param model An NRLB model (required to determine footprint)
#' @param chr.names Chromosome(s) to which profile should be restricted (optional; the default is to include profiles for all chromosomes in the genome)
#'
#' @return A nested list containing forward and reverse NRLB score profiles for each chromosome.
#' 
#' @examples
#' 
#' m = NRLBtools::hox.models()$ExdScr
#' gp = fetch.genomic.profiles("http://bussemakerlab.org/NRLB/dm3/ExdScr", 
#'                             BSgenome.Dmelanogaster.UCSC.dm3, m, 
#'                             c("chrXHet", "chrYHet"))
#' names(gp)
#' names(gp$chrXHet)
#' 
#' 
#' head(gp$chrXHet$fwd)
#' head(gp$chrXHet$rev)
#' 
#' v = Views(BSgenome.Dmelanogaster.UCSC.dm3$chrXHet, start=1:6, width=footprint.size(m))
#' sapply(v, function(x) score.genome(x, model=m))
#'
#' @export
#'
fetch.genomic.profiles = function(url.prefix, genome, model, chr.names = NULL) {
  if (is.null(chr.names)) { chr.names = seqnames(genome) }
  result = lapply(chr.names, function(chr) {
    cat(sprintf("fetching %s [%s] ...", url.prefix, chr), "\n")
    list(fwd = readBin(sprintf("%s/_%s_F.dat", url.prefix, chr),
                       what = "double", endian = "big",
                       n = seqlengths(genome)[chr] - footprint.size(model) ),
         rev = readBin(sprintf("%s/_%s_R.dat", url.prefix, chr),
                       what = "double", endian = "big",
                       n = seqlengths(genome)[chr] - footprint.size(model) ), 
         flag = readBin(sprintf("%s/_%s_Flag.dat", url.prefix, chr),
                        what = "double", endian = "big",
                        n = seqlengths(genome)[chr] - footprint.size(model) ),
         footprint.size = NRLBtools::footprint.size(model),
         chr.length = seqlengths(genome)[chr]
    )
    
  })
  names(result) = chr.names
  return(result)
}


#' Read in binary files containing genomic scoring results
#' 
#' @param genome TO_BE_ADDED
#' @param score.dir TO_BE_ADDED
#' @param prot.name TO_BE_ADDED
#' @param chr.name TO_BE_ADDED
#' @param k TO_BE_ADDED
#'
#' @return TO_BE_ADDED
#' 
#' @examples
#'
#' @export
#' 
score.parser = function(genome, score.dir, prot.name, chr.name, k) {
    s = length(genome[[chr.name]])
    fwd.score = readBin(con = paste0(score.dir,prot.name,"_",chr.name,"_F.dat"), what = "double", n=s, endian = "big")
    rev.score = readBin(con = paste0(score.dir,prot.name,"_",chr.name,"_R.dat"), what = "double", n=s, endian = "big")
    flag = readBin(con = paste0(score.dir,prot.name,"_",chr.name,"_Flag.dat"), what = "integer", n=s, size=1, endian = "big")
    return(list(score=rbind(fwd.score, rev.score), flag=flag, chr=chr.name, k=k))
}