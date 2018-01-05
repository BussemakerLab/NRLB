#' Fetch pre-computed genomic score profiles as GRanges object
#' 
#' @param url.prefix Prefix of URL containing precomputed profiles
#' @param genome Complete genome sequence in BSgenome format
#' @param chr.names Chromosome(s) to which profile should be restricted (optional; the default is to include profiles for all chromosomes in the genome)
#' @param model An NRLB model (required to determine footprint)
#' @param rc Whether to swap forward and reverse profiles
#'
#' @return A nested list containing forward and reverse NRLB score profiles for each chromosome.
#' 
#' @examples
#' 
#' m = NRLBtools::hox.models()$ExdScr
#' gr = fetch.as.track("http://bussemakerlab.org/NRLB/dm3/ExdScr", 
#'                     genome = BSgenome.Dmelanogaster.UCSC.dm3, 
#'                     model = m, normalize = "none",
#'                     chr.names = c("chrXHet", "chrYHet")
#'                     )
#' gr
#' 
#' v = Views(BSgenome.Dmelanogaster.UCSC.dm3$chrXHet, start=1:6, width=footprint.size(m))
#' sapply(v, function(x) score.genome(x, model=m))
#'
#' @export
#'
fetch.as.track = function(url.prefix, genome, chr.names = NULL, normalize = "none", model=NULL) {
 
  if (is.null(chr.names)) { chr.names = seqnames(genome) }
  
  fs = NRLBtools::footprint.size(model)
  
  gr.list = lapply(chr.names, function(chr) {
    
    cat(sprintf("fetching %s [%s] ...", url.prefix, chr), "\n")
    cl = seqlengths(genome)[chr]
    
    url = sprintf("%s/_%s_F.dat", url.prefix, chr)
    score.fwd = readBin(url, what = "double", endian = "big", n = cl-fs+1)
    url = sprintf("%s/_%s_R.dat", url.prefix, chr)
    score.rev = readBin(url, what = "double", endian = "big", n = cl-fs+1)
    #url = sprintf("%s/_%s_Flag.dat", url.prefix, chr)
    #flag = readBin(url, what = "double", endian = "big", n = cl-fs+1)

    if (model$rc == TRUE) { fwd.strand = "-"; rev.strand = "+" } else {fwd.strand = "+"; rev.strand = "-"}
    
    gr.fwd = GRanges(seqnames = rep(chr, length(score.fwd)), 
                     ranges = IRanges(start=1:length(score.fwd), width=fs),
                     strand = fwd.strand,
                     score = score.fwd)
    gr.rev = GRanges(seqnames = rep(chr, length(score.rev)), 
                     ranges = IRanges(start=1:length(score.rev), width=fs),
                     strand = rev.strand,
                     score = score.rev)
    
    c(gr.fwd, gr.rev)
  })
  gr = suppressWarnings(do.call("c", gr.list)) # suppress warning about non-overlapping seqnames
  seqinfo(gr) = seqinfo(genome)[chr.names]
  
  if (normalize == "genome") {
    reference.score = max(gr$score)
    cat("normalizing scores by genomic maximum\n")
    gr$score = gr$score / reference.score
  } else if (normalize == "model") {
    if (is.null(model)) error("model needs to be defined")
    reference.score = NRLBtools::optimal.site(model$fits, model$index, model$mode)$score
    cat("normalizing scores by theoretical maximum\n")
    gr$score = gr$score / reference.score
  } else if (normalize == "none") {
    cat("using raw scores (not recommended)\n")
  } else {
    error("normalize argument needs to be 'model' or 'genome'")
  }
  
  return(gr)
}


#' Fetch pre-computed genomic score profiles as nested list object
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


#' Maximum NRLB score in genomic profile
#' 
#' @param profiles Raw genomic score profiles, created using \code{fetch.genomic.profiles}
#'
#' @return Maximumover all chromosomes and both binding directions
#' 
#' @examples
#'
#' @export
#' 
genomic.max = function(profiles) {
  max(sapply(profiles, function(x) max(x$fwd, x$rev)))
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