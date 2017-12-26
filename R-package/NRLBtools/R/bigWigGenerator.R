#' Read in binary files containing genomic scoring results
#' 
#' @param profiles Raw genomic score profiles, created using \code{fetch.genomic.profiles}
#' @param footprint Footprint size of the NRLB model that was used to generate the genomic profile.
#' @param file Name of the bigwig file that is to be created.
#' @param truncate.tol Cutoff below which affinities are assumed to be zero.
#'
#' @return None
#' 
#' @examples
#' 
#' m = NRLBtools::hox.models()$ExdScr
#' gp = fetch.genomic.profiles("http://bussemakerlab.org/NRLB/dm3/ExdScr", 
#'                             genome = BSgenome.Dmelanogaster.UCSC.dm3, 
#'                             model = m, 
#'                             chr.names = c("chrXHet", "chrYHet"))
#' 
#' create.bigwig(profiles = gp, file = tempfile(fileext=".bw"), model = m)
#'
#' @export
#' 
create.bigwig = function(profiles, file, model, truncate.tol = 1E-4, normalize = "model") {

  if (normalize == "genome") {
    reference.score = NRLBtools::genomic.max(profiles)
    cat("normalizing scores by genomic maximum\n")
  } else if (normalize == "model") {
    reference.score = NRLBtools::optimal.site(model$fits, model$index, model$mode)$score
  } else {
    error("normalize argument needs to be 'model' or 'genome'")
  }
  
  k = NRLBtools::footprint.size(model)
  
  #Filter, create GRanges object, and store bigWig file
  currOut = GRanges()
  for (currChr in names(profiles)) {
    scores = rbind(profiles[[currChr]]$fwd, profiles[[currChr]]$rev)
    currScores = colSums(scores) # add FWD and REV scores
    currScores = currScores / reference.score # normalize
    currScores = signif(currScores, 5)
    idx = currScores<=truncate.tol #Truncate
    currScores[idx] = 0
    chr.length = length(currScores) + k - 1
    processedScores = numeric(length = chr.length) #Create windows
    processedScores[1:k] = currScores[1]
    for (i in 2:length(currScores)) {
      if (processedScores[i] < currScores[i]) {
        processedScores[i:(i+k-1)] = currScores[i]
      }
    }
    idx = processedScores==0 #Remove NAs
    processedScores[idx] = NA
    pos = 1:length(processedScores)
    processedScores = cbind(pos, processedScores)
    processedScores = processedScores[complete.cases(processedScores),]
    #Minimize and find uniques through a flag array
    flag = (processedScores[,2] - c(processedScores[2:nrow(processedScores),2],0))!=0
    fStrand = matrix(data=NA, nrow = sum(flag), ncol = 3)     #output matrix
    idx = flag==1                                             #index array
    fStrand[,1] = processedScores[c(TRUE,idx[1:(nrow(processedScores)-1)]), 1]
    fStrand[,2] = processedScores[idx, 1]
    fStrand[,3] = processedScores[idx, 2]
    seqLen = chr.length
    names(seqLen) = currChr
    currOut = append(currOut, GRanges(seqnames=currChr, ranges=IRanges(start=fStrand[,1], end=fStrand[,2]), strand="*", seqlengths=seqLen, score=fStrand[,3]))
  }
  export(currOut, file)
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