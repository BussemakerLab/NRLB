#' Create bigwig file from raw genomic profile
#' 
#' @param profiles Raw genomic score profiles, created using \code{fetch.genomic.profiles}
#' @param footprint Footprint size of the NRLB model that was used to generate the genomic profile.
#' @param file Name of the bigwig file that is to be created.
#' @param truncate.tol Cutoff below which affinities are assumed to be zero.
#' @param normalize Whether to use genomic maximum ('genome') or optimal site for model ('model') as reference for relative binding affinity. Note that model needs to be specified when the latter option is used.
#' @param model Model used for normalization based on theoretical optimal site (see 'optimize').
#'
#' @return None
#' 
#' @note First, the affinity scores for the forward and reverse strand are added up for each window. Next, the value assigned to any particular base pair is the maximum over all windows containing that base pair.
#' 
#' @examples
#' 
#' m = NRLBtools::hox.models()$ExdScr
#' gp = fetch.genomic.profiles("http://bussemakerlab.org/NRLB/dm3/ExdScr", 
#'                             genome = BSgenome.Dmelanogaster.UCSC.dm3, 
#'                             model = m, 
#'                             chr.names = c("chrXHet", "chrYHet"))
#' 
#' my.bwfile = tempfile(fileext=".bw")
#' create.bigwig.chaitanya(profiles = gp, file = my.bwfile, model = m)
#' 
#' bw = import(my.bwfile)
#' bw
#'
#' @export
#' 
create.bigwig.chaitanya = function(profiles, file, model=NULL, truncate.tol = 1E-4, normalize = "genome") {

  if (normalize == "genome") {
    reference.score = NRLBtools::genomic.max(profiles)
    cat("normalizing scores by genomic maximum\n")
  } else if (normalize == "model") {
    if (is.null) error("model needs to be defined")
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