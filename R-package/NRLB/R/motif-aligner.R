.rc = function(fits, index=NULL) {
  #Dinuc reverse çomplement indicies 
  dinuc.idx = c(16,12,8,4,15,11,7,3,14,10,6,2,13,9,5,1)
  if (is.null(index)) {
    fit.output = fits
  } else {
    fit.output = fits$Values[[index]]
  }
  #Handle multi-round fits
  if (class(fit.output[[1]])=="list" || length(fit.output)==9) {
    #Indices for simple reverse complement
    norm.rev.idx  = c(1, 3, 4, 6, 7, 9)
    #Indices for dinuc reverse complement
    dinuc.rev.idx = c(2, 5, 8)
  } else {
    #Indices for simple reverse complement
    norm.rev.idx  = c(1, 3, 5, 7, 9, 11)
    #Indices for dinuc reverse complement
    dinuc.rev.idx = c(2, 6, 10)
  }
  if (class(fit.output[[1]])=="list") {
    nModes = length(fit.output)
    if ("NSBinding" %in% names(fit.output)) {
      nModes = nModes-1
    }
    for (currMode in 1:nModes) {
      for (i in norm.rev.idx) {
        if (!is.null(fit.output[[currMode]][[i]])) {
          fit.output[[currMode]][[i]] = rev(fit.output[[currMode]][[i]])
        }
      }
      for (i in dinuc.rev.idx) {
        if (!is.null(fit.output[[currMode]][[i]])) {
          temp = matrix(fit.output[[currMode]][[i]], nrow=16)
          temp = temp[dinuc.idx,rev(1:ncol(temp))]
          fit.output[[currMode]][[i]] = as.numeric(temp)
        }
      }
    }
  } else {
    for (i in norm.rev.idx) {
      if (!is.null(fit.output[[i]])) {
        fit.output[[i]] = rev(fit.output[[i]])
      }
    }
    for (i in dinuc.rev.idx) {
      if (!is.null(fit.output[[i]])) {
        temp = matrix(fit.output[[i]], nrow=16)
        temp = temp[dinuc.idx,rev(1:ncol(temp))]
        fit.output[[i]] = as.numeric(temp)
      }
    }
  }
  return(fit.output)
}

#Requires a core fit list
.alignment.score = function(fitA, fitB, offset) {
  #Nucleotide Features
  motifA = fitA$NB
  motifB = fitB$NB
  if (offset<0) {
    motifB = motifB[(abs(offset)*4+1):length(motifB)]
  } else if (offset>0) {
    motifA = motifA[(offset*4+1):length(motifA)]
  }
  minLen = min(length(motifA), length(motifB))
  betasA = motifA[1:minLen]
  betasB = motifB[1:minLen]
  #Dinuc Features, if both exist
  motifA = fitA$DB
  motifB = fitB$DB
  if (!is.null(motifA) && !is.null(motifB)) {
    if (offset<0) {
      motifB = motifB[(abs(offset)*16+1):length(motifB)]
    } else if (offset>0) {
      motifA = motifA[(offset*16+1):length(motifA)]
    }
    minLen = min(length(motifA), length(motifB))
    betasA = c(betasA, motifA[1:minLen])
    betasB = c(betasB, motifB[1:minLen])
  }
  #Shape Features, if both exist
  motifA = fitA$SB
  motifB = fitB$SB
  if (!is.null(motifA) && !is.null(motifB)) {
    block.sizeA = 4*length(fitA$SB)/length(fitA$NB)
    block.sizeB = 4*length(fitB$SB)/length(fitB$NB)
    if (block.sizeA!=block.sizeB) {
      stop("The two fits do not have the same number of shape features")
    }
    if (offset<0) {
      motifB = motifB[(abs(offset)*block.size+1):length(motifB)]
    } else if (offset>0) {
      motifA = motifA[(offset*block.size+1):length(motifA)]
    }
    minLen = min(length(motifA), length(motifB))
    betasA = c(betasA, motifA[1:minLen])
    betasB = c(betasB, motifB[1:minLen])
  }
  summary(lm(betasA~betasB))$adj.r.squared
}

#Requires a core fit list


#' Align models
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param fitA
#' @param fitB
#' @param range
#' @return 
#' @note Requires a core fit list
#' 
#' @examples
#'
#' @export
#' 
motif.aligner = function(fitA, fitB, range=NULL) {
  if (is.null(range)) {
    range = floor(min(length(fitA$NB)/4,length(fitB$NB)/4)/2)
  }
  range = (-range):range
  rcB = .rc(fitB)
  A.B = numeric(length(range))
  A.rcB = A.B
  for (i in 1:length(range)) {
    A.B[i] = .alignment.score(fitA, fitB, range[i])
    A.rcB[i] = .alignment.score(fitA, rcB, range[i])
  }
  max1 = which.max(A.B)
  max2 = which.max(A.rcB)
  glob.max = which.max(c(A.B[max1], A.rcB[max2]))
  if (glob.max==1) {
    return( data.frame(B.RC=FALSE, Offset=range[max1], R2=A.B[max1]))
  } else {
    return( data.frame(B.RC=TRUE, Offset=range[max2], R2=A.rcB[max2]))
  }
}
