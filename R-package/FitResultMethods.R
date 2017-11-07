library(Biostrings)
library(ggplot2)
library(grid)

file.parser = function(fileName, fitInfoOnly=FALSE) {
  nCols = max(count.fields(fileName, sep=","))
  header = read.csv(fileName, header=FALSE, nrows = 1, stringsAsFactors=FALSE, strip.white=TRUE)
  colNames = as.character(as.vector(header[1,]));
  colNames = c(colNames, (length(colNames)+1):nCols)
  output = read.csv(fileName, header=FALSE, fill=TRUE, col.names=colNames, stringsAsFactors=FALSE)
  output = output[2:nrow(output),2:ncol(output)]
  row.names(output) = as.character(1:nrow(output))
  for (i in 7:12) {
    output[,i] = as.logical(output[,i])
  }
  for (i in c(2:6, 15:23)) {
    output[,i] = as.numeric(output[,i])
  }
  output[output==""] <- NA
  output$k[is.na(output$k)] = "Multi"
  if (fitInfoOnly) {
    output = output[,1:24]
    return(output)
  } else {    #Parse Betas, Seed, etc.
    info = output[,1:24]
    values = vector(mode="list", length=nrow(output))
    for (i in 1:nrow(output)) {
      nb = db = sb = nsb = ne = de = se = nsbe = ns = ds = ss = nsbs = NULL
      #check to see if multi-mode fit
      if (info$k[i]=="Multi") {
        #Multi-mode fit; find the number of modes and parse NSB first
        fit = output[i,]
        modeStartIdx = grep("NB>", fit)
        nModes = length(modeStartIdx)
        modeStartIdx = c(modeStartIdx, grep("<EOL>", fit))
        if (info$NS[i]) {
          modes = vector(mode="list", length=(nModes+1))
          names(modes) = c(paste0("Mode", 1:nModes),"NSBinding")
          nsb = as.numeric(fit[grep("NSB>", fit)+1])
          nsbe= as.numeric(fit[grep("NSBE>", fit)+1])
          nsbs= as.numeric(fit[grep("NSBS>", fit)+1])
          modes$NSBinding = list(NSB=nsb, NSBE=nsbe, NSBS=nsbs)
        } else {
          modes = vector(mode="list", length=nModes)
          names(modes) = paste0("Mode", 1:nModes)
        }
        #Now parse the parameters for the modes
        for (currMode in 1:nModes) {
          nb = db = sb = ne = de = se = ns = ds = ss = NULL
          currModeValues = c(fit[modeStartIdx[currMode]:(modeStartIdx[currMode+1]-1)],">")
          delimiters = grep(">", currModeValues)
          for (delimIdx in 1:(length(delimiters)-1)) {
            currValue = as.numeric(currModeValues[(delimiters[delimIdx]+1):(delimiters[delimIdx+1]-1)])
            currType = currModeValues[delimiters[delimIdx]]
            if      (currType=="NB>")   {nb = currValue}
            else if (currType=="DB>")   { db  = currValue }
            else if (currType=="SB>")   { sb  = currValue }
            else if (currType=="NE>")   { ne  = currValue }
            else if (currType=="DE>")   { de  = currValue }
            else if (currType=="SE>")   { se  = currValue }
            else if (currType=="NS>")   { ns  = currValue }
            else if (currType=="DS>")   { ds  = currValue }
            else if (currType=="SS>")   { ss  = currValue }
          }
          modes[[currMode]] = list(NB=nb, DB=db, SB=sb, NE=ne, DE=de, SE=se, NS=ns, DS=ds, SS=ss)
        }
        values[[i]] = modes
      } else {
        fit = output[i,]
        delimiters = grep(">", fit)
        for (delimIdx in 1:(length(delimiters)-1)) {
          currValue = as.numeric(fit[(delimiters[delimIdx]+1):(delimiters[delimIdx+1]-1)])
          currType = fit[delimiters[delimIdx]]
          if      (currType=="NB>")   {nb = currValue}
          else if (currType=="DB>")   { db  = currValue }
          else if (currType=="SB>")   { sb  = currValue }
          else if (currType=="NSB>")  { nsb = currValue }
          else if (currType=="NE>")   { ne  = currValue }
          else if (currType=="DE>")   { de  = currValue }
          else if (currType=="SE>")   { se  = currValue }
          else if (currType=="NSBE>") { nsbe= currValue }
          else if (currType=="NS>")   { ns  = currValue }
          else if (currType=="DS>")   { ds  = currValue }
          else if (currType=="SS>")   { ss  = currValue }
          else if (currType=="NSBS>") { nsbs= currValue }
        }
        values[[i]] = list(NB=nb, DB=db, SB=sb, NSB=nsb, NE=ne, DE=de, SE=se, NSBE=nsbe, NS=ns, DS=ds, SS=ss, NSBS=nsbs)
      }
    }
    return(list(Information=info, Values=values))
  }
}

score.parser = function(genome, score.dir, prot.name, chr.name, k) {
  s = length(genome[[chr.name]])
  fwd.score = readBin(con = paste0(score.dir,prot.name,"_",chr.name,"_F.dat"), what = "double", n=s, endian = "big")
  rev.score = readBin(con = paste0(score.dir,prot.name,"_",chr.name,"_R.dat"), what = "double", n=s, endian = "big")
  flag = readBin(con = paste0(score.dir,prot.name,"_",chr.name,"_Flag.dat"), what = "integer", n=s, size=1, endian = "big")
  return(list(score=rbind(fwd.score, rev.score), flag=flag, chr=chr.name, k=k))
}

print.fit = function(fits, index=NULL) {
  #See if only the info has been provided
  if (class(fits)=="data.frame") {
    #Print all fits if no index is provided
    if (is.null(index)) {
      index = 1:nrow(fits)
    }
    for (i in index) {
      fit = fits[i,]
      cat(paste0("k: ", fit$k, ", f: ", fit$Flank,", shift: ",fit$Shift,"\t"))
      if (fit$NS) { cat("NSBinding  ") } 
      if (fit$Di) { cat("Dinuc  ") }
      if (fit$MG) { cat("MGW  ") }
      if (fit$PT) { cat("ProT  ") }
      if (fit$HT) { cat("HelT  ") }
      if (fit$RO) { cat("Roll  ") }
      cat("\n")
      cat("Fit Steps:",fit$FitSteps,"\t\tFunction Calls:",fit$FncCalls,"\tFitting Time:", fit$FitTime, "\n")
      if (!is.na(fit$NuSym) || !is.na(fit$DiSym)) {
        if (!is.na(fit$NuSym)) {cat(fit$NuSym,"\t")}
        if (!is.na(fit$DiSym)) {cat(fit$DiSym)}
        cat("\n")
      }
      cat("Train -LL:",fit$TrainLPerRead,"\tTest -LL:",fit$TestLPerRead,"\n")
      cat("Nucleotide PSAM:\t",fit$PSAM,"\n")
      cat("------------------------------------------------------------------------------------------------------------\n")
    }
  } else {
    if (is.null(index)) {
      index = 1:nrow(fits$Information)
    }
    #Print the raw values as well here
    for (i in index) {
      fit = fits$Information[i,]
      values = fits$Values[[i]]
      isMulti = fit$k=="Multi"
      if (isMulti) {
        k = "Multi"
      } else {
        k = as.numeric(fit$k)
      }
      cat(paste0("k: ", k, ", f: ", fit$Flank,", shift: ",fit$Shift,"\t"))
      if (fit$NS) { cat("NSBinding  ") } 
      if (fit$Di) { cat("Dinuc  ") }
      if (fit$MG) { cat("MGW  ") }
      if (fit$PT) { cat("ProT  ") }
      if (fit$HT) { cat("HelT  ") }
      if (fit$RO) { cat("Roll  ") }
      cat("\n")
      cat("Fit Steps:",fit$FitSteps,"\t\tFunction Calls:",fit$FncCalls,"\tFitting Time:", fit$FitTime, "\n")
      if (!is.na(fit$NuSym) || !is.na(fit$DiSym)) {
        if (!is.na(fit$NuSym)) {cat(fit$NuSym,"\t")}
        if (!is.na(fit$DiSym)) {cat(fit$DiSym)}
        cat("\n")
      }
      cat("Train -LL:",fit$TrainLPerRead,"\tTest -LL:",fit$TestLPerRead,"\n\n")
      #check to see if fit is a multi-mode fit
      if (isMulti) {
        #First print NS binding, if it exists:
        if (fit$NS) {
          if (!is.null(values$NSBinding$NSBE)) {
            cat(sprintf("NS Binding Beta:\t%7.4f \u00B1 %-7.4f\n",values$NSBinding$NSB, values$NSBinding$NSBE))
          } else {
            cat(sprintf("Nucleotide Betas:\t%8.5f\n",values$NSBinding$NSB))
          }
          nModes = length(values)-1
        } else {
          nModes = length(values)
        }
        #Now loop over modes
        for (currMode in 1:nModes) {
          k = length(values[[currMode]]$NB)/4
          #Next print nucleotide PSAMs:
          cat(paste0("Mode ", currMode, " Nucleotide Betas:\n"))
          .print.psam(values[[currMode]]$NB, values[[currMode]]$NE, k, NULL)
          #Next print dinucleotide PSAMs:
          if (!is.null(values[[currMode]]$DB)) {
            cat(paste0("\nMode ", currMode, " Dinucleotide Betas:\n"))
            .print.psam(values[[currMode]]$DB, values[[currMode]]$DE, k-1, NULL)
          }
          #Next print shape PSAM:
          if (!is.null(values[[currMode]]$SB)) {
            cat(paste0("\nMode ", currMode, " Shape Betas:\n"))
            .print.psam(values[[currMode]]$SB, values[[currMode]]$SE, k, NULL)
          }
          
          #Next print NS Seed, if it exists:
          if (!is.null(values$NSBinding$NSBS)) {
            cat(sprintf("\nNS Binding Seed:\t%8.5f\n",values$NSBinding$NSBS))
          } else {
            cat("\n")
          }
          #Next print nucleotide seed, if it exists:
          if (!is.null(values[[currMode]]$NS)) {
            cat("Nucleotide Seed:\n")
            .print.psam(values[[currMode]]$NS, NULL, k, NULL)
          }
          #Next print dinucleotide seed, if it exists:
          if (!is.null(values[[currMode]]$DS)) {
            cat("\nDinucleotide Seed:\n")
            .print.psam(values[[currMode]]$DS, NULL, k-1, NULL)
          }
          #Next print shape seed, if it exists:
          if (!is.null(values[[currMode]]$SS)) {
            cat("\nShape Seed:\n")
            .print.psam(values[[currMode]]$SS, NULL, k, NULL)
          }
        }
      } else {
        #First print NS binding, if it exists:
        if (fit$NS) {
          if (!is.null(values$NSBE)) {
            cat(sprintf("NS Binding Beta:\t%7.4f \u00B1 %-7.4f\n",values$NSB, values$NSBE))
          } else {
            cat(sprintf("Nucleotide Betas:\t%8.5f\n",values$NSB))
          }
        }
        #Next print nucleotide PSAMs:
        cat("Nucleotide Betas:\n")
        .print.psam(values$NB, values$NE, k, fit$PSAM)
        #Next print dinucleotide PSAMs:
        if (!is.null(values$DB)) {
          cat("\nDinucleotide Betas:\n")
          .print.psam(values$DB, values$DE, k-1, NULL)
        }
        #Next print shape PSAM:
        if (!is.null(values$SB)) {
          cat("\nShape Betas:\n")
          .print.psam(values$SB, values$SE, k, NULL)
        }
        #Next print NS Seed, if it exists:
        if (!is.null(values$NSBS)) {
          cat(sprintf("\nNS Binding Seed:\t%8.5f\n",values$NSBS))
        }
        #Next print nucleotide seed, if it exists:
        if (!is.null(values$NS)) {
          cat("Nucleotide Seed:\n")
          .print.psam(values$NS, NULL, k, NULL)
        }
        #Next print dinucleotide seed, if it exists:
        if (!is.null(values$DS)) {
          cat("\nDinucleotide Seed:\n")
          .print.psam(values$DS, NULL, k-1, NULL)
        }
        #Next print shape seed, if it exists:
        if (!is.null(values$SS)) {
          cat("\nShape Seed:\n")
          .print.psam(values$SS, NULL, k, NULL)
        }
      }
      cat("------------------------------------------------------------------------------------------------------------\n")
    }
  }
}

.print.psam = function(value1, value2 = NULL, positions, idString = NULL) {
  if (is.null(idString)) {
    nucArray = as.character(1:positions)
  } else {
    nucArray = substring(idString, seq(1,nchar(idString),1), seq(1,nchar(idString),1))
  }
  matrix1 = t(matrix(value1, ncol=positions))
  if (is.null(value2)) {
    for (row in 1:positions) {
      cat(sprintf("%s   ",nucArray[row]))
      for (col in 1:ncol(matrix1)) {
        cat(sprintf("%8.5f   ",matrix1[row, col]))
      }
      cat(sprintf("%s\n",nucArray[row]))
    }
  } else {
    matrix2 = t(matrix(value2, ncol=positions))
    for (row in 1:positions) {
      cat(sprintf("%s   ",nucArray[row]))
      for (col in 1:ncol(matrix1)) {
        cat(sprintf("%7.4f \u00B1 %-7.4f  ",matrix1[row, col], matrix2[row,col]))
      }
      cat(sprintf("%s\n",nucArray[row]))
    }
  }
}

#does not work with multi-mode fits (notion of window is unknown)
plot.likelihoods = function(fits, l, test=FALSE) {
  old = list(fits$Information)
  for (colIdx in c(5, 7:12)) {
    new = NULL
    for (currList in 1:length(old)) {
      new = c(new, by(old[[currList]], old[[currList]][,colIdx], function(x) x))
    }
    old = new
  }
  tot.divs = length(old)
  likelihoods = vector(mode="list", length=tot.divs)
  fit.types = NULL
  for (i in 1:tot.divs) {
    topRow = old[[i]][1,]
    dat.string = c(topRow$k,sapply(7:12, function (x) if (as.logical(topRow[x])) {names(topRow)[x]} else {NA}))
    dat.string = dat.string[!is.na(dat.string)]
    fit.types = c(fit.types, paste(dat.string,collapse="-"))
    output = matrix(nrow=nrow(old[[i]]), ncol=3)
    output[,1] = l+2*old[[i]]$Flank-as.numeric(old[[i]]$k)+1
    output[,2] = old[[i]]$Shift
    if (test) {
      output[,3] = old[[i]]$TestLPerRead
    } else {
      output[,3] = old[[i]]$TrainLPerRead
    }
    likelihoods[[i]] = output
  }
  window.range = range(l-as.numeric(fits$Information$k)+2*fits$Information$Flank+1)
  if (test) {
    lik.range = range(fits$Information$TestLPerRead)
    plot(x = 0, y=0, type="n", ylim=lik.range, xlim=window.range, xlab="Windows", ylab="Testing -Log Likelihood")
  } else {
    lik.range = range(fits$Information$TrainLPerRead)
    plot(x = 0, y=0, type="n", ylim=lik.range, xlim=window.range, xlab="Windows", ylab="Training -Log Likelihood")
  }
  col.type = NULL
  lty.type = NULL
  for (i in 1:tot.divs) {
    col.type = 1 #c(col.type, (i %% 3)+1)
    lty.type = c(lty.type, i) #c(lty.type, (1+((i-1) %/% 3)))
    curr.break = by(likelihoods[[i]], likelihoods[[i]][,2], function(x) x)
    for (b in 1:length(curr.break)) {
      lines(x=curr.break[[b]][,1], y=curr.break[[b]][,3], col="grey", lty=4)
      points(x=curr.break[[b]][,1], y=curr.break[[b]][,3], col=1, pch=i, cex=.75)
    }
  }
  legend("topright", legend=fit.types, col=col.type, pch=lty.type, bty="n", ncol = 2)
}

#Works on either a single fit or list of fits + index
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

plot.compare = function(fitA, fitB, indexA=NULL, indexB=NULL, modeIdxA=NULL, modeIdxB=NULL, rcA=FALSE, rcB=FALSE, offset=NULL, labA=NULL, labB=NULL) {
  #See if fit labels have been provided
  if (is.null(labA)) {
    fitA.name = deparse(substitute(fitA))
  } else {
    fitA.name = labA
  }
  if (is.null(labB)) {
    fitB.name = deparse(substitute(fitB))
  } else {
    fitB.name = labB
  }
  #See if a specific fit has been provided or if an index is to be used
  if (is.null(indexA)) {
    #If a specific fit has been provided, handle multiple modes
    if (length(grep("Mode", names(fitA)))>0) {
      #Multi-mode fit provided; ensure that mode index is provided
      if (is.null(modeIdxA)) {
        stop("Multi-Mode Fit Detected: Mode Index Required for Fit A")
      } else {
        fitA = fitA[[modeIdxA]]
      }
    }
    fit.infoA = NULL
    motif.seq = NULL
  } else {
    #Check to see if fitA is a multi-mode fit
    fit.infoA = fitA$Information[indexA,]
    if (fit.infoA$k=="Multi") {
      A.isMulti = TRUE
      #Require mode index
      if (is.null(modeIdxA)) {
        stop("Multi-Mode Fit Detected: Mode Index Required for Fit A")
      } else {
        fitA = fitA$Values[[indexA]][[modeIdxA]]
        motif.seq = NULL
      }
    } else {
      fitA = fitA$Values[[indexA]]
      motif.seq = substring(fit.infoA$PSAM, 1:as.numeric(fit.infoA$k), 1:as.numeric(fit.infoA$k))
    }
  } 
  if (is.null(indexB)) {
    #If a specific fit has been provided, handle multiple modes
    if (length(grep("Mode", names(fitB)))>0) {
      #Multi-mode fit provided; ensure that mode index is provided
      if (is.null(modeIdxB)) {
        stop("Multi-Mode Fit Detected: Mode Index Required for Fit B")
      } else {
        fitB = fitB[[modeIdxB]]
      }
    }
    fit.infoB=NULL
  } else {
    #Check to see if fitB is a multi-mode fit
    fit.infoB = fitB$Information[indexB,]
    if (fit.infoB$k=="Multi") {
      B.isMulti = TRUE
      #Require mode index
      if (is.null(modeIdxB)) {
        stop("Multi-Mode Fit Detected: Mode Index Required for Fit B")
      } else {
        fitB = fitB$Values[[indexB]][[modeIdxB]]
        motif.seq = NULL
      }
    } else {
      fitB = fitB$Values[[indexB]]
    }
  }
  if (is.null(offset)) {
    a = motif.aligner(fitA, fitB)
    offset = a$Offset
    rcB = a$B.RC
  }
  if (rcA) {
    fitA = .rc(fitA)
  }
  if (rcB) {
    fitB = .rc(fitB)
  }
  #First Process & Align Nuc, Dinuc and Shape motifs
  nucA = fitA$NB
  nucB = fitB$NB
  if (offset<0) {
    nucA = c(rep(NA, 4*abs(offset)), nucA)
  } else if (offset>=0) {
    nucB = c(rep(NA, 4*offset), nucB)
  }
  if (which.min(c(length(nucA), length(nucB)))==1) {
    nucA = c(nucA, rep(NA, (length(nucB)-length(nucA))))
  } else {
    nucB = c(nucB, rep(NA, (length(nucA)-length(nucB))))
  }
  if (!is.null(fitA$NE) && !is.null(fitB$NE)) {
    nucA.SE = fitA$NE
    nucB.SE = fitB$NE
    if (offset<0) {
      nucA.SE = c(rep(NA, 4*abs(offset)), nucA.SE)
    } else if (offset>=0) {
      nucB.SE = c(rep(NA, 4*offset), nucB.SE)
    }
    if (which.min(c(length(nucA.SE), length(nucB.SE)))==1) {
      nucA.SE = c(nucA.SE, rep(NA, (length(nucB.SE)-length(nucA.SE))))
    } else {
      nucB.SE = c(nucB.SE, rep(NA, (length(nucA.SE)-length(nucB.SE))))
    }
  } 
  l = length(nucA)/4
  nuc = data.frame(position=as.character(rep(1:l, each=4)), 
                   base=as.character(rep(1:4, l)), nucA=nucA, nucB=nucB)
  nuc = nuc[complete.cases(nuc),]
  is.dinuc = FALSE
  dinucA = fitA$DB
  dinucB = fitB$DB
  if (!is.null(dinucA) && !is.null(dinucB)) {
    is.dinuc = TRUE
    if (offset<0) {
      dinucA = c(rep(NA, 16*abs(offset)), dinucA)
    } else if (offset>=0) {
      dinucB = c(rep(NA, 16*offset), dinucB)
    }
    if (which.min(c(length(dinucA), length(dinucB)))==1) {
      dinucA = c(dinucA, rep(NA, (length(dinucB)-length(dinucA))))
    } else {
      dinucB = c(dinucB, rep(NA, (length(dinucA)-length(dinucB))))
    }
    dinuc = data.frame(position=as.character(rep(1:(l-1), each=16)),
                       base=as.character(rep(1:16, (l-1))), dinucA=dinucA, dinucB=dinucB)
    dinuc = dinuc[complete.cases(dinuc),]
  }
  l = nrow(nuc)/4
  p1 = ggplot(nuc, aes(x=nucA, y=nucB, color=position, shape=base)) +
    geom_point(size=3) +
    coord_fixed(ratio=1, xlim=range(nuc$nucA, nuc$nucB), ylim=range(nuc$nucA, nuc$nucB)) +
    geom_abline(slope=1) +
    scale_colour_hue(name="Position", breaks=as.character(1:l), labels=as.character(1:l)) +
    scale_shape_discrete(name="", labels=c("A","C","G","T"))+
    guides(col=guide_legend(nrow=4)) + 
    labs(title="Nucleotide Features", x=fitA.name, y=fitB.name) + 
    theme(text=element_text(size=17, family="Helvetica"))
  if (is.dinuc) {
    p2 = ggplot(dinuc, aes(x=dinucA, y=dinucB, color=position, shape=base)) +
      geom_point(size=3) +
      coord_fixed(ratio=1, xlim=range(dinuc$dinucA, dinuc$dinucB), ylim=range(dinuc$dinucA, dinuc$dinucB)) +
      geom_abline(slope=1) +
      scale_colour_discrete(name="Position", breaks=as.character(1:(l-1)), labels=as.character(1:(l-1)), guide=FALSE) +
      scale_shape_manual(name="", labels=c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"), 
                         values=c(1:16), guide=guide_legend(nrow = 4))+
      labs(title="Dinucleotide Features", x=fitA.name, y=fitB.name) + 
      theme(text=element_text(size=17, family="Helvetica"))
    .multiplot(p1, p2, cols=1)
  } else {
    print(p1)
  }
  is.shape = FALSE
  shapeA = fitA$SB
  shapeB = fitB$SB
  if (!is.null(shapeA) && !is.null(shapeB)) {
    block.sizeA = 4*length(fitA$SB)/length(fitA$NB)
    block.sizeB = 4*length(fitB$SB)/length(fitB$NB)
    if (block.sizeA!=block.sizeB) {         #CURRENTLY DESIGNED FOR ONLY 1 SHAPE FEATURE!!
      stop("The two fits do not have the same number of shape features")
    }
    is.shape = TRUE
    if (offset<0) {
      shapeA = c(rep(NA, block.sizeA*abs(offset)), shapeA)
    } else if (offset>=0) {
      shapeB = c(rep(NA, block.sizeA*offset), shapeB)
    }
    if (which.min(c(length(shapeA), length(shapeB)))==1) {
      shapeA = c(shapeA, rep(NA, (length(shapeB)-length(shapeA))))
    } else {
      shapeB = c(shapeB, rep(NA, (length(shapeA)-length(shapeB))))
    }
    #Extract error bars
    shapeA.SE = rep(NA, length(shapeA))
    shapeB.SE = rep(NA, length(shapeA))
    if(!is.null(fitA$SE) && !is.null(fitB$SE)) {
      shapeA.SE = fitA$SE
      shapeB.SE = fitB$SE
      if (offset<0) {
        shapeA.SE = c(rep(NA, block.sizeA*abs(offset)), shapeA.SE)
      } else if (offset>=0) {
        shapeB.SE = c(rep(NA, block.sizeA*offset), shapeB.SE)
      }
      if (which.min(c(length(shapeA.SE), length(shapeB.SE)))==1) {
        shapeA.SE = c(shapeA.SE, rep(NA, (length(shapeB.SE)-length(shapeA.SE))))
      } else {
        shapeB.SE = c(shapeB.SE, rep(NA, (length(shapeA.SE)-length(shapeB.SE))))
      }
    }
    if (is.null(motif.seq)) {
      labels = rep("", length(shapeA))
    } else {
      labels = c(motif.seq, rep("", (length(shapeA)-length(motif.seq))))
    }
    shape = data.frame(group=as.character(rep(1:2, each=length(shapeA))), 
                       index=rep(1:length(shapeA), 2), shape=c(shapeA,shapeB), 
                       error=1.95*c(shapeA.SE, shapeB.SE))

    ggplot(shape, aes(x=index, y=shape, group=group, color=group)) +
      geom_line() +
      geom_point(size=2) +
      expand_limits(y=0) +
      theme_bw(base_size=15) +
      ylab(expression(paste("Betas (", Delta, Delta, "G/RT)"))) +
      xlab("Position") +
      geom_errorbar(aes(ymin=shape-error, ymax=shape+error), colour="black", width=.4) + 
      theme(legend.justification=c(1,0), legend.position=c(1,0), legend.title=element_blank()) + 
      scale_colour_hue(name="group", breaks=c("1","2"), labels=c(fitA.name, fitB.name)) + 
      ggtitle(paste("Shape Profile Comparison:", fitA.name, "vs.", fitB.name)) + 
      scale_x_discrete(labels=labels, limits=c(1:length(labels)))
  }
}

logo = function(fits, index=NULL, mode=NULL, rc=FALSE, display=TRUE, betas=NULL, save.path=NULL, isPDF=TRUE, title=NULL, ylim=NULL, fit.name=NULL,
                l.pad=0, r.pad=0, l.del=0, r.del=0) {
  if (is.null(fit.name)) {
    fit.name = deparse(substitute(fits))
  }
  #Check to see if the input is a core fit
  if (is.null(index)) {
    #Check to see if the input is just a vector of betas
    if (is.null(betas)) {
      fit.output = fits
      fit.info = NULL
    } else {
      if (rc) {
        fit.output = list(NB=rev(betas))
      } else {
        fit.output = list(NB=betas)
      }
      info.string = "manual-input"
      fit.info = NULL
    }
  } else {
    fit.output = fits$Values[[index]]
    fit.info = fits$Information[index,]
    info.string = c(fit.name, "k", fit.info$k,"f", fit.info$Flank, 
                    sapply(7:12, function (x) if (as.logical(fit.info[x])) {names(fit.info)[x]} else {NA}))
    info.string = info.string[!is.na(info.string)]
    if (rc) {
      fit.output = .rc(fit.output)
    }
  }
  #Handle Multi-Mode Fits
  if (class(fit.output[[1]])=="list" || length(fit.output)==9) {
    if (is.null(mode)) {
      stop("Multi-Mode Fit Detected: Mode Index Required for Fit")
    } else {
      fit.output = fit.output[[mode]]
      k = length(fit.output$NB)/4
      isMulti = TRUE
      info.string = c(fit.name, "k", fit.info$k,"f", fit.info$Flank, "m", mode, 
                      sapply(7:12, function (x) if (as.logical(fit.info[x])) {names(fit.info)[x]} else {NA}))
      info.string = info.string[!is.na(info.string)]
    }
  } else {
    isMulti = FALSE
  }
  if (isPDF) {
    f.name=paste0(paste(info.string,collapse="-"),".pdf")
  } else {
    f.name=paste0(paste(info.string,collapse="-"),".png")
  }
  if (is.null(fit.output$DB)) {
    motif = exp(as.numeric(fit.output$NB))
  } else {
    max.values = max.seq(fits, index, mode)
    rescale = as.numeric(max.values[1])
    top.string = strsplit(as.character(max.values[2]), split="")[[1]]
    PSAM = matrix(data = NA, nrow=4, ncol=length(top.string))
    row.names(PSAM) = c("A", "C", "G", "T")
    
    for (i in 1:length(top.string)) {
      #Loop over all characters at each position
      for (currChar in c("A", "C", "G", "T")) {
        curr.string = top.string
        curr.string[i] = currChar
        curr.string = DNAString(paste0(curr.string, collapse=""))
        #Evaluate model on this
        out = score.genome(curr.string, fits, index, mode)
        #Store
        PSAM[currChar, i] = out[1,1]/rescale
      }
    }
    if (rc) {
      motif = rev(as.numeric(PSAM))
    } else {
      motif = as.numeric(PSAM)
    }
  }
  motif = c(rep(0, l.pad*4), motif, rep(0, r.pad*4))
  if (l.del>0) {
    motif = motif[-c(1:(l.del*4))]
  }
  if (r.del>0) {
    motif = motif[-c((length(motif)-r.del*4+1):length(motif))]
  }
  k = length(motif)/4
  if (is.null(index) || isMulti) {
    motif.seq = as.character(1:k)
  } else {
    if (rc) {
      motif.seq = substring(as.character(reverseComplement(DNAString(fit.info$PSAM))), 1:k, 1:k)
    } else {
      motif.seq = substring(fit.info$PSAM, 1:k, 1:k)
    }
  }
  dim(motif) = c(4, k)
  motif = t(motif)
  system("mkdir ~/logo_tmp")
  write(paste0("<matrix_reduce>\n<psam_length>", k, "</psam_length>\n<psam>\n")
        , file="~/logo_tmp/data.xml", append=TRUE)
  for (i in 1:k) {
    write(paste0(motif[i,1],"\t\t",motif[i,2],"\t\t",motif[i,3],"\t\t",motif[i,4],"\n")
          , file="~/logo_tmp/data.xml", append=TRUE)
  }
  write("</psam>\n</matrix_reduce>"
        , file="~/logo_tmp/data.xml", append=TRUE)
  bash.string = paste0("~/REDUCE-Suite-v2.2/bin/LogoGenerator -output=~/logo_tmp/ -logo=", f.name)
  bash.string = paste(bash.string, "-style=ddG -file=~/logo_tmp/data.xml") 
  if (isPDF) {
    bash.string = paste(bash.string, "-format=pdf")
  } else {
    bash.string = paste(bash.string, "-format=png")
  }
  if (!is.null(title)) {
    bash.string = paste0(bash.string, " -title=",title)
  }
  if (!is.null(ylim)) {
    bash.string = paste0(bash.string, " -ymin=", ylim[1], " -ymax=", ylim[2])
  }
  system(bash.string)
  if (!is.null(save.path)) {
    system(paste0("mv ~/logo_tmp/",f.name, " ", save.path))
  }
  if (display && isPDF) {
    if (!is.null(save.path)) {
      if (dir.exists(save.path)) {
        .openPDF(paste0(save.path, f.name))
      } else {
        .openPDF(save.path)
      }
    } else {
      .openPDF(paste0("~/logo_tmp/", f.name))
    }
  }

  #Now make shape profile, if it exists
  if (!is.null(fit.output$SB)) {
    shape.motif = fit.output$SB
    shape.motif = t(matrix(data=shape.motif, ncol=k))
    n.features = ncol(shape.motif)
    plot.df = data.frame(matrix(data=NA, ncol=4, nrow=k*n.features))
    names(plot.df) = c("Feature","Position","Value","SD")
    for (i in 1:ncol(shape.motif)) {
      for (j in 1:k) {
        plot.df[((i-1)*k+j),] = c(as.character(i), j, shape.motif[j,i], NA)
      }
    }
    if (!is.null(fit.output$SE)) {
      shape.motif = fit.output$SE
      shape.motif = t(matrix(data=shape.motif, ncol=k))
      for (i in 1:ncol(shape.motif)) {
        for (j in 1:k) {
          plot.df[((i-1)*k+j),4] = as.numeric(1.95*shape.motif[j, i])
        }
      }
    }
    plot.df[,2] = as.numeric(plot.df[,2])
    plot.df[,3] = as.numeric(plot.df[,3])
    plot.df[,4] = as.numeric(plot.df[,4])
    
    f.name=paste0(paste(info.string,collapse="-"),"-shape.pdf")
    pd = position_dodge(0.1)
    ggplot(plot.df, aes(x=Position, y=Value, colour=Feature, group=Feature)) + 
      geom_line(position=pd) +
      geom_point(size=2, position=pd) +
      expand_limits(y=0) +
      theme_bw(base_size = 15) +
      ylab(expression(paste("Betas (", Delta, Delta, "G/RT)"))) +
      geom_errorbar(aes(ymin=Value-SD, ymax=Value+SD), colour="black", width=.4, position=pd) + 
      theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
      scale_colour_hue(name="Feature", breaks=as.character(1:n.features), labels=as.character(1:n.features)) +
      ggtitle(paste("Shape Profile:", paste(info.string, collapse=" "))) + 
      scale_x_discrete(labels=motif.seq, limits=c(1:k))
    ggsave(filename=f.name, path="~/logo_tmp/")
    
    if (!is.null(save.path)) {
      if (dir.exists(save.path)) {
        save.path = paste0(save.path, f.name)
      } else {
        new.dir.path = dirname(save.path)
        new.f.path = paste0(strsplit(basename(save.path), ".pdf")[[1]], "-shape.pdf")
        save.path = paste0(new.dir.path, "/",new.f.path)
      }
      system(paste0("mv ~/logo_tmp/",f.name, " ", save.path))
    }
    if (display) {
      if (!is.null(save.path)) {
        if (dir.exists(save.path)) {
          .openPDF(paste0(save.path, f.name))
        } else {
          .openPDF(save.path)
        }
      } else {
        .openPDF(paste0("~/logo_tmp/", f.name))
      }
    }
  }
  Sys.sleep(1)
  system("rm -rf ~/logo_tmp")
}

animate.traj = function(trajPath, outPath, k) {
  traj = read.table(trajPath, header=FALSE)
  dir.create(path="~/tmp_animate_dir")
  for (i in 1:nrow(traj)) {
    logo(betas=traj[i,1:(k*4)], save.path = paste0("~/tmp_animate_dir/",formatC(i, width=5, format="d", flag="0"),".png"), display = FALSE, isPDF = FALSE, title=i)
  }
  system(paste0("convert -delay 80 ~/tmp_animate_dir/*.png ", outPath))
  system("rm -rf ~/tmp_animate_dir")  
}

.openPDF <- function(f) {
  os <- .Platform$OS.type
  if (os=="windows")
    shell.exec(normalizePath(f))
  else {
    pdf <- getOption("pdfviewer", default='')
    if (nchar(pdf)==0)
      stop("The 'pdfviewer' option is not set. Use options(pdfviewer=...)")
    system2(pdf, args=c(f))
  }
}

plot.score.genome = function(genomicSequence, fits, index, mode=NULL, rc=FALSE, nPeaks=NULL, annotate=FALSE, rescale=NULL, 
                             genomicSequence2=NULL) {
  seq.string = deparse(substitute(genomicSequence))
  fit.string = deparse(substitute(fits))
  score = score.genome(genomicSequence, fits, index, mode, rc)
  if (fits[[1]]$k[index]=="Multi") {
    k = length(fits[[2]][[index]][[mode]]$NB)/4
  } else {
    k = as.numeric(fits[[1]]$k[index])
  }
  if (is.null(rescale)) {
    rescale = max.seq(fits, index, mode)$MaxAffinity
  }
  score = score/rescale
  max.score = max(score)
  #compare two tracks
  if (!is.null(genomicSequence2)) {
    if (length(genomicSequence2)!=length(genomicSequence)) {
      stop("Input genomic sequences are not of the same length.")
    }
    score2 = score.genome(genomicSequence2, fits, index, mode, rc)/rescale
    max.score = max(score, score2)
  }
  #create proper scale
  multiplier = 1
  while(max.score*multiplier<10) {
    multiplier = multiplier*10
  }
  #round up to nearest integer
  upper.bound = ceiling(max.score*multiplier)
  #find optimal step size
  divs = signif(seq(-upper.bound, upper.bound, length.out = 11)/multiplier, digits=2)
  #find top sites and rank them
  if (!is.null(nPeaks)) {
    idx = c(score[1,], score[2,])
    idx = cbind(idx, c(1:ncol(score), -(1:ncol(score))))
    idx = idx[order(-idx[,1]),]
    if (nPeaks>nrow(idx)) {
      nPeaks = nrow(idx)
    }
    Sequence = character(nPeaks)
    for (i in 1:nPeaks) {
      Sequence[i] = as.character(genomicSequence[abs(idx[i,2]):(abs(idx[i,2])+k-1)])
    }
    peaks = data.frame(Affinity=idx[1:nPeaks,1], Position=idx[1:nPeaks,2], Sequence)
    cat(paste0(fit.string," scores in ",seq.string,"\n"))
    print(peaks)
  }
  
  df = data.frame(Group=as.character(rep(1:2,each=ncol(score))), Position=c(1:ncol(score),1:ncol(score)), 
                  Affinity=c(score[1,], -score[2,]))
  if(!is.null(genomicSequence2)) {
    df = rbind(df, data.frame(Group=as.character(rep(3:4,each=ncol(score2))), Position=c(1:ncol(score2),1:ncol(score2)), 
                              Affinity=c(score2[1,], -score2[2,])))
    seq.string = paste0(seq.string, " and ", deparse(substitute(genomicSequence2)))
  }
  
  p = ggplot(df, aes(x=Position, y=Affinity, colour=Group)) +
    theme_bw() +
    geom_line() +
    ylab("Relative Affinity") +
    coord_fixed(ylim=c(-max.score,max.score)) +
    scale_y_continuous(breaks=divs, labels=format(abs(divs), nsmall=2))+
    theme(aspect.ratio=1, text=element_text(size=17, family="Helvetica")) +
    theme(legend.title=element_blank(), legend.position=c(0,0), legend.justification=c(0,0)) +
    theme(axis.title.x = element_text(vjust=0), axis.title.y = element_text(vjust=0)) +
    labs(title = paste0(fit.string," scores in ",seq.string))
  if (is.null(genomicSequence2)) {
    p = p+scale_color_manual(values=c("#000000", "#FF0000"), labels=c("Forward", "Reverse"))
  } else {
    p = p+scale_color_manual(values=c("#000000", "#FF0000", "blue", "green"), 
                             labels=c(paste0(deparse(substitute(genomicSequence)), " Forward"), 
                                      paste0(deparse(substitute(genomicSequence)), " Reverse"),
                                      paste0(deparse(substitute(genomicSequence2)), " Forward"), 
                                      paste0(deparse(substitute(genomicSequence2)), " Reverse"))) 
  }
    if (annotate && !is.null(nPeaks)) {
      xpos = abs(peaks$Position)+ncol(score)*.015
      ypos = sign(peaks$Position)*peaks$Affinity+max.score*.015
      p = p + annotate("text", x=xpos, y=ypos, label=as.character(1:nPeaks))
    }
  return(p)
}

score.genome = function(genomicSequence, fits, index, mode=NULL, rc=FALSE) {
  #Create sequences for rapid scoring
  fSeq = abs(toComplex(genomicSequence, c(A=1, C=2, G=3, T=4)))
  rSeq = abs(toComplex(reverseComplement(genomicSequence), c(A=1, C=2, G=3, T=4)))
  charFSeq = as.character(genomicSequence)
  charRSeq = as.character(reverseComplement(genomicSequence))
  l = length(genomicSequence)
  fit = fits[[2]][[index]]
  #Check to see if input is a multi-mode model
  if (fits[[1]]$k[index]=="Multi") {
    if (is.null(mode)) {
      stop("Multi-Mode Fit Detected: Mode Index Required")
    } else {
      fit = fit[[mode]]
      k = length(fit$NB)/4
    }
  } else {
    k = as.numeric(fits[[1]]$k[index])
  }
  adjK = k-1
  adjL = l+1
  nuc = fit$NB
  dim(nuc) = c(4,k)
  if (is.null(fit$DB)) {
    isDinuc = FALSE
    dinuc = NULL
  } else {
    isDinuc = TRUE
    dinuc = fit$DB
    dim(dinuc) = c(16, k-1)
    rownames(dinuc) = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
  }
  score = sapply(1:(l-adjK), FUN=function(x) .fastScore(substr(charFSeq, x, x+adjK), substr(charRSeq, adjL-x-adjK, adjL-x), 
                                                        fSeq[x:(x+adjK)], rSeq[(adjL-x-adjK):(adjL-x)], k, isDinuc, nuc, dinuc))
  if(rc) {
    score = score[c(2, 1),]
  }
  return(score)
}

#fast sequence scorer, optimized to work with score.genome
.fastScore = function(charFSeq, charRSeq, fSeq, rSeq, k, isDinuc, nuc, dinuc) {
  fTotal = 0
  rTotal = 0
  if (isDinuc) {
    for (j in 1:k) {
      fTotal = fTotal+nuc[fSeq[j],j]
      rTotal = rTotal+nuc[rSeq[j],j]
      if (j<k) {
        fTotal = fTotal+as.numeric(dinuc[substr(charFSeq, j, j+1),j])
        rTotal = rTotal+as.numeric(dinuc[substr(charRSeq, j, j+1),j])
      }
    }
  } else {
    for (j in 1:k) {
      fTotal = fTotal+nuc[fSeq[j],j]
      rTotal = rTotal+nuc[rSeq[j],j]
    }
  }
  output = c(exp(fTotal), exp(rTotal))
  return(output)
}

#Sequence Scorer, requires explicit nuc, dinuc formatting (mode independent)
score.seq = function(sequence, rSequence, k, isDinuc, nuc, dinuc) {
  if (isDinuc) {
    fSeq = abs(toComplex(sequence, c(A=1, C=2, G=3, T=4)))
    rSeq = abs(toComplex(rSequence, c(A=1, C=2, G=3, T=4)))
    fTotal = 0
    rTotal = 0
    for (j in 1:k) {
      fTotal = fTotal+nuc[fSeq[j],j]
      rTotal = rTotal+nuc[rSeq[j],j]
      if (j<k) {
        fTotal = fTotal+as.numeric(dinuc[as.character(sequence[j:(j+1)]),j])
        rTotal = rTotal+as.numeric(dinuc[as.character(rSequence[j:(j+1)]),j])
      }
    }
  } else {
    fSeq = abs(toComplex(sequence, c(A=1, C=2, G=3, T=4)))
    rSeq = abs(toComplex(rSequence, c(A=1, C=2, G=3, T=4)))
    fTotal = 0
    rTotal = 0
    for (j in 1:k) {
      fTotal = fTotal+nuc[fSeq[j],j]
      rTotal = rTotal+nuc[rSeq[j],j]
    }
  }
  output = c(exp(fTotal), exp(rTotal))
  return(output)
}

#Max-seq function (dynamic programming)
max.seq = function(fits, index, mode = NULL) {
  #Get betas and transform them into a matrix
  fit = fits[[2]][[index]]
  k = fits[[1]]$k[index]
  isMulti = (k=="Multi")
  if (isMulti) {
    output = NULL
    #Loop over all modes
    if (is.null(mode)) {
      modes = 1:length(fit)
    } else {
      modes = mode
    }
    for (currMode in modes) {
      if (names(fit)[currMode]=="NSBinding") {
        next
      }
      nuc = fit[[currMode]]$NB
      k = length(nuc)/4
      dim(nuc) = c(4, k)
      if (is.null(fit[[currMode]]$DB)) {
        isDinuc = FALSE
        dinuc = NULL
      } else {
        isDinuc = TRUE
        dinuc = fit[[currMode]]$DB
        dim(dinuc) = c(16, k-1)
        rownames(dinuc) = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
      }
      output = rbind(output, .maxSeqHelper(isDinuc, k, nuc, dinuc))
    }
    return(output)
  } else {
    nuc = fit$NB
    k = as.numeric(k)
    dim(nuc) = c(4,k)
    if (is.null(fit$DB)) {
      isDinuc = FALSE
      dinuc = NULL
    } else {
      isDinuc = TRUE
      dinuc = fit$DB
      dim(dinuc) = c(16, k-1)
      rownames(dinuc) = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
    }
    return(.maxSeqHelper(isDinuc, k, nuc, dinuc))
  }
}

.maxSeqHelper = function(isDinuc, k, nuc, dinuc) {
  #Initialize loop (position 1)
  max.list = nuc[,1] #(A C G T)
  char.list= c("A", "C", "G", "T")
  temp.list= char.list
  curr.list = matrix(data=0, 4, 4)
  #Loop over all positions
  for (currPos in 2:k) {
    #Loop over all previous bases
    for (prevBase in 1:4) {
      curr.list[,prevBase] = max.list[prevBase] + nuc[, currPos]
      if (isDinuc) {
        curr.list[,prevBase] = curr.list[,prevBase] + dinuc[((prevBase-1)*4+1):(prevBase*4), (currPos-1)]
      }
    }
    for (currBase in 1:4) {
      max.list[currBase] = max(curr.list[currBase,])
      temp.list[currBase] = paste0(char.list[which.max(curr.list[currBase,])], c("A", "C", "G", "T")[currBase])
    }
    char.list = temp.list
  }
  return(data.frame(MaxAffinity=exp(max(max.list)), BestSeq=char.list[which.max(max.list)], stringsAsFactors = FALSE))
}

.multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plot.R0.qq = function(fileName, lowerBound, upperBound, divSize, completeCoverage=TRUE, compare=NULL, ylim=NULL, xlim=NULL) {
  counts = read.table(fileName, header=FALSE)
  totCount = colSums(counts)
  nComparisons = nrow(compare)
  if (is.null(nComparisons)) {
    nComparisons = 1
    compare = t(as.matrix(compare))
  }
  if (completeCoverage) {
    totCount = sum(totCount*c(0:(ncol(counts)-1)))
  } else {
    totCount = sum(totCount*c(1:ncol(counts)))
  }
  x = seq(lowerBound, upperBound, by=divSize)
  #truncate zero rows
  all.zeros = apply(counts, 1, function(row) all(row==0))
  #find first and last non-zero elements
  for (i in 1:length(all.zeros)) {
    if (!all.zeros[i]){
      low = i
      break
    }
  }
  for (i in length(all.zeros):1) {
    if (!all.zeros[i]){
      high = i
      break
    }
  }
  x = x[low:high]
  #plot all combinations
  if (is.null(compare)) {
    compare = t(combn(c(1:ncol(counts)), m=2))
  } else if (completeCoverage) {
    compare = compare+1
  }
  #build comparisons
  ratios = NULL
  for (i in 1:nComparisons) {
    ratios = cbind(ratios, log10(counts[,compare[i,2]]/counts[,compare[i,1]]))
  }
  ratios = ratios[low:high,]
  #plot ratios
  colors = 1:nComparisons
  if (is.null(ylim)) {
    if (is.null(xlim)) {
      #ylim and xlim null
      matplot(x, ratios, type="l", lty=1, col=colors, xlab="log(w)", ylab="Log Ratio", main=fileName)
    } else {
      #ylim null and xlim not null
      matplot(x, ratios, type="l", lty=1, col=colors, xlab="log(w)", ylab="Log Ratio", main=fileName, xlim=xlim)
    }
  } else if (is.null(xlim)) {
    #ylim not null but xlim is null
    matplot(x, ratios, type="l", lty=1, col=colors, xlab="log(w)", ylab="Log Ratio", main=fileName, ylim=ylim)
  } else {
    #ylim and xlim not null
    matplot(x, ratios, type="l", lty=1, col=colors, xlab="log(w)", ylab="Log Ratio", main=fileName, xlim=xlim, ylim=ylim)
  }
  if (completeCoverage) {
    compare = compare-1
  }
  for (i in 1:nComparisons) {
    lines(x, .poisson.ratio(10^x*totCount, compare[i,1], compare[i,2]), col="grey", lty=2)
  }
  #Legend
  legend.text = NULL
  for (i in 1:nComparisons) {
    legend.text = c(legend.text, paste0(compare[i,], collapse="-"))
  }
  legend("bottomright", legend=legend.text, lty=1, col=colors, bty="n", ncol=ceiling(nrow(compare)/3), title="Rounds", cex=.75)
  locs = par('usr')
  text(x=.025*(locs[2]-locs[1])+locs[1], y=(.975*(locs[4]-locs[3])+locs[3]), labels=paste0("Total Reads: ", totCount), pos=4, cex=.75)
}

.poisson.ratio = function(x, countA, countB) {
  log10(x^(countB-countA)*factorial(countA)/factorial(countB))
}

ql = function(fits, index, mode=NULL, rc=FALSE) {
  logo(fits, index = index, mode = mode, rc = rc, fit.name=deparse(substitute(fits)))
}

aff = function(fits, index) {
  modes = max.seq(fits, index)$MaxAffinity
  nsb   = exp(fits[[1]]$NSBind[index])
  tot   = c(modes, nsb)
  tot   = tot/max(tot)
  return(log(tot))
}

o = function(fits, n=5) {
  info = fits[[1]]
  info = info[,-c(2, 3, 4, 17, 19, 20, 21)]
  info$k[info$k=="Multi"] = "-1"
  info = info[order(info$TrainLPerRead,-as.numeric(info$k), -info$Di, -as.numeric(as.POSIXct(info$Time))),]
  info$k[info$k=="-1"] = "Multi"
  return(info[1:n,])
}