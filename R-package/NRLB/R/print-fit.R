#' Print fit results
#' 
#' @param fits 
#' @param index
#' @return
#' 
#' @examples
#'
#' @export
#' 
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
