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

#' Generate PSAM
#' 
#' @param fits TO_BE_ADDED
#' @param index TO_BE_ADDED
#' @param mode TO_BE_ADDED
#' @param rc TO_BE_ADDED
#' @param display TO_BE_ADDED
#' @param betas TO_BE_ADDED
#' @param save.path TO_BE_ADDED
#' @param isPDF TO_BE_ADDED
#' @param title TO_BE_ADDED
#' @param ylim TO_BE_ADDED
#' @param fit.name TO_BE_ADDED
#' @param l.pad TO_BE_ADDED
#' @param r.pad TO_BE_ADDED
#' @param l.del TO_BE_ADDED
#' @param r.del TO_BE_ADDED
#' @return TO_BE_ADDED
#' 
#' @examples
#'
#' @export
#' 
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
