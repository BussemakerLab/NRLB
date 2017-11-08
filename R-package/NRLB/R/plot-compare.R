#' Compare two models
#' 
#' @param fitA 
#' @param fitB 
#' @param indexA 
#' @param indexB 
#' @param modeIdxA 
#' @param modeIdxB 
#' @param rcA 
#' @param rcB 
#' @param offset
#' @param labA
#' @param labB
#' @return
#'
#' @examples
#'
#' @export
#' 
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
