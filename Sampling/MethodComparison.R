#' ---
#' title: "Various Methods to Compare NRLB and DeepBind Model Performance"
#' author: "Chaitanya Rastogi"
#' output:
#'  html_document:
#'    fig_width: 12
#'    fig_height: 6
#'    df_print: paged
#'    code_folding: hide
#' ---

#'#Introduction
#'In order to compare the ability of NRLB and DeepBind to explain HT-SELEX data, 
#'we need to bin observed probe counts using the scores $S_i$ from either NRLB 
#'or DeepBind. These binned counts can be analyzed using the following methods:
#'
#'#####Enrichment
#'We can compare the ratios of binned R0 and R1 reads in every score bin $j$ by 
#'assuming
#'$$ \frac{R1_j}{R0_j} = \exp^{S_i-c}$$
#'where $c$ is a scaling factor incorporating library size and other unaccounted
#'effects. $c$ is found by minimizing the RMSE.
#'
#'Load information:
#+ load_info, echo=TRUE, message=TRUE, warning=TRUE
library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
setwd("~/Desktop/dbsampling/")
#setwd("~/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/dbsampling/")
info = read.table("dbIDs.tsv", header=TRUE, stringsAsFactors=FALSE)

#'Some helper functions below:
#+ helper_functions, echo=TRUE, message=TRUE, warning=TRUE
integrator = function(input) {
  sum = 0
  x = input[,1]
  y = input[,2]
  for (i in 2:length(y)) {
    sum = sum + (y[i]+y[i-1])/2*(x[i]-x[i-1])
  }
  return(sum)
}

bin.scores = function(counts, scores, breaks) {
  idx = counts!=0
  counts = counts[idx]
  scores = scores[idx]
  
  out.scores = scores
  for (i in 2:max(counts)) {
    for (j in 2:i) {
      out.scores = c(out.scores, scores[counts==i])
    }
  }
  hist.counts = hist(x = out.scores, breaks=breaks, plot=FALSE)$counts
  return(hist.counts)
}

compute.enrichment = function(R0, R1) {
  Enrichment = data.frame(Score=R1$Score,
                           Enrichment=log(R1$BinnedCounts/R0$BinnedCounts),
                           MinCount=do.call(pmin, as.data.frame(cbind(R1$BinnedCounts, R0$BinnedCounts))))
  Enrichment = Enrichment[complete.cases(Enrichment),]
  Enrichment = Enrichment[Enrichment$Enrichment!=Inf,]
  Enrichment = Enrichment[Enrichment$Enrichment!=-Inf,]
  #Adjust to minimize sum of squared errors
  offset = optimize(f = function(x) sum((Enrichment$Score-Enrichment$Enrichment+x)^2), 
                    interval = c(-40,40))
  Enrichment$Enrichment = Enrichment$Enrichment-offset$minimum
  RMSE = sqrt(offset$objective/nrow(Enrichment))
  temp = Enrichment[Enrichment$MinCount>50,]
  RMSE.50 = sqrt(sum((temp$Score-temp$Enrichment)^2)/nrow(temp))
  Range = range(Enrichment[,1:2])
  Range[1] = floor(Range[1])
  Range[2] = ceiling(Range[2])
  return(list(Enrichment=Enrichment, RMSE=RMSE, RMSE50=RMSE.50, Range=Range))
}

compute.means = function(counts, density) {
  #Compute Means
  Means = NULL
  for (i in 2:nrow(counts)) {
    Means = c(Means, counts$BinnedCounts[i]/
                ((density$ProbeCounts[i]+density$ProbeCounts[i-1])*(density$Score[i]-density$Score[i-1])/2))
  }
  Means = log10(Means/sum(counts$BinnedCounts))
  #Compute a Z
  Z = integrator(input = cbind(density$Score, exp(density$Score)*density$ProbeCounts))
  #Compute observation frequency
  Frequency = NULL
  for (i in 2:nrow(counts)) {
    Frequency = c(Frequency, (exp(density$Score[i])+exp(density$Score[i-1]))/(2*Z))
  }
  Frequency = log10(Frequency)
  Means = data.frame(Frequency, Means, MinCount=counts$BinnedCounts[2:nrow(counts)])
  Means = Means[complete.cases(Means),]
  Means = Means[Means$Means!=Inf,]
  Means = Means[Means$Means!=-Inf,]
  RMSE = sqrt(sum((Means$Frequency-Means$Means)^2)/nrow(Means))
  temp = Means[Means$MinCount>50,]
  RMSE.50 = sqrt(sum((temp$Frequency-temp$Means)^2)/nrow(temp))
  Range = range(Means[,1:2])
  Range[1] = floor(Range[1])
  Range[2] = ceiling(Range[2])
  return(list(Means=Means, RMSE=RMSE, RMSE50 = RMSE.50, Range=Range))
}

roc = function(scores, labels) {
  breaks = scores[order(scores)]
  P = sum(labels>0)
  N = sum(labels==0)
  
  FP.subset = scores[labels==0]
  TP.subset = scores[labels>0]
  
  FP = sapply(breaks, function(x) sum(FP.subset>=x))
  TP = sapply(breaks, function(x) sum(TP.subset>=x))
  
  return(data.frame(FalsePositive=FP/N, TruePositive=TP/P, Threshold=breaks))
}

means.plot = function(df) {
  r = df$Range
  out = ggplot(df$Means, aes(x=Frequency, y=Means, color=MinCount)) +
    coord_fixed(xlim = r, ylim = r, ratio=1) +
    theme_bw() +
    geom_abline(size=1, lty=2, col="grey") + 
    geom_point(alpha=.65) +
    scale_color_gradient("Counts", low="red", high="blue", limits=c(0,100), oob=squish) +
    xlab(expression(paste(log[10]," Predicted Sequencing Rate (per probe)"))) +
    ylab(expression(paste(log[10]," Observed Sequencing Rate (per probe)"))) + 
    annotate("text", x=r[1], y=r[2], 
             label=paste0("RMSE: ", round(df$RMSE, digits=3)), hjust=0, vjust=1) + 
    annotate("text", x=r[1], y=.95*r[2]+.05*r[1], 
             label=paste0("RMSE[50]: ", round(df$RMSE50, digits=3)), parse=TRUE, hjust=0, vjust=1) +  
  theme(legend.position=c(.9,.2), legend.title.align=.5)
  return(out)
}

enrichment.plot = function(df) {
  r = df$Range
  out = ggplot(df$Enrichment, aes(x=Score, y=Enrichment, color=MinCount)) +
    coord_fixed(xlim = r, ylim = r, ratio=1) +
    theme_bw() +
    geom_abline(size=1, lty=2, col="grey") + 
    geom_point(alpha=.65) +
    scale_color_gradient("Counts", low="red", high="blue", limits=c(0,100), oob=squish) +
    ylab("Log Enrichment") +
    xlab("Score") +   
    annotate("text", x=r[1], y=r[2], 
             label=paste0("RMSE: ", round(df$RMSE, digits=3)), hjust=0, vjust=1) + 
    annotate("text", x=r[1], y=.95*r[2]+.05*r[1], 
             label=paste0("RMSE[50]: ", round(df$RMSE50, digits=3)), parse=TRUE, hjust=0, vjust=1) +
  theme(legend.position=c(.9,.2), legend.title.align=.5)
  return(out)
}

density.plots = function(DOS, R0, R1) {
  divSize = DOS$Score[2]-DOS$Score[1]
  r = range(c(range(DOS$Score[DOS$ProbeCounts!=0]), range(DOS$Score[DOS$ProbeCounts!=0])))
  df = rbind(data.frame(Score = DOS$Score, Density = DOS$ProbeCounts/integrator(DOS), Type="DOS"),
             data.frame(Score = DOS$Score, Density = R0$BinnedCounts/sum(R0$BinnedCounts)/divSize, Type="R0 Counts"),
             data.frame(Score = DOS$Score, Density = R1$BinnedCounts/sum(R1$BinnedCounts)/divSize, Type="R1 Counts"))
  out = ggplot(df, aes(x = Score, y = Density, color = Type)) +
    coord_fixed(xlim = r) + 
    theme_bw() +
    geom_path(aes(linetype = Type)) + 
    scale_linetype_manual(values=c(1, 2, 2)) + 
    scale_color_manual(values=c("darkgrey", "red", "blue")) + 
    theme(legend.position=c(.85,.7), legend.title.align=.5, aspect.ratio = 1)
  return(out)
}

rescaled.dos = function(dos) {
  sorted.vals = unique(sort(dos)) 
  for (currCutoff in sorted.vals) {
    cloned.dos = dos
    cloned.dos[cloned.dos<currCutoff] = 0
    cloned.dos[cloned.dos>=currCutoff] = exp(cloned.dos[cloned.dos>=currCutoff] - currCutoff + .01)
    if (sum(range(cloned.dos) %in% Inf) == 0) {
      break
    }
  }
  return(cloned.dos)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
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

#'#Analysis over all proteins
#+ loop, echo=TRUE, message=TRUE, warning=TRUE, results = "asis"
AllEnrichments = NULL
AllMeans = NULL
AllROCs = NULL
Names = NULL
for (i in 1:30) { #nrow(info)) {
  protName = info$Name[i]
  count.file = paste0("CountScores/", protName, "_R0R1.tsv")
  score.file = paste0("CountScores/", protName, "_scores.tsv")
  dbdensity.file = paste0("DBDensities/",protName,".tsv")
  nrlbdensity.file = paste0("NRLBDensities/",protName,".tsv")
  #Check to see if count/score files exist
  if(!(file.exists(count.file) && file.exists(score.file) && file.exists(dbdensity.file))) {
    next
  }
  if (file.info(dbdensity.file)$size==0) {
    next
  }
  cat(paste0("  \n####",protName," \n"))
  Names = c(Names, protName)
  #Load all files
  counts = read.table(count.file, header = FALSE, stringsAsFactors = FALSE)
  names(counts) = c("Sequence", "R0Count", "R1Count")
  scores = read.table(score.file, header = TRUE, stringsAsFactors = FALSE)
  DOS.DB  = read.table(dbdensity.file, header = FALSE, stringsAsFactors = FALSE)
  DOS.NRLB= read.table(nrlbdensity.file, header = FALSE, stringsAsFactors = FALSE)
  
  #Compute DB/NRLB distribution and reset energies
  dbMax = max(DOS.DB$V1[DOS.DB$V2!=0])
  nrlbMax = max(DOS.NRLB$V1[DOS.NRLB$V2!=0])
  DOS.DB$V2 = rescaled.dos(DOS.DB$V2)
  DOS.NRLB$V2 = rescaled.dos(DOS.NRLB$V2)
  Density.DB = data.frame(Score=DOS.DB$V1-dbMax, ProbeCounts=(DOS.DB$V2*4^info$VarLen[i])/integrator(DOS.DB))
  Density.NRLB=data.frame(Score=DOS.NRLB$V1-nrlbMax, ProbeCounts=(DOS.NRLB$V2*4^info$VarLen[i])/integrator(DOS.NRLB))
  scores$DeepBind = scores$DeepBind - dbMax
  scores$NRLB = scores$NRLB - nrlbMax
  #TODO: Optimal histogram spacing?
  
  #Bin counts
  R0.DB = data.frame(Score=Density.DB$Score, 
                    BinnedCounts=c(0, bin.scores(counts = counts$R0Count, 
                                                 scores = scores$DeepBind, breaks = Density.DB$Score)))
  R1.DB = data.frame(Score=Density.DB$Score, 
                    BinnedCounts=c(0, bin.scores(counts = counts$R1Count, 
                                                 scores = scores$DeepBind, breaks = Density.DB$Score)))
  R0.NRLB = data.frame(Score=Density.NRLB$Score, 
                    BinnedCounts=c(0, bin.scores(counts = counts$R0Count, 
                                                 scores = scores$NRLB, breaks = Density.NRLB$Score)))
  R1.NRLB = data.frame(Score=Density.NRLB$Score, 
                    BinnedCounts=c(0, bin.scores(counts = counts$R1Count, 
                                                 scores = scores$NRLB, breaks = Density.NRLB$Score)))
  
  #Plot densities
  p1 = density.plots(Density.NRLB, R0.NRLB, R1.NRLB) + ggtitle("NRLB Densities")
  p2 = density.plots(Density.DB, R0.DB, R1.DB) + ggtitle("DeepBind Densities")
  multiplot(p1, p2, cols=2)
  #grid.arrange(p1, p2, ncol = 2)
  
  #Compute log enrichment for NRLB and DB and plot them
  NRLB.enrichment = compute.enrichment(R0.NRLB, R1.NRLB)
  DB.enrichment = compute.enrichment(R0.DB, R1.DB)
  p.nrlb = enrichment.plot(NRLB.enrichment) + ggtitle("NRLB Enrichment")
  p.db = enrichment.plot(DB.enrichment) + ggtitle("DeepBind Enrichment")
  multiplot(p.nrlb, p.db, cols=2)
  #grid.arrange(p.nrlb, p.db, ncol = 2)
  AllEnrichments = rbind(AllEnrichments, c(NRLB.enrichment$RMSE50, DB.enrichment$RMSE50))
  
  #Compute means/poisson rate
  NRLB.means = compute.means(R1.NRLB, Density.NRLB)
  DB.means = compute.means(R1.DB, Density.DB)
  p.nrlb = means.plot(NRLB.means) + ggtitle("NRLB Poisson Means")
  p.db = means.plot(DB.means) + ggtitle("DeepBind Poisson Means")
  multiplot(p.nrlb, p.db, cols=2)
  #grid.arrange(p.nrlb, p.db, ncol = 2)
  
  #Compute optimal means
  NRLB.opt = read.table(paste0("NRLBPoissonDensities/", protName, ".tsv"), header=FALSE, stringsAsFactors = FALSE)
  Means = data.frame(Frequency=NRLB.opt$V1, Means = log10((NRLB.opt$V3/NRLB.opt$V2)/sum(NRLB.opt$V3)), MinCount=NRLB.opt$V3)
  Means = Means[complete.cases(Means),]
  Means = Means[Means$Means!=Inf,]
  Means = Means[Means$Means!=-Inf,]
  RMSE = sqrt(sum((Means$Frequency-Means$Means)^2)/nrow(Means))
  temp = Means[Means$MinCount>50,]
  RMSE.50 = sqrt(sum((temp$Frequency-temp$Means)^2)/nrow(temp))
  if (is.nan(RMSE.50)) {
    RMSE.50 = RMSE
  }
  Range = range(Means[,1:2])
  Range[1] = floor(Range[1])
  Range[2] = ceiling(Range[2])
  NRLB.opt = list(Means=Means, RMSE=RMSE, RMSE50 = RMSE.50, Range=Range)
  AllMeans = rbind(AllMeans, c(NRLB.opt$RMSE50, NRLB.means$RMSE50, DB.means$RMSE50))
  p1 = means.plot(NRLB.opt) + ggtitle("NRLB Optimal Poisson Means")
  
  #Compute ROC Curves
  categories = read.table(paste0("CountScores/",protName,"_ROC.tsv"), header = FALSE, stringsAsFactors = FALSE)
  scores = read.table(paste0("CountScores/",protName,"_ROC_scores.tsv"), header=TRUE, stringsAsFactors = FALSE)
  ROC.DB = roc(scores = scores$DeepBind, labels = categories$V2)
  ROC.NRLB = roc(scores = scores$NRLB, labels = categories$V2)
  AUROC.DB = -integrator(ROC.DB)
  AUROC.NRLB = -integrator(ROC.NRLB)
  AllROCs = rbind(AllROCs, c(AUROC.NRLB, AUROC.DB))
  
  #Plot
  df = rbind(cbind(ROC.DB, Method="DeepBind"), cbind(ROC.NRLB, Method="NRLB"))
  p2 = ggplot(df, aes(x=FalsePositive, y=TruePositive, color=Method)) +
   coord_fixed(ratio = 1, xlim=c(0,1), ylim=c(0,1)) +
   geom_path(size=.5, show.legend = TRUE) +
   theme_bw() +
   expand_limits(y=0) +
   scale_colour_manual(values=c("black", "red"), labels=c(paste0("DeepBind: ", round(AUROC.DB, 3)),
                                                          paste0("NRLB: ", round(AUROC.NRLB, 3)))) +
   geom_abline(slope=1, intercept = 0, linetype=2, alpha=.5, size=1) +
   xlab("False Positive Rate") +
   ylab("True Positive Rate") +
   ggtitle("ROC Curve") +
   theme(legend.position=c(.85,.2), legend.title.align=.5)
  multiplot(p1, p2, cols=2)
  #grid.arrange(p1, p2, ncol = 2)
  cat("  \n")
}

#'##Compiled Results
#+ final, echo=TRUE, message=TRUE, warning=TRUE
AllROCs = data.frame(AllROCs)
names(AllROCs) = c("NRLB", "DeepBind")
AllROCs = cbind(AllROCs, Names)
AllMeans = data.frame(AllMeans)
names(AllMeans) = c("NRLBOpt", "NRLB", "DeepBind")
AllMeans = cbind(AllMeans, Names)
AllEnrichments = data.frame(AllEnrichments)
names(AllEnrichments) = c("NRLB", "DeepBind")
AllEnrichments = cbind(AllEnrichments, Names)

r = range(AllEnrichments[,1:2])
ggplot(AllEnrichments, aes(x=NRLB, y=DeepBind, label=Names)) + 
  coord_fixed(ratio = 1, xlim=r, ylim=r) + 
  geom_abline(slope=1, intercept = 0, linetype=2, alpha=.5, size=.5) + 
  ggtitle("All Enrichments") +
  theme_bw() + 
  geom_point() + 
  geom_text(check_overlap = TRUE, vjust=1, hjust=0, nudge_x = 0.025*(r[2]-r[1])) +
  xlab("NRLB RMSE") +
  ylab("DeepBind RMSE")

r = range(AllROCs[,1:2])
ggplot(AllROCs, aes(y=NRLB, x=DeepBind, label=Names)) + 
  coord_fixed(ratio = 1, xlim=r, ylim=r) + 
  geom_abline(slope=1, intercept = 0, linetype=2, alpha=.5, size=.5) + 
  ggtitle("All ROCs") +
  theme_bw() + 
  geom_point() + 
  geom_text(check_overlap = TRUE, vjust=1, hjust=0, nudge_x = 0.025*(r[2]-r[1])) +
  ylab("NRLB ROC") +
  xlab("DeepBind ROC")

r = range(AllMeans[,1:3])
ggplot(AllMeans, aes(x=NRLB, y=DeepBind, color=NRLBOpt, label=Names)) + 
  coord_fixed(ratio = 1, xlim=r, ylim=r) + 
  geom_abline(slope=1, intercept = 0, linetype=2, alpha=.5, size=.5) + 
  ggtitle("All Poisson Means") +
  theme_bw() + 
  geom_point() + 
  scale_color_gradient("Optimal\nRMSE", low="blue", high="red", limits=r, oob=squish) + 
  geom_text(check_overlap = TRUE, vjust=1, hjust=0, nudge_x = 0.025*(r[2]-r[1]), color="black") +
  xlab("NRLB RMSE") +
  ylab("DeepBind RMSE")