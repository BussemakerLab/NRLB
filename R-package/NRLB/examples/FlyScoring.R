library(NRLB)

library(BSgenome.Dmelanogaster.UCSC.dm3)
genome = Dmelanogaster

# define E3N enhancer region from Crocker 2015 Cell paper
E3N = genome$chrX[4915195:4915486]     

#Define binding models
Ubx.file = system.file("extdata/models", "Exd-UbxIVa.csv", package = "NRLB")
Ubx.model = file.parser(Ubx.file)
#Index 20 corresponds to the best 18bp nucleotide model

#Score regions, create plots and provide sequences
nrlb.plot.score.genome(E3N, Ubx.model, index=20, mode=1, nPeaks = 10, annotate = TRUE)
