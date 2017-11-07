\name{print.fit}
\alias{print.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print fit results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.fit(fits, index = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fits}{
%%     ~~Describe \code{fits} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fits, index = NULL) 
{
    if (class(fits) == "data.frame") {
        if (is.null(index)) {
            index = 1:nrow(fits)
        }
        for (i in index) {
            fit = fits[i, ]
            cat(paste0("k: ", fit$k, ", f: ", fit$Flank, ", shift: ", 
                fit$Shift, "\t"))
            if (fit$NS) {
                cat("NSBinding  ")
            }
            if (fit$Di) {
                cat("Dinuc  ")
            }
            if (fit$MG) {
                cat("MGW  ")
            }
            if (fit$PT) {
                cat("ProT  ")
            }
            if (fit$HT) {
                cat("HelT  ")
            }
            if (fit$RO) {
                cat("Roll  ")
            }
            cat("\n")
            cat("Fit Steps:", fit$FitSteps, "\t\tFunction Calls:", 
                fit$FncCalls, "\tFitting Time:", fit$FitTime, 
                "\n")
            if (!is.na(fit$NuSym) || !is.na(fit$DiSym)) {
                if (!is.na(fit$NuSym)) {
                  cat(fit$NuSym, "\t")
                }
                if (!is.na(fit$DiSym)) {
                  cat(fit$DiSym)
                }
                cat("\n")
            }
            cat("Train -LL:", fit$TrainLPerRead, "\tTest -LL:", 
                fit$TestLPerRead, "\n")
            cat("Nucleotide PSAM:\t", fit$PSAM, "\n")
            cat("------------------------------------------------------------------------------------------------------------\n")
        }
    }
    else {
        if (is.null(index)) {
            index = 1:nrow(fits$Information)
        }
        for (i in index) {
            fit = fits$Information[i, ]
            values = fits$Values[[i]]
            isMulti = fit$k == "Multi"
            if (isMulti) {
                k = "Multi"
            }
            else {
                k = as.numeric(fit$k)
            }
            cat(paste0("k: ", k, ", f: ", fit$Flank, ", shift: ", 
                fit$Shift, "\t"))
            if (fit$NS) {
                cat("NSBinding  ")
            }
            if (fit$Di) {
                cat("Dinuc  ")
            }
            if (fit$MG) {
                cat("MGW  ")
            }
            if (fit$PT) {
                cat("ProT  ")
            }
            if (fit$HT) {
                cat("HelT  ")
            }
            if (fit$RO) {
                cat("Roll  ")
            }
            cat("\n")
            cat("Fit Steps:", fit$FitSteps, "\t\tFunction Calls:", 
                fit$FncCalls, "\tFitting Time:", fit$FitTime, 
                "\n")
            if (!is.na(fit$NuSym) || !is.na(fit$DiSym)) {
                if (!is.na(fit$NuSym)) {
                  cat(fit$NuSym, "\t")
                }
                if (!is.na(fit$DiSym)) {
                  cat(fit$DiSym)
                }
                cat("\n")
            }
            cat("Train -LL:", fit$TrainLPerRead, "\tTest -LL:", 
                fit$TestLPerRead, "\n\n")
            if (isMulti) {
                if (fit$NS) {
                  if (!is.null(values$NSBinding$NSBE)) {
                    cat(sprintf("NS Binding Beta:\t\%7.4f <U+00B1> \%-7.4f\n", 
                      values$NSBinding$NSB, values$NSBinding$NSBE))
                  }
                  else {
                    cat(sprintf("Nucleotide Betas:\t\%8.5f\n", 
                      values$NSBinding$NSB))
                  }
                  nModes = length(values) - 1
                }
                else {
                  nModes = length(values)
                }
                for (currMode in 1:nModes) {
                  k = length(values[[currMode]]$NB)/4
                  cat(paste0("Mode ", currMode, " Nucleotide Betas:\n"))
                  .print.psam(values[[currMode]]$NB, values[[currMode]]$NE, 
                    k, NULL)
                  if (!is.null(values[[currMode]]$DB)) {
                    cat(paste0("\nMode ", currMode, " Dinucleotide Betas:\n"))
                    .print.psam(values[[currMode]]$DB, values[[currMode]]$DE, 
                      k - 1, NULL)
                  }
                  if (!is.null(values[[currMode]]$SB)) {
                    cat(paste0("\nMode ", currMode, " Shape Betas:\n"))
                    .print.psam(values[[currMode]]$SB, values[[currMode]]$SE, 
                      k, NULL)
                  }
                  if (!is.null(values$NSBinding$NSBS)) {
                    cat(sprintf("\nNS Binding Seed:\t\%8.5f\n", 
                      values$NSBinding$NSBS))
                  }
                  else {
                    cat("\n")
                  }
                  if (!is.null(values[[currMode]]$NS)) {
                    cat("Nucleotide Seed:\n")
                    .print.psam(values[[currMode]]$NS, NULL, 
                      k, NULL)
                  }
                  if (!is.null(values[[currMode]]$DS)) {
                    cat("\nDinucleotide Seed:\n")
                    .print.psam(values[[currMode]]$DS, NULL, 
                      k - 1, NULL)
                  }
                  if (!is.null(values[[currMode]]$SS)) {
                    cat("\nShape Seed:\n")
                    .print.psam(values[[currMode]]$SS, NULL, 
                      k, NULL)
                  }
                }
            }
            else {
                if (fit$NS) {
                  if (!is.null(values$NSBE)) {
                    cat(sprintf("NS Binding Beta:\t\%7.4f <U+00B1> \%-7.4f\n", 
                      values$NSB, values$NSBE))
                  }
                  else {
                    cat(sprintf("Nucleotide Betas:\t\%8.5f\n", 
                      values$NSB))
                  }
                }
                cat("Nucleotide Betas:\n")
                .print.psam(values$NB, values$NE, k, fit$PSAM)
                if (!is.null(values$DB)) {
                  cat("\nDinucleotide Betas:\n")
                  .print.psam(values$DB, values$DE, k - 1, NULL)
                }
                if (!is.null(values$SB)) {
                  cat("\nShape Betas:\n")
                  .print.psam(values$SB, values$SE, k, NULL)
                }
                if (!is.null(values$NSBS)) {
                  cat(sprintf("\nNS Binding Seed:\t\%8.5f\n", 
                    values$NSBS))
                }
                if (!is.null(values$NS)) {
                  cat("Nucleotide Seed:\n")
                  .print.psam(values$NS, NULL, k, NULL)
                }
                if (!is.null(values$DS)) {
                  cat("\nDinucleotide Seed:\n")
                  .print.psam(values$DS, NULL, k - 1, NULL)
                }
                if (!is.null(values$SS)) {
                  cat("\nShape Seed:\n")
                  .print.psam(values$SS, NULL, k, NULL)
                }
            }
            cat("------------------------------------------------------------------------------------------------------------\n")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
