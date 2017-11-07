\name{score.genome}
\alias{score.genome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Score genome with a fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
score.genome(genomicSequence, fits, index, mode = NULL, rc = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genomicSequence}{
%%     ~~Describe \code{genomicSequence} here~~
}
  \item{fits}{
%%     ~~Describe \code{fits} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
  \item{mode}{
%%     ~~Describe \code{mode} here~~
}
  \item{rc}{
%%     ~~Describe \code{rc} here~~
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
function (genomicSequence, fits, index, mode = NULL, rc = FALSE) 
{
    fSeq = abs(toComplex(genomicSequence, c(A = 1, C = 2, G = 3, 
        T = 4)))
    rSeq = abs(toComplex(reverseComplement(genomicSequence), 
        c(A = 1, C = 2, G = 3, T = 4)))
    charFSeq = as.character(genomicSequence)
    charRSeq = as.character(reverseComplement(genomicSequence))
    l = length(genomicSequence)
    fit = fits[[2]][[index]]
    if (fits[[1]]$k[index] == "Multi") {
        if (is.null(mode)) {
            stop("Multi-Mode Fit Detected: Mode Index Required")
        }
        else {
            fit = fit[[mode]]
            k = length(fit$NB)/4
        }
    }
    else {
        k = as.numeric(fits[[1]]$k[index])
    }
    adjK = k - 1
    adjL = l + 1
    nuc = fit$NB
    dim(nuc) = c(4, k)
    if (is.null(fit$DB)) {
        isDinuc = FALSE
        dinuc = NULL
    }
    else {
        isDinuc = TRUE
        dinuc = fit$DB
        dim(dinuc) = c(16, k - 1)
        rownames(dinuc) = c("AA", "AC", "AG", "AT", "CA", "CC", 
            "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", 
            "TT")
    }
    score = sapply(1:(l - adjK), FUN = function(x) .fastScore(substr(charFSeq, 
        x, x + adjK), substr(charRSeq, adjL - x - adjK, adjL - 
        x), fSeq[x:(x + adjK)], rSeq[(adjL - x - adjK):(adjL - 
        x)], k, isDinuc, nuc, dinuc))
    if (rc) {
        score = score[c(2, 1), ]
    }
    return(score)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
