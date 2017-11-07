\name{max.seq}
\alias{max.seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the highest affinity sequence for a model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
max.seq(fits, index, mode = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fits}{
%%     ~~Describe \code{fits} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
  \item{mode}{
%%     ~~Describe \code{mode} here~~
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
function (fits, index, mode = NULL) 
{
    fit = fits[[2]][[index]]
    k = fits[[1]]$k[index]
    isMulti = (k == "Multi")
    if (isMulti) {
        output = NULL
        if (is.null(mode)) {
            modes = 1:length(fit)
        }
        else {
            modes = mode
        }
        for (currMode in modes) {
            if (names(fit)[currMode] == "NSBinding") {
                next
            }
            nuc = fit[[currMode]]$NB
            k = length(nuc)/4
            dim(nuc) = c(4, k)
            if (is.null(fit[[currMode]]$DB)) {
                isDinuc = FALSE
                dinuc = NULL
            }
            else {
                isDinuc = TRUE
                dinuc = fit[[currMode]]$DB
                dim(dinuc) = c(16, k - 1)
                rownames(dinuc) = c("AA", "AC", "AG", "AT", "CA", 
                  "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", 
                  "TC", "TG", "TT")
            }
            output = rbind(output, .maxSeqHelper(isDinuc, k, 
                nuc, dinuc))
        }
        return(output)
    }
    else {
        nuc = fit$NB
        k = as.numeric(k)
        dim(nuc) = c(4, k)
        if (is.null(fit$DB)) {
            isDinuc = FALSE
            dinuc = NULL
        }
        else {
            isDinuc = TRUE
            dinuc = fit$DB
            dim(dinuc) = c(16, k - 1)
            rownames(dinuc) = c("AA", "AC", "AG", "AT", "CA", 
                "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", 
                "TC", "TG", "TT")
        }
        return(.maxSeqHelper(isDinuc, k, nuc, dinuc))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
