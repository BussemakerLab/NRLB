\name{score.seq}
\alias{score.seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Score DNA sequence using a model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
score.seq(sequence, rSequence, k, isDinuc, nuc, dinuc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{
%%     ~~Describe \code{sequence} here~~
}
  \item{rSequence}{
%%     ~~Describe \code{rSequence} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{isDinuc}{
%%     ~~Describe \code{isDinuc} here~~
}
  \item{nuc}{
%%     ~~Describe \code{nuc} here~~
}
  \item{dinuc}{
%%     ~~Describe \code{dinuc} here~~
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
function (sequence, rSequence, k, isDinuc, nuc, dinuc) 
{
    if (isDinuc) {
        fSeq = abs(toComplex(sequence, c(A = 1, C = 2, G = 3, 
            T = 4)))
        rSeq = abs(toComplex(rSequence, c(A = 1, C = 2, G = 3, 
            T = 4)))
        fTotal = 0
        rTotal = 0
        for (j in 1:k) {
            fTotal = fTotal + nuc[fSeq[j], j]
            rTotal = rTotal + nuc[rSeq[j], j]
            if (j < k) {
                fTotal = fTotal + as.numeric(dinuc[as.character(sequence[j:(j + 
                  1)]), j])
                rTotal = rTotal + as.numeric(dinuc[as.character(rSequence[j:(j + 
                  1)]), j])
            }
        }
    }
    else {
        fSeq = abs(toComplex(sequence, c(A = 1, C = 2, G = 3, 
            T = 4)))
        rSeq = abs(toComplex(rSequence, c(A = 1, C = 2, G = 3, 
            T = 4)))
        fTotal = 0
        rTotal = 0
        for (j in 1:k) {
            fTotal = fTotal + nuc[fSeq[j], j]
            rTotal = rTotal + nuc[rSeq[j], j]
        }
    }
    output = c(exp(fTotal), exp(rTotal))
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
