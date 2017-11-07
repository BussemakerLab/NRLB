\name{motif.aligner}
\alias{motif.aligner}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Align models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
motif.aligner(fitA, fitB, range = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitA}{
%%     ~~Describe \code{fitA} here~~
}
  \item{fitB}{
%%     ~~Describe \code{fitB} here~~
}
  \item{range}{
%%     ~~Describe \code{range} here~~
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
function (fitA, fitB, range = NULL) 
{
    if (is.null(range)) {
        range = floor(min(length(fitA$NB)/4, length(fitB$NB)/4)/2)
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
    if (glob.max == 1) {
        return(data.frame(B.RC = FALSE, Offset = range[max1], 
            R2 = A.B[max1]))
    }
    else {
        return(data.frame(B.RC = TRUE, Offset = range[max2], 
            R2 = A.rcB[max2]))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
