\name{aff}
\alias{aff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate relative affinities for a fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
aff(fits, index)
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
function (fits, index) 
{
    modes = max.seq(fits, index)$MaxAffinity
    nsb = exp(fits[[1]]$NSBind[index])
    tot = c(modes, nsb)
    tot = tot/max(tot)
    return(log(tot))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
