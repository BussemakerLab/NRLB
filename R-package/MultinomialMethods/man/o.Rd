\name{o}
\alias{o}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Order fit results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
o(fits, n = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fits}{
%%     ~~Describe \code{fits} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (fits, n = 5) 
{
    info = fits[[1]]
    info = info[, -c(2, 3, 4, 17, 19, 20, 21)]
    info$k[info$k == "Multi"] = "-1"
    info = info[order(info$TrainLPerRead, -as.numeric(info$k), 
        -info$Di, -as.numeric(as.POSIXct(info$Time))), ]
    info$k[info$k == "-1"] = "Multi"
    return(info[1:n, ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
