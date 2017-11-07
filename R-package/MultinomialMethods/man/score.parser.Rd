\name{score.parser}
\alias{score.parser}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read in binary files containing genomic scoring results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
score.parser(genome, score.dir, prot.name, chr.name, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genome}{
%%     ~~Describe \code{genome} here~~
}
  \item{score.dir}{
%%     ~~Describe \code{score.dir} here~~
}
  \item{prot.name}{
%%     ~~Describe \code{prot.name} here~~
}
  \item{chr.name}{
%%     ~~Describe \code{chr.name} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
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
function (genome, score.dir, prot.name, chr.name, k) 
{
    s = length(genome[[chr.name]])
    fwd.score = readBin(con = paste0(score.dir, prot.name, "_", 
        chr.name, "_F.dat"), what = "double", n = s, endian = "big")
    rev.score = readBin(con = paste0(score.dir, prot.name, "_", 
        chr.name, "_R.dat"), what = "double", n = s, endian = "big")
    flag = readBin(con = paste0(score.dir, prot.name, "_", chr.name, 
        "_Flag.dat"), what = "integer", n = s, size = 1, endian = "big")
    return(list(score = rbind(fwd.score, rev.score), flag = flag, 
        chr = chr.name, k = k))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
