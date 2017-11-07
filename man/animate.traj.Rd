\name{animate.traj}
\alias{animate.traj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Animate solution trajectory with PSAMs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
animate.traj(trajPath, outPath, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trajPath}{
%%     ~~Describe \code{trajPath} here~~
}
  \item{outPath}{
%%     ~~Describe \code{outPath} here~~
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
function (trajPath, outPath, k) 
{
    traj = read.table(trajPath, header = FALSE)
    dir.create(path = "~/tmp_animate_dir")
    for (i in 1:nrow(traj)) {
        logo(betas = traj[i, 1:(k * 4)], save.path = paste0("~/tmp_animate_dir/", 
            formatC(i, width = 5, format = "d", flag = "0"), 
            ".png"), display = FALSE, isPDF = FALSE, title = i)
    }
    system(paste0("convert -delay 80 ~/tmp_animate_dir/*.png ", 
        outPath))
    system("rm -rf ~/tmp_animate_dir")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
