\name{plot.likelihoods}
\alias{plot.likelihoods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot likelihoods
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.likelihoods(fits, l, test = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fits}{
%%     ~~Describe \code{fits} here~~
}
  \item{l}{
%%     ~~Describe \code{l} here~~
}
  \item{test}{
%%     ~~Describe \code{test} here~~
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
function (fits, l, test = FALSE) 
{
    old = list(fits$Information)
    for (colIdx in c(5, 7:12)) {
        new = NULL
        for (currList in 1:length(old)) {
            new = c(new, by(old[[currList]], old[[currList]][, 
                colIdx], function(x) x))
        }
        old = new
    }
    tot.divs = length(old)
    likelihoods = vector(mode = "list", length = tot.divs)
    fit.types = NULL
    for (i in 1:tot.divs) {
        topRow = old[[i]][1, ]
        dat.string = c(topRow$k, sapply(7:12, function(x) if (as.logical(topRow[x])) {
            names(topRow)[x]
        } else {
            NA
        }))
        dat.string = dat.string[!is.na(dat.string)]
        fit.types = c(fit.types, paste(dat.string, collapse = "-"))
        output = matrix(nrow = nrow(old[[i]]), ncol = 3)
        output[, 1] = l + 2 * old[[i]]$Flank - as.numeric(old[[i]]$k) + 
            1
        output[, 2] = old[[i]]$Shift
        if (test) {
            output[, 3] = old[[i]]$TestLPerRead
        }
        else {
            output[, 3] = old[[i]]$TrainLPerRead
        }
        likelihoods[[i]] = output
    }
    window.range = range(l - as.numeric(fits$Information$k) + 
        2 * fits$Information$Flank + 1)
    if (test) {
        lik.range = range(fits$Information$TestLPerRead)
        plot(x = 0, y = 0, type = "n", ylim = lik.range, xlim = window.range, 
            xlab = "Windows", ylab = "Testing -Log Likelihood")
    }
    else {
        lik.range = range(fits$Information$TrainLPerRead)
        plot(x = 0, y = 0, type = "n", ylim = lik.range, xlim = window.range, 
            xlab = "Windows", ylab = "Training -Log Likelihood")
    }
    col.type = NULL
    lty.type = NULL
    for (i in 1:tot.divs) {
        col.type = 1
        lty.type = c(lty.type, i)
        curr.break = by(likelihoods[[i]], likelihoods[[i]][, 
            2], function(x) x)
        for (b in 1:length(curr.break)) {
            lines(x = curr.break[[b]][, 1], y = curr.break[[b]][, 
                3], col = "grey", lty = 4)
            points(x = curr.break[[b]][, 1], y = curr.break[[b]][, 
                3], col = 1, pch = i, cex = 0.75)
        }
    }
    legend("topright", legend = fit.types, col = col.type, pch = lty.type, 
        bty = "n", ncol = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
