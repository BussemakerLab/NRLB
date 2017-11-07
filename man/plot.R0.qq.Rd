\name{plot.R0.qq}
\alias{plot.R0.qq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
R0 QQ Plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.R0.qq(fileName, lowerBound, upperBound, divSize, completeCoverage = TRUE, compare = NULL, ylim = NULL, xlim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileName}{
%%     ~~Describe \code{fileName} here~~
}
  \item{lowerBound}{
%%     ~~Describe \code{lowerBound} here~~
}
  \item{upperBound}{
%%     ~~Describe \code{upperBound} here~~
}
  \item{divSize}{
%%     ~~Describe \code{divSize} here~~
}
  \item{completeCoverage}{
%%     ~~Describe \code{completeCoverage} here~~
}
  \item{compare}{
%%     ~~Describe \code{compare} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
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
function (fileName, lowerBound, upperBound, divSize, completeCoverage = TRUE, 
    compare = NULL, ylim = NULL, xlim = NULL) 
{
    counts = read.table(fileName, header = FALSE)
    totCount = colSums(counts)
    nComparisons = nrow(compare)
    if (is.null(nComparisons)) {
        nComparisons = 1
        compare = t(as.matrix(compare))
    }
    if (completeCoverage) {
        totCount = sum(totCount * c(0:(ncol(counts) - 1)))
    }
    else {
        totCount = sum(totCount * c(1:ncol(counts)))
    }
    x = seq(lowerBound, upperBound, by = divSize)
    all.zeros = apply(counts, 1, function(row) all(row == 0))
    for (i in 1:length(all.zeros)) {
        if (!all.zeros[i]) {
            low = i
            break
        }
    }
    for (i in length(all.zeros):1) {
        if (!all.zeros[i]) {
            high = i
            break
        }
    }
    x = x[low:high]
    if (is.null(compare)) {
        compare = t(combn(c(1:ncol(counts)), m = 2))
    }
    else if (completeCoverage) {
        compare = compare + 1
    }
    ratios = NULL
    for (i in 1:nComparisons) {
        ratios = cbind(ratios, log10(counts[, compare[i, 2]]/counts[, 
            compare[i, 1]]))
    }
    ratios = ratios[low:high, ]
    colors = 1:nComparisons
    if (is.null(ylim)) {
        if (is.null(xlim)) {
            matplot(x, ratios, type = "l", lty = 1, col = colors, 
                xlab = "log(w)", ylab = "Log Ratio", main = fileName)
        }
        else {
            matplot(x, ratios, type = "l", lty = 1, col = colors, 
                xlab = "log(w)", ylab = "Log Ratio", main = fileName, 
                xlim = xlim)
        }
    }
    else if (is.null(xlim)) {
        matplot(x, ratios, type = "l", lty = 1, col = colors, 
            xlab = "log(w)", ylab = "Log Ratio", main = fileName, 
            ylim = ylim)
    }
    else {
        matplot(x, ratios, type = "l", lty = 1, col = colors, 
            xlab = "log(w)", ylab = "Log Ratio", main = fileName, 
            xlim = xlim, ylim = ylim)
    }
    if (completeCoverage) {
        compare = compare - 1
    }
    for (i in 1:nComparisons) {
        lines(x, .poisson.ratio(10^x * totCount, compare[i, 1], 
            compare[i, 2]), col = "grey", lty = 2)
    }
    legend.text = NULL
    for (i in 1:nComparisons) {
        legend.text = c(legend.text, paste0(compare[i, ], collapse = "-"))
    }
    legend("bottomright", legend = legend.text, lty = 1, col = colors, 
        bty = "n", ncol = ceiling(nrow(compare)/3), title = "Rounds", 
        cex = 0.75)
    locs = par("usr")
    text(x = 0.025 * (locs[2] - locs[1]) + locs[1], y = (0.975 * 
        (locs[4] - locs[3]) + locs[3]), labels = paste0("Total Reads: ", 
        totCount), pos = 4, cex = 0.75)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
