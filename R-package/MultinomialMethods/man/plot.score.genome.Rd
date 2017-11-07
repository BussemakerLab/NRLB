\name{plot.score.genome}
\alias{plot.score.genome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot affinities of a sequence
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.score.genome(genomicSequence, fits, index, mode = NULL, rc = FALSE, nPeaks = NULL, annotate = FALSE, rescale = NULL, genomicSequence2 = NULL)
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
  \item{nPeaks}{
%%     ~~Describe \code{nPeaks} here~~
}
  \item{annotate}{
%%     ~~Describe \code{annotate} here~~
}
  \item{rescale}{
%%     ~~Describe \code{rescale} here~~
}
  \item{genomicSequence2}{
%%     ~~Describe \code{genomicSequence2} here~~
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
function (genomicSequence, fits, index, mode = NULL, rc = FALSE, 
    nPeaks = NULL, annotate = FALSE, rescale = NULL, genomicSequence2 = NULL) 
{
    seq.string = deparse(substitute(genomicSequence))
    fit.string = deparse(substitute(fits))
    score = score.genome(genomicSequence, fits, index, mode, 
        rc)
    if (fits[[1]]$k[index] == "Multi") {
        k = length(fits[[2]][[index]][[mode]]$NB)/4
    }
    else {
        k = as.numeric(fits[[1]]$k[index])
    }
    if (is.null(rescale)) {
        rescale = max.seq(fits, index, mode)$MaxAffinity
    }
    score = score/rescale
    max.score = max(score)
    if (!is.null(genomicSequence2)) {
        if (length(genomicSequence2) != length(genomicSequence)) {
            stop("Input genomic sequences are not of the same length.")
        }
        score2 = score.genome(genomicSequence2, fits, index, 
            mode, rc)/rescale
        max.score = max(score, score2)
    }
    multiplier = 1
    while (max.score * multiplier < 10) {
        multiplier = multiplier * 10
    }
    upper.bound = ceiling(max.score * multiplier)
    divs = signif(seq(-upper.bound, upper.bound, length.out = 11)/multiplier, 
        digits = 2)
    if (!is.null(nPeaks)) {
        idx = c(score[1, ], score[2, ])
        idx = cbind(idx, c(1:ncol(score), -(1:ncol(score))))
        idx = idx[order(-idx[, 1]), ]
        if (nPeaks > nrow(idx)) {
            nPeaks = nrow(idx)
        }
        Sequence = character(nPeaks)
        for (i in 1:nPeaks) {
            Sequence[i] = as.character(genomicSequence[abs(idx[i, 
                2]):(abs(idx[i, 2]) + k - 1)])
        }
        peaks = data.frame(Affinity = idx[1:nPeaks, 1], Position = idx[1:nPeaks, 
            2], Sequence)
        cat(paste0(fit.string, " scores in ", seq.string, "\n"))
        print(peaks)
    }
    df = data.frame(Group = as.character(rep(1:2, each = ncol(score))), 
        Position = c(1:ncol(score), 1:ncol(score)), Affinity = c(score[1, 
            ], -score[2, ]))
    if (!is.null(genomicSequence2)) {
        df = rbind(df, data.frame(Group = as.character(rep(3:4, 
            each = ncol(score2))), Position = c(1:ncol(score2), 
            1:ncol(score2)), Affinity = c(score2[1, ], -score2[2, 
            ])))
        seq.string = paste0(seq.string, " and ", deparse(substitute(genomicSequence2)))
    }
    p = ggplot(df, aes(x = Position, y = Affinity, colour = Group)) + 
        theme_bw() + geom_line() + ylab("Relative Affinity") + 
        coord_fixed(ylim = c(-max.score, max.score)) + scale_y_continuous(breaks = divs, 
        labels = format(abs(divs), nsmall = 2)) + theme(aspect.ratio = 1, 
        text = element_text(size = 17, family = "Helvetica")) + 
        theme(legend.title = element_blank(), legend.position = c(0, 
            0), legend.justification = c(0, 0)) + theme(axis.title.x = element_text(vjust = 0), 
        axis.title.y = element_text(vjust = 0)) + labs(title = paste0(fit.string, 
        " scores in ", seq.string))
    if (is.null(genomicSequence2)) {
        p = p + scale_color_manual(values = c("#000000", "#FF0000"), 
            labels = c("Forward", "Reverse"))
    }
    else {
        p = p + scale_color_manual(values = c("#000000", "#FF0000", 
            "blue", "green"), labels = c(paste0(deparse(substitute(genomicSequence)), 
            " Forward"), paste0(deparse(substitute(genomicSequence)), 
            " Reverse"), paste0(deparse(substitute(genomicSequence2)), 
            " Forward"), paste0(deparse(substitute(genomicSequence2)), 
            " Reverse")))
    }
    if (annotate && !is.null(nPeaks)) {
        xpos = abs(peaks$Position) + ncol(score) * 0.015
        ypos = sign(peaks$Position) * peaks$Affinity + max.score * 
            0.015
        p = p + annotate("text", x = xpos, y = ypos, label = as.character(1:nPeaks))
    }
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
