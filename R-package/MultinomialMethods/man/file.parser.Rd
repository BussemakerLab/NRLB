\name{file.parser}
\alias{file.parser}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parse fit output
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
file.parser(fileName, fitInfoOnly = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileName}{
%%     ~~Describe \code{fileName} here~~
}
  \item{fitInfoOnly}{
%%     ~~Describe \code{fitInfoOnly} here~~
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
function (fileName, fitInfoOnly = FALSE) 
{
    nCols = max(count.fields(fileName, sep = ","))
    header = read.csv(fileName, header = FALSE, nrows = 1, stringsAsFactors = FALSE, 
        strip.white = TRUE)
    colNames = as.character(as.vector(header[1, ]))
    colNames = c(colNames, (length(colNames) + 1):nCols)
    output = read.csv(fileName, header = FALSE, fill = TRUE, 
        col.names = colNames, stringsAsFactors = FALSE)
    output = output[2:nrow(output), 2:ncol(output)]
    row.names(output) = as.character(1:nrow(output))
    for (i in 7:12) {
        output[, i] = as.logical(output[, i])
    }
    for (i in c(2:6, 15:23)) {
        output[, i] = as.numeric(output[, i])
    }
    output[output == ""] <- NA
    output$k[is.na(output$k)] = "Multi"
    if (fitInfoOnly) {
        output = output[, 1:24]
        return(output)
    }
    else {
        info = output[, 1:24]
        values = vector(mode = "list", length = nrow(output))
        for (i in 1:nrow(output)) {
            nb = db = sb = nsb = ne = de = se = nsbe = ns = ds = ss = nsbs = NULL
            if (info$k[i] == "Multi") {
                fit = output[i, ]
                modeStartIdx = grep("NB>", fit)
                nModes = length(modeStartIdx)
                modeStartIdx = c(modeStartIdx, grep("<EOL>", 
                  fit))
                if (info$NS[i]) {
                  modes = vector(mode = "list", length = (nModes + 
                    1))
                  names(modes) = c(paste0("Mode", 1:nModes), 
                    "NSBinding")
                  nsb = as.numeric(fit[grep("NSB>", fit) + 1])
                  nsbe = as.numeric(fit[grep("NSBE>", fit) + 
                    1])
                  nsbs = as.numeric(fit[grep("NSBS>", fit) + 
                    1])
                  modes$NSBinding = list(NSB = nsb, NSBE = nsbe, 
                    NSBS = nsbs)
                }
                else {
                  modes = vector(mode = "list", length = nModes)
                  names(modes) = paste0("Mode", 1:nModes)
                }
                for (currMode in 1:nModes) {
                  nb = db = sb = ne = de = se = ns = ds = ss = NULL
                  currModeValues = c(fit[modeStartIdx[currMode]:(modeStartIdx[currMode + 
                    1] - 1)], ">")
                  delimiters = grep(">", currModeValues)
                  for (delimIdx in 1:(length(delimiters) - 1)) {
                    currValue = as.numeric(currModeValues[(delimiters[delimIdx] + 
                      1):(delimiters[delimIdx + 1] - 1)])
                    currType = currModeValues[delimiters[delimIdx]]
                    if (currType == "NB>") {
                      nb = currValue
                    }
                    else if (currType == "DB>") {
                      db = currValue
                    }
                    else if (currType == "SB>") {
                      sb = currValue
                    }
                    else if (currType == "NE>") {
                      ne = currValue
                    }
                    else if (currType == "DE>") {
                      de = currValue
                    }
                    else if (currType == "SE>") {
                      se = currValue
                    }
                    else if (currType == "NS>") {
                      ns = currValue
                    }
                    else if (currType == "DS>") {
                      ds = currValue
                    }
                    else if (currType == "SS>") {
                      ss = currValue
                    }
                  }
                  modes[[currMode]] = list(NB = nb, DB = db, 
                    SB = sb, NE = ne, DE = de, SE = se, NS = ns, 
                    DS = ds, SS = ss)
                }
                values[[i]] = modes
            }
            else {
                fit = output[i, ]
                delimiters = grep(">", fit)
                for (delimIdx in 1:(length(delimiters) - 1)) {
                  currValue = as.numeric(fit[(delimiters[delimIdx] + 
                    1):(delimiters[delimIdx + 1] - 1)])
                  currType = fit[delimiters[delimIdx]]
                  if (currType == "NB>") {
                    nb = currValue
                  }
                  else if (currType == "DB>") {
                    db = currValue
                  }
                  else if (currType == "SB>") {
                    sb = currValue
                  }
                  else if (currType == "NSB>") {
                    nsb = currValue
                  }
                  else if (currType == "NE>") {
                    ne = currValue
                  }
                  else if (currType == "DE>") {
                    de = currValue
                  }
                  else if (currType == "SE>") {
                    se = currValue
                  }
                  else if (currType == "NSBE>") {
                    nsbe = currValue
                  }
                  else if (currType == "NS>") {
                    ns = currValue
                  }
                  else if (currType == "DS>") {
                    ds = currValue
                  }
                  else if (currType == "SS>") {
                    ss = currValue
                  }
                  else if (currType == "NSBS>") {
                    nsbs = currValue
                  }
                }
                values[[i]] = list(NB = nb, DB = db, SB = sb, 
                  NSB = nsb, NE = ne, DE = de, SE = se, NSBE = nsbe, 
                  NS = ns, DS = ds, SS = ss, NSBS = nsbs)
            }
        }
        return(list(Information = info, Values = values))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
