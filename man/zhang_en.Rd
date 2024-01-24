\name{zhang_en}
\alias{zhang_en}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
zhang_en epigenetic clock
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
zhang_en(betas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betas}{
%%     ~~Describe \code{betas} here~~
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
function (betas)
{
    Zhang_en_CpGs <- read.table("data/Zhang_en_CpGs.coef", stringsAsFactor = F,
        header = T)
    intercept <- Zhang_en_CpGs[1, 2]
    betas.norm <- apply(t(betas), 1, scale)
    rownames(betas.norm) <- rownames(betas)
    coefs <- setNames(Zhang_en_CpGs$CoefficientTraining, Zhang_en_CpGs$CpGmarker)
    CpGs <- intersect(names(coefs), rownames(betas))
    tt <- betas.norm[CpGs, ] * coefs[CpGs]
    colSums(tt, na.rm = TRUE) + intercept
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
