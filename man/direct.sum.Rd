\name{direct.sum}
\alias{direct.sum}
\title{
Form the direct sum of designs.
}
\description{
  Constructs the direct sum of two or more designs.  Each plot of one design
  is matched with every plot of the other.  (This might also be called the
  Cartesian product of two designs).
}
\usage{
direct.sum(D1, \dots, tiebreak=letters)
}
\arguments{
  \item{D1}{
    First component design.
  }
  \item{...}{
    Additional component designs, if any.
  }
  \item{tiebreak}{
    Series of characters or digits to be used for breaking ties (or repeats) in
    the variable names in the component designs.  Augmented if necessary.
}}
\value{
  The direct sum of all component designs.
}
\section{Side Effects}{
  None.
}
\details{
  Each plot of one design is matched with every plot of the next, (if
  any), and so on recursively.
}
\seealso{
  \code{\link{conf.design}}
}
\examples{
### Generate a half replicate of a 2^3 x 3^2 experiment.  The factors are
### to be A, B, C, D, E.  The fractional relation is to be I = ABC and the
### DE effect is to be confounded with blocks.

### First construct the 2^3 design, confounded in two blocks:
d1 <- conf.design(cbind(A = 1, B = 1, C = 1), p = 2)

### Next the 3^2 design, with DE partially confounded in blocks:
d2 <- conf.design(cbind(D = 1, E = 1), p = 3)

### Now extract the principal block from the 2^3 design and form the direct
### sum withthe 3^2 design
dsn <- direct.sum(subset(d1, Blocks == "0"), d2)

### combine block factors into one
dsn <- within(dsn, {
  Blocks <- join(Blocks, Blocksa)
  Blocksa <- NULL
})
head(dsn)
}
\keyword{design}

