\name{conf.design}
\alias{conf.design}
\title{
Construct symmetric confounded factorial designs.
}
\description{
Construct designs with specified treatment contrasts confounded with
blocks.  All treatment factors must have the sampe (prime) number of
levels.
}
\usage{
conf.design(G, p, block.name = "Blocks",  treatment.names = NULL)
}
\arguments{
\item{G}{

  Matrix whose rows define the contrasts to be confounded.

  The number of columns of \code{G} defines the number of
  factors.

}
\item{p}{
  The common number of levels for each factor.  Must be a prime number.
}
\item{block.name}{
  Name to be given to the factor defining the blocks of the design.
}
\item{treatment.names}{
  Name to be given to the treatment factors of the design.  If
  \code{NULL} and if \code{G} has a \code{dimnames} attribute, then
  \code{dimnames[[2]]} is the default, otherwise T1, T2, \dots
}}
\value{
  A design with a \code{Blocks} factor defining the blocks and treatment
  factors defining the way treatments are allocated to each plot.  (Not
  in randomised order!)
}
\section{Side Effects}{
  None.
}
\details{
  For example in a \code{3^4} experiment with \code{AB^2C} and
  \code{BCD} confounded with blocks (together with their generalized
  interactions), the matrix \code{G} could be given by

  \code{rbind(c(A = 1, B = 2, C = 1, D = 0), c(A = 0, B = 1, C = 1, D =
  1))}

  For this example, \code{p = 3}

  Having column names for the \code{G} matrix implicitly supplies the
  treatment factor names.

  For a single replicate of treatments, blocks are calculated using the
  confounded contrasts in the standard textbook way.  The method is
  related to that of Collings (1989).
}
\references{
  Collings, B. J. (1989) Quick confounding.  Technometrics, v31, pp107-110.
}
\seealso{
  \code{\link{conf.set}}, \code{\link{direct.sum}}
}
\examples{
###
### Generate a 3^4 factorial with A B^2 C and B C D confounded with blocks.
###
d34 <- conf.design(rbind(c(A = 1, B = 2, C = 1, D = 0),
                         c(A = 0, B = 1, C = 1, D = 1)), p = 3)
head(d34)
###   Blocks A B C D
### 1     00 0 0 0 0
### 2     00 1 2 1 0
### 3     00 2 1 2 0
### 4     00 2 2 0 1
### 5     00 0 1 1 1
### 6     00 1 0 2 1

as.matrix(replications(~ .^2, d34))
###           [,1]
###  Blocks      9
###  A          27
###  B          27
###  C          27
###  D          27
###  Blocks:A    3
###  Blocks:B    3
###  Blocks:C    3
###  Blocks:D    3
###  A:B         9
###  A:C         9
###  A:D         9
###  B:C         9
###  B:D         9
###  C:D         9
}
\keyword{design}

