\name{conf.design}
\alias{conf.design}
\title{
Construct symmetric confounded factorial designs.
}
\description{
Construct designs with specified treatment contrasts confounded with blocks.
}
\usage{
conf.design(G, p, block.name="Blocks", treatment.names=<<see below>>)
}
\arguments{
\item{G}{
Matrix whose rows define the contrasts to be confounded.  For example in a
\code{3^4} experiment with \code{A B^2 C} and \code{B C D} confounded with blocks (together with
their generalized interactions), the matrix G would be rbind(c(1,2,1,0)),
c(0,1,1,1)).  The number of columns of G is the number of factors.
}
\item{p}{
The common number of levels for each factor.  Must be a prime number.
}
\item{block.name}{
Name to be given to the factor defining the blocks of the design.  
}
\item{treatment.names}{
Name to be given to the treatment factors of the design.  If G has a
dimnames attribute, then dimnames[[2]] is the default, otherwise T1, T2,
}}
\value{
A design with a Blocks factor defining the blocks and Treatment factors
defining the way treatments are allocated to each plot.  Not in random
order. 
}
\section{Side Effects}{
None.
}
\details{
For a single replicate of treatments, blocks are calculated using the
confounded contrasts in the standard textbook way.  The method is related
to that of Collings (1989).
}
\references{
Collings, B. J. (1989) Quick confounding.  Technometrics, v31, pp107-110.
}
\seealso{
conf.set, direct.sum, fac.design, fractionate
}
\examples{
# Generate a 3^4 factorial with A B^2 C and B C D confounded with blocks.  
 d34 <- conf.design(rbind(c(1,2,1,0), c(0,1,1,1)), p=3, treatment.names
= LETTERS[1:4])
}
\keyword{design}
% Converted by Sd2Rd version 1.10.
