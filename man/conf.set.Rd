\name{conf.set}
\alias{conf.set}
\title{
Find confounded effects.
}
\description{
Find minimal complete sets of confounded effects from a defining set for
symmetric confounded factorial designs. Useful for checking if a low order
interaction will be unintentionally confounded with blocks. As in the
usual convention, only effects whose leading factor has an index of one are
listed.
}
\usage{
conf.set(G, p)
}
\arguments{
\item{G}{
Matrix whose rows define the effects to be confounded with blocks, in the
same way as for conf.design().
}
\item{p}{
Number of levels for each factor.  Must be a prime number.
}}
\value{
A matrix like G with a minimal set of confounded with blocks defined in the
rows.
}
\section{Side Effects}{
None
}
\details{
The function constructs all linear functions of the rows of G (over GF(p)),
and removes those rows whose leading non-zero component is not one.
}
\seealso{
conf.design
}
\examples{
 G <- rbind(c(1,2,1,0), c(0,1,1,1))
 dimnames(G) <- list(NULL, LETTERS[1:4])
 conf.set(G, 3)
 #     A B C D 
 # [1,] 1 2 1 0
 # [2,] 0 1 1 1
 # [3,] 1 0 2 1
 # [4,] 1 1 0 2
 # If A B^2 C and B C D are confounded with blocks, then so are A C^2 D
 # and A B D^2.  Only three-factor interactions are confounded, so the
 # design is presumably useful.
}
\keyword{design}
% Converted by Sd2Rd version 1.10.
