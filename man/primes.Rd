\name{primes}
\alias{primes}
\title{
Prime numbers
}
\description{
Generate a table of prime numbers.
}
\usage{
primes(n)
}
\arguments{
\item{n}{
Positive integer value.
}}
\value{
A vector of all prime numbers less than the argument.  1 is not a prime.
}
\section{Side Effects}{
None
}
\details{
Uses an elementary sieve method, not suitable for very large n.  This
function is no longer used as part of the design library and is included
for curiosity value only.
}
\seealso{
factorize
}
\examples{
 primes(50)
 # [1]  2  3  5  7 11 13 17 19 23 29 31 37 41 43 47
}
\keyword{design}
% Converted by Sd2Rd version 1.10.
