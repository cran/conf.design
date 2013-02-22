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
  A positive integer value, or vector or such values.
}}
\value{
  A vector of all prime numbers less than the \code{max(n)}.

  NB: 1 is not a prime number!
}
\section{Side Effects}{
  None
}
\details{
  Uses an elementary sieve method.
}
\seealso{
\code{\link{factorize}}
}
\examples{
primes(1:50)
### [1]  2  3  5  7 11 13 17 19 23 29 31 37 41 43 47
}
\keyword{design}

