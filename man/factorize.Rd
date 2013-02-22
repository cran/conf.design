\name{factorize}
\alias{factorize}
\alias{factorize.default}
\alias{factorize.factor}
\title{
  S3 generic function and methods for factorization.
}
\description{
 The \code{default} method factorizes positive numeric integer arguments,
  returning a vector of prime factors.  The \code{factor} method can be used to
  generate pseudo-factors.  It accepts a factor, \code{f}, as principal
  argument and returns a data frame with factors \code{fa}, \code{fb},
  \dots{} each with a prime number of levels such that \code{f} is model
  equivalent to \code{join(fa, fb, \dots)}.
}

\usage{
\method{factorize}{default}(x, divisors = primes(max(x)), \dots)
\method{factorize}{factor}(x, name = deparse(substitute(x)), extension =
letters, \dots)
}
\arguments{
  \item{x}{
    Principal argument.

    The \code{default} method factorizes (smallish) positive integers;

    The \code{factor} method generates prime pseudo-factors from a
    factor with a composite number of levels (as required for partial
    confounding).

  }
  \item{divisors}{
    Candidate prime divisors for all numers to be factorized.
  }
  \item{name}{
    Stem of the name to be given to component pseudo-factors.
  }
  \item{extension}{
    Distinguishing strings to be added to the stem to nominate the
    pseudo-factors.
  }
\item{\dots}{
  Additional arguments, if any.  (Presently ignored.)}
}
\value{
  For the \code{default} method a vector, or list of vectors, of prime
  integer divisors of the components of \code{x}, (including repeats).

  For the \code{factor} method, a design with factors having prime
  numbers of levels for factor arguments.

}
\section{Side Effects}{
  None.
}
\details{
  Primarily intended to split a factor with a non-prime number of levels
  into a series of pseudo-factors, each with a prime number of levels
  and which jointly define the same classes as the factor itself.

  The main reason to do this would be to confound one or more of the
  pseudo-factors, or their interactions, with blocks using constructions
  that only apply for prime numbers of levels.  In this way the
  experiment can be made smaller, at the cost of some treatment
  contrasts being confounded with blocks.

  The default method factorizes integers by a clumsy, though effective
  enough way for small integers.  The function is vectorized in the
  sense that if a vector of integers to factorize is specified, the
  object returned is a list of numeric vectors, giving the prime
  divisors (including repeats) of the given integers respectively.

  As with any method of factorizing integers, it may become very slow if
  the prime factors are large.
}
\seealso{
conf.design, join
}
\examples{
factorize(12321)
### [1]  3  3 37 37

f <- factor(1:6)
data.frame(f, factorize(f))
###   f fa fb
### 1 1  0  0
### 2 2  1  0
### 3 3  0  1
### 4 4  1  1
### 5 5  0  2
### 6 6  1  2

des <- with(list(f = factor(rep(6:1, 1:6))),
            data.frame(f, factorize(f)))
head(des, 7)
##   f fa fb
## 1 6  1  2
## 2 5  0  2
## 3 5  0  2
## 4 4  1  1
## 5 4  1  1
## 6 4  1  1
## 7 3  0  1
}
\keyword{design}

