\name{factorize}
\alias{factorize}
\alias{factorize.default}
\alias{factorize.factor}
\title{
Generic function.
}
\description{
The default method factorizes positive numeric integer arguments, returning
a vector of prime factors.  The factor method can be used to generate
pseudo-factors.  It accepts a factor, f, as principal argument and returns
a design with factors fa, fb, \dots{} each with a prime number of levels such
that f is model equivalent to join(fa, fb, \dots{}).
}
\usage{
factorize(x, \dots)
}
\arguments{
  \item{x}{Principal argument.  At this stage, it may be a numeric
    vector to elicit the default method, or a factor to elicit the
    factor method.}
  \item{\dots}{Additional arguments, if any.}
}
\value{
  A vector of (numeric) factors for numeric arguments, or a design with
  (S-PLUS) factors with prime numbers of levels for factor arguments.
}
\section{Side Effects}{
  None.
}
\details{
  Factorizes by a clumsy though effective enough way for small integers.
  May become very slow if some prime factors are large.  For the factor
  method it generates pseudo factors in the usual way.
}
\seealso{
conf.design, join
}
\examples{
 factorize(18)
 # [1] 2 3 3
 f <- factor(rep(0,5), rep(6,5))
 fd <- factorize(f)
}
\keyword{design}
% Converted by Sd2Rd version 1.10.
