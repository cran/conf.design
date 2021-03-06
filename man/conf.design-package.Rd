\name{conf.design-package}
\alias{conf.design-package}
% \alias{conf.design}
\docType{package}
\title{
Confounded Factorial Block Design
}
\description{
Construct confounded designs with specific contrasts confounded with
blocks.  The package only directly handles the \code{p^k} case, that is,
all treatment factors having the same (prime) number of levels.  Some
simple facilities are provided for combining component designs into
larger ones, thus providing some facilities for generating interesting
designs for the more general case.  Some fractional replication is also
possible with the tools proveded.
}
\details{
\tabular{ll}{
Package: \tab conf.design\cr
Type: \tab Package\cr
Version: \tab 2.0.0\cr
Date: \tab 2013-02-23\cr
License: \tab GPL-2\cr
}
The key functions are \code{conf.design}, \code{rjoin} and
\code{direct.sum}.  The help information for these functions contain
some fairly detailed examples.

Some ancillary functions may be of independent interest, for example
\code{primes} for generating prime numbers and \code{factorize}, which
has an experimental design application, but a default method can be used
to factorize (not too large) integers into prime factors.
}
\author{
Bill Venables

Maintainer: Bill Venables <Bill.Venables@gmail.com>
}
\references{
  Collings, B. J. (1989) Quick confounding.  Technometrics, v31, pp107-110.
}
\keyword{ design }
\seealso{
The CRAN task view on Design of Experiments
}
\examples{

### Generate a half replicate of a 2^3 x 3^2 experiment.  The factors
### are to be A, B, C, D, E.  The fractional relation is to be I = ABC
### and the DE effect is to be confounded with blocks.

### First construct the 2^3 design, confounded in two blocks:
d1 <- conf.design(c(A = 1, B = 1, C = 1), p=2)

### Next the 3^2 design, with DE confounded in blocks:
d2 <- conf.design(c(D = 1, E = 1), p=3)

### Now extract the principal block from the 2^3 design and form the direct
### sum withthe 3^2 design
dsn <- direct.sum(subset(d1, Blocks == "0"), d2)
head(dsn)
###    Blocks A B C Blocksa D E
###  1      0 0 0 0       0 0 0
###  2      0 0 0 0       0 2 1
###  3      0 0 0 0       0 1 2
###  4      0 0 0 0       1 1 0
###  5      0 0 0 0       1 0 1
###  6      0 0 0 0       1 2 2
###
### Combine the two "Blocks" factors into a single block factor:
dsn <- within(dsn, {
  Blocks <- join(Blocks, Blocksa)
  Blocksa <- NULL
})
### Now to do some checks.
as.matrix(replications( ~ .^2, dsn))

### Blocks     12
### A          18
### B          18
### C          18
### D          12
### E          12
### Blocks:A    6
### Blocks:B    6
### Blocks:C    6
### Blocks:D    4
### Blocks:E    4
### A:B         9
### A:C         9
### A:D         6
### A:E         6
### B:C         9
### B:D         6
### B:E         6
### C:D         6
### C:E         6
### D:E         4

### We can check the confounding by analysing some dummy data:

dsn$y <- rnorm(nrow(dsn))
dummyAov <- aov(y ~ A*B*C*D*E + Error(Blocks), data=dsn)
summary(dummyAov)

###  Error: Blocks
###      Df Sum Sq Mean Sq
###  D:E  2  8.915   4.458
###
###  Error: Within
###        Df Sum Sq Mean Sq
###  A      1  2.077   2.077
###  B      1  1.111   1.111
###  C      1  3.311   3.311
###  D      2  1.929   0.964
###  E      2  0.848   0.424
###  A:D    2  3.421   1.711
###  B:D    2  3.231   1.615
###  C:D    2  2.484   1.242
###  A:E    2  0.214   0.107
###  B:E    2  0.006   0.003
###  C:E    2  0.349   0.174
###  D:E    2  1.442   0.721
###  A:D:E  4  2.560   0.640
###  B:D:E  4  4.454   1.114
###  C:D:E  4  7.942   1.986

### Two of the D:E degrees of freedom are confounded with Blocks, as desired.
}
