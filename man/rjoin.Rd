\name{rjoin}
\alias{rjoin}
\title{
  Concatenate designs by rows.
}
\description{
  Combine two or more designs with the same names into a single design by row
  concatenation.
}
\usage{
rjoin(\dots, part.name="Part")
}
\arguments{
\item{...}{
  Two or more designs with identical component names.
}
\item{part.name}{
  Name for an additional factor to identify the original components in the result.
}}
\value{
  A single design with the arguments stacked above each other (in a similar
  manner to \code{rbind}), together with an additional factor whose levels
  identify the original component designs, or \code{parts}.
}
\section{Side Effects}{
  None.
}
\details{
Almost the same as \code{rbind}, but an additional factor in the result
separates the original components.
}
\seealso{
\code{\link{rbind}}
}
\examples{
### A two replicate partially confounded factorial design.
d1 <- conf.design(c(A = 1, B = 1, C = 1), 2)
d2 <- conf.design(c(A = 0, B = 1, C = 1), 2)
dsn <- within(rjoin(d1, d2), {
   Blocks <- join(Part, Blocks)
   Part <- NULL
})
as.matrix(replications(~ .^2, dsn))
###          [,1]
### Blocks      4
### A           8
### B           8
### C           8
### Blocks:A    2
### Blocks:B    2
### Blocks:C    2
### A:B         4
### A:C         4
### B:C         4
}
\keyword{design}


