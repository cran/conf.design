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
manner to rbind()), together with an additional factor whose levels
identify the original component designs, or \code{parts}.
}
\section{Side Effects}{
None.
}
\details{
Almost the same as rbind(), but an additional factor in the result
separates the original components.
}
\seealso{
conf.design, join, direct.sum
}
\examples{
 # A two replicate partially confounded factorial design.
 d1 <- conf.design(c(1,1,1), 2, treatment.names=LETTERS[1:3])
 d2 <- conf.design(c(0,1,1), 2, treatment.names=LETTERS[1:3])
 dsn <- rjoin(d1, d2)
}
\keyword{design}
% Converted by Sd2Rd version 1.10.

