\name{join}
\alias{join}
\title{
Amalgamate two or more factors.
}
\description{
Joins two or more factors together into a single composite factor defining
the subclasses.  In a model formula join(f1, f2, f3) is equivalent to
f1:f2:f3.
}
\usage{
join(\dots)
}
\arguments{
\item{...}{
Two or more factors or numeric vectors, or objects of mode list containing
these kinds of component.
}}
\value{
A single composite factor with levels made up of the distinct combinations
of levels or values of the arguments which occur.
}
\section{Side Effects}{
None.
}
\details{
Similar in effect to paste(), which it uses.
}
\seealso{
paste, rjoin, direct.sum
}
\examples{
 d1 <- conf.design(c(1,1,1), 2, treatment.names=LETTERS[1:3])
 d2 <- conf.design(c(0,1,1), 2, treatment.names=LETTERS[1:3])
 fd1d2 <- join(d1,d2)
}
\keyword{design}
% Converted by Sd2Rd version 1.10.
