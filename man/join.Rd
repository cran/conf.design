\name{join}
\alias{join}
\title{
Amalgamate two or more factors.
}
\description{
  Joins two or more factors together into a single composite factor defining
  the subclasses.  In a model formula \code{join(f1, f2, f3)} is equivalent to
  \code{f1:f2:f3}.
}
\usage{
join(\dots)
}
\arguments{
\item{...}{
  Two or more factors or numeric vectors, or lists containing
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
  Similar in effect to \code{paste}, which it uses.
}
\seealso{
  \code{\link{:}}, \code{\link{paste}}, \code{\link{rjoin}},
  \code{\link{direct.sum}}
}
\examples{
within(data.frame(f = gl(2, 3)), {
  g <- gl(3,2,length(f))
  fg <- join(f, g)
})
###   f  fg g
### 1 1 1:1 1
### 2 1 1:1 1
### 3 1 1:2 2
### 4 2 2:2 2
### 5 2 2:3 3
### 6 2 2:3 3
}
\keyword{design}

