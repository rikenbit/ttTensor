\name{dtensor}
\alias{dtensor}
\title{Dense Tensor Creation}
\description{
  Creates a dense tensor representation. This is a compatibility function
  that simply returns the input as-is.
}
\usage{
dtensor(x)
}
\arguments{
  \item{x}{An array or matrix}
}
\value{
  The input array or matrix unchanged
}
\details{
  This function is provided for compatibility with code that previously
  used the tensorr package. It simply returns the input without modification.
}
\examples{
# Create a 3D array
x <- array(rnorm(24), dim = c(2, 3, 4))

# Create dense tensor (returns x unchanged)
dense_x <- dtensor(x)
}
\keyword{methods}