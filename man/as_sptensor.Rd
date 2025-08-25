\name{as_sptensor}
\alias{as_sptensor}
\title{Convert to Simple Sparse Tensor}
\description{
  Converts an array or matrix to a simple sparse tensor format.
  This is a minimal implementation to replace the tensorr dependency.
}
\usage{
as_sptensor(x)
}
\arguments{
  \item{x}{An array or matrix to convert}
}
\value{
  A simple_sparse_tensor object
}
\details{
  This function provides a minimal sparse tensor implementation to support
  the TTCross function without requiring the archived tensorr package.
  For production use with actual sparse data, consider using specialized
  sparse tensor packages.
}
\examples{
# Create a 3D array
x <- array(rnorm(24), dim = c(2, 3, 4))

# Convert to sparse tensor format
sparse_x <- as_sptensor(x)
}
\keyword{methods}