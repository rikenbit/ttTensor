\name{unfold}
\alias{unfold}
\title{Unfold a Tensor}
\description{
  Unfolds a tensor along a specified mode into a matrix representation.
}
\usage{
unfold(x, mode)
}
\arguments{
  \item{x}{A tensor object (simple_sparse_tensor, Tensor, array, or matrix)}
  \item{mode}{The mode along which to unfold the tensor}
}
\value{
  A matrix representation of the unfolded tensor
}
\details{
  This function unfolds a tensor along the specified mode into a matrix.
  It supports simple_sparse_tensor objects, rTensor Tensor objects, and
  regular arrays/matrices. The function uses rTensor's rs_unfold internally.
}
\examples{
library(rTensor)

# Create a 3D tensor
x <- array(rnorm(24), dim = c(2, 3, 4))
tensor_x <- as.tensor(x)

# Unfold along mode 1 (using ttTensor's unfold function)
unfolded <- ttTensor::unfold(tensor_x, mode = 1)
}
\keyword{methods}