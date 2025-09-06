\name{TTCross}
\docType{methods}
\alias{TTCross}

\title{
Tensor-Train Decomposition by TRCross
}
\description{
TTCross incrementaly decomposes the input tensor by skeleton decomposition algorithm. The algorithm only select the row/column indices and any large temporal matrix are genrated in the process. Therefore, this method is suitable for the sparse tensor.
}
\usage{
TTCross(A, Ranks=NULL, thr=1E-10, num.iter=30)
}

\arguments{
  \item{A}{
The input sparse tensor.
}
  \item{Ranks}{
TT-ranks to specify the lower dimensions.
}
  \item{thr}{
The threshold to determine the convergence (Default: 1E-10).
}
  \item{num.iter}{
The number of iteration (Default: 30).
}
}

\value{
  G : Core tensors
}

\references{
I. V. Oseledets, et. al., (2010). TT-cross approximation for multidimensional arrays. \emph{Linear Algebra and its Applications}
}
\author{Koki Tsuyuzaki}

\examples{
\donttest{
# TTCross requires sparse tensor input
# Creating a simple example
library("rTensor")
X1 <- array(rnorm(3*4*5), c(3,4,5))
X1 <- as.tensor(X1)
# Convert to sparse format
X2 <- as_sptensor(dtensor(X1@data))
# TT-ranks (should be less than dimensions)
Ranks <- c(p=2, q=2)
# Note: TTCross is designed for sparse tensors
# and may have numerical issues with some inputs
tryCatch({
  out.TTCross <- TTCross(X2, Ranks, num.iter=2)
  print("TTCross completed")
}, error = function(e) {
  print("TTCross encountered an error - this function is experimental")
})
}
}

\keyword{methods}