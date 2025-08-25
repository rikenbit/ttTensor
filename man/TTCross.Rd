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
\dontrun{
library("rTensor")
# Sparse Tensor data
X1 <- array(rnorm(3*5*7), c(3,5,7))
X1 <- as.tensor(X1)
X2 <- as_sptensor(dtensor(X1@data))
# TT-ranks
Ranks <- c(p=2, q=4)
# TT-Cross
out.TTCross <- TTCross(X2, Ranks, num.iter=2)
}
}

\keyword{methods}