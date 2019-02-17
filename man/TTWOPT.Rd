\name{TTWOPT}
\docType{methods}
\alias{TTWOPT}

\title{
Tensor-Train Decomposition by Tensor-train Weighted OPTimization
}
\description{
TTWOPT incrementaly decomposes the input tensor by gradient desecent. The tensor with missing entries is also specified with weight tensor W.
}
\usage{
TTWOPT(X, Ranks, W=NULL, eta=1E-7, thr=1E-10, num.iter=100)
}

\arguments{
  \item{X}{
The input tensor.
}
  \item{Ranks}{
TT-ranks to specify the lower dimensions.
}
  \item{W}{
The weight tensor to specify the missing entries (0: missing, 1: existing). The size must be same as that of X.
}
  \item{eta}{
The learning rate parameter of the gradient descent algorithm (Default : 1E-10).
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
  RelChange : The relative change of the error
  f : The values of the object function
  RecError : The reconstruction error between data tensor and reconstructed tensor from C, U, and R
}

\references{
Yuan, Longhao, et. al., (2017). Completion of high order tensor data with missing entries via tensor-train decomposition. \emph{International Conference on Neural Information Processing}
}
\author{Koki Tsuyuzaki}

\examples{
library("rTensor")
# Tensor data
X1 <- array(rnorm(3*5*7*9*11), c(3,5,7,9,11))
dimnames(X1) <- list(
    I=paste0("i", 1:3),
    J=paste0("j", 1:5),
    K=paste0("k", 1:7),
    L=paste0("l", 1:9),
    M=paste0("m", 1:11)
    )
X1 <- as.tensor(X1)
# TT-ranks
Ranks <- c(p=2, q=4, r=6, s=8)
# TTWOPT
out.TTWOPT <- TTWOPT(X1, Ranks, eta=1E-7)
}

\keyword{methods}