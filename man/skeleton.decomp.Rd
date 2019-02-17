\name{skeleton.decomp}
\docType{methods}
\alias{skeleton.decomp}

\title{
Skeleton Decomposition
}
\description{
skeleton.decomp decomposes the input sparse matrix (n*m) and return the three matrices C (n*r), U (r*r), and R (r*m).
Only sparse matrix defined by the Matrix package is acceptable as the input.
}
\usage{
skeleton.decomp(A, r, thr=1E-10, num.iter=30)
}

\arguments{
  \item{A}{
The input sparse matrix.
}
  \item{r}{
Rank parameter to specify the lower dimension (r <= min(A)).
}
  \item{thr}{
The threshold to determine the convergence (Default: 1E-10).
}
  \item{num.iter}{
The number of iteration (Default: 30).
}
}

\value{
  C : A[I, :]
  U : inverse(A[I, J])
  R : A[:, J]
  rowidx :The indices of rows
  colidx : The indices of columns
  RecError : The reconstruction error between data matrix and reconstructed matrix from C, U, and R
  RelChange : The relative change of the error
}

\references{
I. V. Oseledets, et. al., (2010). TT-cross approximation for multidimensional arrays. \emph{Linear Algebra and its Applications}
}
\author{Koki Tsuyuzaki}

\seealso{
\code{\link{maxvol}}
}
\examples{
library("Matrix")
# Matrix data
X3 <- matrix(runif(10*20), nrow=10)
X3 <- as(X3, "sparseMatrix")
# Skeleton Decomposition
out.SKD <- skeleton.decomp(X3, r=3, num.iter=2, thr=1E-5)
}

\keyword{methods}