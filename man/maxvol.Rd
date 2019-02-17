\name{maxvol}
\docType{methods}
\alias{maxvol}

\title{
maxvol algorithm
}
\description{
maxvol finds the r*r submatrix of maximal volume in C (n*r) by greedily searching the vector of max norm, and subtractction of its projection from the rest of rows.
See also http://tensorly.org/stable/_modules/tensorly/contrib/decomposition/mps_decomposition_cross.html#matrix_product_state_cross
}
\usage{
maxvol(C)
}

\arguments{
  \item{C}{
The input sparse matrix.
}
}

\value{
  row_idx : The indices of rows, which make the determinant as large
}

\references{
Ali Civril, et. al., (2009). On selecting a maximum volume sub-matrix of a matrix and related problems. \emph{Theoretical Computer Science}
}
\author{Koki Tsuyuzaki}

\seealso{
\code{\link{skeleton.decomp}}
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