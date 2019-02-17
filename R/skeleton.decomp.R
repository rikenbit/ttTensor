skeleton.decomp <-
function (A, r, thr = 1e-10, num.iter = 30) 
{
    .is.sparse(A)
    if (min(dim(A)) < r) {
        stop("Please specify the rank as min(dim(A)) > rank")
    }
    J = seq_len(r)
    iter <- 1
    RecError <- c()
    RelChange <- c()
    RecError[1] <- 1e+10
    RelChange[1] <- thr * 10
    while ((RelChange[iter] > thr) && (iter <= num.iter)) {
        A_bar <- A
        C <- A[, J]
        C <- as(C, "sparseMatrix")
        Q <- qr.Q(qr(C))
        Q <- as(Q, "sparseMatrix")
        I <- maxvol(Q)
        R <- A[I, ]
        R <- t(R)
        R <- as(R, "sparseMatrix")
        Q <- qr.Q(qr(R))
        Q <- as(Q, "sparseMatrix")
        J <- maxvol(Q)
        iter <- iter + 1
        Q_hat <- Q[J, ]
        A <- A[, J] %*% t(Q %*% solve(Q_hat))
        RecError[iter] <- sqrt(sum((A - A_bar)^2))
        RelChange[iter] <- abs(RecError[iter - 1] - RecError[iter])/RecError[iter]
    }
    U <- solve(A[I, J])
    U <- as(U, "sparseMatrix")
    R <- t(R)
    R <- as(R, "sparseMatrix")
    .is.sparse(C)
    .is.sparse(U)
    .is.sparse(R)
    list(C = C, U = U, R = R, rowidx = I, colidx = J, RecError = RecError, 
        RelChange = RelChange)
}
