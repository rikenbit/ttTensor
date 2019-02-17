TTSVD <-
function (A, Ranks = NULL, accuracy = NULL) 
{
    checkDense <- "Tensor" %in% is(A)
    if (!checkDense) {
        stop("Please specify A as Tensor (cf. rTensor)")
    }
    if (is.null(Ranks) && is.null(accuracy)) {
        stop("Please specify the Ranks or accuracy parameter")
    }
    if (!is.null(accuracy)) {
        l <- length(dim(A)) - 1
        Ranks <- rep(NA, length = l)
        names(Ranks) <- letters[16:26][seq_len(l)]
    }
    A <- A@data
    D <- dim(A)
    nModes <- length(D)
    G <- rep(list(NULL), nModes)
    if (!is.null(accuracy)) {
        delta <- accuracy/sqrt(nModes) * sqrt(sum(A^2))
    }
    C <- A
    for (k in 1:(nModes - 1)) {
        if (k == 1) {
            rk1 <- 1
        }
        else {
            rk1 <- Ranks[k - 1]
        }
        dim(C) <- c(rk1 * D[k], prod(dim(C))/(rk1 * D[k]))
        res.svd <- svd(C)
        if (is.null(accuracy)) {
            rk <- Ranks[k]
        }
        else {
            rk <- .rankDelta(C, res.svd$u, res.svd$d, t(res.svd$v), 
                delta)
            Ranks[k] <- rk
        }
        G[[k]] <- res.svd$u[, 1:rk]
        if (k == 1) {
            dim(G[[1]]) <- c(D[1], prod(dim(G[[1]]))/D[1])
        }
        else {
            dim(G[[k]]) <- c(rk1, D[k], rk)
        }
        Sigma <- res.svd$d[seq_len(rk)]
        if (length(Sigma) >= 2) {
            Sigma <- diag(res.svd$d[seq_len(rk)])
        }
        V <- as.matrix(res.svd$v[, seq_len(rk)])
        C <- Sigma %*% t(V)
    }
    G[[nModes]] <- C
    .setDimNames(G, A, Ranks)
}
