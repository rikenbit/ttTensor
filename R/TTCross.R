TTCross <-
function (A, Ranks = NULL, thr = 1e-10, num.iter = 30) 
{
    .is.sparse(A)
    if (is.null(Ranks)) {
        stop("Please specify the Ranks")
    }
    D <- dim(A)
    nModes <- length(D)
    G <- rep(list(NULL), nModes)
    C <- A
    for (k in 1:(nModes - 1)) {
        if (k == 1) {
            rk1 <- 1
            C <- tensorr::unfold(C, 1)
            C <- C@mat
        }
        else {
            rk1 <- Ranks[k - 1]
            dim(C) <- c(rk1 * D[k], prod(dim(C))/(rk1 * D[k]))
        }
        rk <- Ranks[k]
        res.skd <- skeleton.decomp(C, rk, thr, num.iter)
        G[[k]] <- as.matrix(res.skd$C)
        if (k == 1) {
            dim(G[[1]]) <- c(D[1], prod(dim(G[[1]]))/D[1])
        }
        else {
            dim(G[[k]]) <- c(rk1, D[k], rk)
        }
        C <- res.skd$U[seq_len(rk), seq_len(rk)] %*% res.skd$R[seq_len(rk), 
            ]
    }
    G[[nModes]] <- as.matrix(C)
    .setDimNames(G, A, Ranks)
}
