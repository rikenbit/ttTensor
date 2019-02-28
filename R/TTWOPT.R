TTWOPT <-
function (X, Ranks, W = NULL, eta = 1e-07, thr = 1e-10, num.iter = 100) 
{
    checkDense <- "Tensor" %in% is(X)
    if (!checkDense) {
        stop("Please specify X as Tensor (cf. rTensor)")
    }
    if (is.null(W)) {
        W <- X
        W@data[] <- 1
    }
    else {
        checkW <- identical(dim(X), dim(W))
        if (!checkW) {
            stop("The size of W must be same as that of X")
        }
    }
    nModes <- length(dim(X))
    Y <- X * W
    G <- .genCores(X, Ranks)
    Z <- W * .recTensor(G)
    f = c()
    RecError = c()
    RelChange = c()
    f[1] <- 0.5 * .tnorm(Y) - sum(Y@data * Z@data) + 0.5 * .tnorm(Z)
    X_bar <- .recTensor(G)
    RecError[1] <- .recError(X@data, X_bar)
    RelChange[1] <- thr * 1e+10
    iter <- 1
    while ((RelChange[iter] > thr) && (iter < num.iter)) {
        Z <- W * .recTensor(G)
        for (n in seq_len(nModes)) {
            Zn <- cs_unfold(Z, n)
            Yn <- cs_unfold(Y, n)
            if (n == 1) {
                Glgn <- .recTensor(G[2:nModes])
                Glen <- 1
            }
            else if (n == nModes) {
                Glgn <- 1
                Glen <- .recTensor(G[1:(nModes - 1)])
            }
            else {
                Glgn <- .recTensor(G[(n + 1):nModes])
                Glen <- .recTensor(G[1:(n - 1)])
            }
            if (length(dim(Glgn)) > 2) {
                uGlgn <- cs_unfold(as.tensor(Glgn), 1)@data
            }
            else {
                uGlgn <- Glgn
            }
            if (length(dim(Glen)) > 2) {
                uGlen <- cs_unfold(as.tensor(Glen), length(dim(Glen)))@data
            }
            else {
                uGlen <- Glen
            }
            if (n >= nModes - 1) {
                grad <- t(Zn@data - Yn@data) %*% t(kronecker(uGlgn, 
                  t(uGlen)))
            }
            else {
                grad <- t(Zn@data - Yn@data) %*% kronecker(uGlgn, 
                  uGlen)
            }
            if (n == 1) {
                grad <- t(grad)
            }
            dimGn <- dim(G[[n]])
            tmpGn <- rs_unfold(as.tensor(G[[n]]), 2)@data
            tmpGn = tmpGn - eta * grad
            G[[n]] <- fold(tmpGn, row_idx = 2, col_idx = setdiff(seq_len(length(dim(G[[n]]))), 
                2), modes = dimGn)@data
        }
        iter <- iter + 1
        if (iter > 1) {
            f[iter] <- 0.5 * .tnorm(Y) - sum(Y@data * Z@data) + 
                0.5 * .tnorm(Z)
            X_bar <- .recTensor(G)
            RecError[iter] <- .recError(X@data, X_bar)
            RelChange[iter] <- abs(f[iter - 1] - f[iter])/f[iter]
        }
    }
    list(G = G, RelChange = RelChange, f = f, RecError = RecError)
}
