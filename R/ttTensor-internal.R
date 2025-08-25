.rankDelta <-
function (C, left, mid, right, delta) 
{
    error <- lapply(seq_along(mid), function(x) {
        U <- as.matrix(left[, seq_len(x)])
        Sigma <- mid[seq_len(x)]
        V <- as.matrix(right[seq_len(x), ])
        if (length(Sigma) >= 2) {
            Sigma <- diag(mid[seq_len(x)])
        }
        else {
            V <- t(right[seq_len(x), ])
        }
        sqrt(sum((C - U %*% Sigma %*% V)^2))
    })
    min(which(error <= delta))
}
.dimList <-
function (R, len) 
{
    eval(parse(text = paste0("list(", R, " = paste0(tolower(R), seq_len(len)))")))
}
.setDimNames <-
function (G, A, Ranks) 
{
    names(G)[1] <- paste0(names(dimnames(A)[1]), "-", names(Ranks[1]))
    dimnames(G[[1]]) <- c(dimnames(A)[1], .dimList(names(Ranks[1]), 
        Ranks[1]))
    nModes <- length(G)
    for (k in 2:(nModes - 1)) {
        names(G)[k] <- paste0(names(Ranks[k - 1]), "-", names(dimnames(A)[k]), 
            "-", names(Ranks[k]))
        dimnames(G[[k]]) <- c(.dimList(names(Ranks[k - 1]), Ranks[k - 
            1]), dimnames(A)[k], .dimList(names(Ranks[k]), Ranks[k]))
    }
    names(G)[nModes] <- paste0(names(Ranks[nModes - 1]), "-", 
        names(dimnames(A)[nModes]))
    dimnames(G[[nModes]]) <- c(.dimList(names(Ranks[nModes - 
        1]), Ranks[nModes - 1]), dimnames(A)[nModes])
    G
}
.recError <-
function (X, X_bar) 
{
    v <- unlist(lapply(seq_len(length(X)), function(x) {
        X[[x]] - X_bar[[x]]
    }))
    sqrt(sum(v * v))
}
.recTensor <-
function (G) 
{
    nCores <- length(G)
    out <- G[[1]]
    if (nCores >= 2) {
        for (k in 2:(nCores)) {
            out <- CONTRACTION(X = out, z = G[[k]], Xwiz = length(dim(out)), 
                zwiX = 1)
        }
    }
    out
}
.tnorm <-
function (X) 
{
    if ("Tensor" %in% is(X)) {
        out <- sqrt(sum((X * X)@data))
    } else if ("simple_sparse_tensor" %in% class(X)) {
        out <- sqrt(sum(X$data * X$data))
    } else {
        # For regular arrays/matrices
        out <- sqrt(sum(X * X))
    }
    out
}
.genCores <-
function (X, Ranks) 
{
    nCores <- length(dim(X))
    D <- dim(X)
    G <- rep(list(NULL), length = nCores)
    G[[1]] <- matrix(runif(D[1] * Ranks[1]), nrow = D[1], ncol = Ranks[1])
    for (i in 2:(nCores - 1)) {
        G[[i]] <- rand_tensor(c(Ranks[i - 1], D[i], Ranks[i]))@data
    }
    G[[nCores]] <- matrix(runif(Ranks[nCores - 1] * D[nCores]), 
        nrow = Ranks[nCores - 1], ncol = D[nCores])
    .setDimNames(G, X@data, Ranks)
}
.is.sparse <-
function (A) 
{
    checkSimpleSparseTensor <- "simple_sparse_tensor" %in% class(A)
    checkSparseMatrix <- "sparseMatrix" %in% is(A)
    checkTensor <- "Tensor" %in% is(A)
    checkArray <- is.array(A) || is.matrix(A)
    
    if (!checkSimpleSparseTensor && !checkSparseMatrix && !checkTensor && !checkArray) {
        stop("Please specify A as a tensor, array, matrix, or sparseMatrix")
    }
}
