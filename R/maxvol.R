maxvol <-
function (C) 
{
    .is.sparse(C)
    nr <- dim(C)
    n <- nr[1]
    r <- nr[2]
    if (n < r) {
        stop("Please confirm that the number of rows are larger than the number of columns of matrix.")
    }
    row_idx <- rep(0, r)
    rest_of_rows <- seq_len(n)
    i = 1
    C_new = C
    while (i <= r) {
        if (is.null(dim(C_new))) {
            mask = 1
            rows_norms = sum(C_new^2)
            l <- 1
        }
        else {
            mask = seq_len(nrow(C_new))
            rows_norms = apply(C_new^2, 1, sum)
            l <- length(which(rows_norms != 0))
        }
        if (l == 1) {
            row_idx[i] <- rest_of_rows
            break
        }
        if (l != length(rows_norms)) {
            zero_idx <- which(rows_norms == min(rows_norms))
            mask <- mask[setdiff(seq_len(length(mask)), zero_idx)]
            rest_of_rows <- rest_of_rows[mask]
            C_new = C_new[mask, ]
            next
        }
        max_row_idx = which(rows_norms == max(rows_norms))[1]
        max_row = C_new[max_row_idx, ]
        projection = C_new %*% max_row
        normalization = sqrt(rows_norms[max_row_idx] * rows_norms)
        projection = projection/normalization
        for (j in 1:ncol(C_new)) {
            C_new[, j] = C_new[, j] - C_new[, j] * projection
        }
        mask = mask[setdiff(seq_along(mask), max_row_idx)]
        C_new = C_new[mask, ]
        row_idx[i] = rest_of_rows[max_row_idx]
        rest_of_rows = rest_of_rows[mask]
        i = i + 1
    }
    row_idx
}
