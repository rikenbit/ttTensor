# Simple sparse tensor implementation to replace tensorr dependency
# This is a minimal implementation to support TTCross function

# Simple sparse tensor class (minimal implementation)
as_sptensor <- function(x) {
    if (is.array(x) || is.matrix(x)) {
        # For simplicity, just wrap the array/matrix
        # In a real sparse implementation, we'd store only non-zero values
        structure(list(data = x, dims = dim(x)), class = "simple_sparse_tensor")
    } else {
        stop("Input must be an array or matrix")
    }
}

# Convert dense tensor to simple sparse format
dtensor <- function(x) {
    # Just return the input as-is for dense representation
    # This maintains compatibility with the test code
    x
}

# Unfold for our simple sparse tensor
unfold <- function(x, mode) {
    if ("simple_sparse_tensor" %in% class(x)) {
        # Use rTensor's unfold functions
        tensor_obj <- as.tensor(x$data)
        result <- rs_unfold(tensor_obj, m = mode)
        return(result)
    } else if ("Tensor" %in% class(x)) {
        # For regular tensors, use rTensor's unfold
        return(rs_unfold(x, m = mode))
    } else if (is.array(x) || is.matrix(x)) {
        # Convert to tensor first
        tensor_obj <- as.tensor(x)
        return(rs_unfold(tensor_obj, m = mode))
    } else {
        stop("Unsupported type for unfold")
    }
}

# Add dim method for simple_sparse_tensor
dim.simple_sparse_tensor <- function(x) {
    x$dims
}

# Add indexing support for simple_sparse_tensor
`[.simple_sparse_tensor` <- function(x, ...) {
    x$data[...]
}

# Add dimnames support
`dimnames<-.simple_sparse_tensor` <- function(x, value) {
    dimnames(x$data) <- value
    x
}

dimnames.simple_sparse_tensor <- function(x) {
    dimnames(x$data)
}