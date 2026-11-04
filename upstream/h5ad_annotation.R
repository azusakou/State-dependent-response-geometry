read_h5ad_annotation <- function(h5, group, field) {
  node <- h5[[group]][[field]]
  if (inherits(node, "H5Group")) {
    categories <- node[["categories"]][]
    codes <- as.integer(node[["codes"]][])
    out <- rep(NA, length(codes))
    keep <- codes >= 0L
    out[keep] <- categories[codes[keep] + 1L]
    out
  } else {
    node[]
  }
}

read_h5ad_obs <- function(h5, field) {
  read_h5ad_annotation(h5, "obs", field)
}

read_h5ad_var <- function(h5, field) {
  read_h5ad_annotation(h5, "var", field)
}

