suppressPackageStartupMessages(library(data.table))
safe_cosine <- function(x,y,eps=1e-10){den<-sqrt(sum(x*x))*sqrt(sum(y*y));if(!is.finite(den)||den<eps)return(NA_real_);sum(x*y)/den}
safe_norm <- function(x)sqrt(sum(x*x))
standard_error <- function(x){x<-x[is.finite(x)];if(length(x)<2)return(NA_real_);sd(x)/sqrt(length(x))}

# Read losslessly partitioned condition/GEM response profiles.
read_norman_input <- function(input_dir = "inputs") {
  x <- readRDS(file.path(input_dir, "norman_summary.rds"))
  layout <- attr(x, "gem_storage", exact = TRUE)
  if (is.null(layout)) stop("Missing Norman GEM storage layout.")
  covered <- unlist(lapply(layout$parts, `[[`, "indices"), use.names = FALSE)
  if (!identical(covered, seq_len(layout$dim[3L]))) stop("Invalid Norman GEM partition order.")
  response <- array(NA_real_, dim = layout$dim, dimnames = layout$dimnames)
  for (part in layout$parts) {
    path <- file.path(input_dir, part$file)
    if (!file.exists(path)) stop("Missing Norman input: ", part$file)
    if (unname(tools::md5sum(path)) != part$md5) stop("Norman input checksum mismatch: ", part$file)
    block <- readRDS(path)
    expected_names <- layout$dimnames
    if (!is.null(expected_names)) expected_names[[3L]] <- expected_names[[3L]][part$indices]
    if (!identical(dim(block), c(layout$dim[1:2], length(part$indices))) ||
        !identical(dimnames(block), expected_names)) stop("Norman GEM layout mismatch: ", part$file)
    response[, , part$indices] <- block
  }
  attr(x, "gem_storage") <- NULL
  x$delta_gem <- response
  x
}
