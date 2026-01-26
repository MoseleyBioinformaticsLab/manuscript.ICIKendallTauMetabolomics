# Copied from correlation package, matrix_inverse.R, which is
# licensed under the GPL3 license.
# 
#' Matrix Inversion
#'
#' Performs a Moore-Penrose generalized inverse (also called the Pseudoinverse).
#'
#' @inheritParams cor_to_pcor
#' @examples
#' m <- cor(iris[1:4])
#' matrix_inverse(m)
#' @param m Matrix for which the inverse is required.
#'
#' @return An inversed matrix.
#' @seealso pinv from the pracma package
#' @export
matrix_inverse <- function(m, tol = .Machine$double.eps^(2 / 3)) {
  # valid matrix checks
  # valid matrix checks
  if (!isSquare(m)) {
    stop("The matrix should be a square matrix.", call. = FALSE)
  }
  
  stopifnot(is.numeric(m), length(dim(m)) == 2, is.matrix(m))
  
  s <- svd(m)
  
  p <- (s$d > max(tol * s$d[1], 0))
  if (all(p)) {
    mp <- s$v %*% (1 / s$d * t(s$u))
  } else if (any(p)) {
    mp <- s$v[, p, drop = FALSE] %*% (1 / s$d[p] * t(s$u[, p, drop = FALSE]))
  } else {
    mp <- matrix(0, nrow = ncol(m), ncol = nrow(m))
  }
  
  colnames(mp) <- colnames(m)
  row.names(mp) <- row.names(m)
  mp
}


invert_matrix <- function(m, tol = .Machine$double.eps^(2 / 3)) {
  if (det(m) < tol) {
    # The inverse of variance-covariance matrix is calculated using
    # Moore-Penrose generalized matrix inverse due to its determinant of zero.
    out <- matrix_inverse(m, tol)
    colnames(out) <- colnames(m)
    row.names(out) <- row.names(m)
  } else {
    out <- solve(m)
  }
  out
}

# copied from correlation package, cor_to_pcor.R
# licensed under the GPL3 license
cor_to_pcor <- function(cor, tol = .Machine$double.eps^(2 / 3)) {
  
  inverted <- invert_matrix(cor, tol = tol)
  out <- -stats::cov2cor(inverted)
  
  diag(out) <- 1
  out
}

# copied from correlation package, utils.R
# licensed under the GPL3 license
isSquare <- function(m) {
  if (dim(m)[1] != dim(m)[2]) {
    FALSE
  } else {
    TRUE
  }
}
