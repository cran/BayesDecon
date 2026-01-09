#' Generate a Penalty Matrix
#'
#' @param K An integer dimensino value.
#'
#' @return A first-order difference penalty matrix.
#' @keywords internal
P1.mat <- function(K) {
  P <- diff(diag(K))
  t(P) %*% P
}

#' Generate a Penalty Matrix
#'
#' @param K An integer dimensino value.
#'
#' @return A second-order difference penalty matrix.
#' @keywords internal
P2.mat <- function(K) {
  P <- diff(diff(diag(K)))
  t(P) %*% P
}
