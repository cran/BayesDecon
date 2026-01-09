#' Printing summary of MCMC samples for `bdeconv_result` object
#' @param x A `bdeconv_result` object.
#' @param ... Unused.
#' @return Prints the five number summary for the relevant variables.
#' @details
#' For episodic components, this gives five number summary for the observed surrogates (y), the latent surrogates corresponding to the consumption days (w), the true intake patterns (x) and the errors (u). For regular components, this gives five number summary for y, x and u only.
#'
#' @exportS3Method
print.bdeconv_result <- function(x, ...) {
  obj <- x
  if (!is.null(obj$call)) {
    cat('Call: ')
    print(obj$call)
  }
  cat('Number of samples: ')
  if (x$res_type == 'mean') {
    cat('1 mean of', obj$control$simsize - obj$control$burnin, '\n')
  } else {
    cat(length(obj), '\n')
  }
  cat('Number of variables:', obj$q + obj$p, '\n\n')
  NextMethod()
}

#' @exportS3Method
print.zi_res <- function(x, ...) {
  obj <- x
  qntl <- data.frame(y = stats::quantile(obj$y), w = stats::quantile(obj$w), x = stats::quantile(obj$x),
                     u = stats::quantile(obj$u))
  colnames(qntl) <- paste0(c('y','w', 'x', 'u'))
  rownames(qntl) <- c('Min.', '1Q', 'Med.', '3Q', 'Max.')
  cat('Variable:', obj$names, '\n')
  print(t(qntl))
}

#' @exportS3Method
print.nc_res <- function(x, ...) {
  obj <- x
  qntl <- data.frame(stats::quantile(obj$w), x = stats::quantile(obj$x),
                     u = stats::quantile(obj$u))
  colnames(qntl) <- paste0(c('y', 'x', 'u'))
  rownames(qntl) <- c('Min.', '1Q', 'Med.', '3Q', 'Max.')
  cat('Variable:', obj$names, '\n')
  print(t(qntl))
}

#' @exportS3Method
print.mv_res <- function(x, ...) {
  obj <- x
  if (obj$q > 0) {
    for (qq in 1:obj$q) {
      print.zi_res(obj[, qq])
      cat('\n')
    }
  }
  if (obj$p > 0) {
    for (pp in 1:obj$p) {
      print.nc_res(obj[, pp + obj$q])
      if (pp != obj$p) {
        cat('\n')
      }
    }
  }
}
