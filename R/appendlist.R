#' Append Two Lists Recursively
#'
#' This utility function appends values from one list to another.
#'
#' @param x A list to which values will be appended.
#' @param val A list of values to append.
#'
#' @return A merged list with values from both \code{x} and \code{val}.
#'
#' @examples
#' list1 <- list(a = 1, b = list(c = 2))
#' list2 <- list(b = list(d = 3), e = 4)
#' combined <- appendList(list1, list2)
#' combined
#'
#' @export
appendList <- function(x, val) {
  stopifnot(is.list(x), is.list(val))
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
      appendList(x[[v]], val[[v]])
    else
      c(x[[v]], val[[v]])
  }
  x
}
