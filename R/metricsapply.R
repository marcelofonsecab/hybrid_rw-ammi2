#' Apply a Metric to a List of Data and Optionally Test Hypotheses
#'
#' This function applies a specified metric to each element in a list of data,
#' and can optionally perform hypothesis testing on the metric values.
#'
#' @param list A list of data elements.
#' @param metric A character string specifying the metric to extract from each element.
#' @param hypothesis.testing Logical; if \code{TRUE}, performs a t-test for each metric against \code{h0}. Default is \code{FALSE}.
#' @param h0 A numeric value for the null hypothesis mean in t-tests. Default is 100.
#' @param side A character string specifying the alternative hypothesis for t-tests (\code{"greater"} or \code{"less"}). Default is \code{"greater"}.
#'
#' @return A list containing the calculated metric values and, if requested, the corresponding p-values.
#'
#' @examples
#' # Assuming a list of simulated metric outputs:
#' metric_list <- list(model1 = list(Metric = c(101, 102)),
#'                     model2 = list(Metric = c(98, 99)))
#' result <- Metricsapply(metric_list, metric = "Metric", hypothesis.testing = TRUE, h0 = 100)
#' result
#'
#' @export
Metricsapply <- function(list, metric, hypothesis.testing = FALSE, h0 = 100, side = "greater") {
  len <- length(list)
  tmp <- tmp_aux <- test_aux <- NULL
  for (i in 1:len) {
    tmp_aux <- mean(list[[i]][[metric]])
    tmp <- c(tmp, tmp_aux)

    if (hypothesis.testing) {
      test_aux <- c(test_aux, t.test(x = list[[i]][[metric]],
                                     alternative = side,
                                     mu = h0)$p.value)
    }
  }

  if (hypothesis.testing) {
    Values <- list('Metric' = tmp, 'pvalues' = test_aux)
  } else {
    Values <- list('Metric' = tmp)
  }

  return(Values)
}
