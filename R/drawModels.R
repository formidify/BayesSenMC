#' Plot Model
#'
#' Plot the posterior distribution of adjusted odds ratio given the stanfit object. 
#' @param model A stanfit object.
#' @param a # of exposed subjects in the case group. Along with N1, c, N0, se and sp, they are used to plot probability density with no misclassification and constant misclassification as a comparison.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param se sensitivity. Default to 1. If no other values are specified for either se or sp, then only the density curve of corrected model will be drawn.
#' @param sp specificity. Default to 1.
#' @param binwidth default to \code{0.25}
#' @param fill default to \code{"gray"}
#' @param ... optional additional arguments passed to \code{geom_histogram}
#' @return It returns a \code{\href{https://www.rdocumentation.org/packages/ggplot2/versions/0.9.0/topics/ggplot}{ggplot}} that can be further 
#' customized using the ggplot2 package.
#' @import rstan
#' @import Rcpp
#' @import ggplot2
#' @export
#' @examples
#' my_model <- crudeOR(a = 126, N1 = 218, c = 71, N0 = 295, se = 0.94, sp = 0.97, chains = 5, iter = 20000, warmup = 2000)
#' my_hist <- OR_hist(my_model, a = 126, N1 = 218, c = 71, N0 = 295, se = 0.94, sp = 0.97)

OR_hist <- function(model, a, N1, c, N0, se = 1, sp = 1, binwidth = 0.25, fill = "gray", ...) {
  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0) & (se >= 0) & (se <= 1) & (sp >= 0) & (sp <= 1))) {
    stop("The value(s) for a/N0/c/N1/se/sp is not valid.")
  }

  b <- N1 - a
  d <- N0 - c

  # when se and sp are both equal to 1 (corrected OR)
  ORadj <- ((a + (1 - 1) * N1) * (1 * N0 - c))/((c + (1 - 1) * N0) * (1 * N1 - a))
  term1 <- (a * b * (a + b))/(((-b + (a + b) * 1)^2) * ((1 * (a + b) - a)^2))
  term2 <- (c * d * (c + d))/(((-d + (c + d) * 1)^2) * ((1 * (c + d) - c)^2))
  se.log.or <- (1 + 1 - 1) * sqrt(term1 + term2)

  xx <- seq(0, 10, 0.1)
  yy1 <- dlnorm(xx, log(ORadj), se.log.or)

  y1 <- data.frame(x = xx, y = yy1)

  hist <- stan_hist(model, binwidth = binwidth, fill = fill, ...) +
            geom_line(mapping = aes(x = x, y = y), data = y1, linetype = 1, colour = "black") 

  if (!(se == 1 & sp == 1)) {
    # when se and sp are equal to user input values (crude OR)
    ORadj <- ((a + (sp - 1) * N1) * (se * N0 - c))/((c + (sp - 1) * N0) * (se * N1 - a))
    term1 <- (a * b * (a + b))/(((-b + (a + b) * sp)^2) * ((se * (a + b) - a)^2))
    term2 <- (c * d * (c + d))/(((-d + (c + d) * sp)^2) * ((se * (c + d) - c)^2))
    se.log.or <- (se + sp - 1) * sqrt(term1 + term2)

    yy2 <- dlnorm(xx, log(ORadj), se.log.or)

    y2 <- data.frame(x = xx, y = yy2)

    hist <- hist + geom_line(mapping = aes(x = x, y = y), data = y2, linetype = 2, colour = "black") +
              scale_linetype_manual(name = "", values = c(2, 1), labels = c("crude", "corrected"))
  }

  else {
    hist <- hist + scale_linetype_manual(name = "", values = 1, labels = "corrected")
  }
            
  print(hist)
  return(hist)
}