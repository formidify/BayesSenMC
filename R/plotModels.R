#' Plot Model
#'
#' Plot the posterior distribution of adjusted odds ratio given the stanfit object. It also plots the density lines of crude/uncorrected odds ratio and corrected odds ratio with constant misclassification, assuming log-normality is true. If both Se and Sp are set to 1 (i.e., no misclassification), then only the density line of crude OR will be plotted.
#' @param model A stanfit object.
#' @param a number of exposed subjects in the case group. Along with N1, c, N0, se and sp, they are used to plot probability density with no misclassification and constant misclassification as a comparison.
#' @param N1 number of total subjects in the case group.
#' @param c number of exposed subjects in the control group.
#' @param N0 number of total subjects in the control group.
#' @param se sensitivity. Default to 1. If no other values are specified for either se or sp, then only the density curve of corrected model will be drawn.
#' @param sp specificity. Default to 1.
#' @param x.min shows only samples with corrected odds ratio larger or equal to \code{x.min}. Default to 0.
#' @param x.max shows only samples with corrected odds ratio smaller or equal to \code{x.max}. Default to the largest OR in the posterior samples.
#' @param y.max shows only samples or density line within the range of (0, \code{y.max}).
#' @param binwidth default to \code{0.25}
#' @param fill default to \code{"gray"}
#' @param ... optional additional arguments passed to \code{geom_histogram}
#' @return It returns a \link[ggplot2]{ggplot} that can be further
#' customized using the ggplot2 package.
#' @import rstan
#' @import ggplot2
#' @importFrom stats dlnorm
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#'
#' library(ggplot2)
#' my.mod <- randCorrOR(a = 66, N1 = 11782, c = 243, N0 = 57973, m.lg.se = 1.069,
#' m.lg.sp = 1.126, s.lg.se = 0.893, s.lg.sp = 0.712, m.z = -0.399, s.z = 0.139,
#' seed = 0)
#'
#' my.plot <- plotOR(my.mod, a = 66, N1 = 11782, c = 243, N0 = 57973, se = 0.744,
#' sp = 0.755, x.max = 3, y.max = 5, binwidth = 0.1) + ggtitle("Model with random correlation")
#'
#' # the user can also directly extract the data from a stanfit object using the following
#' my.data <- as.data.frame(my.mod)

plotOR <- function(model, a, N1, c, N0, se = 1, sp = 1, x.min = 0, x.max = NULL, y.max = NULL, binwidth = 0.25, fill = "gray", ...) {
  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0) & (se >= 0) & (se <= 1) & (sp >= 0) & (sp <= 1))) {
    stop("The value(s) for a/N0/c/N1/se/sp is not valid.")
  }

  or_type <- NULL
  b <- N1 - a
  d <- N0 - c

  if (is.null(x.max)) {
    # locate the maximum log odds ratio from the input model posterior samples
    x.max <- summary(model, probs = c(1))$summary['ORadj','100%']
  }

  # when se and sp are both equal to 1 (crude OR)
  ORadj <- ((a + (1 - 1) * N1) * (1 * N0 - c))/((c + (1 - 1) * N0) * (1 * N1 - a))
  term1 <- (a * b * (a + b))/(((-b + (a + b) * 1)^2) * ((1 * (a + b) - a)^2))
  term2 <- (c * d * (c + d))/(((-d + (c + d) * 1)^2) * ((1 * (c + d) - c)^2))
  se.log.or <- (1 + 1 - 1) * sqrt(term1 + term2)

  xx <- seq(x.min, x.max, 0.01)
  yy1 <- dlnorm(xx, log(ORadj), se.log.or)

  y1 <- data.frame(xx, yy1)
  y1$or_type <- 'crude'

  hist <- stan_hist(model, pars = "ORadj", binwidth = binwidth, fill = fill, ...) +
            geom_line(mapping = aes(x = xx, y = yy1, linetype = or_type), data = y1, colour = "black")

  if (!(se == 1 & sp == 1)) {
    # when se and sp are equal to user input values (corrected OR)
    ORadj <- ((a + (sp - 1) * N1) * (se * N0 - c))/((c + (sp - 1) * N0) * (se * N1 - a))
    term1 <- (a * b * N1)/(((-b + N1 * sp)^2) * ((se * N1 - a)^2))
    term2 <- (c * d * N0)/(((-d + N0 * sp)^2) * ((se * N0 - c)^2))
    se.log.or <- (se + sp - 1) * sqrt(term1 + term2)

    yy2 <- dlnorm(xx, log(ORadj), se.log.or)

    y2 <- data.frame(xx, yy2)
    y2$or_type <- 'corrected'

    hist <- hist + geom_line(mapping = aes(x = xx, y = yy2, linetype = or_type), data = y2, colour = "black")
  }

  # update minimum and maximum x scale values
  if (is.null(y.max)) {
    hist <- hist + coord_cartesian(xlim = c(x.min, x.max))
  }
  else {
    hist <- hist + coord_cartesian(xlim = c(x.min, x.max), ylim = c(0, y.max), expand = FALSE)
  }

  hist <- hist + theme_classic() + theme(legend.justification = "top")

  return(hist)
}
