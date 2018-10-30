#' Meta-analysis data on Bipolar Disorder diagnosis accuracy
#'
#' Records the true positive, true negative, false positive and false negative
#' of each diagnosis accuracy study. Also includes the type of screening instruments
#' (Bipolar Spectrum diagnostic scale / HCL-21 / Mood disorder questionnaire), the cut-off
#' value for diagnostics, and the percentage of bipolar cases that were of bipolar disorder
#' type II or not specified.
#'
#' @docType data
#'
#' @usage data(bd_meta)
#'
#' @keywords dataset
#'
#' @references Carvalho et al. (2015) "Screening for bipolar spectrum disorders: A comprehensive meta-analysis of accuracy studies". Journal of Affective Disorders 172: 337 - 346.
#' (\href{http://www.sciencedirect.com/science/article/pii/S0165032714006466}{ScienceDirect})
#'
#' @source \url{https://www.sciencedirect.com/science/article/pii/S0165032714006466}
#' @examples
#' data(bd_meta)
#' \donttest{iplotCurves(phe, times)}
"bd_meta"
