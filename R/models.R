#' Model without misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of uncorrected odds ratio given no misclassification.
#' @param a number of exposed subjects in the case group.
#' @param N1 number of total subjects in the case group.
#' @param c number of exposed subjects in the control group.
#' @param N0 number of total subjects in the control group.
#' @param logitpi0_prior mean and sd of the prior normal distribution of \code{logit(pi0)}. Default to \code{c(0,10)}.
#' @param lor_prior mean and sd of the prior normal distribution of corrected log odds ratio. Default to \code{c(0,2)}.
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the \link[rstan]{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot.
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}.
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \link[rstan]{stan} for more details.
#' @param ... optional parameters passed to \link[rstan]{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See \link[rstan]{rstan} for more details.
#' @import rstan
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#'
#' # 3 MCMC chains with 10000 iterations each
#' crudeOR(a = 66, N1 = 11782, c = 243, N0 = 57973, chains = 3, iter = 10000)

crudeOR <- function(a, N1, c, N0, logitpi0_prior = c(0,10), lor_prior = c(0,2), chains = 2, traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = NA, ...) {
  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0))) {
    stop("The value(s) for a/N0/c/N1 is not valid.")
  }

  if (length(logitpi0_prior) != 2 | length(lor_prior) != 2 |
    !all(sapply(logitpi0_prior, is.numeric)) | !all(sapply(lor_prior, is.numeric))) {
    stop("The value(s) for logitpi0_prior/lor_prior is not valid.")
  }

  model <- rstan::sampling(stanmodels$crude, data = list(a = a, N1 = N1, c = c, N0 = N0, mLogit_pi0 = logitpi0_prior[1], sLogit_pi0 = logitpi0_prior[2],
                mLOR_c = lor_prior[1], sLOR_c = lor_prior[2]), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh, seed = seed, ...)

  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return(model)
}

#' Model with constant nondifferential misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification with Se and Sp (i.e., both are constant and at least one of Se or Sp is lower than 1).
#' @param a number of exposed subjects in the case group.
#' @param N1 number of total subjects in the case group.
#' @param c number of exposed subjects in the control group.
#' @param N0 number of total subjects in the control group.
#' @param prior_list list of priors. Can be replaced by the function call to \code{paramEst}, or a list of prior parameters (\code{se}, \code{sp}).
#' If \code{prior_list} is specified, the values for the function parameters \code{se} and \code{sp} will be disregarded.
#' @param se sensitivity. Do not have to specify this if \code{prior_list} is given - this will be disregarded.
#' @param sp specificity. Do not have to specify this if \code{prior_list} is given - this will be disregarded.
#' @param logitpi0_prior mean and sd of the prior normal distribution of \code{logit(pi0)}. Default to \code{c(0,10)}.
#' @param lor_prior mean and sd of the prior normal distribution of corrected log odds ratio. Default to \code{c(0,2)}.
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the \link[rstan]{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot.
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}.
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \link[rstan]{stan} for more details.
#' @param ... optional parameters passed to \link[rstan]{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See \link[rstan]{rstan} for more details.
#' @import rstan
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}\
#'
#' mod <- nlmeNDiff(bd_meta, lower = 0) # see \code{nlmeNDiff()} for detailed example.
#' prior_list <- paramEst(mod)
#' correctedOR(a = 66, N1 = 11782, c = 243, N0 = 57973, prior_list = prior_list,
#' chains = 3, iter = 10000)


correctedOR <- function(a, N1, c, N0, prior_list=NULL, se=NULL, sp=NULL, logitpi0_prior = c(0,10), lor_prior = c(0,2), chains = 2, traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = NA, ...) {

  if (!is.null(prior_list)) {
    if (length(prior_list) < 2 | !is.numeric(prior_list[['se']]) | !is.numeric(prior_list[['sp']])){
      message('Parameter "prior_list" does not contain all required prior parameters. Trying individual input instead.')
    }
    else {
      se <- prior_list[['se']]
      sp <- prior_list[['sp']]
    }
  }

  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0) & (se >= 0) & (se <= 1) & (sp >= 0) & (sp <= 1))) {
    stop("The value(s) for a/N0/c/N1/se/sp is not valid.")
  }

  if (length(logitpi0_prior) != 2 | length(lor_prior) != 2 |
    !all(sapply(logitpi0_prior, is.numeric)) | !all(sapply(lor_prior, is.numeric))) {
      stop("The value(s) for logitpi0_prior/lor_prior is not valid.")
    }

  model <- rstan::sampling(stanmodels$corrected, data = list(a = a, N1 = N1, c = c, N0 = N0, Se = se, Sp = sp, mLogit_pi0 = logitpi0_prior[1], sLogit_pi0 = logitpi0_prior[2],
                mLOR_c = lor_prior[1], sLOR_c = lor_prior[2]), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh, seed = seed, ...)

  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return(model)
}

#' Model with nondifferential, logit normal-distributed misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification under a logit-transformed scaled bivariate normal distribution.
#
#' @param a number of exposed subjects in the case group.
#' @param N1 number of total subjects in the case group.
#' @param c number of exposed subjects in the control group.
#' @param N0 number of total subjects in the control group.
#' @param prior_list list of priors. Can be replaced by the function call to \code{paramEst}, or a list of prior parameters
#' (\code{m.lg.se}, \code{s.lg.se}, \code{m.lg.sp}, \code{s.lg.sp}).
#' If \code{prior_list} is specified, the values for the corresponding function parameters will be disregarded.
#' @param m.lg.se normal distribution of logit Se with (mean = m.lg.se, sd = s.lg.se). Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param m.lg.sp normal distribution of logit Sp with (m.lg.sp, s.lg.sp). Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param s.lg.se standard deviation of logit Se. Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param s.lg.sp standard deviation of logit Sp. Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param lg.se used as an initial value for logit Se. Default to \code{m.lg.se}
#' @param lg.sp used as an initial value for logit Sp. Default to \code{m.lg.sp}
#' @param logitpi0_prior mean and sd of the prior normal distribution of \code{logit(pi0)}. Default to \code{c(0,10)}.
#' @param lor_prior mean and sd of the prior normal distribution of corrected log odds ratio. Default to \code{c(0,2)}.
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the \link[rstan]{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot.
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}.
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \link[rstan]{stan} for more details.
#' @param ... optional parameters passed to \link[rstan]{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See \link[rstan]{rstan} for more details.
#' @import rstan
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#'
#' mod <- nlmeNDiff(bd_meta, lower = 0) # see \code{nlmeNDiff()} for detailed example.
#' prior_list <- paramEst(mod)
#' logitOR(a = 66, N1 = 11782, c = 243, N0 = 57973, prior_list = prior_list,
#' chains = 3, iter = 10000)

logitOR <- function(a, N1, c, N0, prior_list = NULL, m.lg.se = NULL, m.lg.sp = NULL, s.lg.se = NULL, s.lg.sp = NULL, lg.se = NULL, lg.sp = NULL,
                    logitpi0_prior = c(0,10), lor_prior = c(0,2), chains = 2,
                    traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = NA, ...) {

  if (!is.null(prior_list)) {
    if (length(prior_list) < 4 | !is.numeric(prior_list[['m.lg.se']]) | !is.numeric(prior_list[['m.lg.sp']]) |
        !is.numeric(prior_list[['s.lg.se']]) | !is.numeric(prior_list[['s.lg.sp']])){
      message('Parameter "prior_list" does not contain all required prior parameters. Trying individual input instead.')
    }
    else {
      m.lg.se <- prior_list[['m.lg.se']]
      m.lg.sp <- prior_list[['m.lg.sp']]
      s.lg.se <- prior_list[['s.lg.se']]
      s.lg.sp <- prior_list[['s.lg.sp']]
    }
  }

  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0)
        & is.numeric(m.lg.se) & is.numeric(m.lg.sp) & is.numeric(s.lg.se) & is.numeric(s.lg.sp))) {
    stop("The value(s) for a/N0/c/N1/m.lg.se/m.lg.sp/s.lg.se/s.lg.sp is not valid.")
  }

  # default to mean if not specified
  if (is.null(lg.se)) {lg.se <- m.lg.se}
  if (is.null(lg.sp)) {lg.sp <- m.lg.sp}

  if (!is.numeric(lg.se) | !is.numeric(lg.sp)) {
    stop("The value(s) for lg.se/lg.sp is not numeric.")
  }

  if (length(logitpi0_prior) != 2 | length(lor_prior) != 2 |
    !all(sapply(logitpi0_prior, is.numeric)) | !all(sapply(lor_prior, is.numeric))) {
      stop("The value(s) for logitpi0_prior/lor_prior is not valid.")
    }

  model <- rstan::sampling(stanmodels$logit, data=list(a = a, N1 = N1, c = c, N0 = N0, mX = m.lg.se, mY = m.lg.sp,
                sdX = s.lg.se, sdY = s.lg.sp, mLogit_pi0 = logitpi0_prior[1], sLogit_pi0 = logitpi0_prior[2],
                mLOR_c = lor_prior[1], sLOR_c = lor_prior[2]), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh,
                init = rep(list(list(X = lg.se, Y = lg.sp)), chains), seed = seed, ...)


  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return(model)
}

#' Model with nondifferential, correlated misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification that extends from the logit model but allows there to be a fixed correlation between sentivity and specificity.
#' @param a number of exposed subjects in the case group.
#' @param N1 number of total subjects in the case group.
#' @param c number of exposed subjects in the control group.
#' @param N0 number of total subjects in the control group.
#' @param prior_list list of priors. Can be replaced by the function call to \code{paramEst}, or a list of prior parameters
#' (\code{m.lg.se}, \code{s.lg.se}, \code{m.lg.sp}, \code{s.lg.sp}, \code{rho}).
#' If \code{prior_list} is specified, the values for the corresponding function parameters will be disregarded.
#' @param m.lg.se normal distribution of logit Se with (mean = m.lg.se, sd = s.lg.se). Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param m.lg.sp conditional normal distribution of logit Sp given Se with (m.lg.sp, s.lg.sp). Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param s.lg.se standard deviation of logit Se. Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param s.lg.sp standard deviation of logit Sp. Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param lg.se used as an initial value for logit Se. Default to \code{m.lg.se}
#' @param lg.sp used as an initial value for logit Sp. Default to \code{m.lg.sp}
#' @param rho correlation between Se and Sp. Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param logitpi0_prior mean and sd of the prior normal distribution of \code{logit(pi0)}. Default to \code{c(0,10)}.
#' @param lor_prior mean and sd of the prior normal distribution of corrected log odds ratio. Default to \code{c(0,2)}.
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the \link[rstan]{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot.
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}.
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \link[rstan]{stan} for more details.
#' @param ... optional parameters passed to \link[rstan]{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See \link[rstan]{rstan} for more details.
#' @import rstan
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#'
#' mod <- nlmeNDiff(bd_meta, lower = 0) # see \code{nlmeNDiff()} for detailed example.
#' prior_list <- paramEst(mod)
#' fixedCorrOR(a = 66, N1 = 11782, c = 243, N0 = 57973, prior_list = prior_list,
#' chains = 3, iter = 10000)

fixedCorrOR <- function(a, N1, c, N0, prior_list = NULL, m.lg.se = NULL, m.lg.sp = NULL, s.lg.se = NULL, s.lg.sp = NULL, lg.se = NULL, lg.sp = NULL,
                        rho = NULL, logitpi0_prior = c(0,10), lor_prior = c(0,2), chains = 2,
                        traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = NA, ...) {

  if (!is.null(prior_list)) {
    if (length(prior_list) < 5 | !is.numeric(prior_list[['m.lg.se']]) | !is.numeric(prior_list[['m.lg.sp']]) |
       !is.numeric(prior_list[['s.lg.se']]) | !is.numeric(prior_list[['s.lg.sp']]) | !is.numeric(prior_list[['rho']])){
         message('Parameter "prior_list" does not contain all required prior parameters. Trying individual input instead.')
       }
     else {
       m.lg.se <- prior_list[['m.lg.se']]
       m.lg.sp <- prior_list[['m.lg.sp']]
       s.lg.se <- prior_list[['s.lg.se']]
       s.lg.sp <- prior_list[['s.lg.sp']]
       rho <- prior_list[['rho']]
     }
  }

  if (!(((a <= N1)) & (a >= 0) & (c <= N0) & (c >= 0) & (-1 <= rho) & (rho <= 1)
    & is.numeric(m.lg.se) & is.numeric(m.lg.sp) & is.numeric(s.lg.se) & is.numeric(s.lg.sp))) {
    stop("The value(s) for a/N0/c/N1/rho/m.lg.se/m.lg.sp/s.lg.se/s.lg.sp is not valid.")
  }

  # default to mean if not specified
  if (is.null(lg.se)) {lg.se <- m.lg.se}
  if (is.null(lg.sp)) {lg.sp <- m.lg.sp}

  if (!is.numeric(lg.se) | !is.numeric(lg.sp)) {
    stop("The value(s) for lg.se/lg.sp is not numeric.")
  }

  if (length(logitpi0_prior) != 2 | length(lor_prior) != 2 |
    !all(sapply(logitpi0_prior, is.numeric)) | !all(sapply(lor_prior, is.numeric))) {
      stop("The value(s) for logitpi0_prior/lor_prior is not valid.")
    }

  model <- rstan::sampling(stanmodels$fixedCorr, data = list(a = a, N1 = N1, c = c, N0 = N0, mX0 = m.lg.se, mX1 = m.lg.sp,
                precX0 = 1 / (s.lg.se)^2, precX1 = 1 / (s.lg.sp)^2, rhoSe = rho, mLogit_pi0 = logitpi0_prior[1], sLogit_pi0 = logitpi0_prior[2],
                mLOR_c = lor_prior[1], sLOR_c = lor_prior[2]), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh,
                init = rep(list(list(X0 = lg.se, X1 = lg.sp)), chains), seed = seed, ...)

  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return(model)
}

#' Model with nondifferential, randomly correlated misclassification

#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification that extends from the logit model but allows a random correlation between Sensitivity and Specificity.
#' @param a number of exposed subjects in the case group.
#' @param N1 number of total subjects in the case group.
#' @param c number of exposed subjects in the control group.
#' @param N0 number of total subjects in the control group.
#' @param prior_list list of priors. Can be replaced by the function call to \code{paramEst}, or a list of prior parameters
#' (\code{m.lg.se}, \code{s.lg.se}, \code{m.lg.sp}, \code{s.lg.sp}, \code{m.z}, \code{s.z}).
#' If \code{prior_list} is specified, the values for the corresponding function parameters will be disregarded.
#' @param m.lg.se normal distribution of logit Se with (mean = m.lg.se, sd = s.lg.se). Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param m.lg.sp conditional normal distribution of logit Sp given Se with (m.lg.sp, s.lg.sp). Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param s.lg.se standard deviation of logit Se. Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param s.lg.sp standard deviation of logit Sp. Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param lg.se used as an initial value for logit Se. Default to m.lg.se. Do not have to specify this if \code{prior_list} is given - it will be disregarded. Default to \code{m.lg.se}
#' @param lg.sp used as an initial value for logit Sp. Default to m.lg.sp. Do not have to specify this if \code{prior_list} is given - it will be disregarded. Default to \code{m.lg.sp}
#' @param z used as an initial value of Fisher's Z transformed of rho, where correlation rho = (exp(2*z)-1)/(1+exp(2*z))). Do not have to specify this if \code{prior_list} is given - it will be disregarded. Default to \code{m.z}
#' @param m.z normal distribution of Z with (mean = m.z, sd = s.z). Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param s.z normal distribution of Z with (mean = m.z, sd = s.z). Do not have to specify this if \code{prior_list} is given - it will be disregarded.
#' @param logitpi0_prior mean and sd of the prior normal distribution of \code{logit(pi0)}. Default to \code{c(0,10)}.
#' @param lor_prior mean and sd of the prior normal distribution of corrected log odds ratio. Default to \code{c(0,2)}.
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the \link[rstan]{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot.
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}.
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \link[rstan]{stan} for more details.
#' @param ... optional parameters passed to \link[rstan]{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See \link[rstan]{rstan} for more details.
#' @import rstan
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#'
#' mod <- nlmeNDiff(bd_meta, lower = 0) # see \code{nlmeNDiff()} for detailed example.
#' prior_list <- paramEst(mod)
#' randCorrOR(a = 66, N1 = 11782, c = 243, N0 = 57973, prior_list = prior_list,
#' chains = 3, iter = 10000)

randCorrOR <- function(a, N1, c, N0, prior_list = NULL, m.lg.se = NULL, m.lg.sp = NULL, s.lg.se = NULL, s.lg.sp = NULL, lg.se = NULL, lg.sp = NULL, m.z = NULL, s.z = NULL, z = NULL,
                       logitpi0_prior = c(0,10), lor_prior = c(0,2), chains = 2,
                       traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = NA, ...) {

  if (!is.null(prior_list)) {
    if (length(prior_list) < 6 | !is.numeric(prior_list[['m.lg.se']]) | !is.numeric(prior_list[['m.lg.sp']]) |
     !is.numeric(prior_list[['s.lg.se']]) | !is.numeric(prior_list[['s.lg.sp']]) | !is.numeric(prior_list[['m.z']]) |
      !is.numeric(prior_list[['s.z']])){
       message('Parameter "prior_list" does not contain all required prior parameters. Trying individual input instead.')
     }
     else {
       m.lg.se <- prior_list[['m.lg.se']]
       m.lg.sp <- prior_list[['m.lg.sp']]
       s.lg.se <- prior_list[['s.lg.se']]
       s.lg.sp <- prior_list[['s.lg.sp']]
       m.z <- prior_list[['m.z']]
       s.z <- prior_list[['s.z']]
     }
  }

  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0) & (-1 <= m.z) & (m.z <= 1) & is.numeric(s.z)
     & is.numeric(m.lg.se) & is.numeric(m.lg.sp) & is.numeric(s.lg.se) & is.numeric(s.lg.sp))) {
       stop("The value(s) for a/N0/c/N1/rho/m.lg.se/m.lg.sp/s.lg.se/s.lg.sp is not valid.")
    }

  # default to mean if not specified
  if (is.null(lg.se)) {lg.se <- m.lg.se}
  if (is.null(lg.sp)) {lg.sp <- m.lg.sp}
  if (is.null(z)) {z <- m.z}

  if (!is.numeric(lg.se) | (!is.numeric(lg.sp)) | (!is.numeric(m.z))) {
    stop("The value(s) for lg.se/lg.sp/m.z is not numeric.")
  }

  if (length(logitpi0_prior) != 2 | length(lor_prior) != 2 |
    !all(sapply(logitpi0_prior, is.numeric)) | !all(sapply(lor_prior, is.numeric))) {
      stop("The value(s) for logitpi0_prior/lor_prior is not valid.")
    }

  model <- rstan::sampling(stanmodels$randCorr, data = list(a = a, N1 = N1, c = c, N0 = N0, mX0 = m.lg.se, mX1 = m.lg.sp,
                precX0 = 1 / (s.lg.se)^2, precX1 = 1 / (s.lg.sp)^2, mZ = m.z, sZ = s.z, mLogit_pi0 = logitpi0_prior[1], sLogit_pi0 = logitpi0_prior[2],
                mLOR_c = lor_prior[1], sLOR_c = lor_prior[2]), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh,
                init = rep(list(list(X0 = lg.se, X1 = lg.sp, Z = z)), chains), seed = seed, ...)

  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return (model)
}

#' Model with differential misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given a four-variate differential misclassification.
#' @param a number of exposed subjects in the case group.
#' @param N1 number of total subjects in the case group.
#' @param c number of exposed subjects in the control group.
#' @param N0 number of total subjects in the control group.
#' @param mu vector of length 4; multivariate normal distribution of \eqn{z \sim (mu, varz)}, where each \eqn{\mu} corresponds to the logit mean of \eqn{Se_0}, \eqn{Se_1}, \eqn{Sp_0} and \eqn{Sp_1} (0 for controls, 1 for cases group).
#' @param s.lg.se0 standard deviation of logit Se in the control group.
#' @param s.lg.se1 standard deviation of logit Se in the case group.
#' @param s.lg.sp0 standard deviation of logit Sp in the control group.
#' @param s.lg.sp1 standard deviation of logit Sp in the case group.
#' @param corr.sesp0 correlation between Se_0 and Sp_0.
#' @param corr.sesp1 correlation between Se_1 and Sp_1.
#' @param corr.group correlation between Se_0 and Se_1, Sp_0 and Sp_1. Default to 0.
#' @param z vector of length 4; used as an initial value for \eqn{z \sim (mu, varz)}. Default to mu.
#' @param logitpi0_prior mean and sd of the prior normal distribution of \code{logit(pi0)}. Default to \code{c(0,10)}.
#' @param lor_prior mean and sd of the prior normal distribution of corrected log odds ratio. Default to \code{c(0,2)}.
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the \link[rstan]{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot.
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}.
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \link[rstan]{stan} for more details.
#' @param ... optional parameters passed to \link[rstan]{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See \link[rstan]{rstan} for more details.
#' @import rstan
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#'
#' diffOR(a = 66, N1 = 11782, c = 243, N0 = 57973, mu = c(1.069, 1.069, 1.126, 1.126),
#'   s.lg.se0 = 0.712, s.lg.se1 = 0.712, s.lg.sp0 = 0.893, s.lg.sp1 = 0.893, corr.sesp0 = -0.377,
#'   corr.sesp1 = -0.377, corr.group = 0, chains = 3, iter = 10000)


diffOR <- function(a, N1, c, N0, mu, s.lg.se0, s.lg.se1, s.lg.sp0, s.lg.sp1, corr.sesp0, corr.sesp1,
                      corr.group = 0, z = NULL, logitpi0_prior = c(0,10), lor_prior = c(0,2), chains = 2,
                      traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = 0, ...) {

  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0))) {
    stop("The value(s) for a/N0/c/N1 is not valid.")
  }

  # default to mean if not specified
  if (is.null(z)) {z <- mu}

  if (length(logitpi0_prior) != 2 | length(lor_prior) != 2 |
    !all(sapply(logitpi0_prior, is.numeric)) | !all(sapply(lor_prior, is.numeric))) {
      stop("The value(s) for logitpi0_prior/lor_prior are not valid.")
    }

  varz = matrix(c(s.lg.se0^2, corr.group*s.lg.se0*s.lg.se1, corr.sesp0*s.lg.se0*s.lg.sp0, corr.group*s.lg.se0*s.lg.sp1,
                  corr.group*s.lg.se0*s.lg.se1, s.lg.se1^2, corr.group*s.lg.se1*s.lg.sp0, corr.sesp1*s.lg.se1*s.lg.sp1,
                  corr.sesp0*s.lg.se0*s.lg.sp0, corr.group*s.lg.se1*s.lg.sp0, s.lg.sp0^2, corr.group*s.lg.sp0*s.lg.sp1,
                  corr.group*s.lg.se0*s.lg.sp1, corr.sesp1*s.lg.se1*s.lg.sp1, corr.group*s.lg.sp0*s.lg.sp1, s.lg.sp1^2), 4, 4)

  model <- rstan::sampling(stanmodels$diff, data = list(a = a, N1 = N1, c = c, N0 = N0, Mu = mu,
    varZ = varz, mLogit_pi0 = logitpi0_prior[1], sLogit_pi0 = logitpi0_prior[2], mLOR_c = lor_prior[1], sLOR_c = lor_prior[2]),
    pars = c("LOR_c", "ORadj"), chains = chains, init = rep(list(list(Z = z)), chains), refresh = refresh, seed = seed, ...)

  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return(model)
}
