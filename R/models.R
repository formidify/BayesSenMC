#' Model without misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of uncorrected odds ratio given no misclassification.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param name a string of the name of the model. Default to "Corrected Model".
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot} corresponding to one or more Markov chains. 
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan} for more details.
#' @param ... optional parameters passed to \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#' 
#' # 3 MCMC chains with 10000 iterations each
#' \dontrun{
#' correctedOR(a = 66, N1 = 11782, c = 243, N0 = 57973, chains = 3, iter = 10000, seed = 0)
#' }
#' correctedOR(a = 66, N1 = 11782, c = 243, N0 = 57973, traceplot = TRUE)

correctedOR <- function(a, N1, c, N0, name = "Corrected Model", chains = 2, traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = NA, ...) {
  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0))) {
    stop("The value(s) for a/N0/c/N1 is not valid.")
  }
  
  options <- list(...)
  
  code <- "
  data {
    int<lower=0> a;
    int<lower=0> N1;
    int<lower=0> c;
    int<lower=0> N0;
  }
  parameters {
    real logit_pi0;
    real LOR_c;
  }
  transformed parameters {
    real<lower=0, upper=1> pi1;
    real<lower=0, upper=1> pi0;
    real ORadj;
    pi0 = exp(logit_pi0) / (exp(logit_pi0) + 1);
    pi1 = exp(logit_pi0 + LOR_c) / (exp(logit_pi0 + LOR_c) + 1);
    ORadj = exp(LOR_c);
  }
  model {
    a ~ binomial(N1, pi1);
    c ~ binomial(N0, pi0);
    logit_pi0 ~ normal(0, 10);
    LOR_c ~ normal(0, 2);
  }"
  
  # if user does not specify control parameters
  # default set to smaller step size to improve divergence in some cases
  if ('control' %in% names(options)) {
    model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh, seed = seed, ...)
  }
  else {
    model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0), 
                  pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh, seed = seed,
                  control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 50), ...)
  }
  
  print(summary(model)$summary)
  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return(model)
}

#' Model with constant nondifferential misclassification 
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification with constant Se and Sp.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param se sensitivity
#' @param sp specificity
#' @param name a string of the name of the model. Default to "Constant Misclassification Model".
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan} for more details.
#' @param ... optional parameters passed to \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}\
#' 
#' \dontrun{
#' crudeOR(a = 66, N1 = 11782, c = 243, N0 = 57973, se = 0.744, sp = 0.755, chains = 3, 
#' iter = 10000, seed = 0)
#' }
#' crudeOR(a = 66, N1 = 11782, c = 243, N0 = 57973, se = 0.744, sp = 0.755, traceplot = TRUE)

crudeOR <- function(a, N1, c, N0, se, sp, name = "Constant Misclassification Model", chains = 2, traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = NA, ...) {
  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0) & (se >= 0) & (se <= 1) & (sp >= 0) & (sp <= 1))) {
    stop("The value(s) for a/N0/c/N1/se/sp is not valid.")
  }
  
  options <- list(...)
  
  code <- "
  data {
    int<lower=0> a;
    int<lower=0> N1;
    int<lower=0> c;
    int<lower=0> N0;
    real<lower=0, upper=1> Se;
    real<lower=0, upper=1> Sp;
  }
  parameters {
    real logit_pi0;
    real LOR_c;
  }
  transformed parameters {
    real<lower=0, upper=1> pi1;
    real<lower=0, upper=1> pi0;
    real ORadj;
    real<lower=0, upper=1> p1;
    real<lower=0, upper=1> p0;
    pi0 = exp(logit_pi0) / (exp(logit_pi0) + 1);
    pi1 = exp(logit_pi0 + LOR_c) / (exp(logit_pi0 + LOR_c) + 1);
    p1 = pi1 * Se + (1 - pi1) * (1 - Sp);
    p0 = pi0 * Se + (1 - pi0) * (1 - Sp);
    ORadj = exp(LOR_c);
  }
  model {
    a ~ binomial(N1, p1);
    c ~ binomial(N0, p0);
    logit_pi0 ~ normal(0, 10);
    LOR_c ~ normal(0, 2);
  }"
  
  # if user does not specify control parameters
  # default set to smaller step size to improve divergence in some cases
  if ('control' %in% names(options)) {
    model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, Se = se, Sp = sp), 
                  pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh, seed = seed, ...)
  }
  else {
    model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, Se = se, Sp = sp), 
                  pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh, seed = seed,
                  control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 50), ...)
  }
  print(summary(model)$summary)
  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return(model)
}

#' Model with nondifferential, logit normal-distributed misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification under a logit-transformed scaled bivariate normal distribution.
#
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param m.lg.se normal distribution of logit Se with (mean = m.lg.se, sd = s.lg.se).
#' @param m.lg.sp normal distribution of logit Sp with (m.lg.sp, s.lg.sp).
#' @param s.lg.se standard deviation of logit Se
#' @param s.lg.sp standard deviation of logit Sp
#' @param lg.se used as an initial value for logit Se. Default to m.lg.se
#' @param lg.sp used as an initial value for logit Sp. Default to m.lg.sp
#' @param name a string of the name of the model. Default to "Logit Normal Misclassification Model".
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan} for more details.
#' @param ... optional parameters passed to \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#' 
#' \dontrun{
#' logitOR(a = 66, N1 = 11782, c = 243, N0 = 57973, m.lg.se = 1.069, m.lg.sp = 1.126,
#'   s.lg.se = 0.893, s.lg.sp = 0.712, chains = 3, iter = 10000, seed = 0)
#' }
#' logitOR(a = 66, N1 = 11782, c = 243, N0 = 57973, m.lg.se = 1.069, m.lg.sp = 1.126,
#'   s.lg.se = 0.893, s.lg.sp = 0.712, lg.se = 2.197, lg.sp = 2.197, traceplot = TRUE)

logitOR <- function(a, N1, c, N0, m.lg.se, m.lg.sp, s.lg.se, s.lg.sp, lg.se = NULL, lg.sp = NULL, 
                    name = "Logit Normal Misclassification Model", chains = 2,
                    traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = NA, ...) {
  
  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0))) {
    stop("The value(s) for a/N0/c/N1 is not valid.")
  }
  
  # default to mean if not specified
  if (is.null(lg.se)) {lg.se <- m.lg.se}
  if (is.null(lg.sp)) {lg.sp <- m.lg.sp}
  
  options <- list(...)
  
  code <- "
  data {
    int<lower=0> a;
    int<lower=0> N1;
    int<lower=0> c;
    int<lower=0> N0;
    real mX;
    real mY;
    real sdX;
    real sdY;
  }
  parameters {
    real logit_pi0;
    real LOR_c;
    real X;
    real Y;
  }
  transformed parameters {
    real<lower=0, upper=1> pi1;
    real<lower=0, upper=1> pi0;
    real ORadj;
    real p1;
    real p0;
    real Se;
    real Sp;
    Se=(1 + exp(X)/(1 + exp(X))) / 2;
    Sp=(1 + exp(Y)/(1 + exp(Y))) / 2;
    pi0 = exp(logit_pi0) / (exp(logit_pi0) + 1);
    pi1 = exp(logit_pi0 + LOR_c) / (exp(logit_pi0 + LOR_c) + 1);
    p1 = pi1 * Se + (1 - pi1) * (1 - Sp);
    p0 = pi0 * Se + (1 - pi0) * (1 - Sp);
    ORadj = exp(LOR_c);
  }
  model {
    X ~ normal(mX, sdX);
    Y ~ normal(mY, sdY);
    a ~ binomial(N1, p1);
    c ~ binomial(N0, p0);
    logit_pi0 ~ normal(0, 10);
    LOR_c ~ normal(0, 2);
  }
  "
  
  # if user does not specify control parameters
  # default set to smaller step size to improve divergence in some cases
  if ('control' %in% names(options)) {
    model <- stan(model_code = code, model_name = name, data=list(a = a, N1 = N1, c = c, N0 = N0, mX = m.lg.se, mY = m.lg.sp, 
                  sdX = s.lg.se, sdY = s.lg.sp), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh, 
                  init = rep(list(list(X = lg.se, Y = lg.sp)), chains), seed = seed, ...)
  }
  else {
    model <- stan(model_code = code, model_name = name, data=list(a = a, N1 = N1, c = c, N0 = N0, mX = m.lg.se, mY = m.lg.sp, 
                  sdX = s.lg.se, sdY = s.lg.sp), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh, 
                  init = rep(list(list(X = lg.se, Y = lg.sp)), chains), seed = seed,
                  control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 50), ...)
  }
  
  print(summary(model)$summary)
  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return(model)
}

#' Model with nondifferential, correlated misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification that extends from the logit model but allows there to be a fixed correlation between sentivity and specificity.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param m.lg.se normal distribution of logit Se with (mean = m.lg.se, sd = s.lg.se).
#' @param m.lg.sp conditional normal distribution of logit Sp given Se with (m.lg.sp, s.lg.sp).
#' @param s.lg.se standard deviation of logit Se
#' @param s.lg.sp standard deviation of logit Sp
#' @param lg.se used as an initial value for logit Se. Default to m.lg.se
#' @param lg.sp used as an initial value for logit Sp. Default to m.lg.sp
#' @param rho correlation between Se and Sp
#' @param name a string of the name of the model. Default to "Logit Model with Fixed Correlation".
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan} for more details.
#' @param ... optional parameters passed to \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#' \dontrun{
#' fixedCorrOR(a = 66, N1 = 11782, c = 243, N0 = 57973, m.lg.se = 1.069, m.lg.sp = 1.126,
#'   s.lg.se = 0.893, s.lg.sp = 0.712, rho = -0.379, chains = 3, iter = 10000, seed = 0)
#' }
#' fixedCorrOR(a = 66, N1 = 11782, c = 243, N0 = 57973, m.lg.se = 1.069, m.lg.sp = 1.126,
#'   s.lg.se = 0.893, s.lg.sp = 0.712, lg.se = 2.197, lg.sp = 0.744, rho = -0.379, 
#'   traceplot = TRUE)

fixedCorrOR <- function(a, N1, c, N0, m.lg.se, m.lg.sp, s.lg.se, s.lg.sp, lg.se = NULL, lg.sp = NULL, 
                        rho, name = "Logit Model with Fixed Correlation", chains = 2,
                        traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = NA, ...) {
  
  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0) & (-1 <= rho) & (rho <= 1))) {
    stop("The value(s) for a/N0/c/N1/rho is not valid.")
  }
  
  # default to mean if not specified
  if (is.null(lg.se)) {lg.se <- m.lg.se}
  if (is.null(lg.sp)) {lg.sp <- m.lg.sp}
  
  options <- list(...)
  
  code <- "
  data {
    int<lower=0> a;
    int<lower=0> N1;
    int<lower=0> c;
    int<lower=0> N0;
    real mX0;
    real precX0;
    real mX1;
    real precX1;
    real rhoSe;
  }
  parameters {
    real logit_pi0;
    real LOR_c;
    real X0;
    real X1;
  }
  transformed parameters {
    real<lower=0, upper=1> pi1;
    real<lower=0, upper=1> pi0;
    real ORadj;
    real<lower=0, upper=1> Se;
    real<lower=0, upper=1> Sp;
    real mcx1;
    real preccx1;
    real p1;
    real p0;
    mcx1 = mX1 + rhoSe * (X0 - mX0) * ((precX1 / precX0)^0.5);
    preccx1 = precX1 / (1 - rhoSe^2);
    pi0 = exp(logit_pi0) / (exp(logit_pi0) + 1);
    pi1 = exp(logit_pi0 + LOR_c) / (exp(logit_pi0 + LOR_c) + 1);
    Se = (1 + exp(X0) / (1 + exp(X0))) / 2;
    Sp = (1 + exp(X1) / (1 + exp(X1))) / 2;
    p1 = pi1 * Se + (1 - pi1) * (1 - Sp);
    p0 = pi0 * Se + (1 - pi0) * (1 - Sp);
    ORadj = exp(LOR_c);
  }
  model {
    a ~ binomial(N1, p1);
    c ~ binomial(N0, p0);
    X0 ~ normal(mX0, (1 / precX0)^0.5);
    X1 ~ normal(mcx1, (1 / preccx1)^0.5);
    logit_pi0~normal(0, 10);
    LOR_c~normal(0, 2);
  }
  "
  
  # if user does not specify control parameters
  # default set to smaller step size to improve divergence in some cases
  if ('control' %in% names(options)) {
    model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, mX0 = m.lg.se, mX1 = m.lg.sp, 
                  precX0 = 1 / (s.lg.se)^2, precX1 = 1 / (s.lg.sp)^2, rhoSe = rho), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh, 
                  init = rep(list(list(X0 = lg.se, X1 = lg.sp)), chains), seed = seed, ...)
  }
  else {
    model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, mX0 = m.lg.se, mX1 = m.lg.sp, 
                  precX0 = 1 / (s.lg.se)^2, precX1 = 1 / (s.lg.sp)^2, rhoSe = rho), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh, 
                  init = rep(list(list(X0 = lg.se, X1 = lg.sp)), chains), seed = seed,
                  control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 50), ...)
  }
  
  print(summary(model)$summary)
  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return(model)
}

#' Model with nondifferential, randomly correlated misclassification

#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification that extends from the logit model but allows a random correlation between Sensitivity and Specificity.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param m.lg.se normal distribution of logit Se with (mean = m.lg.se, sd = s.lg.se).
#' @param m.lg.sp conditional normal distribution of logit Sp given Se with (m.lg.sp, s.lg.sp).
#' @param s.lg.se standard deviation of logit Se
#' @param s.lg.sp standard deviation of logit Sp
#' @param lg.se used as an initial value for logit Se. Default to m.lg.se
#' @param lg.sp used as an initial value for logit Sp. Default to m.lg.sp
#' @param z used as an initial value of Fisher's Z transformed of rho, where correlation rho = (exp(2*z)-1)/(1+exp(2*z))).
#' @param m.z normal distribution of Z with (mean = m.z, sd = s.z).
#' @param s.z normal distribution of Z with (mean = m.z, sd = s.z).
#' @param name a string of the name of the model. Default to "Logit Model with Random Correlation".
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan} for more details.
#' @param ... optional parameters passed to \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#' \dontrun{
#' randCorrOR(a = 66, N1 = 11782, c = 243, N0 = 57973, m.lg.se = 1.069, m.lg.sp = 1.126,
#'   s.lg.se = 0.893, s.lg.sp = 0.712, m.z = -0.399, s.z = 0.139, chains = 3, 
#'   iter = 10000, seed = 0)
#' }
#' randCorrOR(a = 66, N1 = 11782, c = 243, N0 = 57973, m.lg.se = 1.069, m.lg.sp = 1.126,
#'   s.lg.se = 0.893, s.lg.sp = 0.712, lg.se = 2.197, lg.sp = 0.744, m.z = -0.399, 
#'   s.z = 0.139, traceplot = TRUE)

randCorrOR <- function(a, N1, c, N0, m.lg.se, m.lg.sp, s.lg.se, s.lg.sp, lg.se = NULL, lg.sp = NULL, m.z, s.z, z = NULL, 
                       name = "Logit Model with Random Correlation", chains = 2,
                       traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = NA, ...) {
  
  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0))) {
    stop("The value(s) for a/N0/c/N1 is not valid.")
  }
  
  # default to mean if not specified
  if (is.null(lg.se)) {lg.se <- m.lg.se}
  if (is.null(lg.sp)) {lg.sp <- m.lg.sp}
  if (is.null(z)) {z <- m.z}

  options <- list(...)
  
  code <- "
  data {
    int<lower=0> a;
    int<lower=0> N1;
    int<lower=0> c;
    int<lower=0> N0;
    real mX0;
    real precX0;
    real mX1;
    real precX1;
    real mZ;
    real sZ;
  }
  parameters {
    real logit_pi0;
    real LOR_c;
    real Z;
    real X0;
    real X1;
  }
  transformed parameters {
    real<lower=0, upper=1> pi1;
    real<lower=0, upper=1> pi0;
    real ORadj;
    real<lower=0, upper=1> Se;
    real<lower=0, upper=1> Sp;
    real rhoSe;
    real mcx1;
    real preccx1;
    real p1;
    real p0;
    rhoSe = (exp(2 * Z) - 1)/(1 + exp(2 * Z));
    mcx1 = mX1 + rhoSe * (X0 - mX0) * pow(precX1 / precX0, 0.5);
    preccx1 = precX1 / (1 - pow(rhoSe, 2));
    Se = (1 + exp(X0) / (1 + exp(X0))) / 2;
    Sp = (1 + exp(X1) / (1 + exp(X1))) / 2;
    pi0 = exp(logit_pi0) / (exp(logit_pi0) + 1);
    pi1 = exp(logit_pi0 + LOR_c) / (exp(logit_pi0 + LOR_c) + 1);
    p1 = pi1 * Se + (1 - pi1) * (1 - Sp);
    p0 = pi0 * Se + (1 - pi0) * (1 - Sp);
    ORadj = exp(LOR_c);
  }
  model {
    a ~ binomial(N1, p1);
    c ~ binomial(N0, p0);
    Z ~ normal(mZ, sZ);
    X0 ~ normal(mX0, (1/precX0)^0.5);
    X1 ~ normal(mcx1, (1/preccx1)^0.5);
    logit_pi0 ~ normal(0, 10);
    LOR_c ~ normal(0, 2);
  }
  "
  # if user does not specify control parameters
  # default set to smaller step size to improve divergence in some cases
  
  if ('control' %in% names(options)) {
    model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, mX0 = m.lg.se, mX1 = m.lg.sp, 
                  precX0 = 1 / (s.lg.se)^2, precX1 = 1 / (s.lg.sp)^2, mZ = m.z, sZ = s.z), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh,
                  init = rep(list(list(X0 = lg.se, X1 = lg.sp, Z = z)), chains), seed = seed, ...)
  }
  else {
    model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, mX0 = m.lg.se, mX1 = m.lg.sp, 
                  precX0 = 1 / (s.lg.se)^2, precX1 = 1 / (s.lg.sp)^2, mZ = m.z, sZ = s.z), pars = c("LOR_c", "ORadj"), chains = chains, refresh = refresh,
                  init = rep(list(list(X0 = lg.se, X1 = lg.sp, Z = z)), chains), seed = seed,
                  control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 50), ...)
  }
  
  print(summary(model)$summary)
  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return (model)
}

#' Model with differential misclassification 
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given a four-variate differential misclassification.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param mu vector of length 4; multivariate normal distribution of \eqn{z \sim (mu, varz)}, where each \eqn{\mu} corresponds to the logit mean of \eqn{Se_0$}, \eqn{Se_1}, \eqn{Sp_0} and \eqn{Sp_1} (0 for controls, 1 for cases group).
#' @param s.lg.se0 standard deviation of logit Se in the control group.
#' @param s.lg.se1 standard deviation of logit Se in the case group.
#' @param s.lg.sp0 standard deviation of logit Sp in the control group.
#' @param s.lg.sp1 standard deviation of logit Sp in the case group.
#' @param corr.sesp0 correlation between Se_0 and Sp_0.
#' @param corr.sesp1 correlation between Se_1 and Sp_1.
#' @param corr.group correlation between Se_0 and Se_1, Sp_0 and Sp_1. Default to 0. 
#' @param z vector of length 4; used as an initial value for \eqn{z \sim (mu, varz)}. Default to mu. 
#' @param name a string of the name of the model. Default to "Model with differential misclassification".
#' @param chains number of Markov Chains. Default to 2.
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param refresh an integer value used to control how often the progress of sampling is reported. By default, the progress indicator is turned off, thus refresh <= 0.
#' If on, refresh = max(iter/10, 1) is generally recommended.
#' @param seed the seed for random number generation. See \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan} for more details.
#' @param ... optional parameters passed to \href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' # Case-control study data of Bipolar Disorder with rheumatoid arthritis (Farhi et al. 2016)
#' # Data from \url{https://www.sciencedirect.com/science/article/pii/S0165032715303864#bib13}
#' \dontrun{
#' nonDiffOR(a = 66, N1 = 11782, c = 243, N0 = 57973, chains = 3, mu = c(1.069, 1.069, 1.126, 1.126),
#'   s.lg.se0 = 0.893, s.lg.se1 = 0.893, s.lg.sp0 = 0.712, s.lg.sp1 = 0.712, corr.sesp0 = -0.377, 
#'   corr.sesp1 = -0.377, corr.group = 0, iter = 10000, seed = 0)
#' }
#' nonDiffOR(a = 66, N1 = 11782, c = 243, N0 = 57973, , mu = c(1.069, 1.069, 1.126, 1.126),
#'   s.lg.se0 = 0.893, s.lg.se1 = 0.893, s.lg.sp0 = 0.712, s.lg.sp1 = 0.712, corr.sesp0 = -0.377, 
#'   corr.sesp1 = -0.377, corr.group = 0, traceplot = TRUE)


nonDiffOR <- function(a, N1, c, N0, mu, s.lg.se0, s.lg.se1, s.lg.sp0, s.lg.sp1, corr.sesp0, corr.sesp1,
                      corr.group = 0, z = NULL, name = "Model with differential classification", chains = 2,
                      traceplot = FALSE, inc_warmup = FALSE, window = NULL, refresh = 0, seed = 0, ...) {
  
  if (!((a <= N1) & (a >= 0) & (c <= N0) & (c >= 0))) {
    stop("The value(s) for a/N0/c/N1 is not valid.")
  }
  
  # default to mean if not specified
  if (is.null(z)) {z <- mu}
  
  options <- list(...)
  
  code <- "
  data {
    int<lower=0> a;
    int<lower=0> N1;
    int<lower=0> c;
    int<lower=0> N0;
    vector[4] Mu;
    matrix[4,4] varZ;
  }
  parameters {
    real logit_pi0;
    real LOR_c;
    vector[4] Z;
  }
  transformed parameters {
    real<lower=0, upper=1> pi1;
    real<lower=0, upper=1> pi0;
    real ORadj;
    real Se0;
    real Se1;
    real Sp0;
    real Sp1;
    real<lower=0, upper=1> p1;
    real<lower=0, upper=1> p0;
    Se0 = (1 + exp(Z[1]) / (1 + exp(Z[1]))) / 2;
    Se1 = (1 + exp(Z[2]) / (1 + exp(Z[2]))) / 2;
    Sp0 = (1 + exp(Z[3]) / (1 + exp(Z[3]))) / 2;
    Sp1 = (1 + exp(Z[4]) / (1 + exp(Z[4]))) / 2;
    pi0 = exp(logit_pi0) / (exp(logit_pi0) + 1);
    pi1 = exp(logit_pi0 + LOR_c) / (exp(logit_pi0 + LOR_c) + 1);
    p1 = pi1 * Se1 + (1 - pi1) * (1 - Sp1);
    p0 = pi0 * Se0 + (1 - pi0) * (1 - Sp0);
    ORadj = exp(LOR_c);
  }
  model {
    a ~ binomial(N1, p1);
    c ~ binomial(N0, p0);
    logit_pi0 ~ normal(0, 10);
    LOR_c ~ normal(0, 2);
    Z ~ multi_normal(Mu, varZ);
  }
  "
  
  varz = matrix(c(s.lg.se0^2, corr.group*s.lg.se0*s.lg.se1, corr.sesp0*s.lg.se0*s.lg.sp0, corr.group*s.lg.se0*s.lg.sp1,
                  corr.group*s.lg.se0*s.lg.se1, s.lg.se1^2, corr.group*s.lg.se1*s.lg.sp0, corr.sesp1*s.lg.se1*s.lg.sp1,
                  corr.sesp0*s.lg.se0*s.lg.sp0, corr.group*s.lg.se1*s.lg.sp0, s.lg.sp0^2, corr.group*s.lg.sp0*s.lg.sp1,
                  corr.group*s.lg.se0*s.lg.sp1, corr.sesp1*s.lg.se1*s.lg.sp1, corr.group*s.lg.sp0*s.lg.sp1, s.lg.sp1^2), 4, 4)
  
  # if user does not specify control parameters
  # default set to smaller step size to improve divergence in some cases
  if ('control' %in% names(options)) {
    model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, Mu = mu, 
      varZ = varz), pars = c("LOR_c", "ORadj"), chains = chains, init = rep(list(list(Z = z)), chains), refresh = refresh, seed = seed, ...)
  }
  else {
    model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, Mu = mu, 
                  varZ = varz), pars = c("LOR_c", "ORadj"), chains = chains, init = rep(list(list(Z = z)), chains), refresh = refresh, seed = seed,
                  control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 50), ...)
  }
  print(summary(model)$summary)
  
  if (traceplot) {
    print(traceplot(model, pars = "LOR_c", inc_warmup = inc_warmup, window = window) + xlab("Iterations shown"))
  }
  return(model)
}