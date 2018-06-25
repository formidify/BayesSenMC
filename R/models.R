#' Model without misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of uncorrected odds ratio given no misclassification.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param name a string of the name of the model. Default to "Corrected Model".
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot}} corresponding to one or more Markov chains. 
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param ... optional parameters passed to \code{\href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here}} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' correctedOR(a = 126, N1 = 218, c = 71, N0 = 295, chains = 3, iter = 10000)
#'
#' correctedOR(a = 126, N1 = 218, c = 71, N0 = 295, traceplot = TRUE, window = c(1, 1000))

correctedOR <- function(a, N1, c, N0, name = "Corrected Model", traceplot = FALSE, inc_warmup = FALSE, window = NULL, ...) {
  code <- "
  data {
    int<lower=0> a;
    int<lower=0> N1;
    int<lower=0> c;
    int<lower=0> N0;
  }
  parameters {
    real alpha0;
    real alpha1;
  }
  transformed parameters {
    real<lower=0, upper=1> pi1;
    real<lower=0, upper=1> pi0;
    real ORadj;
    pi0 = exp(alpha0) / (exp(alpha0) + 1);
    pi1 = exp(alpha0 + alpha1) / (exp(alpha0 + alpha1) + 1);
    ORadj = exp(alpha1);
  }
  model {
    a ~ binomial(N1, pi1);
    c ~ binomial(N0, pi0);
    alpha0 ~ normal(0, 10);
    alpha1 ~ normal(0, 2);
  }"
  model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0), pars = c("ORadj"), ...)
  print(summary(model))
  if (traceplot) {
    print(traceplot(model, inc_warmup = inc_warmup, window = window))
  }
  return(model)
}

#' Model with constant nondifferential classification 
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification with constant Se and Sp.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param se sensitivity
#' @param sp specificity
#' @param name a string of the name of the model. Default to "Crude Model".
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot}} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param ... optional parameters passed to \code{\href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here}} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' crudeOR(a = 126, N1 = 218, c = 71, N0 = 295, se = 0.94, sp = 0.97, chains = 3, iter = 10000)
#'
#' crudeOR(a = 126, N1 = 218, c = 71, N0 = 295, se = 0.94, sp = 0.97, traceplot = TRUE, window = c(1, 1000))

crudeOR <- function(a, N1, c, N0, se, sp, name = "Crude Model", traceplot = FALSE, inc_warmup = FALSE, window = NULL, ...) {
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
    real alpha0;
    real alpha1;
  }
  transformed parameters {
    real<lower=0, upper=1> pi1;
    real<lower=0, upper=1> pi0;
    real ORadj;
    real<lower=0, upper=1> p1s;
    real<lower=0, upper=1> p0s;
    pi0 = exp(alpha0) / (exp(alpha0) + 1);
    pi1 = exp(alpha0 + alpha1) / (exp(alpha0 + alpha1) + 1);
    p1s = pi1 * Se + (1 - pi1) * (1 - Sp);
    p0s = pi0 * Se + (1 - pi0) * (1 - Sp);
    ORadj = exp(alpha1);
  }
  model {
    a ~ binomial(N1, p1s);
    c ~ binomial(N0, p0s);
    alpha0 ~ normal(0, 10);
    alpha1 ~ normal(0, 2);
  }"
  model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, Se = se, Sp = sp), pars = c("ORadj"), ...)
  print(summary(model))
  if (traceplot) {
    print(traceplot(model, inc_warmup = inc_warmup, window = window))
  }
  return(model)
}

#' Model with nondifferential distributed misclassification (logit model)
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification under a logit-transformed scaled bivariate normal distribution.
#
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param mx normal distribution of X with (mX, varX).
#' @param my normal distribution of Y with (mY, varY).
#' @param varx normal distribution of X with (mX, varX).
#' @param vary normal distribution of Y with (mY, varY).
#' @param x used as an initial value of X.
#' @param y used as an initial value of Y.
#' @param name a string of the name of the model. Default to "Logit Model".
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot}} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param ... optional parameters passed to \code{\href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here}} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' logitOR(a = 126, N1 = 218, c = 71, N0 = 295, mx = 1.9814, my = 2.8586,
#'   varx = 0.9872, vary = 0.5806, x = 2.76, y = 3.58, chains = 3, iter = 10000)
#'
#' logitOR(a = 126, N1 = 218, c = 71, N0 = 295, mx = 1.9814, my = 2.8586,
#'   varx = 0.9872, vary = 0.5806, x = 2.76, y = 3.58, traceplot = TRUE, window = c(1, 1000))

logitOR <- function(a, N1, c, N0, mx, my, varx, vary, x, y, name = "Logit Model", 
                    traceplot = FALSE, inc_warmup = FALSE, window = NULL, ...) {
  code <- "
  data {
    int<lower=0> a;
    int<lower=0> N1;
    int<lower=0> c;
    int<lower=0> N0;
    real mX;
    real mY;
    real varX;
    real varY;
  }
  parameters {
    real alpha0;
    real alpha1;
    real X;
    real Y;
  }
  transformed parameters {
    real<lower=0, upper=1> pi1;
    real<lower=0, upper=1> pi0;
    real ORadj;
    real p1s;
    real p0s;
    real Se;
    real Sp;
    Se=(1 + exp(X)/(1 + exp(X))) / 2;
    Sp=(1 + exp(Y)/(1 + exp(Y))) / 2;
    pi0 = exp(alpha0) / (exp(alpha0) + 1);
    pi1 = exp(alpha0 + alpha1) / (exp(alpha0 + alpha1) + 1);
    p1s = pi1 * Se + (1 - pi1) * (1 - Sp);
    p0s = pi0 * Se + (1 - pi0) * (1 - Sp);
    ORadj = exp(alpha1);
  }
  model {
    X ~ normal(mX, varX^0.5);
    Y ~ normal(mY, varY^0.5);
    a ~ binomial(N1, p1s);
    c ~ binomial(N0, p0s);
    alpha0 ~ normal(0, 10);
    alpha1 ~ normal(0, 2);
  }
  "
  model <- stan(model_code = code, model_name = name, data=list(a = a, N1 = N1, c = c, N0 = N0, mX = mx, mY = my, 
    varX = varx, varY = vary), pars = c("ORadj"), ...)
  print(summary(model))
  if (traceplot) {
    print(traceplot(model, inc_warmup = inc_warmup, window = window))
  }
  return(model)
}

#' Model with nondifferential, correlated misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given nondifferential misclassification that extends from the logit model but allows there to be a fixed correlation between Sentivity and Specificity.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param mx0 normal distribution of X0 with (mx0, varx0)
#' @param mx1 conditional distribution of X1 with (mx1, varx1)
#' @param varx0 normal distribution of X0 with (mx0, varx0)
#' @param varx1 conditional normal distribution of X1 with (mx1, varx1)
#' @param x0 used as an initial value of X0
#' @param x1 used as an initial value of X1
#' @param rhose correlation between Se and Sp
#' @param name a string of the name of the model. Default to "Model with fixed correlation".
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot}} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param ... optional parameters passed to \code{\href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here}} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' fixedCorrOR(a = 126, N1 = 218, c = 71, N0 = 295, mx0 = 1.9814, mx1 = 2.8586,
#'   varx0 = 0.9872, varx1 = 0.5806, x0 = 1.9814, x1 = 2.8586, rhose = -0.6271, chains = 3, iter = 10000)
#'
#' fixedCorrOR(a = 126, N1 = 218, c = 71, N0 = 295, mx0 = 1.9814, mx1 = 2.8586,
#'   varx0 = 0.9872, varx1 = 0.5806, x0 = 1.9814, x1 = 2.8586, rhose = -0.6271, traceplot = TRUE, window = c(1, 1000))

fixedCorrOR <- function(a, N1, c, N0, mx0, mx1, varx0, varx1, x0, x1, rhose, name = "Model with fixed correlation", 
                        traceplot = FALSE, inc_warmup = FALSE, window = NULL, ...) {
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
    real alpha0;
    real alpha1;
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
    real p1s;
    real p0s;
    mcx1 = mX1 + rhoSe * (X0 - mX0) * ((precX1 / precX0)^0.5);
    preccx1 = precX1 / (1 - rhoSe^2);
    pi0 = exp(alpha0) / (exp(alpha0) + 1);
    pi1 = exp(alpha0 + alpha1) / (exp(alpha0 + alpha1) + 1);
    Se = (1 + exp(X0) / (1 + exp(X0))) / 2;
    Sp = (1 + exp(X1) / (1 + exp(X1))) / 2;
    p1s = pi1 * Se + (1 - pi1) * (1 - Sp);
    p0s = pi0 * Se + (1 - pi0) * (1 - Sp);
    ORadj = exp(alpha1);
  }
  model {
    a ~ binomial(N1, p1s);
    c ~ binomial(N0, p0s);
    X0 ~ normal(mX0, (1 / precX0)^0.5);
    X1 ~ normal(mcx1, (1 / preccx1)^0.5);
    alpha0~normal(0, 10);
    alpha1~normal(0, 2);
  }
  "
  model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, mX0 = mx0, mX1 = mx1, 
    precX0 = 1 / varx0, precX1 = 1 / varx1, rhoSe = rhose), pars = c("ORadj"), ...)
  print(summary(model))
  if (traceplot) {
    print(traceplot(model, inc_warmup = inc_warmup, window = window))
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
#' @param mx0 normal distribution of X0 with (mx0, varx0).
#' @param mx1 conditional distribution of X1 with (mx1, varx1).
#' @param varx0 normal distribution of X0 with (mx0, varx0).
#' @param varx1 conditional normal distribution of X1 with (mx1, varx1).
#' @param x0 used as an initial value of X0.
#' @param x1 used as an initial value of X1.
#' @param z used as an initial value of Z, where correlation (rhoSe) = (exp(2*z)-1)/(1+exp(2*z))).
#' @param mz normal distribution of Z with (mz, varz).
#' @param varz normal distribution of Z with (mz, varz).
#' @param name a string of the name of the model. Default to "Model with random correlation".
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot}} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param ... optional parameters passed to \code{\href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here}} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' randCorrOR(a = 126, N1 = 218, c = 71, N0 = 295, mx0 = 1.9814, mx1 = 2.8586, varx0 = 0.9872,
#'   varx1 = 0.5806, x0 = 1.9814, x1 = 2.8586, z = -0.7366, mz = -0.725, varz = 0.25, chains = 3, iter = 10000)
#'
#' randCorrOR(a = 126, N1 = 218, c = 71, N0 = 295, mx0 = 1.9814, mx1 = 2.8586, varx0 = 0.9872,
#'   varx1 = 0.5806, x0 = 1.9814, x1 = 2.8586, z = -0.7366, mz = -0.725, varz = 0.25, traceplot = TRUE, window = c(1, 1000))

randCorrOR <- function(a, N1, c, N0, mx0, mx1, varx0, varx1, x0, x1, z, mz, varz, name = "Model with random correlation", 
                        traceplot = FALSE, inc_warmup = FALSE, window = NULL, ...) {
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
    real varZ;
  }
  parameters {
    real alpha0;
    real alpha1;
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
    real p1s;
    real p0s;
    rhoSe = (exp(2 * Z) - 1)/(1 + exp(2 * Z));
    mcx1 = mX1 + rhoSe * (X0 - mX0) * pow(precX1 / precX0, 0.5);
    preccx1 = precX1 / (1 - pow(rhoSe, 2));
    Se = (1 + exp(X0) / (1 + exp(X0))) / 2;
    Sp = (1 + exp(X1) / (1 + exp(X1))) / 2;
    pi0 = exp(alpha0) / (exp(alpha0) + 1);
    pi1 = exp(alpha0 + alpha1) / (exp(alpha0 + alpha1) + 1);
    p1s = pi1 * Se + (1 - pi1) * (1 - Sp);
    p0s = pi0 * Se + (1 - pi0) * (1 - Sp);
    ORadj = exp(alpha1);
  }
  model {
    a ~ binomial(N1, p1s);
    c ~ binomial(N0, p0s);
    Z ~ normal(mZ, varZ^0.5);
    X0 ~ normal(mX0, (1/precX0)^0.5);
    X1 ~ normal(mcx1, (1/preccx1)^0.5);
    alpha0 ~ normal(0, 10);
    alpha1 ~ normal(0, 2);
  }
  "
  model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, mX0 = mx0, mX1 = mx1, 
    precX0 = 1 / varx0, precX1 = 1 / varx1, mZ = mz, varZ = varz), pars = c("ORadj"), ...)
  print(summary(model))
  if (traceplot) {
    print(traceplot(model, inc_warmup = inc_warmup, window = window))
  }
  return (model)
}

#' Model with differential misclassification 
#'
#' Generate a stanfit object corresponding to a posterior distribution of corrected odds ratio given a four-variate differential misclassification that extends from the logit model.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param mx0 normal distribution of X0 with (mx0, varx0).
#' @param mx1 conditional distribution of X1 with (mx1, varx1).
#' @param mu vector of length 4; multivariate normal distribution of z ~ (mu, varz).
#' @param varz covariance matrix of length 4*4, passed in as an one dimentional array of length 16; multivariate normal distribution of z ~ (mu, varz).
#' @param z vector of length 4; used as an initial value for Z.
#' @param name a string of the name of the model. Default to "Model with differential misclassification".
#' @param traceplot Logical, defaulting to \code{FALSE}. If \code{TRUE} it will draw the 
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-method-traceplot.html}{traceplot}} corresponding to one or more Markov chains.
#' @param inc_warmup Only evaluated when \code{traceplot = TRUE}. \code{TRUE} or \code{FALSE}, indicating whether or not to include the warmup sample in the
#' traceplot; defaults to \code{FALSE}.
#' @param window Only evaluated when \code{traceplot = TRUE}. A vector of length 2. Iterations between \code{window[1]} and \code{window[2]} will be shown in the plot. 
#' The default shows all iterations if \code{inc_warmup} is \code{TRUE} and all iterations from the sampling period only if \code{inc_warmup} is \code{FALSE}. 
#' If \code{inc_warmup} is \code{FALSE} the iterations specified in \code{window} do not include iterations from the warmup period.
#' The default number of iterations is 2000 unless otherwise specified in the optional \code{iter} argument.
#' @param ... optional parameters passed to \code{\href{https://www.rdocumentation.org/packages/rstan/versions/2.17.3/topics/stan}{stan}}.
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' \code{\href{https://mc-stan.org/rstan/reference/stanfit-class.html}{here}} for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' nonDiffOR(a = 126, N1 = 218, c = 71, N0 = 295, mu = c(2.76, 2.76, 3.58, 3.58),
#'   varz = c(0.9100002, 0.8190002, -0.5053286, -0.5053286, 0.8190002, 0.9100002,
#'   -0.5053286, -0.5053286, -0.5053286, -0.5053286, 0.7300000, 0.6570000,
#'   -0.5053286, -0.5053286, 0.6570000, 0.7300000), z = c(2.76, 2.76, 3.58, 3.58), chains = 3, iter = 10000)
#'
#' nonDiffOR(a = 126, N1 = 218, c = 71, N0 = 295, mu = c(2.76, 2.76, 3.58, 3.58),
#'   varz = c(0.9100002, 0.8190002, -0.5053286, -0.5053286, 0.8190002, 0.9100002,
#'   -0.5053286, -0.5053286, -0.5053286, -0.5053286, 0.7300000, 0.6570000,
#'   -0.5053286, -0.5053286, 0.6570000, 0.7300000), z = c(2.76, 2.76, 3.58, 3.58), traceplot = TRUE, window = c(1, 1000))

nonDiffOR <- function(a, N1, c, N0, mu, varz, z, name = "Model with differential classification", 
                      traceplot = FALSE, inc_warmup = FALSE, window = NULL, ...) {
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
    real alpha0;
    real alpha1;
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
    real<lower=0, upper=1> p1s;
    real<lower=0, upper=1> p0s;
    Se0 = (1 + exp(Z[1]) / (1 + exp(Z[1]))) / 2;
    Se1 = (1 + exp(Z[2]) / (1 + exp(Z[2]))) / 2;
    Sp0 = (1 + exp(Z[3]) / (1 + exp(Z[3]))) / 2;
    Sp1 = (1 + exp(Z[4]) / (1 + exp(Z[4]))) / 2;
    pi0 = exp(alpha0) / (exp(alpha0) + 1);
    pi1 = exp(alpha0 + alpha1) / (exp(alpha0 + alpha1) + 1);
    p1s = pi1 * Se1 + (1 - pi1) * (1 - Sp1);
    p0s = pi0 * Se0 + (1 - pi0) * (1 - Sp0);
    ORadj = exp(alpha1);
  }
  model {
    a ~ binomial(N1, p1s);
    c ~ binomial(N0, p0s);
    alpha0~normal(0, 10);
    alpha1~normal(0, 2);
    Z ~ multi_normal(Mu, varZ);
  }
  "
  model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, Mu = mu, 
    varZ = structure(.Data = varz, .Dim = c(4, 4))), pars = c("ORadj"), ...)
  print(summary(model))
  if (traceplot) {
    print(traceplot(model, inc_warmup = inc_warmup, window = window))
  }
  return(model)
}