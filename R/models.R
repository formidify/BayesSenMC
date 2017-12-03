#' Model with no misclassification
#'
#' This function generates a posterior distribution of uncorrected odds ratio given no misclassification.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param name a string of the name of the model. Default to "Corrected Model".
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See
#' http://mc-stan.org/rstan/reference/stanfit-class.html for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' correctedOR(a = 126, N1 = 218, c = 71, N0 = 295)

correctedOR <- function(a, N1, c, N0, name = "Corrected Model"){
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
  model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0), init = list(list(alpha0 = 0, alpha1 = 0), list(alpha0 = 2, alpha1 = 2)), pars = c("ORadj"), chains = 2, warmup = 10000, iter = 500000, verbose = FALSE)
  print(summary(model))
  return(model)
}

#' Model with nondifferential classification with constant Sensitivity and Specificity
#'
#' This function generates a posterior distribution of corrected odds ratio given nondifferential misclassification with constant Se and Sp.
#' @param a # of exposed subjects in the case group.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param se sensitivity
#' @param sp specificity
#' @param name a string of the name of the model. Default to "Crude Model".
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See http://mc-stan.org/rstan/reference/stanfit-class.html for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' crudeOR(a = 126, N1 = 218, c = 71, N0 = 295, se = 0.94, sp = 0.97)

crudeOR <- function(a, N1, c, N0, se, sp, name = "Crude Model"){
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
  model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, Se = se, Sp = sp), init = list(list(alpha0 = 0, alpha1 = 0), list(alpha0 = 2, alpha1 = 2)), pars = c("ORadj"), chains = 2, warmup = 10000, iter = 500000, verbose = FALSE)
  print(summary(model))
  return(model)
}

#' Model with nondifferential misclassification with Sensitivity and Specificity each under a logit-transformed scaled bivariate normal distribution
#'
#' This function generates a posterior distribution of corrected odds ratio given nondifferential misclassification with Se=0.5+0.5*expit(X) and Sp=0.5+0.5*expit(Y), where X and Y are independently distributed.
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
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See http://mc-stan.org/rstan/reference/stanfit-class.html for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' logitOR(a = 126, N1 = 218, c = 71, N0 = 295, mx = 1.9814, my = 2.8586,
#'   varx = 0.9872, vary = 0.5806, x = 2.76, y = 3.58)

logitOR <- function(a, N1, c, N0, mx, my, varx, vary, x, y, name = "Logit Model"){
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
  model <- stan(model_code = code, model_name = name, data=list(a = a, N1 = N1, c = c, N0 = N0, mX = mx, mY = my, varX = varx, varY = vary), init = list(list(X = x, Y = y, alpha0 = 2, alpha1 = 2), list(X = x, Y = y, alpha0 = 0, alpha1 = 0)), pars = c("ORadj"), chains = 2, warmup = 10000, iter = 500000, verbose = FALSE)
  print(summary(model))
  return(model)
}

#' Model that extends from logitModel but allows there to be a fixed correlation between Sentivity and Specificity
#'
#' This function generates a posterior distribution of corrected odds ratio given nondifferential misclassification with Se=0.5+0.5*expit(X) and Sp=0.5+0.5*expit(Y), where X and Y are independently distributed, and Se and Sp are correlated
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
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See http://mc-stan.org/rstan/reference/stanfit-class.html for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' fixedCorrOR(a = 126, N1 = 218, c = 71, N0 = 295, mx0 = 1.9814, mx1 = 2.8586,
#'   varx0 = 0.9872, varx1 = 0.5806, x0 = 1.9814, x1 = 2.8586, rhose = -0.6271)

fixedCorrOR <- function(a, N1, c, N0, mx0, mx1, varx0, varx1, x0, x1, rhose, name = "Model with fixed correlation"){
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
  model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, mX0 = mx0, mX1 = mx1, precX0 = 1 / varx0, precX1 = 1 / varx1, rhoSe = rhose), init = list(list(X0 = x0, X1 = x1, alpha0 = 0, alpha1 = 0), list(X0 = x0, X1 = x1, alpha0 = 2, alpha1 = 2)), pars = c("ORadj"), chains = 2, warmup = 10000, iter = 500000, verbose = FALSE)
  print(summary(model))
  return(model)
}

#' Model that extends the logit model but allows a random correlation between Sensitivity and Specificity, i.e., a Fisher's Z-transformed normal distribution is allowed for the correlation
#'
#' This function generates a posterior distribution of corrected odds ratio given nondifferential misclassification with Se=0.5+0.5*expit(X) and Sp=0.5+0.5*expit(Y), where X and Y are independently distributed, and Se and Sp are correlated and the correlation is a random variable, which is equal to 0.5*[In(1+Z)-In(1-Z)] with a normally distributed Z.
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
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See http://mc-stan.org/rstan/reference/stanfit-class.html for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' randCorrOR(a = 126, N1 = 218, c = 71, N0 = 295, mx0 = 1.9814, mx1 = 2.8586, varx0 = 0.9872,
#'   varx1 = 0.5806, x0 = 1.9814, x1 = 2.8586, z = -0.7366, mz = -0.725, varz = 0.25)

randCorrOR <- function(a, N1, c, N0, mx0, mx1, varx0, varx1, x0, x1, z, mz, varz, name = "Model with random correlation"){
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
  model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, mX0 = mx0, mX1 = mx1, precX0 = 1 / varx0, precX1 = 1 / varx1, mZ = mz, varZ = varz), init = list(list(X0 = x0, X1 = x1, alpha0 = 0, alpha1 = 0, Z = z), list(X0 = x0, X1 = x1, alpha0 = 2, alpha1 = 2, Z = z)), pars = c("ORadj"), chains = 2, warmup = 10000, iter = 500000, verbose = FALSE)
  print(summary(model))
  return (model)
}

#' Model with differential misclassification with Se0, Se1, Sp0, and Sp1 to be four-variate logit-transformed rescaled normal distributions
#'
#' This function generates a posterior distribution of corrected odds ratio given nondifferential misclassification with Se=0.5+0.5*expit(X) and Sp=0.5+0.5*expit(Y), where X and Y are independently distributed, and Se and Sp are correlated and the correlation is a random variable, which is equal to 0.5*[In(1+Z)-In(1-Z)] with a normally distributed Z.
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
#' @return It returns a stanfit object of this model, which inherits stanfit class methods. See http://mc-stan.org/rstan/reference/stanfit-class.html for more details.
#' @import rstan
#' @import Rcpp
#' @export
#' @examples
#' nonDiffOR(a = 126, N1 = 218, c = 71, N0 = 295, mu = c(2.76, 2.76, 3.58, 3.58),
#'   varz = c(0.9100002, 0.8190002, -0.5053286, -0.5053286, 0.8190002, 0.9100002,
#'   -0.5053286, -0.5053286, -0.5053286, -0.5053286, 0.7300000, 0.6570000,
#'   -0.5053286, -0.5053286, 0.6570000, 0.7300000), z = c(2.76, 2.76, 3.58, 3.58))

nonDiffOR <- function(a, N1, c, N0, mu, varz, z, name = "Model with differential classification"){
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
  model <- stan(model_code = code, model_name = name, data = list(a = a, N1 = N1, c = c, N0 = N0, Mu = mu, varZ = structure(.Data = varz, .Dim = c(4, 4))), init = list(list(Z = z, alpha1 = 0.5, alpha0 = 0.5), list(Z = z, alpha1 = 1.5, alpha0 = 1.5)), pars = c("ORadj"), chains = 2, warmup = 10000, iter = 500000, verbose = FALSE)
  print(summary(model))
  return(model)
}

#' Draw Models
#'
#' This function plots any of the above six models and a solid and dotted line that are the probability density of crude and corrected OR, respectively (more information in the paper of Chu et al. 2006). A file named "draws.pdf" will be returned to your local directory.
#' @param modellist a list of stanfit objects that are returned by stan() command.
#' @param graphNames a vector of strings of names given to the plots, must be of the same length as "modelList".
#' @param a # of exposed subjects in the case group. Along with N1, c, N0, se and sp, they are used to plot the solid and dotted line that each represent a model with no misclassification and misclassification with constant Se and Sp.
#' @param N1 # of total subjects in the case group.
#' @param c # of exposed subjects in the control group.
#' @param N0 # of total subjects in the control group.
#' @param se sensitivity
#' @param sp specificity
#' @param lowerprob the lowest value of all the lowerprob-percentiles of all the models that will be plotted will be taken as the lower limit for all the plots. Default to 0.01.
#' @param upperprob the highest value of all the upperprob-percentiles of all the models that will be plotted will be taken as the upper limit for all the plots. Default to 0.99.
#' @return It creates a models.pdf file under the working directory, that consists of the plots
#' of the models arranged in an appropriate grid.
#' @import rstan
#' @import Rcpp
#' @import ggplot2
#' @import gridExtra
#' @export
#' @examples
#' A <- correctedOR(a = 126, N1 = 218, c = 71, N0 = 295)
#' B <- crudeOR(a = 126, N1 = 218, c = 71, N0 = 295, se = 0.94, sp = 0.97)
#' C <- logitOR(a = 126, N1 = 218, c = 71, N0 = 295, mx = 1.9814, my = 2.8586,
#' varx = 0.9872, vary = 0.5806, x = 2.76, y = 3.58)
#' drawModels(modelList = list(A, B, C), graphNames = c("(A)", "(B)", "(C)"),
#' a = 126, N1 = 218, c = 71, N0 = 295, se = 0.94, sp = 0.97)

drawModels <- function(modelList, graphNames, a, N1, c, N0, se, sp, lowerprob = 0.01, upperprob = 0.99) {
  if (length(modelList) != length(graphNames)) {
    stop("The length of list of models and that of vector of graph names do not match.")
  }
  b <- N1 - a
  d <- N0 - c
  #when se and sp are equal to user input values (crude OR)
  ORadj <- ((a + (sp - 1) * N1) * (se * N0 - c))/((c + (sp - 1) * N0) * (se * N1 - a))
  term1 <- (a * b * (a + b))/(((-b + (a + b) * sp)^2) * ((se * (a + b) - a)^2))
  term2 <- (c * d * (c + d))/(((-d + (c + d) * sp)^2) * ((se * (c + d) - c)^2))
  se.log.or <- (se + sp - 1) * sqrt(term1 + term2)

  xx <- seq(0, 10, 0.1)
  yy <- dlnorm(xx, log(ORadj), se.log.or)

  #when se and sp are both equal to 1 (corrected OR)
  ORadj <- ((a + (1 - 1) * N1) * (1 * N0 - c))/((c + (1 - 1) * N0) * (1 * N1 - a))
  term1 <- (a * b * (a + b))/(((-b + (a + b) * 1)^2) * ((1 * (a + b) - a)^2))
  term2 <- (c * d * (c + d))/(((-d + (c + d) * 1)^2) * ((1 * (c + d) - c)^2))
  se.log.or <- (1 + 1 - 1) * sqrt(term1 + term2)
  yy2 <- dlnorm(xx, log(ORadj), se.log.or)

  y1 <- data.frame(x = xx, y = yy)
  y2 <- data.frame(x = xx, y = yy2)

  #used to extract a stanfit object from a list of stanfit objects
  separateList <- function(modelList, i) {
    return (modelList[[i]])
  }
  l <- length(modelList)
  lowerbounds <- c()
  upperbounds <- c()
  hists <- list()
  updatedhists <- list()
  #create a list of stanhist objects
  for (i in 1:l){
    model <- separateList(modelList, i)
    lowerbounds[i] <- (summary(model, pars = c("ORadj"), probs = c(lowerprob, upperprob))$summary)[4]
    upperbounds[i] <- (summary(model, pars = c("ORadj"), probs = c(lowerprob, upperprob))$summary)[5]
    hists[[i]] <- stan_hist(model, binwidth = 0.25, fill = "gray", color = "white")
  }
  #create PDF of different dimensions depending on the # of graphs to be plotted
  if (sqrt(l) %% 1 == 0) {
    pdf(file = "models.pdf", height = 3 * sqrt(l), width = 3 * sqrt(l), colormodel = "rgb")
  } else {
    pdf(file = "models.pdf", height=9 * ((l - 1) %/% 3 + 1) / 3, width = 9, colormodel = "rgb")
  }
  low <- min(lowerbounds)
  up <- max(upperbounds)
  print(paste("Lower limit of all the plots is ", low, ".", sep = ""))
  print(paste("Upper limit is ", up, ".", sep = ""))
  #update with appropriate x scale and geom lines
  for (i in 1:l) {
    updatedhists[[i]] <- hists[[i]] + scale_x_continuous(name = graphNames[i], limits = c(low,up)) +geom_line(mapping = aes(x = x, y = y), data = y1, linetype = 2, colour = "black") + geom_line(mapping = aes(x = x, y = y), data = y2, linetype = 1, colour = "black")
  }
  #arrange them properly in the PDF
  if (sqrt(l) %% 1 == 0) {
    grid.arrange(grobs = updatedhists, ncol = sqrt(l), nrow = sqrt(l))
  } else {
    grid.arrange(grobs = updatedhists, ncol = 3, nrow = ceiling(l / 3))
  }
  dev.off()
}
