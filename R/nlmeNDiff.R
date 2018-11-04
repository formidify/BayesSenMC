#' Non-differential Generalized Linear Mixed Effects Model
#'
#' Fit a bivariate generalized linear mixed-effects model (GLMM) for non-differential sensitivity and specificity using the \code{glmer} function in \code{lme4}. 
#' Lower and upper bounds for Se and Sp can be specified according to the assumptions of the study. 
#' @param data a data frame containing the 2 by 2 data of the diagnostics table of exposure status for every study in a meta-analysis.
#' It contains at least 4 columns in the data named as following: \code{TP} indicates the true positives, \code{FN} the false negatives,
#' \code{TN} the true negatives and \code{FN} the false negatives. Each column is a vector of same length,
#' which is the number of meta-analysis study results used in the model.
#' @param lower an optional argument specifying the lower bound assumption of Se and Sp. Default to 0.5, which provides the mild assumption that Se and Sp are better than chance.
#' @param upper an optional argument specifying the upper bound assumption of Se and Sp. Default to 1.
#' @param ... optional parameters passed to \href{https://www.rdocumentation.org/packages/lme4/versions/1.1-18-1/topics/glmer}{glmer}.
#' @return It returns an object of class \href{https://www.rdocumentation.org/packages/lme4/versions/1.1-18-1/topics/merMod-class}{mermod}. 
#' Besides generic class methods, \code{paramEst()} is implemented in \code{BayesSenMC} to get the parameter estimates used in the Bayesian misclassification model functions. 
#' @import lme4
#' @import dplyr
#' @importFrom stats binomial
#' @export
#' @examples
#' data(bd_meta)
#'
#' mod <- nlme_nondiff(bd_meta, lower = 0)
#' summary(mod) 

nlme_nondiff <- function(data, lower = 0.5, upper = 1, ...) {
  
  # to get rid of cran global not definied warning
  N1 <- TP <- FN <- N0 <- TN <- FP <- a <- b <- Y <- N <- Se <- NULL
  
  dat <- data %>% mutate(N1 = TP + FN,
                        N0 = TN + FP) %>%
    rename(a = TP, b = TN) %>%
    dplyr::select(a, b, N1, N0)
  
  sid <- seq(1, nrow(dat))
  # Se = 1 represents sensitivity, otherwise specificity
  dat_final <- merge(
    dat %>% mutate(sid, Y = a, N = N1, Se = 1) %>% select(sid, Y, N, Se),
    dat %>% mutate(sid, Y = b, N = N0, Se = 0) %>% select(sid, Y, N, Se), all = TRUE
  )
  
  dat_final <- cbind(dat_final, as.numeric(dat_final$Se=='0'))
  dat_final <- cbind(dat_final, 1 - as.numeric(dat_final$Se=='0'))
  names(dat_final)[c(5, 6)] <- c("is_sp","is_se")

  user_logis <- function(a, b) {
    linkfun <- function(mu) log((mu - a) / (a + b - mu))
    linkinv <- function(eta) a + b / (exp(-1*eta) + 1)
    mu.eta <- function(eta) b * exp(-1*eta) / (exp(-1*eta) + 1)**2
    valideta <- function(eta) TRUE
    link <- paste0("logit(", a, ", ", b, ")")
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
  }
  
  mod <- glmer(cbind(Y, N - Y) ~ ((0 + is_sp + is_se) | sid) + Se, 
               data = dat_final, family = binomial(link = user_logis(lower, upper)), ...)
  
  return(mod)
}

#' Parameter estimates of the GLMM model
#'
#' Get parameter estimates of the GLMM model to plug into modeling functions in \code{BayesSenMC} for Bayesian inference of adjusted odds ratio.
#' @param model a GLMM model built with the \code{nlme_nondiff()} function.
#' @return It returns a list of parameter estimates which can be input into the Bayesian model functions in
#' \code{BayesSenMC}. \code{(mean_logSe, var_logSe)} and \code{(mean_logSp, var_logSp)} are the logit prior distributions for Se and Sp.
#' \code{Se} and \code{Sp} are the corresponding mean values given the logit prior means. \code{rho} is the correlation estimate between Se and
#' Sp. \code{fisher_mean} is the Fisher's mean of the correlation assume a Fisher's distribution.
#' @import lme4
#' @export
#' @examples
#' data(bd_meta)
#' 
#' mod <- nlme_nondiff(bd_meta, lower = 0) # see nlme_nondiff() for detailed example.
#' pList <- paramEst(mod)

paramEst <- function(model) {
  s <- summary(model)
  alpha <- s$coefficients[1, 1] + s$coefficients[2, 1]
  beta <- s$coefficients[1, 1]
  s2Se <- s$varcor$sid[1]
  s2Sp <- s$varcor$sid[4]
  SeSp <- s$varcor$sid[2]
  
  return(list(
    mean_logSe = alpha,
    mean_logSp = beta,
    var_logSe = s2Se,
    var_logSp = s2Sp,
    Se = exp(alpha) / (1 + exp(alpha)),
    Sp = exp(beta) / (1 + exp(beta)),
    rho = SeSp/sqrt(s2Se*s2Sp),
    fisher_mean = 0.5*(log(1+SeSp/sqrt(s2Se*s2Sp))-log(1-SeSp/sqrt(s2Se*s2Sp)))
  ))
}