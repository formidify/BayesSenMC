#' Non-differential Generalized Linear Mixed Effects Model
#'
#' Fit a bivariate generalized linear mixed-effects model (GLMM) for non-differential sensitivity and specificity using the \code{glmer} function in \code{lme4}. 
#' Lower and upper bounds for Se and Sp can be specified according to the assumptions of the study. 
#' @param data a data frame containing the 2 by 2 data of the diagnostics table of exposure status for every study in a meta-analysis.
#' There are 4 columns in the data (\code{n_11} indicates the true positives, \code{n_11} the false negatives,
#' \code{n_00} the true negatives and \code{n_01} the false negatives). Each column is a vector of same length,
#' which is the number of meta-analysis study results used in the model.
#' @param lower an optional argument specifying the lower bound assumption of Se and Sp. Default to 0.5, which provides the mild assumption that it is better than chance.
#' @param upper an optional argument specifying the upper bound assumption of Se and Sp. Default to 1.
#' @param ... optional parameters passed to \code{\href{https://www.rdocumentation.org/packages/lme4/versions/1.1-18-1/topics/glmer}{glmer}}.
#' @return It returns an object of class \code{\href{https://www.rdocumentation.org/packages/lme4/versions/1.1-18-1/topics/merMod-class}{mermod}}. 
#' Besides generic class methods, \code{paramEst()} is implemented in \code{BayesSenMC} to get the parameter estimates used in the Bayesian misclassification model functions. 
#' @import lme4
#' @export
#' @examples
#' raw <- data.frame(n11 = c(59, 312, 177, 149, 213, 214, 126, 115, 1357, 19),
#' n10 = c(5, 24, 4, 8, 47, 0, 11, 7, 185, 2),
#' n00 = c(227, 594, 178, 370, 76, 20, 181, 179, 3322, 96),
#' n01 = c(5, 45, 9, 15, 2, 2, 0, 2, 68, 1))
#'
#' mod <- nlme_nondiff(raw, lower = 0)

nlme_nondiff <- function(data, lower = 0.5, upper = 1, ...) {
  
  dat <- raw %>% mutate(N1 = n11 + n10,
                        N0 = n00 + n01) %>%
    rename(a = n11, b = n00) %>%
    select(a, b, N1, N0)
  
  dat_final <- merge(
    dat %>% mutate(sid = seq(1, 10), Y = a, N = N1, Se = 1) %>% dplyr::select(sid, Y, N, Se),
    dat %>% mutate(sid = seq(1, 10), Y = b, N = N0, Se = 0) %>% dplyr::select(sid, Y, N, Se), all = TRUE
  )
  
  dat_final <- cbind(dat_final, as.numeric(dat_final$Se=='0'))
  dat_final <- cbind(dat_final, 1 - as.numeric(dat_final$Se=='0'))
  names(dat_final)[c(5, 6)] <- c("control","experimental")

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
  
  # supposedly for differential logistic function (does not work)
  # might have to go inside glmer() and change things around a bit
  user_logis_diff <- function(ab1, ab2, se) {
    if (se == 0) {
      a = ab1[1]
      b = ab1[2]
    }
    else {
      a = ab2[1]
      b = ab2[2]
    }
    linkfun <- function(mu) log((mu - a) / (a + b - mu))
    linkinv <- function(eta) a + b / (exp(-1*eta) + 1)
    mu.eta <- function(eta) b * exp(-1*eta) / (exp(-1*eta) + 1)**2
    valideta <- function(eta) TRUE
    link <- paste0("logit(", a, ", ", b, ")")
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
  }
  
  mod <- glmer(cbind(Y, N - Y) ~ ((0 + control + experimental) | sid) + Se, 
               data = dat_final, family = binomial(link = user_logis(lower, upper)), ...)
  
  return(mod)
}

#' Model without misclassification
#'
#' Generate a stanfit object corresponding to a posterior distribution of uncorrected odds ratio given no misclassification.
#' @param model a GLMM model built with the \code{nlme_nondiff()} function.
#' @return It returns a list of parameter estimates which can be input into the Bayesian model functions in
#' \code{BayesSenMC}. \code{(mean_logSe, var_logSe)} and \code{(mean_logSp, var_logSp)} are the logit prior distributions for Se and Sp.
#' \code{Se} and \code{Sp} are the corresponding mean values given the logit prior means. \rho is the correlation estimate between Se and
#' Sp. \code{fisher_mean} is the Fisher's mean of the correlation assume a Fisher's distribution.
#' @import lme4
#' @export
#' @examples
#' mod <- nlme_nondiff(raw) # see nlme_nondiff() for detailed example.
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