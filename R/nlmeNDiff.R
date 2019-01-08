#' Non-differential Generalized Linear Mixed Effects Model
#'
#' Fit a bivariate generalized linear mixed-effects model (GLMM) for non-differential sensitivity and specificity using the \code{glmer} function in \code{lme4}. 
#' Lower and upper bounds for Se and Sp can be specified according to the assumptions of the study. 
#' @param data a data frame containing the 2 by 2 data of the diagnostics table of exposure status for every study in a meta-analysis.
#' It contains at least 4 columns in the data named as following: \code{n11} indicates the true positives, \code{n01} the false positives,
#' \code{n00} the true negatives and \code{n10} the false negatives. Each column is a vector of same length,
#' which is the number of meta-analysis study results used in the model.
#' @param lower an optional argument specifying the lower bound assumption of Se and Sp. Default to 0.5 (or the lowest Se/Sp of all studies, whichever is lower), which provides the mild assumption that Se and Sp are better than chance.
#' @param upper an optional argument specifying the upper bound assumption of Se and Sp. Default to 1.
#' @param id a TRUE of FALSE argument indicating if the supplied data has a \code{sid} column that gives same studies
#' the same subject ID. Default to FALSE, which assumes that all studies have different IDs.
#' @param ... optional parameters passed to \href{https://www.rdocumentation.org/packages/lme4/versions/1.1-18-1/topics/glmer}{glmer}.
#' @return It returns an object of class \href{https://www.rdocumentation.org/packages/lme4/versions/1.1-18-1/topics/merMod-class}{mermod}. 
#' Besides generic class methods, \code{paramEst()} is implemented in \code{BayesSenMC} to get the parameter estimates used in the Bayesian misclassification model functions. 
#' @import lme4
#' @importFrom dplyr mutate rename select %>%
#' @importFrom stats binomial
#' @export
#' @examples
#' data(bd_meta)
#'
#' mod <- nlme_nondiff(bd_meta, lower = 0)

nlme_nondiff <- function(data, lower = 0.5, upper = 1, id = FALSE, ...) {
  
  # to get rid of cran global not definied warning
  N1 <- n11 <- n10 <- N0 <- n01 <- n00 <- a <- b <- Y <- N <- Se <- sid <- NULL
  
  if (id == FALSE) {
    data$sid <- seq(1, nrow(data))
  }
  
  dat <- data %>% mutate(N1 = n11 + n10,
                        N0 = n01 + n00) %>%
    rename(a = n11, b = n00) %>%
    select(a, b, N1, N0, sid)
  
  # Se = 1 represents sensitivity, otherwise specificity
  dat_final <- merge(
    dat %>% mutate(Y = a, N = N1, Se = 1) %>% select(sid, Y, N, Se),
    dat %>% mutate(Y = b, N = N0, Se = 0) %>% select(sid, Y, N, Se), all = TRUE
  )
  
  # make sure it will not crash with significantly high lower bound value
  lower_data <- min(dat_final$Y / dat_final$N)
  if (lower > min(dat_final$Y / dat_final$N)) {
    lower <- min(lower, min(dat_final$Y / dat_final$N))
    print(paste0("The lower bound used in this model is changed to ", lower, 
                 " to prevent NaNs; see model specification for more details."))
  }
  
  dat_final <- cbind(dat_final, as.numeric(dat_final$Se==0))
  names(dat_final)[c(5)] <- c("Sp")

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
  
  mod <- glmer(cbind(Y, N - Y) ~ ((0 + Se + Sp) | sid) + Se, 
               data = dat_final, family = binomial(link = user_logis(lower, upper)), ...)
  
  attributes(mod)$type <- "nondiff"
  attributes(mod)$lower <- lower
  attributes(mod)$upper <- upper
  
  print(mod)
  return(mod)
}

#' Differential Generalized Linear Mixed Effects Model
#'
#' Fit a bivariate generalized linear mixed-effects model (GLMM) for differential sensitivity and specificity using the \code{glmer} function in \code{lme4}. 
#' Thus, meta-study data for both case and control subjects must be provided. Lower and upper bounds for Se and Sp can be specified according to the assumptions of the study. 
#' @param data a data frame containing the 2 by 2 data of the diagnostics table of exposure status for every study in a meta-analysis.
#' It contains at least 4 columns in the data named as following: \code{n11} indicates the true positives, \code{n01} the false positives,
#' \code{n00} the true negatives and \code{n10} the false negatives. Each column is a vector of same length,
#' which is the number of meta-analysis study results used in the model. It should also have a column \code{group} (0 for controls, 1 for cases).
#' @param lower an optional argument specifying the lower bound assumption of Se and Sp. Default to 0.5 (or the lowest Se/Sp of all studies, whichever is lower), which provides the mild assumption that Se and Sp are better than chance.
#' @param upper an optional argument specifying the upper bound assumption of Se and Sp. Default to 1.
#' @param id a TRUE of FALSE argument indicating if the supplied data has a \code{sid} column that gives same studies
#' the same subject ID. Default to FALSE, which assumes that all studies have different IDs.
#' @param ... optional parameters passed to \href{https://www.rdocumentation.org/packages/lme4/versions/1.1-18-1/topics/glmer}{glmer}.
#' @return It returns an object of class \href{https://www.rdocumentation.org/packages/lme4/versions/1.1-18-1/topics/merMod-class}{mermod}. 
#' Besides generic class methods, \code{paramEst()} is implemented in \code{BayesSenMC} to get the parameter estimates used in the Bayesian misclassification model functions. 
#' @import lme4
#' @importFrom dplyr mutate rename select %>%
#' @importFrom stats binomial
#' @export
#' @examples
#' data(bd_meta)
#' library(dplyr)
#' bd_meta <- bd_meta %>% mutate(group = ifelse(methods == "BSDS", 1, 0))
#' mod <- nlme_diff(bd_meta, lower = 0)

nlme_diff <- function(data, lower = 0.5, upper = 1, id = FALSE, ...) {
  
  # to get rid of cran global not definied warning
  N1 <- n11 <- n10 <- N0 <- n01 <- n00 <- a <- b <- Y <- N <- Se <- group <- sid <- NULL
  
  if (id == FALSE) {
    data$sid <- seq(1, nrow(data))
  }
  
  dat <- data %>% mutate(N1 = n11 + n10,
                         N0 = n01 + n00) %>%
    rename(a = n11, b = n00) %>%
    select(a, b, N1, N0, group)
  
  sid <- seq(1, nrow(dat))
  # Se = 1 represents sensitivity, otherwise specificity
  dat_final <- merge(
    dat %>% mutate(sid, Y = a, N = N1, Se = 1) %>% select(sid, Y, N, Se, group),
    dat %>% mutate(sid, Y = b, N = N0, Se = 0) %>% select(sid, Y, N, Se, group), all = TRUE
  )
  
  # make sure it will not crash with significantly high lower bound value
  lower_data <- min(dat_final$Y / dat_final$N)
  if (lower > lower_data) {
    lower <- min(lower,lower_data)
    print(paste0("The lower bound used in this model is changed to ", lower, 
                 " to prevent NaNs; see model specification for more details."))
  }
  
  dat_final <- cbind(dat_final, as.numeric(dat_final$Se==0))
  names(dat_final)[c(6)] <- c("Sp")
  
  
  dat_final <- cbind(dat_final, as.numeric(dat_final$group==0))
  dat_final <- cbind(dat_final, 1 - as.numeric(dat_final$group==0))
  names(dat_final)[c(7, 8)] <- c("controls", "cases")
  
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
  
  mod <- glmer(cbind(Y, N-Y) ~ ((0 + cases + controls) | sid) + ((0 + Se + Sp) | sid) + Se + cases, 
               data = dat_final, family = binomial(link = user_logis(lower, upper)), ...)
  
  attributes(mod)$type <- "diff"
  attributes(mod)$lower <- lower
  attributes(mod)$upper <- upper
  
  print(mod)
  return(mod)
}

#' Parameter estimates of the GLMM model
#'
#' Get parameter estimates of the GLMM model to plug into modeling functions in \code{BayesSenMC} for Bayesian inference of adjusted odds ratio.
#' @param model a GLMM model built with the \code{nlme_nondiff()} function.
#' @param lower an optional argument matching the lower bound assumption of Se and Sp of the input \code{model}. Default to 0.5 as in \code{nlme_nondiff()}.
#' @param upper an optional argument matching the upper bound assumption of Se and Sp. Default to 1 as in \code{nlme_nondiff()}.
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

paramEst <- function(model, lower = 0.5, upper = 1) {
  if (is.null(attributes(model)$type) || is.null(attributes(model)$lower) 
      || is.null(attributes(model)$upper)) {
    stop("Model must be generated from the nlme_nondiff() or the nlme_diff() function.")
  }
  
  lower <- attr(model, "lower")
  upper <- attr(model, "upper")
  
  if (attr(model, "type") == "nondiff") {
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
      Se = (exp(alpha) / (1 + exp(alpha))) * (upper - lower) + lower,
      Sp = (exp(beta) / (1 + exp(beta))) * (upper - lower) + lower,
      rho = SeSp/sqrt(s2Se*s2Sp),
      fisher_mean = 0.5*(log(1+SeSp/sqrt(s2Se*s2Sp))-log(1-SeSp/sqrt(s2Se*s2Sp)))
    ))
  }
  
  else if (attr(model, "type") == "diff") {
    s <- summary(model)
    intercept <- s$coefficients[1,1]
    se <- s$coefficients['Se', 1]
    group <- s$coefficients['group', 1]
    
    
    return(list(
      m.lg.se0 = intercept + se,
      m.lg.se1 = intercept + se + group,
      m.lg.sp0 = intercept,
      m.lg.sp1 = intercept + group,
      s.lg.se = sqrt(s$varcor$sid.1[1]),
      s.lg.sp = sqrt(s$varcor$sid.1[4]),
      corr.group = s$varcor$sid[3] / sqrt(s$varcor$sid[1] * s$varcor$sid[4]),
      corr.sesp = s$varcor$sid.1[3] / sqrt(s$varcor$sid.1[1] * s$varcor$sid.1[4])
    ))
  }
  
  else {stop("Model is not generated from nlme_nondiff() or nlme_diff().")}
}