data {
  int<lower=0> a;
  int<lower=0> N1;
  int<lower=0> c;
  int<lower=0> N0;
  real<lower=0, upper=1> Se;
  real<lower=0, upper=1> Sp;
  real mLogit_pi0;
  real<lower=0> sLogit_pi0;
  real mLOR_c;
  real<lower=0> sLOR_c;
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
  logit_pi0 ~ normal(mLogit_pi0, sLogit_pi0);
  LOR_c ~ normal(mLOR_c, sLOR_c);
}
