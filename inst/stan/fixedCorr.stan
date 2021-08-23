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
  real mLogit_pi0;
  real<lower=0> sLogit_pi0;
  real mLOR_c;
  real<lower=0> sLOR_c;
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
  logit_pi0 ~ normal(mLogit_pi0, sLogit_pi0);
  LOR_c ~ normal(mLOR_c, sLOR_c);
}
