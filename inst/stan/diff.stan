data {
  int<lower=0> a;
  int<lower=0> N1;
  int<lower=0> c;
  int<lower=0> N0;
  vector[4] Mu;
  matrix[4,4] varZ;
  real mLogit_pi0;
  real<lower=0> sLogit_pi0;
  real mLOR_c;
  real<lower=0> sLOR_c;
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
  logit_pi0 ~ normal(mLogit_pi0, sLogit_pi0);
  LOR_c ~ normal(mLOR_c, sLOR_c);
  Z ~ multi_normal(Mu, varZ);
}
