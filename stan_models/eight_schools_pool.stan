
data {
  int <lower = 0> J;
  real y[J];
  real<lower = 0> sigma[J];
}

transformed data {
  real<lower = 0> tau = 2.5;
}

parameters {
  vector[J] theta;
  real mu;
}

model {
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}
