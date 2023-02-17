
data {
}

transformed data {
  real mu = 1.5;
}

parameters {
  real x;
  real y;
}

model {
  x ~ normal(mu, 1);
  x ~ normal(-mu, 1);
  y ~ normal(mu, 1);
  y ~ normal(-mu, 1);
}

