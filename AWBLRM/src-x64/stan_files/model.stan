data {
  int<lower=0> N;   // No. of observations
  int<lower=0, upper=1> y[N];   // Binary outcome for DLT
  real dose[N];   // Dose level assigned to patient 
  real<lower=0, upper=1> W[N];   // Adaptive weight
}

parameters{
  real<lower = 0> beta;   // slope of dose-toxicity curve
  real alpha;   // intercept of dose-toxicity curve
}

model{
  beta ~ normal(0,2);         // prior for slope 
  alpha ~ normal(0,2);        // prior for alpha
  
  for (n in 1:N) {
  y[n] ~ bernoulli(W[n] * (exp(beta * dose[n] + alpha) / (1 + exp(beta * dose[n] + alpha))));      // likelihood
  }
} 

