data {
  int N; // number of observations
  int P; // number of treatment groups
  array[N] int x; // treatment indicator
  vector[N] yc; // outcome response
  real<lower=0> beta_sd; // prior marginal standard deviation scaling for beta
  int<lower=0,upper=1> prior;
}

transformed data{
  real med = quantile(yc, 0.5);
  real mad = quantile(abs(yc - med), 0.5);
  // enforce sum-to-zero contraint on beta
  matrix[P, P] A = diag_matrix(rep_vector(1, P)) - inv(P);
  matrix[P, P - 1] A_qr;
  A_qr = qr_Q(A)[ , 1:(P - 1)];
}

parameters {
  vector[P-1] beta_raw;
  real alpha;
  real<lower=0> sigma;
}

transformed parameters {
  vector[P] beta =  A_qr * (beta_sd * beta_raw);
  vector[P] mu = alpha + beta;
}

model {
  // prior
  // - beta ~ Normal(0, beta_sd)
  // - alpha ~ Normal(median(y), MAD(y))
  // - sigma ~ Exponential(lambda)
  target += std_normal_lpdf(beta_raw) +
            normal_lpdf(alpha | med, mad) +
            exponential_lpdf(sigma | inv(mad));
  // likelihood
  if (!prior) {
    for(n in 1:N)
    // yc ~ normal(mu, sigma) T[0, 100]
      target += normal_lpdf(yc[n] | mu[x[n]], sigma) +
                -log_diff_exp(
                  normal_lcdf(100 | mu[x[n]], sigma),
                  normal_lcdf(0 | mu[x[n]], sigma));
  }
}

generated quantities {
  real c3v1 = mu[3] - mu[1];
  real c3v2 = mu[3] - mu[2];
  vector[P] y_ppc;
  for (i in 1:P)
    y_ppc[i] = normal_rng(mu[i], sigma);
}
