data {
  int N; // number of observations
  int P; // number of design parameters
  int K; // number of outcome levels
  matrix[N, P] X; // design matrix (no intercept)
  array[N] int y; // outcome level
  vector[P] beta_sd; // prior standard deviation for beta parameters
  vector<lower=0> [K] p_par; // dirichlet prior hyper-parameters
}

parameters {
  vector[P] beta_raw;
  real alpha_int;
  vector<lower=0>[K-2] alpha_raw;
  real<lower=1> tau;
}

transformed parameters {
  vector[P] beta = beta_sd .* beta_raw;
  ordered[K-1] alpha;
  alpha[1] = 10 * alpha_int;
  alpha[2:(K-1)] = alpha[1] + tau * cumulative_sum(alpha_raw);
}

model {
  // lp
  vector[N] eta = X * beta;
  // prior
  // - beta ~ Normal(0, beta_sd)
  // - tau ~ HalfT(3, 0, 1)
  // - alpha ~ RW(1)
  target += std_normal_lpdf(beta_raw) +
            std_normal_lpdf(alpha_int) +
            std_normal_lpdf(alpha_raw) - 1 * std_normal_lccdf(0) +
            student_t_lpdf(tau | 3, 0, 1) - 1 * student_t_lccdf(0 | 3, 0, 1);
  // likelihood
  target += ordered_logistic_lpmf(y | eta, alpha);
}
