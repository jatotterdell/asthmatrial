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
  simplex[K] p;
}

transformed parameters {
  vector[P] beta = beta_sd .* beta_raw;
  ordered[K-1] alpha = logit(cumulative_sum(p[1:(K-1)]));
}

model {
  // lp
  vector[N] eta = X * beta;
  // prior
  // - beta ~ Normal(0, beta_sd)
  // - p(eta = 0) ~ Dirichlet(p_par)
  target += std_normal_lpdf(beta_raw) +
            dirichlet_lpdf(p | p_par);
  // likelihood
  target += ordered_logistic_lpmf(y | eta, alpha);
}
