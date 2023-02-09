data {
  int N; // number of observations
  int P; // number of treatment groups
  int K; // number of outcome levels
  array[N] int x; // treatment indicator
  array[N] int y; // outcome level
  real<lower=0> beta_sd; // prior marginal standard deviation scaling for beta
  vector<lower=0> [K] p_par; // dirichlet prior hyper-parameters
  int<lower=0,upper=1> prior;
}

transformed data{
  // enforce sum-to-zero contraint on beta
  matrix[P, P] A = diag_matrix(rep_vector(1, P)) - inv(P);
  matrix[P, P - 1] A_qr;
  A_qr = qr_Q(A)[ , 1:(P - 1)];
}

parameters {
  vector[P-1] beta_raw;
  simplex[K] p;
}

transformed parameters {
  vector[P] beta =  A_qr * (beta_sd * beta_raw);
  ordered[K-1] alpha = logit(cumulative_sum(p[1:(K-1)]));
}

model {
  // prior
  // - beta ~ Normal(0, beta_sd)
  // - p(eta = 0) ~ Dirichlet(p_par)
  target += std_normal_lpdf(beta_raw) +
            dirichlet_lpdf(p | p_par);
  // likelihood
  if (!prior) {
    for(n in 1:N)
      target += ordered_logistic_lpmf(y[n] | beta[x[n]], alpha);
  }

}
