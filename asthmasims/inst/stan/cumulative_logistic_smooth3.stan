functions {
  vector inv_ordered_logit(vector x) {
    int K = num_elements(x);
    vector[K] cp = inv_logit(x);
    vector[K+1] p;
    for (k in 1:K+1) {
      if (k == 1) p[k] = cp[1];
      else if (k == K+1) p[k] = 1 - cp[K];
      else p[k] = cp[k] - cp[k-1];
    }
    return p;
  }
}

data {
  int N; // number of observations
  int P; // number of treatment groups
  int K; // number of outcome levels
  int P_Isp; // dimension of Ispline basis
  array[N] int x; // treatment indicator
  array[N] int y; // outcome level
  vector[K] y_weights;
  matrix[K-1, P_Isp] X_Isp;
  real<lower=0> beta_sd; // prior marginal standard deviation scaling for beta
  real<lower=0> alpha_int_sd; // scale on intercept term
  vector[P_Isp] theta_sd; // prior scale on Ispline coefficients
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
  vector<lower=0>[P_Isp] theta;
  real alpha_int;
}

transformed parameters {
  vector[P] beta =  A_qr * (beta_sd * beta_raw);
  vector[K-1] alpha = alpha_int + X_Isp * (theta_sd .* theta);
}

model {
  // prior
  // - beta ~ Normal(0, beta_sd)
  target += std_normal_lpdf(beta_raw) +
            normal_lpdf(alpha_int | 0, alpha_int_sd) +
            normal_lpdf(theta | 0, 1);
  // likelihood
  if (!prior) {
    for(n in 1:N)
      target += ordered_logistic_lpmf(y[n] | beta[x[n]], alpha);
  }
}

generated quantities {
  matrix[K, P] p;
  vector[P] mu;
  vector[P] y_ppc;
  for (i in 1:P) {
    p[, i] = inv_ordered_logit(alpha - beta[i]);
    mu[i] = dot_product(y_weights, p[, i]);
    y_ppc[i] = y_weights[ordered_logistic_rng(beta[i], alpha)];
  }
}
