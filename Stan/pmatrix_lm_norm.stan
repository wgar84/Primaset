data {
  int<lower = 0> N;      // individuals
  int<lower = 1> J;      // number of fixed effects
  int<lower = 1> K;      // traits
  vector[K] Y[N];        // measurements
  vector[J] X[N];        // model matrix
}

transformed data {
  vector[K] Y_std[N];
  vector[K] Y_sd;
  vector[K] Y_mean;

  for(k in 1:K){
    Y_sd[k] = sd(Y[,k]);
    Y_mean[k] = mean(Y[,k]);
    for(n in 1:N)
      Y_std[n,k] = (Y[n,k] - Y_mean[k]) / Y_sd[k];
  }
}

parameters {
  cholesky_factor_corr[K] Omega_P; // correlation structure
  vector[K] sigma_P; // variances
  matrix[K, J] beta; // fixed effects
}

model {
  vector[K] mu[N];
  matrix[K, K] Sigma_P;
  
  /** priors **/

  to_vector(beta) ~ normal(0, 1);
  sigma_P ~ cauchy(0, 2.5);
  Omega_P ~ lkj_corr_cholesky(4);

  /** model **/
  
  for(n in 1:N)
    mu[n] = beta * X[n];

  Sigma_P = diag_pre_multiply(sigma_P, Omega_P);
  
  Y_std ~ multi_normal_cholesky(mu, Sigma_P);
}

generated quantities {
  cov_matrix[K] P;
  vector[K] sigma_P_rescale;

  sigma_P_rescale = Y_sd .* sigma_P;

  P = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_P_rescale, Omega_P));
}

