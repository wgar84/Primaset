data {
  int k; // traits (PCs retained)
  int m; // taxa
  cov_matrix[2 * (m - 1)] cov_phylo; // inverse cov phylo
  int ni[m]; // sample sizes
  int ni_max; // max ni
  vector[k] Y[m, ni_max]; // data (projections along PCs)
  matrix[m, ni_max] X; // predictor variable (centroid size)
}

parameters {

  /** regression parameters **/
  
  vector[k] as_term[m];
  vector[k] bs_term[m];

  vector[k] as_anc[m - 2];
  vector[k] bs_anc[m - 2];

  vector[k] as_root;
  vector[k] bs_root;

  /** variances **/

  cov_matrix[k] Sigma_e; // cov between projections
  cov_matrix[k] Sigma_a; // cov between intercepts
  cov_matrix[k] Sigma_b; // cov between slopes

  /** no covariances between slopes and intercepts ffs **/
  
}

model {

  matrix[k, k] invSigma_a;
  matrix[k, k] invSigma_b;
  matrix[2 * (m - 1), k] as_center;
  matrix[2 * (m - 1), k] bs_center;
  vector[k] mu[m, ni_max];
  real ldetSigma_a;
  real ldetSigma_b;

  /** priors? **/

  /**
  
  /** intercepts **/
  for(i in 1:m)
    as_center[i] = to_row_vector(as_term[i] - as_root);

  for(i in 1:(m - 2))
    as_center[i + m] = to_row_vector(as_term[i] - as_root);

  invSigma_a = inverse_spd(Sigma_a);

  ldetSigma_a = log_determinant(Sigma_a);
  
  target +=
    - 0.5 * trace_gen_quad_form(invSigma_a, cov_phylo, as_center)
    - (m - 1) * ldetSigma_a;

  /** slopes **/
  for(i in 1:m)
    bs_center[i] = to_row_vector(bs_term[i] - bs_root);

  for(i in 1:(m - 2))
    bs_center[i + m] = to_row_vector(bs_term[i] - bs_root);

  invSigma_b = inverse_spd(Sigma_b);

  ldetSigma_b = log_determinant(Sigma_b);
  
  target +=
    - 0.5 * trace_gen_quad_form(invSigma_b, cov_phylo, bs_center)
    - (m - 1) * ldetSigma_b;
  
  /** pops **/
  for(i in 1:m)
    for(j in 1:(ni[i]))
      {
	mu[i,j] = as_term[i] + bs_term[i] * X[i,j];
	Y[i,j] ~ multi_normal(mu[i,j], Sigma_e);
      }
}
