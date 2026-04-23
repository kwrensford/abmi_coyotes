functions {
  real nmix_site_lpmf(array[] int y_site,
                      int Kmin,
                      int K,
                      real log_lambda_i,
                      vector logit_p_site,
                      real phi) {

    int n_obs = size(y_site);
    vector[K - Kmin + 1] lpN;

    for (n_idx in 1:(K - Kmin + 1)) {
      int N = Kmin + n_idx - 1;

      if (max(y_site) > N) {
        lpN[n_idx] = negative_infinity();
      } else {
        real lp_y = 0;
        for (m in 1:n_obs) {
          lp_y += binomial_logit_lpmf(y_site[m] | N, logit_p_site[m]);
        }

        // Use log-mean parameterization to avoid exp underflow
        lpN[n_idx] =
          neg_binomial_2_log_lpmf(N | log_lambda_i, phi)
          + lp_y;
      }
    }
    return log_sum_exp(lpN);
  }
}


data {
  int<lower=1> I;
  array[I] int<lower=1> J;
  int<lower=1> max_J;
  int<lower=1> max_K;

  array[I, max_J, max_K] int<lower=0> y;
  array[I, max_J, max_K] int<lower=0, upper=1> use_y;

  vector[I] clim_pc1;
  vector[I] clim_pc2;
  vector[I] human_fp;

  array[I, max_J] int<lower=0, upper=1> baited;

  int<lower=0> K;
}

transformed data {
  array[I] int n_obs_site;
  array[I] int Kmin;
  array[I, max_J * max_K] int y_flat;
  array[I] int idx_start;
  array[I] int idx_end;

  {
    int idx;
    for (i in 1:I) {
      idx = 0;
      Kmin[i] = 0;

      for (j in 1:J[i]) {
        for (k in 1:max_K) {
          if (use_y[i, j, k] == 1) {
            idx += 1;
            y_flat[i, idx] = y[i, j, k];
            if (y[i, j, k] > Kmin[i])
              Kmin[i] = y[i, j, k];
          }
        }
      }

      n_obs_site[i] = idx;
      idx_start[i] = 1;
      idx_end[i]   = idx;
    }
  }
}

parameters {
  real alpha_lambda;
  real beta_clim1;
  real beta_clim2;
  real beta_human;

  real beta_interact1;   // NEW
  real beta_interact2;   // NEW

  real alpha_p;
  real beta_bait;

  vector[sum(J)] u_cam_raw;
  real<lower=0> sigma_cam;

  real<lower=1e-6> phi;
  
  real logit_psi;
}

transformed parameters {
  vector[sum(J)] u_cam = sigma_cam * u_cam_raw;
  real psi = inv_logit(logit_psi);
}

model {
  alpha_lambda ~ normal(0, 2);
  beta_clim1   ~ normal(0, 2);
  beta_clim2   ~ normal(0, 2);
  beta_human   ~ normal(0, 2);
  beta_interact1 ~ normal(0, 2);
  beta_interact2 ~ normal(0, 2);


  alpha_p      ~ normal(0, 2);
  beta_bait    ~ normal(0, 2);

  sigma_cam    ~ exponential(1);
  u_cam_raw    ~ normal(0, 1);

  phi          ~ gamma(2, 0.5);
  
  logit_psi ~ normal(0, 2);   // NEW

  {
    int cam_idx = 1;

    for (i in 1:I) {
      real log_lambda_i =
    alpha_lambda
    + beta_clim1 * clim_pc1[i]
    + beta_clim2 * clim_pc2[i]
    + beta_human * human_fp[i]
    + beta_interact1 * clim_pc1[i] * human_fp[i]
    + beta_interact2 * clim_pc2[i] * human_fp[i];


      vector[n_obs_site[i]] logit_p_site;
      int idx = 0;

      for (j in 1:J[i]) {
        real logit_p_cam =
          alpha_p +
          beta_bait * baited[i, j] +
          u_cam[cam_idx];

        for (k in 1:max_K) {
          if (use_y[i, j, k] == 1) {
            idx += 1;
            logit_p_site[idx] = logit_p_cam;
          }
        }

        cam_idx += 1;
      }

      if (Kmin[i] == 0) {
        real lp_struct_zero =
          bernoulli_logit_lpmf(1 | logit_psi);

        real lp_nb_part =
          bernoulli_logit_lpmf(0 | logit_psi) +
          nmix_site_lpmf(
            y_flat[i, idx_start[i]:idx_end[i]] |
            Kmin[i], K, log_lambda_i, logit_p_site, phi);

        target += log_sum_exp(lp_struct_zero, lp_nb_part);

      } else {
        target +=
          bernoulli_logit_lpmf(0 | logit_psi) +
          nmix_site_lpmf(
            y_flat[i, idx_start[i]:idx_end[i]] |
            Kmin[i], K, log_lambda_i, logit_p_site, phi);
      }   
    }
  }
}
