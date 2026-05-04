functions {
  // =========================================================================
  // nmix_site_lpmf()
  //
  // PURPOSE:
  //   Compute the marginal log-likelihood of all observations at a site under
  //   an N‑mixture model with:
  //
  //     • Latent abundance:
  //         N ~ NegBinomial2(log_lambda_i, phi)
  //
  //     • Detection model:
  //         y_m ~ Binomial(N, p_m)
  //
  //   The function marginalizes over all possible N values from:
  //       N = Kmin ... K
  //
  //   where:
  //     Kmin = max observed count (minimum feasible abundance)
  //     K    = user-specified truncation limit
  //
  // ARGUMENTS:
  //   y_site[]       : vector of observed counts at the site
  //   Kmin           : minimum feasible abundance
  //   K              : maximum abundance allowed in truncation
  //   max_y          : precomputed max(y_site)
  //   log_lambda_i   : log-scale mean abundance for the site
  //   logit_p_site[] : detection logits for each observation
  //   phi            : NegBinomial dispersion parameter
  //
  // RETURNS:
  //   log ∑_N [ P(N | λ, φ) * ∏_m Binomial(y_m | N, p_m) ]
  //
  // OPTIMIZATIONS:
  //   • Vectorized binomial_logit_lpmf
  //   • Avoid repeated max(y_site)
  //   • Avoid slicing inside loops
  // =========================================================================
  real nmix_site_lpmf(array[] int y_site,
                      int Kmin,
                      int K,
                      int max_y,
                      real log_lambda_i,
                      vector logit_p_site,
                      real phi) {

    int n_obs = size(y_site);
    vector[K - Kmin + 1] lpN;

    // Loop over all possible latent abundances N
    for (n_idx in 1:(K - Kmin + 1)) {
      int N = Kmin + n_idx - 1;

      // If any observed count exceeds N, this N is impossible
      if (max_y > N) {
        lpN[n_idx] = negative_infinity();
      } else {

        // Vectorized detection likelihood:
        //   ∑_m log Binomial(y_m | N, p_m)
        real lp_y = binomial_logit_lpmf(y_site | N, logit_p_site);

        // Combine abundance prior + detection likelihood
        lpN[n_idx] =
          neg_binomial_2_log_lpmf(N | log_lambda_i, phi)
          + lp_y;
      }
    }

    // Return log-sum-exp over all possible N values
    return log_sum_exp(lpN);
  }
}


data {
  // =========================================================================
  // DATA INPUTS
  //
  // I          : number of sites
  // J[i]       : number of cameras at site i
  // max_J      : maximum number of cameras across sites
  // max_K      : maximum number of sampling occasions per camera
  //
  // y[i,j,k]   : observed counts
  // use_y      : indicator for whether y[i,j,k] is a valid observation
  //
  // covariates : climate PCs, human footprint, baiting indicator
  //
  // K          : maximum latent abundance allowed in N‑mixture truncation
  // =========================================================================

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
  // =========================================================================
  // PREPROCESSING FOR EFFICIENT LIKELIHOOD EVALUATION
  //
  // For each site i:
  //
  //   n_obs_site[i] : number of valid observations
  //   Kmin[i]       : minimum feasible abundance (max observed count)
  //   max_y_site[i] : max observed count (used inside nmix function)
  //
  //   y_flat[i,*]   : flattened vector of all y for site i
  //   start[i], len[i] : indexing for fast segment() extraction
  //
  // This avoids repeated looping and slicing inside the model block.
  // =========================================================================

  array[I] int n_obs_site;
  array[I] int Kmin;
  array[I] int max_y_site;
  array[I, max_J * max_K] int y_flat;
  array[I] int start;
  array[I] int len;

  {
    for (i in 1:I) {
      int idx = 0;
      int max_y = 0;

      // Flatten all valid y[i,j,k] into y_flat[i,*]
      for (j in 1:J[i]) {
        for (k in 1:max_K) {
          if (use_y[i,j,k] == 1) {
            idx += 1;
            y_flat[i, idx] = y[i,j,k];

            if (y[i,j,k] > max_y)
              max_y = y[i,j,k];
          }
        }
      }

      n_obs_site[i] = idx;
      Kmin[i]       = max_y;
      max_y_site[i] = max_y;

      start[i] = 1;
      len[i]   = idx;
    }
  }
}

parameters {
  // =========================================================================
  // PARAMETERS: ABUNDANCE MODEL
  // =========================================================================
  real alpha_lambda;     // intercept for log λ
  real beta_clim1;       // climate PC1 effect
  real beta_clim2;       // climate PC2 effect
  real beta_human;       // human footprint effect
  real beta_interact1;   // climate PC1 × human footprint
  real beta_interact2;   // climate PC2 × human footprint

  // =========================================================================
  // PARAMETERS: DETECTION MODEL
  // =========================================================================
  real alpha_p;          // intercept for logit p
  real beta_bait;        // baiting effect on detection

  // =========================================================================
  // RANDOM EFFECTS: CAMERA-LEVEL DETECTABILITY
  // =========================================================================
  vector[sum(J)] u_cam_raw;  // non-centered random effects
  real<lower=0> sigma_cam;   // SD of camera-level variation

  // =========================================================================
  // OVERDISPERSION IN ABUNDANCE
  // =========================================================================
  real<lower=1e-6> phi;      // NegBinomial dispersion

  // =========================================================================
  // ZERO-INFLATION: STRUCTURAL ABSENCE
  // =========================================================================
  real logit_psi;            // logit probability of structural zero
}

transformed parameters {
  // =========================================================================
  // NON-CENTERED RANDOM EFFECTS
  // =========================================================================
  vector[sum(J)] u_cam = sigma_cam * u_cam_raw;

  // =========================================================================
  // PRECOMPUTE CAMERA-LEVEL DETECTION LOGITS
  //
  // This avoids recomputing:
  //   alpha_p + beta_bait * baited + u_cam
  //
  // inside the likelihood loop.
  // =========================================================================
  vector[sum(J)] logit_p_cam;
  {
    int idx = 1;
    for (i in 1:I)
      for (j in 1:J[i]) {
        logit_p_cam[idx] =
          alpha_p +
          beta_bait * baited[i,j] +
          u_cam[idx];
        idx += 1;
      }
  }
}

model {

  // =========================================================================
  // PRIORS: ABUNDANCE PROCESS
  // =========================================================================
  alpha_lambda ~ normal(0, 2);
  beta_clim1   ~ normal(0, 2);
  beta_clim2   ~ normal(0, 2);
  beta_human   ~ normal(0, 2);
  beta_interact1 ~ normal(0, 2);
  beta_interact2 ~ normal(0, 2);

  // =========================================================================
  // PRIORS: DETECTION PROCESS
  // =========================================================================
  alpha_p ~ normal(0, 2);
  beta_bait ~ normal(0, 2);

  // =========================================================================
  // CAMERA-LEVEL RANDOM EFFECTS
  // =========================================================================
  sigma_cam ~ exponential(1);
  u_cam_raw ~ normal(0, 1);

  // =========================================================================
  // OVERDISPERSION IN ABUNDANCE
  // =========================================================================
  phi ~ gamma(2, 0.5);

  // =========================================================================
  // ZERO-INFLATION
  // =========================================================================
  logit_psi ~ normal(0, 2);

  // Precompute zero-inflation log-probabilities
  real log_psi1 = bernoulli_logit_lpmf(1 | logit_psi);
  real log_psi0 = bernoulli_logit_lpmf(0 | logit_psi);

  // =========================================================================
  // LIKELIHOOD
  //
  // For each site:
  //   1. Compute log λ_i from covariates
  //   2. Build detection logits for all valid observations
  //   3. Evaluate zero-inflated N‑mixture likelihood
  //
  // All heavy computation is now vectorized or precomputed.
  // =========================================================================

  {
    int cam_idx = 1;  // global camera index

    for (i in 1:I) {

      // ----------------------------------------------------------------------
      // SITE-LEVEL ABUNDANCE MODEL
      // ----------------------------------------------------------------------
      real log_lambda_i =
        alpha_lambda
        + beta_clim1 * clim_pc1[i]
        + beta_clim2 * clim_pc2[i]
        + beta_human * human_fp[i]
        + beta_interact1 * clim_pc1[i] * human_fp[i]
        + beta_interact2 * clim_pc2[i] * human_fp[i];

      // ----------------------------------------------------------------------
      // DETECTION LINEAR PREDICTORS
      //
      // Build a vector of detection logits for all valid observations at site i.
      // ----------------------------------------------------------------------
      vector[n_obs_site[i]] logit_p_site;
      {
        int idx = 1;
        for (j in 1:J[i]) {
          real lp_cam = logit_p_cam[cam_idx];
          for (k in 1:max_K) {
            if (use_y[i,j,k] == 1) {
              logit_p_site[idx] = lp_cam;
              idx += 1;
            }
          }
          cam_idx += 1;
        }
      }

      // ----------------------------------------------------------------------
      // ZERO-INFLATED N-MIXTURE LIKELIHOOD
      // ----------------------------------------------------------------------

      if (Kmin[i] == 0) {

        // Case: no detections → mixture of structural zero and ecological zero
        target += log_sum_exp(
          log_psi1,
          log_psi0 +
          nmix_site_lpmf(
            segment(y_flat[i], start[i], len[i]) |
            Kmin[i], K,
            max_y_site[i],
            log_lambda_i,
            logit_p_site,
            phi)
        );

      } else {

        // Case: at least one detection → site cannot be structurally unoccupied
        target +=
          log_psi0 +
          nmix_site_lpmf(
            segment(y_flat[i], start[i], len[i]) |
            Kmin[i], K,
            max_y_site[i],
            log_lambda_i,
            logit_p_site,
            phi);
      }
    }
  }
}

