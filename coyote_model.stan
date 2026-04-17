functions{
real lp_occu(array[] int y, real logit_psi, vector logit_p, int Kmin){
  real out;
  out = log_inv_logit(logit_psi) + bernoulli_logit_lpmf(y | logit_p);
  if(Kmin == 1){
    return out;
  }
  return log_sum_exp(out, log1m_inv_logit(logit_psi));
}
vector get_loglik_occu(array[] int y, int M, array[,] int J, array[,] int si, vector logit_psi,
                  vector logit_p, array[] int Kmin){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_occu(y[si[i,1]:si[i,2]], logit_psi[i],
                     logit_p[si[i,1]:si[i,2]], Kmin[i]);
   }
  return out;
}
real lp_rn(array[] int y, real log_lambda, vector logit_r, int J, int K, int Kmin){
  int numN = K - Kmin + 1;
  vector[J] q = 1 - inv_logit(logit_r);
  vector[numN] lp;
  vector[J] p;
  int N;
  for (i in 1:numN){
    N = K - i + 1;
    for (j in 1:J) p[j] = 1 - q[j]^N;
    lp[i] = poisson_log_lpmf(N | log_lambda) + bernoulli_lpmf(y | p);
  }
  return log_sum_exp(lp);
}
vector get_loglik_rn(array[] int y, int M, array[,] int J, array[,] int si, vector log_lambda,
                     vector logit_p, int K, array[] int Kmin){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_rn(y[si[i,1]:si[i,2]], log_lambda[i], logit_p[si[i,1]:si[i,2]],
                   J[i,1], K, Kmin[i]);
  }
  return out;
}
real lp_pcount_pois(array[] int y, real log_lambda, vector logit_p, int K, int Kmin){
  real fac = 1;
  real ff = exp(log_lambda) * prod(1 - inv_logit(logit_p));
  real N;
  real ky;
  int numN = K - Kmin;
  for (i in 1:numN){
    N = K - i + 1;
    ky = 1;
    for (j in 1:size(y)){
      ky *= N / (N - y[j]);
    }
    fac = 1 + fac * ff * ky / N;
  }
  return  poisson_log_lpmf(Kmin | log_lambda) +
          binomial_logit_lpmf(y | Kmin, logit_p) +
          log(fac);
}
vector get_loglik_pcount(array[] int y, int M, array[,] int J, array[,] int si, vector log_lambda,
                         vector logit_p, int z_dist, real beta_scale, int K, array[] int Kmin){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_pcount_pois(y[si[i,1]:si[i,2]], log_lambda[i],
                            logit_p[si[i,1]:si[i,2]], K, Kmin[i]);
  }
  return out;
}
real int_halfnorm_point(real sigma, real a, real b){
  real s2 = pow(sigma, 2);
  return s2 * ((1 - exp(-b*b / (2*s2))) - (1 - exp(-a*a / (2*s2))));
}
real int_halfnorm_line(real sigma, real a, real b){
  real den = sqrt(2) * sigma;
  return sqrt(pi()/2) * sigma * (erf(b/den) - erf(a/den));
}
real int_halfnorm(real log_sigma, real a, real b, int point){
  real out;
  real sigma = exp(log_sigma);
  if(point){
    out = int_halfnorm_point(sigma, a, b);
  } else{
    out = int_halfnorm_line(sigma, a, b);
  }
  return out;
}
real int_negexp_point(real rate, real a, real b){
  return rate * exp(-a/rate) * (a+rate) -
         rate * exp(-b/rate) * (b+rate);
}
real int_negexp_line(real rate, real a, real b){
  return rate * (exp(-a/rate) - exp(-b/rate));
}
real int_negexp(real log_rate, real a, real b, int point){
  real out;
  real rate = exp(log_rate);
  if(point){
    out = int_negexp_point(rate, a, b);
  } else{
    out = int_negexp_line(rate, a, b);
  }
  return out;
}
//real p_hazard_line(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i){
//  return(1 - exp(-1 * pow((x/theta[1]), (-1*theta[2]))));
//}
real p_hazard_line(real x, array[] real theta){
  return(1 - exp(-1 * pow((x/theta[1]), (-1*theta[2]))));
}
//real p_hazard_point(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i){
//  return((1 - exp(-1 * pow(x/theta[1], -1*theta[2]))) * x);
//}
real p_hazard_point(real x, array[] real theta){
  return((1 - exp(-1 * pow(x/theta[1], -1*theta[2]))) * x);
}
real trap_rule_line(array[] real theta, real a, real b){
  int n = 100;
  real h = (b - a) / n;
  real int_sum = 0;
  for (i in 1:(n-1)){
    int_sum += p_hazard_line(a+i*h, theta);
  }
  return h/2 * (p_hazard_line(a, theta) + 2*int_sum + p_hazard_line(b, theta));
}
real trap_rule_point(array[] real theta, real a, real b){
  int n = 100;
  real h = (b - a) / n;
  real int_sum = 0;
  for (i in 1:(n-1)){
    int_sum += p_hazard_point(a+i*h, theta);
  }
  return( h/2 * (p_hazard_point(a, theta) + 2*int_sum + p_hazard_point(b, theta)) );
}
real int_hazard(real log_shape, real log_scale, real a, real b, int point){
  real out;
  real shape = exp(log_shape);
  real scale = exp(log_scale);
  array[2] real theta;
  theta[1] = shape;
  theta[2] = scale;
  if(point){
    out = trap_rule_point(theta, a, b);
  } else{
    out = trap_rule_line(theta, a, b);
  }
  return out;
}
real prob_dist(real par1, real par2, int keyfun, real a, real b, int point){
  real out;
  if(keyfun == 0){
    out = int_halfnorm(par1, a, b, point);
  } else if(keyfun == 1){
    out = int_negexp(par1, a, b, point);
  } else if(keyfun == 2){
    out = int_hazard(par1, par2, a, b, point);
  }
  return out;
}
real lp_distsamp(array[] int y, vector db, real log_lambda, real par1, real par2,
                 int point, int keyfun, vector conv_const){
  real lam = exp(log_lambda);
  int J = num_elements(db) - 1;
  real loglik = 0.0;
  real cp;
  for (j in 1:J){
    cp = prob_dist(par1, par2, keyfun, db[j], db[j+1], point);
    cp = cp * conv_const[j];
    loglik += poisson_lpmf(y[j] | lam * cp);
  }
  return loglik;
}
vector get_loglik_distsamp(array[] int y, int M, vector db, array[,] int si,
                           vector log_lambda, vector trans_par1, int z_dist,
                           real trans_par2, int point, int keyfun, vector conv_const){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_distsamp(y[si[i,1]:si[i,2]], db, log_lambda[i], trans_par1[i],
                         trans_par2, point, keyfun, conv_const[si[i,1]:si[i,2]]);
  }
  return out;
}
//Removal sampling
//p can be any length > 1
vector pi_removal(vector p){
  int J = num_elements(p);
  vector[J] pi_out;
  pi_out[1] = p[1];
  for (j in 2:J){
    pi_out[j] = pi_out[j-1] / p[j-1] * (1-p[j-1]) * p[j];
  }
  return pi_out;
}
//Double observer
//p must have 2 elements
vector pi_double(vector p){
  vector[3] pi_out;
  pi_out[1] = p[1] * (1 - p[2]);
  pi_out[2] = p[2] * (1 - p[1]);
  pi_out[3] = p[1] * p[2];
  return pi_out;
}
vector pi_fun(int pi_type, vector p, int J){
  vector[J] out;
  if(pi_type == 0){
    out = pi_double(p);
  } else if(pi_type == 1){
    out = pi_removal(p);
  } else {
    reject("Invalid pi function type");
  }
  return out;
}
real lp_multinomPois(array[] int y, real log_lambda, vector logit_p, int pi_type){
  real loglik = 0.0;
  real lam = exp(log_lambda);
  int J = num_elements(y);
  int np = num_elements(logit_p);
  vector[np] p;
  vector[J] cp;
  for (i in 1:np){
    p[i] = inv_logit(logit_p[i]);
  }
  cp = pi_fun(pi_type, p, J);
  for (j in 1:J){
    loglik += poisson_lpmf(y[j] | lam * cp[j]);
  }
  return loglik;
}
vector get_loglik_multinomPois(array[] int y, int M, array[,] int si, vector log_lambda,
                               vector logit_p, int pi_type){
  vector[M] out;
  int J = num_elements(logit_p) %/% M; // use %/% in future version of stan
  int pstart = 1;
  int pend;
  for (i in 1:M){
    pend = pstart + J - 1;
    out[i] = lp_multinomPois(y[si[i,1]:si[i,2]], log_lambda[i],
                             logit_p[pstart:pend], pi_type);
    pstart += J;
  }
  return out;
}
vector ttd_prob_exp(vector y, vector log_lam, array[] int delta){
  int J = num_elements(y);
  vector[J] e_lamt;
  real lam;
  for (j in 1:J){
    lam = exp(log_lam[j]);
    e_lamt[j] = pow(lam, delta[j]) * exp(-lam*y[j]);
  }
  return e_lamt;
}
vector ttd_prob_weib(vector y, vector log_lam, array[] int delta, real log_k){
  int J = num_elements(y);
  vector[J] e_lamt;
  real k = exp(log_k);
  real lam;
  for (j in 1:J){
    lam = exp(log_lam[j]);
    e_lamt[j] = pow(k*lam*pow(lam*y[j], (k-1)), delta[j]) *
                exp(-1*pow(lam*y[j], k));
  }
  return e_lamt;
}
real lp_occuTTD(vector y, real logit_psi, vector log_lam,
                real log_k, array[] int delta, int ydist){
  int J = num_elements(y);
  vector[J] e_lamt;
  real psi = inv_logit(logit_psi);
  real lik;
  if(ydist == 1){ //exponential
    e_lamt = ttd_prob_exp(y, log_lam, delta);
  } else if(ydist == 3){ //weibull
    e_lamt = ttd_prob_weib(y, log_lam, delta, log_k);
  }
  lik = psi * prod(e_lamt) + (1-psi) * (1-max(delta));
  return log(lik);
}
vector get_loglik_occuTTD(vector y, int M, array[,] int si, vector logit_psi,
                          vector log_lam, real log_k, array[] int delta, int ydist){
  vector[M] out;
  for (i in 1:M){
    out[i] = lp_occuTTD(y[si[i,1]:si[i,2]], logit_psi[i], log_lam[si[i,1]:si[i,2]],
                        log_k, delta[si[i,1]:si[i,2]], ydist);
  }
  return out;
}
real lp_single_prior(vector x, int dist, row_vector pars1,
                     row_vector pars2, row_vector pars3){
  real out = 0.0;
  if(dist == 1){
    out += normal_lpdf(x | pars1, pars2);
  } else if(dist == 2){
    out += uniform_lpdf(x | pars1, pars2);
  } else if(dist == 3){
    out += student_t_lpdf(x | pars1, pars2, pars3);
  } else if(dist == 4){
    out += logistic_lpdf(x | pars1, pars2);
  } else if(dist == 5){
    out += gamma_lpdf(x | pars1, pars2);
  } else if(dist == 6){
    out += double_exponential_lpdf(x | pars1, pars2);
  }
  return out;
}
real lp_priors(vector beta, array[] int dist, matrix pars){
  int idx;
  real out = 0.0;
  int nb = num_elements(beta);
  if(nb == 0) return out;
  idx = dist[1] == 0 ? 1 : 2;
  // intercept prior, if intercept exists
  if(dist[1] != 0){
    out += lp_single_prior(beta[1:1], dist[1], pars[1,1:1],
                          pars[2,1:1], pars[3,1:1]);
  }
  // regression coefficients priors, if there are any
  if(dist[2] != 0){
    out += lp_single_prior(beta[idx:nb], dist[2], pars[1,idx:nb],
                          pars[2,idx:nb], pars[3,idx:nb]);
  }
  return out;
}
real lp_random_prior(int has_random, int n_group_vars, vector b,
                     array[] int n_random, vector sigma, int dist, matrix pars){
  int idx = 1;
  real out = 0;
  int par_idx = cols(pars);
  row_vector[n_group_vars] rep_par1 = rep_row_vector(pars[1,par_idx], n_group_vars);
  row_vector[n_group_vars] rep_par2 = rep_row_vector(pars[2,par_idx], n_group_vars);
  row_vector[n_group_vars] rep_par3 = rep_row_vector(pars[3,par_idx], n_group_vars);
  if(has_random){
    out += lp_single_prior(sigma, dist, rep_par1, rep_par2, rep_par3);
    for (i in 1:n_group_vars){
      out += normal_lpdf(b[idx:(n_random[i]+idx-1)] | 0, sigma[i]);
      idx += n_random[i];
    }
  }
  return out;
}
}
data{
//Basic data setup for single-season model
int model_code;
int M;
int T;
int Tsamp_size;
array[Tsamp_size] int Tsamp;
int R;
array[M,T] int J;
array[R] int y;
array[M, 6] int si;
int K;
array[M,T] int Kmin;
int y_dist;
int z_dist;
int n_aux1;
int n_aux2;
int n_aux3;
array[n_aux1] int aux1; //Used for various auxiliary data
vector[n_aux2] aux2;
vector[n_aux3] aux3;
int has_random_state;
int has_random_det;
int n_obs_state;
int n_obs_det;
int n_fixed_state;
int n_fixed_det;
int n_group_vars_state;
int n_group_vars_det;
array[has_random_state ? n_group_vars_state : 1] int n_random_state;
array[has_random_det ? n_group_vars_det: 1] int n_random_det;
matrix[n_obs_state, n_fixed_state] X_state;
matrix[n_obs_det, n_fixed_det] X_det;
vector[n_obs_state] offset_state;
vector[n_obs_det] offset_det;
array[5] int Zdim_state;
vector[Zdim_state[3]] Zw_state;
array[Zdim_state[4]] int Zv_state;
array[Zdim_state[5]] int Zu_state;
array[5] int Zdim_det;
vector[Zdim_det[3]] Zw_det;
array[Zdim_det[4]] int Zv_det;
array[Zdim_det[5]] int Zu_det;
// Stuff for custom priors
array[3] int prior_dist_state;
array[3] int prior_dist_det;
array[3] int prior_dist_shape;
array[3] int prior_dist_scale;
matrix[3, (n_fixed_state+1)] prior_pars_state;
matrix[3, (n_fixed_det+1)] prior_pars_det;
matrix[3, 2] prior_pars_shape;
matrix[3, 2] prior_pars_scale;
}
transformed data{
int include_scale;
int include_shape;
include_scale = prior_dist_scale[1] == 0 ? 0 : 1;
include_shape = prior_dist_shape[1] == 0 ? 0 : 1;
}
parameters{
//Basic parameters for single-season model
vector[n_fixed_state] beta_state;
vector[n_fixed_det] beta_det;
vector[include_scale] beta_scale; //Used in NB models
vector[include_shape] beta_shape; //used in Weibull models
vector<lower=0>[n_group_vars_state] sigma_state;
vector<lower=0>[n_group_vars_det] sigma_det;
vector[sum(n_random_state)] b_state;
vector[sum(n_random_det)] b_det;
}
transformed parameters{
vector[M] lp_state;
vector[n_obs_det] lp_det;
vector[M] log_lik;
real log_scale;
real log_shape;
lp_state = X_state * beta_state + offset_state;
lp_det = X_det * beta_det + offset_det;
if(has_random_state){
  lp_state = lp_state +
              csr_matrix_times_vector(Zdim_state[1], Zdim_state[2], Zw_state,
                                      Zv_state, Zu_state, b_state);
}
if(has_random_det){
  lp_det = lp_det +
            csr_matrix_times_vector(Zdim_det[1], Zdim_det[2], Zw_det,
                                    Zv_det, Zu_det, b_det);
}
log_scale = 0;
if(include_scale){
  log_scale = beta_scale[1];
}
log_shape = 0;
if(include_shape){
  log_shape = beta_shape[1];
}
if(model_code == 0){
  log_lik = get_loglik_occu(y, M, J, si, lp_state, lp_det, Kmin[,1]);
} else if(model_code == 1){
  log_lik = get_loglik_rn(y, M, J, si, lp_state, lp_det, K, Kmin[,1]);
} else if(model_code == 2){
  log_lik = get_loglik_pcount(y, M, J, si, lp_state, lp_det, z_dist,
                              log_scale, K, Kmin[,1]);
} else if(model_code == 4){
  log_lik = get_loglik_distsamp(y, M, aux2, si, lp_state, lp_det, z_dist,
                                log_scale, aux1[1], y_dist, aux3);
} else if(model_code == 5){
  log_lik = get_loglik_multinomPois(y, M, si, lp_state, lp_det, y_dist);
} else if(model_code == 6){
  log_lik = get_loglik_occuTTD(aux2, M, si, lp_state, lp_det, log_shape,
                               aux1, y_dist);
}
}
model{
// Fixed effects priors
target += lp_priors(beta_state, prior_dist_state, prior_pars_state);
target += lp_priors(beta_det, prior_dist_det, prior_pars_det);
target += lp_priors(beta_scale, prior_dist_scale, prior_pars_scale);
target += lp_priors(beta_shape, prior_dist_shape, prior_pars_shape);
// Random effects priors
target += lp_random_prior(has_random_state, n_group_vars_state, b_state,
                          n_random_state, sigma_state, prior_dist_state[3],
                          prior_pars_state);
target += lp_random_prior(has_random_det, n_group_vars_det, b_det,
                          n_random_det, sigma_det, prior_dist_det[3],
                          prior_pars_det);
target += sum(log_lik);
}
