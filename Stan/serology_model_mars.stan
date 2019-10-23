
data{
    int<lower=0> N;
    int<lower=0> N_rcva;
    int<lower=1, upper=N> i_rcva[N_rcva];
    int<lower=0> y_rcva[N_rcva];
    int<lower=0> n_rcva[N_rcva];
    int<lower=0> N_rhd1;
    int<lower=1, upper=N> i_rhd1[N_rhd1];
    int<lower=0> y_rhd1[N_rhd1];
    int<lower=0> n_rhd1[N_rhd1];
    int<lower=0> N_rhd2;
    int<lower=1, upper=N> i_rhd2[N_rhd2];
    int<lower=0> y_rhd2[N_rhd2];
    int<lower=0> n_rhd2[N_rhd2];
    int<lower=0> nsites;
    int<lower=0> season[N];
    int<lower=0> start[nsites];
    int<lower=0> end[nsites];
}


parameters {
  matrix[N, 3] logit_prev;
  matrix[3, 3] beta;
  real eta[3];
  matrix[3, nsites] site;
  cholesky_factor_corr[3] Lcorr;
  vector<lower=0>[3] sigma_proc;
  
  vector[nsites] ns_rcva;
  vector[nsites] ns_rhd1;
  vector[nsites] ns_rhd2;
  
}

transformed parameters {

  matrix[N, 3] logit_prev_mu;
 
  for(k in 1:nsites) {

    logit_prev_mu[start[k],1] = ns_rcva[k];
    logit_prev_mu[start[k],2] = ns_rhd1[k];
    logit_prev_mu[start[k],3] = ns_rhd2[k];
    
      for(t in start[k]:(end[k]-1)) {
        for(m in 1:3){
          logit_prev_mu[t+1, m] = site[m, k] + logit_prev[t] * beta[m]' + eta[m]*season[t]; 
        } //m
      } //t
    } //k
    
    
}

    
model{

for(k in 1:nsites){
    for(t in start[k]:end[k]) {
      logit_prev[t] ~ multi_normal_cholesky(logit_prev_mu[t], diag_pre_multiply(sigma_proc, Lcorr)); 
    }
  }
  
  ns_rcva ~ normal(0, 5);
  ns_rhd1 ~ normal(0, 5);
  ns_rhd2 ~ normal(0, 5);
 
 for(i in 1:3) 
  site[i] ~ normal(0, 5);

  sigma_proc ~ student_t(4, 0, 1);

  for(i in 1:3) 
    beta[i] ~ normal(0, 5);
   
  eta ~ normal(0, 5);
 
  Lcorr ~ lkj_corr_cholesky(1); 
 
   // Likelihood
    for (i in 1:N_rcva) 
      y_rcva[i] ~ binomial_logit(n_rcva[i], logit_prev[i_rcva[i],1]);
    for(i in 1:N_rhd1)
      y_rhd1[i] ~ binomial_logit(n_rhd1[i], logit_prev[i_rhd1[i],2]);
    for(i in 1:N_rhd2)
      y_rhd2[i] ~ binomial_logit(n_rhd2[i], logit_prev[i_rhd2[i],3]);
       
}

generated quantities {
  vector[N] prev_rcva;
  vector[N] prev_rhd1;
  vector[N] prev_rhd2;
  matrix[3, 3] Omega;
  matrix[3, 3] Sigma;
  
  prev_rcva = inv_logit(logit_prev[,1]);
  prev_rhd1 = inv_logit(logit_prev[,2]);
  prev_rhd2 = inv_logit(logit_prev[,3]);
  
  Omega = multiply_lower_tri_self_transpose(Lcorr);
  Sigma = quad_form_diag(Omega, sigma_proc); 
}

