
data{
    int<lower=0> N;
    int<lower=-1> y[N, 3];
    int<lower=0> ncounts[N];
    int<lower=0> nsites;
    int<lower=0, upper=1> season[N];
    int<lower=0, upper=1> rhdv2[N];
    int<lower=0, upper=1> Rel[N];
    vector<lower=0>[N] TL;
    int<lower=0> start[nsites];
    int<lower=0> end[nsites];
    int<lower=0> site[N];
    int<lower=-1> max_y[N];
    int<lower=1> K5[nsites];
}
    
transformed data {
  int<lower=0> K[N]; 
  for (i in 1:N){
    K[i] = max_y[i]*3 + 100; //upper limit of integration
  }

}

parameters {

  vector[nsites] beta_raw;
  vector[nsites] gam_raw;
  real eta;
  real<lower=0> sigma_k5;
  real r_mean;
  real<lower=0> sigma_site;
  real<lower=0> sigma_proc;
  vector[N] p_raw;
  vector[nsites] site_r_raw;
  vector[nsites] nstart;
  real mu_gam;
  real<lower=0> sigma_gam;
  vector[2] mu_k5;
  real mu_b;
  real<lower=0> sigma_b;
  vector[nsites] delta_raw;
  vector[nsites] mu_p;
  real<lower=0> sigma_p;
  vector[N] r_rabbits;
}

transformed parameters {
  
  vector[N] mu_rabbits;
  vector[N] mu_counts;
  vector[N] roi;
  vector[N] mu_p_site;
  vector[nsites] site_r;
  vector[nsites] gam;
  vector[nsites] beta;
  vector[nsites] mu_k5_site;
  vector[nsites] delta;
  vector[N] logit_p;
  
  
  mu_k5_site = mu_k5[K5];
  mu_p_site = mu_p[site];
 
 // non-centered parameterisation for the following
  delta = mu_k5_site + sigma_k5 * delta_raw;
  site_r = r_mean + sigma_site * site_r_raw;
  gam = mu_gam + sigma_gam * gam_raw;
  beta = mu_b + sigma_b * beta_raw; 
  logit_p = mu_p_site + sigma_p * p_raw;

  
  for(k in 1:nsites) {
   mu_rabbits[start[k]] = nstart[k];
   roi[start[k]] = 0;
   
      for(t in (start[k]+1):end[k]) {
        roi[t] =  mu_rabbits[t-1]*beta[k] + eta*season[t] +
                                    gam[k]*rhdv2[t] +
                                    delta[k]*Rel[t] + 
                                    site_r[k];
                                 
        mu_rabbits[t] = mu_rabbits[t-1] + r_rabbits[t]; 
      }
  }
  
  mu_counts = mu_rabbits + TL;
}

    
model{

  r_rabbits ~ normal(roi, sigma_proc);
  site_r_raw ~ std_normal();
  delta_raw ~ std_normal();
  beta_raw ~ std_normal();
  gam_raw ~ std_normal();
  p_raw ~ std_normal();
  
  nstart ~ normal(0, 5);
  
  r_mean ~ normal(0, 5);
  sigma_proc ~ student_t(4, 0, 1);
  sigma_site ~ student_t(4, 0, 1);
  
  eta ~ normal(0, 5);
  mu_gam ~ normal(0, 5);
  sigma_gam ~ student_t(4, 0, 1);
  mu_b ~ normal(0, 5);
  sigma_b ~ student_t(4, 0, 1);
  
  mu_k5 ~ normal(0, 5);
  sigma_k5 ~ student_t(4, 0, 1);

  mu_p ~ normal(0, 1.6);
  sigma_p ~ student_t(4, 0, 1);
  
  // Likelihood
    for (i in 1:N) {
      if(ncounts[i] > 0) {
      vector[K[i] - max_y[i] + 1] lp;
        for (j in 1:(K[i] - max_y[i] + 1))
          lp[j] = poisson_log_lpmf(max_y[i] + j - 1 | mu_counts[i])
            + binomial_logit_lpmf(y[i, 1:ncounts[i]] | max_y[i] + j - 1, logit_p[i]);
        target += log_sum_exp(lp);
      }
    }
}

generated quantities {
  vector[nsites] EQB;
  vector[nsites] EQA;
  vector[nsites] DE;
  vector<lower=0, upper=1>[N] p;
   
  EQB = -(site_r + eta) ./ beta;
  EQA = -(site_r + eta + gam) ./ beta;
  DE = exp(EQA - EQB) - 1;
  p = inv_logit(logit_p);
}

