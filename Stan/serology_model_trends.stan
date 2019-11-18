
data{
    int<lower=0> N;
    int<lower=0> NT;
    int<lower=-1> y_rcva_j[N];
    int<lower=-1> n_rcva_j[N];
    int<lower=-1> y_rcva_a[N];
    int<lower=-1> n_rcva_a[N];
    int<lower=-1> y_rhd1_j[N];
    int<lower=-1> n_rhd1_j[N];
    int<lower=-1> y_rhd1_a[N];
    int<lower=-1> n_rhd1_a[N];
    int<lower=-1> y_rhd2_j[N];
    int<lower=-1> n_rhd2_j[N];
    int<lower=-1> y_rhd2_a[N];
    int<lower=-1> n_rhd2_a[N];
    int<lower=0> time_ind[N];
    int<lower=0> nsites;
    int<lower=0> site_ind[N];
}


parameters {
  matrix[6, NT] prev_raw;
  
  vector<lower=0>[6] sigma_proc;
  vector[6] nstart;
  vector[6] A;
  row_vector<lower=0>[nsites] sigma_site;
  matrix[6, N] sigma_raw;
}

transformed parameters {
  matrix[6, NT] logit_prev;
  matrix[6, NT] logit_prev_mu;
  matrix[6, N] eps;

    for(j in 1:6) 
      logit_prev_mu[j,1] = nstart[j];
      
    for(m in 1:6){
      for(t in 1:(NT-1)) {
          logit_prev_mu[m, t+1] = logit_prev_mu[m, t] +  A[m]; 
        } //m
      } //t
    
      for(k in 1:6) {
        eps[k] = sigma_site[site_ind] .* sigma_raw[k];
        logit_prev[k] = logit_prev_mu[k] + sigma_proc[k] * prev_raw[k];
      }
}

    
model{

  A ~ normal(0, 5);
  nstart ~ normal(0, 1.6);
  sigma_proc ~ student_t(4, 0, 1);
  sigma_site ~ student_t(4, 0, 1);
  to_vector(sigma_raw) ~ std_normal();
  to_vector(prev_raw) ~ std_normal();
  
   // Likelihood
   for(i in 1:N){
     if(n_rcva_j[i] > 0)
        y_rcva_j[i] ~ binomial_logit(n_rcva_j[i], logit_prev[1,time_ind[i]] + eps[1,i]);
     if(n_rcva_a[i] > 0)  
        y_rcva_a[i] ~ binomial_logit(n_rcva_a[i], logit_prev[2,time_ind[i]] + eps[2,i]);
     if(n_rhd1_j[i] > 0)
        y_rhd1_j[i] ~ binomial_logit(n_rhd1_j[i], logit_prev[3,time_ind[i]] + eps[3,i]);
     if(n_rhd1_a[i] > 0)  
        y_rhd1_a[i] ~ binomial_logit(n_rhd1_a[i], logit_prev[4,time_ind[i]] + eps[4,i]);
     if(n_rhd2_j[i] > 0)
        y_rhd2_j[i] ~ binomial_logit(n_rhd2_j[i], logit_prev[5,time_ind[i]] + eps[5,i]);
     if(n_rhd2_a[i] > 0)  
        y_rhd2_a[i] ~ binomial_logit(n_rhd2_a[i], logit_prev[6,time_ind[i]] + eps[6,i]);
   }
       
}

generated quantities {
  row_vector[NT] mean_rcva_j;
  row_vector[NT] mean_rcva_a;
  row_vector[NT] mean_rhd1_j;
  row_vector[NT] mean_rhd1_a;
  row_vector[NT] mean_rhd2_j;
  row_vector[NT] mean_rhd2_a;
  
  row_vector[N] prev_rcva_j;
  row_vector[N] prev_rcva_a;
  row_vector[N] prev_rhd1_j;
  row_vector[N] prev_rhd1_a;
  row_vector[N] prev_rhd2_j;
  row_vector[N] prev_rhd2_a;
  
  mean_rcva_j = inv_logit(logit_prev[1]);
  mean_rcva_a = inv_logit(logit_prev[2]);
  mean_rhd1_j = inv_logit(logit_prev[3]);
  mean_rhd1_a = inv_logit(logit_prev[4]);
  mean_rhd2_j = inv_logit(logit_prev[5]);
  mean_rhd2_a = inv_logit(logit_prev[6]);
  
  prev_rcva_j = inv_logit(logit_prev[1,time_ind] + eps[1,]);
  prev_rcva_a = inv_logit(logit_prev[2,time_ind] + eps[2,]);
  prev_rhd1_j = inv_logit(logit_prev[3,time_ind] + eps[3,]);
  prev_rhd1_a = inv_logit(logit_prev[4,time_ind] + eps[4,]);
  prev_rhd2_j = inv_logit(logit_prev[5,time_ind] + eps[5,]);
  prev_rhd2_a = inv_logit(logit_prev[6,time_ind] + eps[6,]);
}

