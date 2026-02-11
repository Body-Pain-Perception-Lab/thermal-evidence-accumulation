
// based on codes/comments by Guido Biele, Joseph Burling, Andrew Ellis, and potentially others @ Stan mailing lists

// defining variables 
  data {
  
  int<lower=0> trials;
  int<lower=0> S;
  array[trials] int S_id;
  matrix[S,2] minRT_cold;                        // minimum RT of the observed data
  matrix[S,2] minRT_warm;
  
  int N_alpha;
  int N_delta;
  int N_tau;
  int N_beta;
  
  
  matrix[trials, N_alpha] X_alpha;
  matrix[trials, N_delta] X_delta;
  matrix[trials, N_tau] X_tau;
  matrix[trials, N_beta] X_beta;
  
  
  vector[trials] RT;
  array[trials] int resp;
  array[trials] int quality;
  array[trials] int task;
  
  
  
}
transformed data{
  
  int N = N_alpha+N_delta+N_tau+N_beta; // number of parameters total 
}


parameters {
  
  // hierarchical group level means 
  vector [N] gm;
  // hierarchical group level deviations
  vector<lower = 0>[N]  tau_u;
  // Subject-level estimate matrix 
  matrix[N, S] z_expo;
  // for the cholesky decomposition
  cholesky_factor_corr[N] L_u;
  
  

  
}

transformed parameters{
  
  vector<lower=0>[trials] alpha; // defining empty vectors of of a, z, t, d yhe len of trial 
  vector<lower = 0, upper = 1> [trials] beta;  
  vector[trials] delta;
  vector[trials] tau;

  // 
  // vector<lower=0, upper = minRT_warm>[S] tau_warm;
  // vector<lower=0, upper = minRT_cold>[S] tau_cold;
  // 
  
  matrix[S, N] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';
  
  matrix[S, N] param;
  
  for(n in 1:N){
    param[,n]= gm[n] + indi_dif[,n];
  }
  
  
  
  matrix[S,N_alpha] alpha_p = param[,1:N_alpha];
  
  matrix[S,N_beta] beta_p = param[,(N_alpha+1):(N_alpha+N_beta)];
  
  matrix[S,N_delta] delta_p = param[,(N_alpha+N_beta+1):(N_alpha+N_beta+N_delta)];
  
  matrix[S,N_tau] tau_p = param[,(N_alpha+N_beta+N_delta+1):N];
    
  
  for(n in 1:trials){

    alpha[n] = exp(dot_product(X_alpha[n,], alpha_p[S_id[n],])); //S_id is the participant id 
    
    beta[n] = inv_logit(dot_product(X_beta[n,], beta_p[S_id[n],]));
    
    delta[n] = dot_product(X_delta[n,], delta_p[S_id[n],]);
    
    if(quality[n] == 0 && task[n] == 0) {
      tau[n] = inv_logit(dot_product(X_tau[n,], tau_p[S_id[n],]))*(minRT_warm[S_id[n], 1]);
    }else if(quality[n] == 1 && task[n] == 0){
      tau[n] = inv_logit(dot_product(X_tau[n,], tau_p[S_id[n],]))*(minRT_cold[S_id[n], 1]);
    }else if(quality[n] == 1 && task[n] == 1){
      tau[n] = inv_logit(dot_product(X_tau[n,], tau_p[S_id[n],]))*(minRT_cold[S_id[n], 2]);
    }else{
      tau[n] = inv_logit(dot_product(X_tau[n,], tau_p[S_id[n],]))*(minRT_warm[S_id[n], 2]);  
    }
  }
}

model {
  
  int c;

// priors 
  target += normal_lpdf(gm[1:N_alpha] | 0, 3); //global mean of alpha
  target += normal_lpdf(gm[(N_alpha+1):(N_alpha+N_beta)] | 0, 1); //global mean of beta
  target += normal_lpdf(gm[(N_alpha+N_beta+1):(N_alpha+N_beta+N_delta)] | 0, 3); //global mean of delta
  target += normal_lpdf(gm[(N_alpha+N_beta+N_delta+1):N] | 0, 2); //global mean of tau


  target += std_normal_lpdf(to_vector(z_expo));
  target += normal_lpdf(tau_u | 0, 3);
  
  target += lkj_corr_cholesky_lpdf(L_u | 2);

  
  // trial loop 
    for(n in 1:trials){
      c = resp[n];
      
      if(c == 1){
        target += wiener_lpdf(RT[n] | alpha[n], tau[n], beta[n], delta[n]); 
      } else {
        target += wiener_lpdf(RT[n] | alpha[n], tau[n], 1-beta[n], -delta[n]);
        }
    }
}


generated quantities{

  vector[trials] log_lik;
  log_lik = rep_vector(0.0, trials);

  for(n in 1:trials){

    if(resp[n] == 1){
      log_lik[n] += wiener_lpdf(RT[n] | alpha[n], tau[n], beta[n], delta[n]);
    } else {
      log_lik[n] += wiener_lpdf(RT[n] | alpha[n], tau[n], 1-beta[n], -delta[n]);
      }
  }


}
