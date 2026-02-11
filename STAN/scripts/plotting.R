
# funciton used to extract parameters for the predictive checks

extract_params = function(df1,model,dq_draw){
  
  bound = as_draws_df(model$draws(c(c(paste0("alpha[",df1$trial_n,"]"))))) %>% 
    select(-contains(".")) %>% 
    mutate(draw = 1:n()) %>% 
    filter(draw %in% dq_draw) %>% 
    pivot_longer(-draw) %>% 
    mutate(trial = rep(df1$trial_n,length(dq_draw)),
           name = "alpha")
  
  drift = as_draws_df(model$draws(c(c(paste0("delta[",df1$trial_n,"]"))))) %>% 
    select(-contains(".")) %>% 
    mutate(draw = 1:n()) %>% 
    filter(draw %in% dq_draw) %>% 
    pivot_longer(-draw) %>% 
    mutate(trial = rep(df1$trial_n,length(dq_draw)),
           name = "drift")
  
  
  beta = as_draws_df(model$draws(c(c(paste0("beta[",df1$trial_n,"]"))))) %>% 
    select(-contains(".")) %>% 
    mutate(draw = 1:n()) %>% 
    filter(draw %in% dq_draw) %>% 
    pivot_longer(-draw) %>% 
    mutate(trial = rep(df1$trial_n,length(dq_draw)),
           name = "beta")
  
  ndt = as_draws_df(model$draws(c(c(paste0("tau[",df1$trial_n,"]"))))) %>% 
    select(-contains(".")) %>% 
    mutate(draw = 1:n()) %>% 
    filter(draw %in% dq_draw) %>% 
    pivot_longer(-draw) %>% 
    mutate(trial = rep(df1$trial_n,length(dq_draw)),
           name = "ndt")
  
  
  
  preds = rbind(bound,drift,ndt,beta) %>% 
    inner_join(.,df1 %>% select(response_cw,trial_n,task) %>% rename(trial = trial_n)) %>% 
    pivot_wider()
  
  
  preds = preds %>% rowwise() %>% 
    mutate(pred_rt = RWiener::rwiener(1,alpha,ndt,beta,drift)$q)
  
  return(preds)
}



