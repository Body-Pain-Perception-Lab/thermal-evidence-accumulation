

make_param_recov = function(){
  
  # Parameter recovery for fit_SS Stan model
  # This script loads the fitted model, extracts group-level parameters, and simulates new subjects
  
  library(posterior)
  library(bayesplot)
  library(dplyr)
  library(tidyverse)
  # Load fitted model
  fit_SS <- readRDS(here::here("STAN","taskmodels","fit_SS.rds"))
  
  

  # Extract posterior samples
  group_mean_posterior <- as_draws_df(fit_SS$draws("gm")) %>% select(-contains("."))
  tau_u_posterior <- as_draws_df(fit_SS$draws("tau_u")) %>% select(-contains("."))
  L_u_posterior <- as_draws_df(fit_SS$draws("L_u")) %>% select(-contains("."))
  
  # Number of subjects to simulate per posterior sample
  
  n_trials <- 300
  n_subjects <- 16
  
  # Sample one posterior draw
  draw_idx <- sample(1:4000, 1)
  gms <- group_mean_posterior[draw_idx, ]
  tau_u <- tau_u_posterior[draw_idx, ]
  L_u_vec <- L_u_posterior[draw_idx, ]
  
  # Number of parameters
  N <- ncol(gms)
  
  # Reconstruct the Cholesky factor matrix from the vector
  L_u <- matrix(0, N, N)
  idx <- 1
  for (i in 1:N) {
    for (j in 1:i) {
      L_u[i, j] <- as.numeric(L_u_vec[idx])
      idx <- idx + 1
    }
  }
  
  # Simulate subject-level parameters properly accounting for correlation
  parameters <- data.frame(subjects = 1:n_subjects)
  for (s in 1:n_subjects) {
    # Sample standard normal deviates
    z <- rnorm(N)
    # Apply Cholesky factor and scale by tau_u
    scaled_deviates <- diag(as.numeric(tau_u)) %*% L_u %*% z
    # Add to group means
    subject_params <- as.numeric(gms) + scaled_deviates
    parameters[s, 2:(N+1)] <- subject_params
  }
  
  
  t_data <- do.call(rbind, lapply(1:n_subjects, function(s) {
    data.frame(
      subject     = s,
      rating_dif  = rnorm(n_trials, 0, 1),
      task        = c(rep("FastSS", n_trials), rep("SlowSS", n_trials)),
      response_cw = sample(c("warm","cold"), n_trials, replace = TRUE),
      quality     = sample(c("warm","cold"), n_trials, replace = TRUE),
      area_size   = sample(c(2,3,4,5), n_trials, replace = TRUE)
    )
  }))
  
  trial_ddm = data.frame()
  
  for(s in unique(t_data$subject)){
  
  df = t_data %>% filter(subject == s) %>% 
    mutate(task_num = ifelse(task == "FastSS",0,1),
                                quality_num = ifelse(quality == "warm",0,1)) %>%
    mutate(boundary = exp(parameters[s,2] + parameters[s,3] * task_num),
           
           bias = plogis(
             parameters[s,4]  + parameters[s,5]  * area_size + parameters[s,6]  * quality_num +
               parameters[s,7]  * task_num + parameters[s,8]  * rating_dif +
               parameters[s,9]  * area_size * quality_num + parameters[s,10] * area_size * task_num +
               parameters[s,11] * task_num * quality_num + parameters[s,12] * quality_num * rating_dif +
               parameters[s,13] * task_num * rating_dif +
               parameters[s,14] * quality_num * task_num * area_size +
               parameters[s,15] * quality_num * rating_dif * task_num
           ),
           
           drift =
             parameters[s,16] + parameters[s,17] * area_size + parameters[s,18] * quality_num +
             parameters[s,19] * task_num + parameters[s,20] * rating_dif +
             parameters[s,21] * area_size * quality_num + parameters[s,22] * area_size * task_num +
             parameters[s,23] * task_num * quality_num + parameters[s,24] * quality_num * rating_dif +
             parameters[s,25] * task_num * rating_dif +
             parameters[s,26] * quality_num * task_num * area_size +
             parameters[s,27] * quality_num * rating_dif * task_num
           ,
           
           ndt =
             ifelse(task == "FastSS" & quality == "warm", brms::inv_logit_scaled(parameters[s,28]) * 0.3,
                    ifelse(task == "FastSS" & quality == "cold", brms::inv_logit_scaled(parameters[s,29]) * 0.2,
                           ifelse(task == "SlowSS" & quality == "warm", brms::inv_logit_scaled(parameters[s,30]) * 2.3,
                                  ifelse(task == "SlowSS" & quality == "cold", brms::inv_logit_scaled(parameters[s,31]) * 2.2,
                                         NA))))
    )  # enforce positivity
  
  trial_ddm = rbind(trial_ddm,df)
  }
  
  # get the trial wise repsonses and response times
  trialwisedata = trial_ddm %>% rowwise() %>% mutate(rt_resp = list(RWiener::rwiener(1,boundary,ndt,bias,drift))) %>% unnest("rt_resp")
  
  
  # just for ease having warm the reference
  trialwisedata = trialwisedata %>% mutate(quality = ifelse(quality == "cold","zcold","warm"))
  
  
  
  # boundary separation
  X_alpha = model.matrix(~1+task, data = trialwisedata)
  
  
  # non-decision time
  X_tau = model.matrix(~0+quality:task, data = trialwisedata)
  
  
  X_deltas = c("~ 1 + area_size * quality * task + quality * rating_dif * task")
  
  
  X_betas = c("~ 1 + area_size * quality * task + quality * rating_dif * task")
  
  
  X_delta <- model.matrix(as.formula(X_deltas[1]), data = trialwisedata)
  X_beta <- model.matrix(as.formula(X_betas[1]), data = trialwisedata)
  
  

  
  DDM = cmdstanr::cmdstan_model(here::here("STAN","stanmodels", "stanmodel_task.stan"))
  
  
  data_stan = list(
    trials = nrow(trialwisedata), # How many trials in total 
    S = length(unique(trialwisedata$subject )), # How many subjects in total 
    S_id = as.numeric(as.factor(as.numeric(trialwisedata$subject))), # Index of each ID 
    minRT_cold = as.matrix(trialwisedata %>% filter(quality == "zcold") %>% group_by(subject ,task) %>% summarize(minrt = min(q))  %>% ungroup() %>% select(task,minrt) %>% pivot_wider(names_from = "task",values_from = "minrt") %>% unnest()),
    minRT_warm = as.matrix(trialwisedata %>% filter(quality == "warm") %>% group_by(subject ,task) %>% summarize(minrt = min(q))  %>% ungroup() %>% select(task,minrt) %>% pivot_wider(names_from = "task",values_from = "minrt") %>% unnest()),
    RT = trialwisedata$q,
    quality = trialwisedata %>% mutate(quality = ifelse(quality == "warm", 0, 1)) %>% .$quality,
    task = trialwisedata %>% mutate(task = ifelse(task == "FastSS", 0, 1)) %>% .$task,
    X_alpha = X_alpha,
    X_delta = X_delta,
    X_beta = X_beta,
    X_tau = X_tau,
    N_alpha = ncol(X_alpha), # number of cols = number of parameters 
    N_beta = ncol(X_beta),
    N_delta = ncol(X_delta),
    N_tau = ncol(X_tau),
    resp = trialwisedata %>% mutate(resp = ifelse(resp ==  "lower",0, 1)) %>% .$resp
  )
  
  # Run the model
  fit_models_simp <- DDM$sample(
    data = data_stan,
    chains = 4,
    iter_sampling = 500, 
    iter_warmup = 500, 
    parallel_chains = 4,
    adapt_delta = 0.9,
    refresh = 50,
    init = 0,
    max_treedepth = 12
  )
  
  # extract group level parameters
  
  estimated_group_means = fit_models_simp$summary("gm")
  estimated_between_subject = fit_models_simp$summary("tau_u")
  
  # extract subject level parameters
  
  subj_parameters = fit_models_simp$summary(c("alpha_p","beta_p","delta_p","tau_p")) %>% 
    mutate(
      subject = as.integer(str_extract(variable, "(?<=\\[)\\d+")),                     # first number
      parameter_number = as.integer(str_extract(variable, "(?<=,)\\d+(?=\\])")),    # second number
      variable = str_remove(variable, "\\[.*\\]")                                    # remove brackets and numbers
    ) %>% select(variable,mean,subject,parameter_number) %>%   
    pivot_wider(
    id_cols = subject,                  # rows = subjects
    names_from = c(variable, parameter_number),  # column names = variable + parameter_number
    values_from = mean
  )
  
  # put these in lists to store
  
  simulated_params = list(groupmeans = gms, between_subj_var = tau_u, subj_params = parameters)
  estimated_params = list(groupmeans = estimated_group_means, between_subj_var = estimated_between_subject, subj_params = subj_parameters)
  
  div = fit_models_simp$diagnostic_summary()
  
  return(list(simulated_params,estimated_params,div))
  
}

# use parallel processing to speed up the parameter recovery analysis.

library(furrr)
library(purrr)
library(tidyverse)

plan(multisession, workers = 15)  # Windows-friendly

# test it works

qq = make_param_recov()
plot(as.numeric(estimated_params$groupmeans$mean),as.numeric(simulated_params$groupmeans))

# if the simulating function or the stan models runs into an error just return "Error" for that list

safe_function <- possibly(make_param_recov, otherwise = "Error")

# Run X times in parallel across the 15 workers
results_list <- future_map(1:15, ~ safe_function())

# Save the results

# save.image(here::here("parameter_recovery_result.RData"))



sim <- lapply(seq_along(results_list[1:100]), function(i) {
  df <- as.data.frame(results_list[[i]][[3]])
  df$id <- i                 # add identifier
  df
}) %>%
  bind_rows()

# check convergence fof the models:

hist(sim$num_divergent)
hist(sim$num_max_treedepth)


# Simulated values
sim <- lapply(seq_along(results_list[1:100]), function(i) {
  df <- as.data.frame(results_list[[i]][[1]][[1]])
  df$id <- i                 # add identifier
  df
}) %>%
  bind_rows() %>%
  pivot_longer(
    cols = -id,
    names_to = "variable",
    values_to = "simulated"
  ) 

# Estimated values
esti <- lapply(seq_along(results_list[1:100]), function(i) {
  df <- as.data.frame(results_list[[i]][[2]][[1]])
  df$id <- i
  df
}) %>%
  bind_rows()

# Join by id (and variable if needed)
groupmean = left_join(sim, esti)


# plot simulated vs estimated parameters and save plots for each group level parameter across the DDM parameters

groupmean %>% 
  ggplot(aes(x = simulated, y = mean, ymin = q5,ymax = q95))+
  geom_pointrange(alpha = 0.5)+
  facet_wrap(~variable, scales = "free")+
  geom_abline(col = "red")

groupmean_alpha = groupmean %>% 
  filter(variable %in% c("gm[1]","gm[2]")) %>% 
  mutate(parameter = ifelse(variable == "gm[1]","Intercept","Task"),
         variable = "Alpha") %>% 
  ggplot(aes(x = simulated, y = mean, ymin = q5,ymax = q95))+
  geom_pointrange(alpha = 0.5)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  labs(x = "Simulated",y = "Estimated [95HDI]")+
  facet_wrap(~parameter, scales = "free")+
  geom_abline(col = "red")+
  theme_classic()

groupmean_alpha


ggsave(here::here("alpha.png"), groupmean_alpha, dpi = 300, height = 3, width = 8)


groupmean_beta = groupmean %>% 
  filter(variable %in% c(paste0("gm[",3:(3+11),"]"))) %>% 
  mutate(parameter = ifelse(variable == "gm[3]","Intercept",
                            ifelse(variable == "gm[4]","Area",
                                   ifelse(variable == "gm[5]","Modality",
                                          ifelse(variable == "gm[6]","Task",
                                                 ifelse(variable == "gm[7]","Rating",
                                                        ifelse(variable == "gm[8]","Area:Modality",
                                                               ifelse(variable == "gm[9]","Area:Task",
                                                                      ifelse(variable == "gm[10]","Modality:Task",
                                                                             ifelse(variable == "gm[11]","Modality:Rating",
                                                                                    ifelse(variable == "gm[12]","Task:Rating",
                                                                                           ifelse(variable == "gm[13]","Area:Modality:Task",
                                                                                                  ifelse(variable == "gm[14]","Modality:Task:Rating",NA)))))))))))),
         variable = "Beta",
         parameter = factor(parameter, levels = c(
           "Intercept", "Area", "Modality", "Task", "Rating",
           "Area:Modality", "Area:Task", "Modality:Task", 
           "Modality:Rating", "Task:Rating", "Area:Modality:Task", "Modality:Task:Rating"
         ))) %>% 
  ggplot(aes(x = simulated, y = mean, ymin = q5,ymax = q95))+
  geom_pointrange(alpha = 0.5)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  labs(x = "Simulated",y = "Estimated [95HDI]")+
  facet_wrap(~parameter, scales = "free")+
  geom_abline(col = "red")+
  theme_classic()

groupmean_beta

ggsave(here::here("beta.png"), groupmean_beta, dpi = 300, height = 3*3, width = 8)


groupmean_delta = groupmean %>% 
  filter(variable %in% c(paste0("gm[",(3+12):(15+11),"]"))) %>% 
  mutate(parameter = ifelse(variable == "gm[15]","Intercept",
                            ifelse(variable == "gm[16]","Area",
                                   ifelse(variable == "gm[17]","Modality",
                                          ifelse(variable == "gm[18]","Task",
                                                 ifelse(variable == "gm[19]","Rating",
                                                        ifelse(variable == "gm[20]","Area:Modality",
                                                               ifelse(variable == "gm[21]","Area:Task",
                                                                      ifelse(variable == "gm[22]","Modality:Task",
                                                                             ifelse(variable == "gm[23]","Modality:Rating",
                                                                                    ifelse(variable == "gm[24]","Task:Rating",
                                                                                           ifelse(variable == "gm[25]","Area:Modality:Task",
                                                                                                  ifelse(variable == "gm[26]","Modality:Task:Rating",NA)))))))))))),
         variable = "Delta",
         parameter = factor(parameter, levels = c(
           "Intercept", "Area", "Modality", "Task", "Rating",
           "Area:Modality", "Area:Task", "Modality:Task", 
           "Modality:Rating", "Task:Rating", "Area:Modality:Task", "Modality:Task:Rating"
         ))) %>% 
  ggplot(aes(x = simulated, y = mean, ymin = q5,ymax = q95))+
  geom_pointrange(alpha = 0.5)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  labs(x = "Simulated",y = "Estimated [95HDI]")+
  facet_wrap(~parameter, scales = "free")+
  geom_abline(col = "red")+
  theme_classic()

groupmean_delta

ggsave(here::here("delta.png"), groupmean_delta, dpi = 300, height = 3*3, width = 8)


groupmean_tau = groupmean %>% 
  filter(variable %in% c(paste0("gm[",(27):(30),"]"))) %>% 
  mutate(parameter = ifelse(variable == "gm[27]","Warm_(Task0)",
                            ifelse(variable == "gm[28]","Cold_(Task0)",
                                   ifelse(variable == "gm[29]","Warm_(Task1)",
                                          ifelse(variable == "gm[30]","Cold_(Task1)",NA)))),
         variable = "Delta") %>% 
  ggplot(aes(x = simulated, y = mean, ymin = q5,ymax = q95))+
  geom_pointrange(alpha = 0.5)+
  facet_wrap(~parameter, scales = "free", ncol = 4)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  labs(x = "Simulated",y = "Estimated [95HDI]")+
  geom_abline(col = "red")+
  theme_classic()

groupmean_tau

ggsave(here::here("tau.png"), groupmean_tau, dpi = 300, height = 3, width = 8)
