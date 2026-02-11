
## This function in extracts useful comparisons from the two STAN-models. some of the comparisons are not returned by the model but live within the function

get_reporting  = function(df){

  
  #########################
  ### spatial summation ###
  #########################
  
  # extract data for NDT
  
  df_fastSS = df %>% 
    filter(ID != 1) %>% 
    group_by(ID, quality) %>%
    filter(corrected_RT > 0.2) %>%        # ascending: fastest first
    ungroup() %>% 
    dplyr::select(ID,config,area_size,corrected_RT,quality,response_cw,rating,task,trial_index) %>% 
    arrange(trial_index) %>% 
    select(-trial_index) %>% 
    drop_na()
  
  df_fastSS = df_fastSS %>% group_by(ID,quality,area_size) %>% mutate(rating_dif = rating-mean(rating)) %>% ungroup()
  df_fastSS = df_fastSS %>% mutate(quality = ifelse(quality == "cold","zcold","warm"))
  
  #minimum response times
  min_cold = as.matrix(df_fastSS %>% filter(quality == "zcold") %>% group_by(ID,task) %>% summarize(minrt = min(corrected_RT))  %>% ungroup() %>% select(task,minrt) %>% pivot_wider(names_from = "task",values_from = "minrt") %>% unnest())
  min_warm = as.matrix(df_fastSS %>% filter(quality == "warm") %>% group_by(ID,task) %>% summarize(minrt = min(corrected_RT))  %>% ungroup() %>% select(task,minrt) %>% pivot_wider(names_from = "task",values_from = "minrt") %>% unnest())
  
  #model
  fit_SS <- readRDS(here::here("STAN","taskmodels","fit_SS.rds"))

  # boundary separation
  X_alpha = model.matrix(~1+task, data = df_fastSS)
  # non-decision time
  X_tau = model.matrix(~0+quality*task, data = df_fastSS)
  #drit rate
  X_deltas = c("~ 1 + area_size * quality * task + quality * rating_dif * task")
  #biases
  X_betas = c("~ 1 + area_size * quality * task + quality * rating_dif * task")
  
  
  #names of the parametesr
  names_delta = colnames(model.matrix(as.formula(X_deltas[1]), data = df_fastSS))
  names_beta = colnames(model.matrix(as.formula(X_betas[1]), data = df_fastSS))
  names_ndt = colnames(X_tau)
  names_alpha = colnames(X_alpha)
  
  # extract the pure parameters from the model for the drift (not used as we need the magnitude)
  drift = as_draws_df(fit_SS$draws("gm")) %>% 
    select(-contains(".")) %>% 
    rename_with(~c(paste0("alpha_",names_alpha),paste0("beta_",names_delta),paste0("delta_",names_delta), paste0("ndt_",names_ndt))) %>% 
    select(contains("delta")) %>% 
    mutate(draw = 1:n()) %>% 
    pivot_longer(-draw) %>% group_by(name) %>% 
    summarize(mean = mean(value), 
              q5 = HDInterval::hdi(value)[1],
              q95 = HDInterval::hdi(value)[2],
              pseudo_p_val(value))
    
  
  # Extract the draws from the drift
  delta_draws <- as_draws_df(fit_SS$draws("gm")) %>% 
    select(-contains(".")) %>% 
    rename_with(~c(
      paste0("alpha_", colnames(X_alpha)),
      paste0("beta_", names_delta),
      paste0("delta_", names_delta),
      paste0("ndt_", names_ndt)
    )) %>% 
    select(starts_with("delta_")) %>% 
    as.matrix()  # dim: [4000 x K]
  

  # Extract the design matrix
  X_delta = data.frame(model.matrix(as.formula(X_deltas[1]), data = df_fastSS %>% 
                                      select(area_size,quality,task) %>% distinct() %>% 
                                      mutate(rating_dif = list(c(0,1))) %>% unnest()))
    

  # back calculate trialwise drift rates
  b = as.matrix(X_delta) %*% t(delta_draws) 
  
  #combine back with the data to have the conditions of the task
  b_long <- as.data.frame(b) %>%
    mutate(obs_id = row_number()) %>%
    pivot_longer(
      cols = -obs_id,
      names_to = "draw",
      values_to = "prediction"
    )
  
  
  b_long$draw <- as.integer(gsub("V", "", b_long$draw))
  
  conditions <- df_fastSS %>%
    select(area_size, quality, task) %>%
    distinct() %>%
    mutate(rating_dif = list(c(0,1))) %>% unnest() %>% 
    mutate(obs_id = row_number())
  
  # Here we combine and flip the warm predictions in order to look at the magnitudes
  b_combined <- b_long %>%
    left_join(conditions, by = "obs_id") %>% 
    mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
  
  
  # three way interaction
  SS_SS_threeway_drift = b_combined %>% filter(rating_dif == 0 & (area_size == 1 | area_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
    pivot_wider(names_from = c(quality,area_size,task), values_from = prediction) %>% 
    mutate(three_way = ((warm_2_FastSS - warm_1_FastSS) - (zcold_2_FastSS - zcold_1_FastSS)) - ((warm_2_SlowSS - warm_1_SlowSS) - (zcold_2_SlowSS - zcold_1_SlowSS))) %>% 
    summarize(mean = mean(three_way), 
              q5 = HDInterval::hdi(three_way)[1],
              q95 = HDInterval::hdi(three_way)[2],
              pseudo_p_val(three_way))
  
  # two way
  SS_SS_toway_drift =b_combined %>% filter(rating_dif == 0 & (area_size == 1 | area_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
    pivot_wider(names_from = c(quality,area_size), values_from = prediction) %>% 
    group_by(task) %>%
    summarize(mean = mean((warm_2 - warm_1) - (zcold_2 - zcold_1)), 
              q5 = HDInterval::hdi((warm_2 - warm_1) - (zcold_2 - zcold_1))[1],
              q95 = HDInterval::hdi((warm_2 - warm_1) - (zcold_2 - zcold_1))[2],
              pseudo_p_val((warm_2 - warm_1) - (zcold_2 - zcold_1)))
  
  # Now we look at the directional effects for each main effect of area
  b_combined <- b_long %>%
    left_join(conditions, by = "obs_id")
  
  SS_SS_main_drift_directional = b_combined %>% filter(rating_dif == 0 & (area_size == 1 | area_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
    pivot_wider(names_from = c(area_size), values_from = prediction) %>% 
    group_by(task,quality) %>%
    summarize(mean = mean(`2` - `1`), 
              q5 = HDInterval::hdi(`2` - `1`)[1],
              q95 = HDInterval::hdi(`2` - `1`)[2],
              pseudo_p_val(`2` - `1`))
  
  
  # now the interaction between quality and task: (also magnitude)
  b_combined <- b_long %>%
    left_join(conditions, by = "obs_id") %>% 
    mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
  
  
  SS_task_quality_toway_drift = b_combined %>% filter(rating_dif == 0 & (area_size == 1))%>% select(-c(obs_id,rating_dif,area_size)) %>% 
    pivot_wider(names_from = c(task,quality), values_from = prediction) %>%
    # group_by(quality,task) %>% 
    mutate(prediction = (FastSS_warm - SlowSS_warm) - (FastSS_zcold - SlowSS_zcold)) %>%
    summarize(mean = mean(prediction), 
              q5 = HDInterval::hdi(prediction)[1],
              q95 = HDInterval::hdi(prediction)[2],
              pseudo_p_val(prediction))
  
  # and the main effects of quality (here also magnitude)
  SS_quality_main_drift =  b_combined %>% filter(rating_dif == 0 & (area_size == 1))%>% select(-c(obs_id,rating_dif,area_size)) %>% 
    pivot_wider(names_from = c(task,quality), values_from = prediction) %>%
    # group_by(quality,task) %>%
    mutate(fast = (FastSS_warm - FastSS_zcold),
           slow = (SlowSS_warm - SlowSS_zcold)) %>%
    pivot_longer(cols = c(fast,slow), values_to = "prediction") %>% group_by(name) %>%
    summarize(mean = mean(prediction), 
              q5 = HDInterval::hdi(prediction)[1],
              q95 = HDInterval::hdi(prediction)[2],
              pseudo_p_val(prediction))
  
  
  ## Now we look at the ratings part of the model on the drift rate in the same manner:
  
  ## ratings:
  b_combined <- b_long %>%
    left_join(conditions, by = "obs_id") %>% 
    mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
  
  
  # three wat magnitude:
  SS_rat_threeway_drift = b_combined %>% filter(area_size == 1) %>% select(-c(obs_id, area_size)) %>% 
    pivot_wider(names_from = c(quality,rating_dif,task), values_from = prediction) %>% 
    mutate(three_way = ((warm_1_FastSS - warm_0_FastSS) - (zcold_1_FastSS - zcold_0_FastSS)) - ((warm_1_SlowSS - warm_0_SlowSS) - (zcold_1_SlowSS - zcold_0_SlowSS))) %>% 
    summarize(mean = mean(three_way), 
              q5 = HDInterval::hdi(three_way)[1],
              q95 = HDInterval::hdi(three_way)[2],
              pseudo_p_val(three_way))
  
  
  SS_rat_twoway_drift = b_combined %>% filter(area_size == 1) %>% select(-c(obs_id, area_size)) %>% 
    pivot_wider(names_from = c(quality,rating_dif), values_from = prediction) %>% 
    group_by(task) %>%
    summarize(mean = mean((warm_1 - warm_0) - (zcold_1 - zcold_0)), 
              q5 = HDInterval::hdi((warm_1 - warm_0) - (zcold_1 - zcold_0))[1],
              q95 = HDInterval::hdi((warm_1 - warm_0) - (zcold_1 - zcold_0))[2],
              pseudo_p_val((warm_1- warm_0) - (zcold_1 - zcold_0)))
  
  b_combined <- b_long %>%
    left_join(conditions, by = "obs_id")
  
  SS_rat_main_drift_directional = b_combined %>% filter(area_size == 1) %>% select(-c(obs_id, area_size)) %>% 
    pivot_wider(names_from = c(rating_dif), values_from = prediction) %>% 
    group_by(task,quality) %>%
    summarize(mean = mean(`1` - `0`), 
              q5 = HDInterval::hdi(`1` - `0`)[1],
              q95 = HDInterval::hdi(`1` - `0`)[2],
              pseudo_p_val(`1` - `0`))
  
  
  
  
  Spatial_summation_drift = list(Spatial_summation_SS = list(SS_SS_threeway_drift = SS_SS_threeway_drift,
                                                             SS_SS_toway_drift = SS_SS_toway_drift,
                                                             SS_SS_main_drift_directional = SS_SS_main_drift_directional),
                                 
       Spatial_summation_task_qual = list(SS_task_quality_toway_drift = SS_task_quality_toway_drift,
                                          SS_quality_main_drift = SS_quality_main_drift),
       
       Spatial_summation_rat = list(SS_rat_threeway_drift,SS_rat_threeway_drift,
                                    SS_rat_twoway_drift = SS_rat_twoway_drift,
                                    SS_rat_main_drift_directional = SS_rat_main_drift_directional)
       )
  
  
  # now for the boudary seperation
  
  
  bound = as_draws_df(fit_SS$draws("gm")) %>% 
    select(-contains(".")) %>% 
    rename_with(~c(paste0("alpha_",colnames(X_alpha)),paste0("beta_",names_delta),paste0("delta_",names_delta), paste0("ndt_",names_ndt))) %>% 
    select(contains("alpha")) %>% 
    mutate(draw = 1:n()) %>% 
    pivot_longer(-draw) %>% group_by(name) %>% 
    summarize(mean = mean(value), 
              q5 = HDInterval::hdi(value)[1],
              q95 = HDInterval::hdi(value)[2],
              pseudo_p_val(value))
  
  
  alpha_draws <- as_draws_df(fit_SS$draws("gm")) %>% 
    select(-contains(".")) %>% 
    rename_with(~c(
      paste0("alpha_", colnames(X_alpha)),
      paste0("beta_", names_delta),
      paste0("delta_", names_delta),
      paste0("ndt_", names_ndt)
    )) %>% 
    select(starts_with("alpha_")) %>% 
    as.matrix()  # dim: [4000 x K]
  
  
  
  X_alpha = data.frame(model.matrix(~1+task, data = df_fastSS %>% 
                                      select(task) %>% distinct() %>% 
                                      unnest()))
  
  
  
  b = as.matrix(X_alpha) %*% t(alpha_draws) 
  
  b_long <- as.data.frame(b) %>%
    mutate(obs_id = row_number()) %>%
    pivot_longer(
      cols = -obs_id,
      names_to = "draw",
      values_to = "prediction"
    )
  
  
  b_long$draw <- as.integer(gsub("V", "", b_long$draw))
  
  
  conditions <- df_fastSS %>%
    select(task) %>%
    distinct() %>%
    mutate(obs_id = row_number())
  
  
  b_combined <- b_long %>%
    left_join(conditions, by = "obs_id")
  
  
  SS_task_main_bound = b_combined %>% select(-c(obs_id)) %>% 
    pivot_wider(names_from = c(task), values_from = prediction) %>% 
    mutate(dif = (FastSS)-(SlowSS)) %>% 
    summarize(mean = mean(dif), 
              q5 = HDInterval::hdi(dif)[1],
              q95 = HDInterval::hdi(dif)[2],
              pseudo_p_val(dif))
  
  
  
  Spatial_summation_boundary = list(SS_task_main_bound = SS_task_main_bound)
  
  
  # And the Non-decision time
  
  ## NDT
  ndt = as_draws_df(fit_SS$draws("gm")) %>% 
    select(-contains(".")) %>% 
    rename_with(~c(paste0("alpha_",colnames(X_alpha)),paste0("beta_",names_delta),paste0("delta_",names_delta), paste0("ndt_",names_ndt))) %>% 
    select(contains("ndt")) %>% 
    mutate(draw = 1:n()) %>% 
    pivot_longer(-draw) %>% group_by(name) %>% 
    summarize(mean = mean(value), 
              q5 = HDInterval::hdi(value)[1],
              q95 = HDInterval::hdi(value)[2],
              pseudo_p_val(value))
  
    
    ndt_draws <- as_draws_df(fit_SS$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(
        paste0("alpha_", colnames(X_alpha)),
        paste0("beta_", names_delta),
        paste0("delta_", names_delta),
        paste0("ndt_", names_ndt)
      )) %>% 
      select(starts_with("ndt_")) %>% 
      as.matrix()  # dim: [4000 x K]
    

    
    X_tau = data.frame(model.matrix(~0+quality*task, data = df_fastSS %>% 
                                        select(task,quality) %>% distinct() %>% 
                                        unnest()))
    
    
    
    b = as.matrix(X_tau) %*% t(ndt_draws) 
    
    b_long <- as.data.frame(b) %>%
      mutate(obs_id = row_number()) %>%
      pivot_longer(
        cols = -obs_id,
        names_to = "draw",
        values_to = "prediction"
      )
    
    
    b_long$draw <- as.integer(gsub("V", "", b_long$draw))
    
    
    conditions <- df_fastSS %>%
      select(task,quality) %>%
      distinct() %>%
      mutate(obs_id = row_number())
    
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id")
    
    
    
    all_ndt = b_combined %>% select(-c(obs_id)) %>% 
      pivot_wider(names_from = c(task,quality), values_from = prediction) %>% 
      mutate(warm_fast = (brms::inv_logit_scaled(FastSS_warm) * mean(min_warm[,1])),
             cold_fast = (brms::inv_logit_scaled(FastSS_zcold) * mean(min_cold[,1])),
             warm_slow = (brms::inv_logit_scaled(SlowSS_warm) * mean(min_warm[,2])),
             cold_slow = (brms::inv_logit_scaled(SlowSS_zcold) * mean(min_cold[,2])),
             dif_fast = warm_fast-cold_fast,
             dif_slow = warm_slow-cold_slow,
             avg_fast = (warm_fast + cold_fast) / 2,
             avg_slow = (warm_slow + cold_slow) / 2,
             main_task = avg_fast - avg_slow             
      ) %>% select(dif_fast,dif_slow,main_task) %>% pivot_longer(cols = c(dif_fast, dif_slow,main_task)) %>% 
      group_by(name) %>% 
      summarize(mean = mean(value), 
                q5 = HDInterval::hdi(value)[1],
                q95 = HDInterval::hdi(value)[2],
                pseudo_p_val(value))
    
    Spatial_summation_ndt = list(all_ndt = all_ndt)
    
    
    #and lastly the bias
    
    # bias
    
    bias = as_draws_df(fit_SS$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(paste0("alpha_",colnames(X_alpha)),paste0("beta_",names_beta),paste0("delta_",names_delta), paste0("ndt_",names_ndt))) %>% 
      select(contains("beta")) %>% 
      mutate(draw = 1:n()) %>% 
      pivot_longer(-draw) %>% group_by(name) %>% 
      summarize(mean = mean(value), 
                q5 = HDInterval::hdi(value)[1],
                q95 = HDInterval::hdi(value)[2],
                pseudo_p_val(value))
    
    
    beta_draws <- as_draws_df(fit_SS$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(
        paste0("alpha_", colnames(X_alpha)),
        paste0("beta_", names_beta),
        paste0("delta_", names_delta),
        paste0("ndt_", names_ndt)
      )) %>% 
      select(starts_with("beta_")) %>% 
      as.matrix()  # dim: [4000 x K]
    
    
    
    X_beta = data.frame(model.matrix(as.formula(X_betas[1]), data = df_fastSS %>% 
                                        select(area_size,quality,task) %>% distinct() %>% 
                                        mutate(rating_dif = list(c(0,1))) %>% unnest()))
    
    
    
    b = as.matrix(X_beta) %*% t(beta_draws) 
    
    b_long <- as.data.frame(b) %>%
      mutate(obs_id = row_number()) %>%
      pivot_longer(
        cols = -obs_id,
        names_to = "draw",
        values_to = "prediction"
      )
    
    
    b_long$draw <- as.integer(gsub("V", "", b_long$draw))
    
    
    
    conditions <- df_fastSS %>%
      select(area_size, quality, task) %>%
      distinct() %>%
      mutate(rating_dif = list(c(0,1))) %>% unnest() %>% 
      mutate(obs_id = row_number())
    
    
    # here again we flip the predictions of the warm quality to estimate the magnitude differences:
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") %>% 
      mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
    
    
    SS_SS_threeway_bias = b_combined %>% filter(rating_dif == 0 & (area_size == 1 | area_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
      pivot_wider(names_from = c(quality,area_size,task), values_from = prediction) %>% 
      mutate(three_way = ((warm_2_FastSS - warm_1_FastSS) - (zcold_2_FastSS - zcold_1_FastSS)) - ((warm_2_SlowSS - warm_1_SlowSS) - (zcold_2_SlowSS - zcold_1_SlowSS))) %>% 
      summarize(mean = mean(three_way), 
                q5 = HDInterval::hdi(three_way)[1],
                q95 = HDInterval::hdi(three_way)[2],
                pseudo_p_val(three_way))
    
    
    SS_SS_toway_bias = b_combined %>% filter(rating_dif == 0 & (area_size == 1 | area_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
      pivot_wider(names_from = c(quality,area_size), values_from = prediction) %>% 
      group_by(task) %>%
      summarize(mean = mean((warm_2 - warm_1) - (zcold_2 - zcold_1)), 
                q5 = HDInterval::hdi((warm_2 - warm_1) - (zcold_2 - zcold_1))[1],
                q95 = HDInterval::hdi((warm_2 - warm_1) - (zcold_2 - zcold_1))[2],
                pseudo_p_val((warm_2 - warm_1) - (zcold_2 - zcold_1)))
    
    # directional effects
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id")
    
    SS_main_directional_bias = b_combined %>% filter(rating_dif == 0 & (area_size == 1 | area_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
      pivot_wider(names_from = c(area_size), values_from = prediction) %>% 
      group_by(task,quality) %>%
      summarize(mean = mean(`2` - `1`), 
                q5 = HDInterval::hdi(`2` - `1`)[1],
                q95 = HDInterval::hdi(`2` - `1`)[2],
                pseudo_p_val(`2` - `1`))
    
    # now the task and quality interactions
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") %>% 
      mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
    

    SS_task_quality_toway_bias = b_combined %>% filter(rating_dif == 0 & (area_size == 1))%>% select(-c(obs_id,rating_dif,area_size)) %>% 
      pivot_wider(names_from = c(task,quality), values_from = prediction) %>% 
      mutate(prediction = (FastSS_warm - SlowSS_warm) - (FastSS_zcold - SlowSS_zcold)) %>% 
      summarize(mean = mean(prediction), 
                q5 = HDInterval::hdi(prediction)[1],
                q95 = HDInterval::hdi(prediction)[2],
                pseudo_p_val(prediction))
    
    # here we backflip the predictions to see the directionality
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id")
    
    
    SS_quality_main_bias = b_combined %>% filter(rating_dif == 0 & (area_size == 1))%>% select(-c(obs_id,rating_dif,area_size)) %>% 
      group_by(task,quality) %>%
      mutate(prediction = (prediction)) %>% 
      summarize(mean = mean(prediction), 
                q5 = HDInterval::hdi(prediction)[1],
                q95 = HDInterval::hdi(prediction)[2],
                pseudo_p_val(prediction))
    
    # and lastly the ratings
    
    ## ratings:
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") %>% 
      mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
    
    
    SS_rat_threeway_bias = b_combined %>% filter(area_size == 1) %>% select(-c(obs_id, area_size)) %>% 
      pivot_wider(names_from = c(quality,rating_dif,task), values_from = prediction) %>% 
      mutate(three_way = ((warm_1_FastSS - warm_0_FastSS) - (zcold_1_FastSS - zcold_0_FastSS)) - ((warm_1_SlowSS - warm_0_SlowSS) - (zcold_1_SlowSS - zcold_0_SlowSS))) %>% 
      summarize(mean = mean(three_way), 
                q5 = HDInterval::hdi(three_way)[1],
                q95 = HDInterval::hdi(three_way)[2],
                pseudo_p_val(three_way))
    
    
    SS_rat_twoway_bias = b_combined %>% filter(area_size == 1) %>% select(-c(obs_id, area_size)) %>% 
      pivot_wider(names_from = c(quality,rating_dif), values_from = prediction) %>% 
      group_by(task) %>%
      summarize(mean = mean((warm_1 - warm_0) - (zcold_1 - zcold_0)), 
                q5 = HDInterval::hdi((warm_1 - warm_0) - (zcold_1 - zcold_0))[1],
                q95 = HDInterval::hdi((warm_1 - warm_0) - (zcold_1 - zcold_0))[2],
                pseudo_p_val((warm_1- warm_0) - (zcold_1 - zcold_0)))
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id")
    
    SS_rat_main_bias_directional = b_combined %>% filter(area_size == 1) %>% select(-c(obs_id, area_size)) %>% 
      pivot_wider(names_from = c(rating_dif), values_from = prediction) %>% 
      group_by(task,quality) %>%
      summarize(mean = mean(`1` - `0`), 
                q5 = HDInterval::hdi(`1` - `0`)[1],
                q95 = HDInterval::hdi(`1` - `0`)[2],
                pseudo_p_val(`1` - `0`))
    
    
    
    
    Spatial_summation_bias = list(Spatial_summation_SS = list(SS_SS_threeway_bias = SS_SS_threeway_bias,
                                                              SS_SS_toway_bias = SS_SS_toway_bias,
                                                              SS_main_directional_bias = SS_main_directional_bias),
                                   
                                   Spatial_summation_task_qual = list(SS_task_quality_toway_bias = SS_task_quality_toway_bias,
                                                                      SS_quality_main_bias = SS_quality_main_bias),
                                   
                                   Spatial_summation_rat = list(SS_rat_threeway_bias,SS_rat_threeway_bias,
                                                                SS_rat_twoway_bias = SS_rat_twoway_bias,
                                                                SS_rat_main_bias_directional = SS_rat_main_bias_directional)
    )
    
    
    SS = list(Spatial_summation_drift = Spatial_summation_drift,
                Spatial_summation_boundary = Spatial_summation_boundary,
                Spatial_summation_ndt = Spatial_summation_ndt,
                Spatial_summation_bias = Spatial_summation_bias)
    
    ###########################################
    ################# Lateral inhibition ######
    ###########################################
    
    # Now the extract same just for lateral inhibition
    
    df_fastSS = df %>% 
      filter(ID != 1) %>% 
      group_by(ID, quality) %>%
      filter(corrected_RT > 0.2) %>%        # ascending: fastest first
      ungroup() %>% 
      dplyr::select(ID,config,distance_size,corrected_RT,quality,response_cw,rating,task,trial_index) %>% 
      arrange(trial_index) %>% 
      select(-trial_index) %>% 
      drop_na()
    
    df_fastSS = df_fastSS %>% group_by(ID,quality,distance_size) %>% mutate(rating_dif = rating-mean(rating)) %>% ungroup()
    df_fastSS = df_fastSS %>% mutate(quality = ifelse(quality == "cold","zcold","warm"))
    
    
    min_cold = as.matrix(df_fastSS %>% filter(quality == "zcold") %>% group_by(ID,task) %>% summarize(minrt = min(corrected_RT))  %>% ungroup() %>% select(task,minrt) %>% pivot_wider(names_from = "task",values_from = "minrt") %>% unnest())
    min_warm = as.matrix(df_fastSS %>% filter(quality == "warm") %>% group_by(ID,task) %>% summarize(minrt = min(corrected_RT))  %>% ungroup() %>% select(task,minrt) %>% pivot_wider(names_from = "task",values_from = "minrt") %>% unnest())
    
    
    fit_LI <- readRDS(here::here("STAN","taskmodels","fit_LI.rds"))
    
    # boundary separation
    X_alpha = model.matrix(~1+task, data = df_fastSS)
    
    
    # non-decision time
    X_tau = model.matrix(~0+quality*task, data = df_fastSS)
    
    
    X_deltas = c("~ 1 + distance_size * quality * task + quality * rating_dif * task")
    
    
    X_betas = c("~ 1 + distance_size * quality * task + quality * rating_dif * task")
    
    
    
    names_delta = colnames(model.matrix(as.formula(X_deltas[1]), data = df_fastSS))
    names_beta = colnames(model.matrix(as.formula(X_betas[1]), data = df_fastSS))
    names_ndt = colnames(X_tau)
    
    
    
    drift = as_draws_df(fit_LI$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(paste0("alpha_",colnames(X_alpha)),paste0("beta_",names_delta),paste0("delta_",names_delta), paste0("ndt_",names_ndt))) %>% 
      select(contains("delta")) %>% 
      mutate(draw = 1:n()) %>% 
      pivot_longer(-draw) %>% group_by(name) %>% 
      summarize(mean = mean(value), 
                q5 = HDInterval::hdi(value)[1],
                q95 = HDInterval::hdi(value)[2],
                pseudo_p_val(value))
    
    
    delta_draws <- as_draws_df(fit_LI$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(
        paste0("alpha_", colnames(X_alpha)),
        paste0("beta_", names_delta),
        paste0("delta_", names_delta),
        paste0("ndt_", names_ndt)
      )) %>% 
      select(starts_with("delta_")) %>% 
      as.matrix()  # dim: [4000 x K]
    
    
    
    X_delta = data.frame(model.matrix(as.formula(X_deltas[1]), data = df_fastSS %>% 
                                        select(distance_size,quality,task) %>% distinct() %>% 
                                        mutate(rating_dif = list(c(0,1))) %>% unnest()))
    
    
    
    b = as.matrix(X_delta) %*% t(delta_draws) 
    
    b_long <- as.data.frame(b) %>%
      mutate(obs_id = row_number()) %>%
      pivot_longer(
        cols = -obs_id,
        names_to = "draw",
        values_to = "prediction"
      )
    
    
    b_long$draw <- as.integer(gsub("V", "", b_long$draw))
    
    
    conditions <- df_fastSS %>%
      select(distance_size, quality, task) %>%
      distinct() %>%
      mutate(rating_dif = list(c(0,1))) %>% unnest() %>% 
      mutate(obs_id = row_number())
    
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") %>% 
      mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
    

  
    
    
    
    LI_LI_threeway_drift = b_combined %>% filter(rating_dif == 0 & (distance_size == 1 | distance_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
      pivot_wider(names_from = c(quality,distance_size,task), values_from = prediction) %>% 
      mutate(three_way = ((warm_2_FastLI - warm_1_FastLI) - (zcold_2_FastLI - zcold_1_FastLI)) - ((warm_2_SlowLI - warm_1_SlowLI) - (zcold_2_SlowLI - zcold_1_SlowLI))) %>% 
      summarize(mean = mean(three_way), 
                q5 = HDInterval::hdi(three_way)[1],
                q95 = HDInterval::hdi(three_way)[2],
                pseudo_p_val(three_way))
    
    
    LI_LI_toway_drift = b_combined %>% filter(rating_dif == 0 & (distance_size == 1 | distance_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
      pivot_wider(names_from = c(quality,distance_size), values_from = prediction) %>% 
      group_by(task) %>%
      summarize(mean = mean((warm_2 - warm_1) - (zcold_2 - zcold_1)), 
                q5 = HDInterval::hdi((warm_2 - warm_1) - (zcold_2 - zcold_1))[1],
                q95 = HDInterval::hdi((warm_2 - warm_1) - (zcold_2 - zcold_1))[2],
                pseudo_p_val((warm_2 - warm_1) - (zcold_2 - zcold_1)))
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id")
    
    LI_main_directional_drift = b_combined %>% filter(rating_dif == 0 & (distance_size == 1 | distance_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
      pivot_wider(names_from = c(distance_size), values_from = prediction) %>% 
      group_by(task,quality) %>%
      summarize(mean = mean(`2` - `1`), 
                q5 = HDInterval::hdi(`2` - `1`)[1],
                q95 = HDInterval::hdi(`2` - `1`)[2],
                pseudo_p_val(`2` - `1`))
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") %>% 
      mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
    
    
    LI_task_quality_toway_drift = b_combined %>% filter(rating_dif == 0 & (distance_size == 1))%>% select(-c(obs_id,rating_dif,distance_size)) %>% 
      pivot_wider(names_from = c(task,quality), values_from = prediction) %>%
      # group_by(quality,task) %>% 
      mutate(prediction = (FastLI_warm - SlowLI_warm) - (FastLI_zcold - SlowLI_zcold)) %>%
      summarize(mean = mean(prediction), 
                q5 = HDInterval::hdi(prediction)[1],
                q95 = HDInterval::hdi(prediction)[2],
                pseudo_p_val(prediction))
    
    LI_quality_main_drift = b_combined %>% filter(rating_dif == 0 & (distance_size == 1))%>% select(-c(obs_id,rating_dif,distance_size)) %>% 
      pivot_wider(names_from = c(task,quality), values_from = prediction) %>%
      # group_by(quality,task) %>%
      mutate(fast = (FastLI_warm - FastLI_zcold),
             slow = (SlowLI_warm - SlowLI_zcold)) %>%
      pivot_longer(cols = c(fast,slow), values_to = "prediction") %>% group_by(name) %>%
      summarize(mean = mean(prediction), 
                q5 = HDInterval::hdi(prediction)[1],
                q95 = HDInterval::hdi(prediction)[2],
                pseudo_p_val(prediction))
    
    
    

    ## ratings:
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") %>% 
      mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
    
    
    
  
    
    LI_rat_threeway_drift = b_combined %>% filter(distance_size == 1) %>% select(-c(obs_id, distance_size)) %>% 
      pivot_wider(names_from = c(quality,rating_dif,task), values_from = prediction) %>% 
      mutate(three_way = ((warm_1_FastLI - warm_0_FastLI) - (zcold_1_FastLI - zcold_0_FastLI)) - ((warm_1_SlowLI - warm_0_SlowLI) - (zcold_1_SlowLI - zcold_0_SlowLI))) %>% 
      summarize(mean = mean(three_way), 
                q5 = HDInterval::hdi(three_way)[1],
                q95 = HDInterval::hdi(three_way)[2],
                pseudo_p_val(three_way))
    
    
    LI_rat_twoway_drift = b_combined %>% filter(distance_size == 1) %>% select(-c(obs_id, distance_size)) %>% 
      pivot_wider(names_from = c(quality,rating_dif), values_from = prediction) %>% 
      group_by(task) %>%
      summarize(mean = mean((warm_1 - warm_0) - (zcold_1 - zcold_0)), 
                q5 = HDInterval::hdi((warm_1 - warm_0) - (zcold_1 - zcold_0))[1],
                q95 = HDInterval::hdi((warm_1 - warm_0) - (zcold_1 - zcold_0))[2],
                pseudo_p_val((warm_1- warm_0) - (zcold_1 - zcold_0)))
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") 
    
    
    LI_rat_main_drift_directional = b_combined %>% filter(distance_size == 1) %>% select(-c(obs_id, distance_size)) %>% 
      pivot_wider(names_from = c(rating_dif), values_from = prediction) %>% 
      group_by(task,quality) %>%
      summarize(mean = mean(`1` - `0`), 
                q5 = HDInterval::hdi(`1` - `0`)[1],
                q95 = HDInterval::hdi(`1` - `0`)[2],
                pseudo_p_val(`1` - `0`))
    
    
    
    Lateral_inhibition_drift = list(Spatial_summation_SS = list(LI_LI_threeway_drift = LI_LI_threeway_drift,
                                                              LI_LI_toway_drift = LI_LI_toway_drift,
                                                              LI_main_directional_drift = LI_main_directional_drift),
                                  
                                  Spatial_summation_task_qual = list(LI_task_quality_toway_drift = LI_task_quality_toway_drift,
                                                                     LI_quality_main_drift = LI_quality_main_drift),
                                  
                                  Spatial_summation_rat = list(LI_rat_threeway_drift,LI_rat_threeway_drift,
                                                               LI_rat_twoway_drift = LI_rat_twoway_drift,
                                                               LI_rat_main_drift_directional = LI_rat_main_drift_directional)
    )
    
    
    
    bound = as_draws_df(fit_LI$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(paste0("alpha_",colnames(X_alpha)),paste0("beta_",names_delta),paste0("delta_",names_delta), paste0("ndt_",names_ndt))) %>% 
      select(contains("alpha")) %>% 
      mutate(draw = 1:n()) %>% 
      pivot_longer(-draw) %>% group_by(name) %>% 
      summarize(mean = mean(value), 
                q5 = HDInterval::hdi(value)[1],
                q95 = HDInterval::hdi(value)[2],
                pseudo_p_val(value))
    
    
    alpha_draws <- as_draws_df(fit_LI$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(
        paste0("alpha_", colnames(X_alpha)),
        paste0("beta_", names_delta),
        paste0("delta_", names_delta),
        paste0("ndt_", names_ndt)
      )) %>% 
      select(starts_with("alpha_")) %>% 
      as.matrix()  # dim: [4000 x K]
    
    
    
    X_alpha = data.frame(model.matrix(~1+task, data = df_fastSS %>% 
                                        select(task) %>% distinct() %>% 
                                        unnest()))
    
    
    
    b = as.matrix(X_alpha) %*% t(alpha_draws) 
    
    b_long <- as.data.frame(b) %>%
      mutate(obs_id = row_number()) %>%
      pivot_longer(
        cols = -obs_id,
        names_to = "draw",
        values_to = "prediction"
      )
    
    
    b_long$draw <- as.integer(gsub("V", "", b_long$draw))
    
    
    conditions <- df_fastSS %>%
      select(task) %>%
      distinct() %>%
      mutate(obs_id = row_number())
    
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id")
    
    
    LI_task_main_bound = b_combined %>% select(-c(obs_id)) %>% 
      pivot_wider(names_from = c(task), values_from = prediction) %>% 
      mutate(dif = (FastLI)-(SlowLI)) %>% 
      summarize(mean = mean(dif), 
                q5 = HDInterval::hdi(dif)[1],
                q95 = HDInterval::hdi(dif)[2],
                pseudo_p_val(dif))
    
    
    
    Lateral_inhibition_boundary = list(LI_task_main_bound = LI_task_main_bound)
    
    ## NDT
    
    ndt = as_draws_df(fit_LI$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(paste0("alpha_",colnames(X_alpha)),paste0("beta_",names_delta),paste0("delta_",names_delta), paste0("ndt_",names_ndt))) %>% 
      select(contains("ndt")) %>% 
      mutate(draw = 1:n()) %>% 
      pivot_longer(-draw) %>% group_by(name) %>% 
      summarize(mean = mean(value), 
                q5 = HDInterval::hdi(value)[1],
                q95 = HDInterval::hdi(value)[2],
                pseudo_p_val(value))
    
    
    ndt_draws <- as_draws_df(fit_LI$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(
        paste0("alpha_", colnames(X_alpha)),
        paste0("beta_", names_delta),
        paste0("delta_", names_delta),
        paste0("ndt_", names_ndt)
      )) %>% 
      select(starts_with("ndt_")) %>% 
      as.matrix()  # dim: [4000 x K]
    
    
    
    X_tau = data.frame(model.matrix(~0+quality*task, data = df_fastSS %>% 
                                      select(task,quality) %>% distinct() %>% 
                                      unnest()))
    
    
    
    b = as.matrix(X_tau) %*% t(ndt_draws) 
    
    b_long <- as.data.frame(b) %>%
      mutate(obs_id = row_number()) %>%
      pivot_longer(
        cols = -obs_id,
        names_to = "draw",
        values_to = "prediction"
      )
    
    
    b_long$draw <- as.integer(gsub("V", "", b_long$draw))
    
    
    conditions <- df_fastSS %>%
      select(task,quality) %>%
      distinct() %>%
      mutate(obs_id = row_number())
    
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id")
    
    
    ndt_all_li = b_combined %>% select(-c(obs_id)) %>% 
      pivot_wider(names_from = c(task,quality), values_from = prediction) %>% 
      mutate(warm_fast = (brms::inv_logit_scaled(FastLI_warm) * mean(min_warm[,1])),
             cold_fast = (brms::inv_logit_scaled(FastLI_zcold) * mean(min_cold[,1])),
             warm_slow = (brms::inv_logit_scaled(SlowLI_warm) * mean(min_warm[,2])),
             cold_slow = (brms::inv_logit_scaled(SlowLI_zcold) * mean(min_cold[,2])),
             dif_fast = warm_fast-cold_fast,
             dif_slow = warm_slow-cold_slow,
             avg_fast = (warm_fast + cold_fast) / 2,
             avg_slow = (warm_slow + cold_slow) / 2,
             main_task = avg_fast - avg_slow             
      ) %>% 
      select(dif_fast,dif_slow,main_task) %>%
      pivot_longer(cols = c(dif_fast, dif_slow,main_task)) %>%
      # pivot_longer(everything()) %>% 
      group_by(name) %>%
      summarize(mean = mean(value), 
                q5 = HDInterval::hdi(value)[1],
                q95 = HDInterval::hdi(value)[2],
                pseudo_p_val(value))
    
    
    Lateral_inhibition_ndt = list(ndt_all_li = ndt_all_li)
    
    
    # bias
    
    
    bias = as_draws_df(fit_LI$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(paste0("alpha_",colnames(X_alpha)),paste0("beta_",names_beta),paste0("delta_",names_delta), paste0("ndt_",names_ndt))) %>% 
      select(contains("beta")) %>% 
      mutate(draw = 1:n()) %>% 
      pivot_longer(-draw) %>% group_by(name) %>% 
      summarize(mean = mean(value), 
                q5 = HDInterval::hdi(value)[1],
                q95 = HDInterval::hdi(value)[2],
                pseudo_p_val(value))
    
    
    beta_draws <- as_draws_df(fit_LI$draws("gm")) %>% 
      select(-contains(".")) %>% 
      rename_with(~c(
        paste0("alpha_", colnames(X_alpha)),
        paste0("beta_", names_beta),
        paste0("delta_", names_delta),
        paste0("ndt_", names_ndt)
      )) %>% 
      select(starts_with("beta_")) %>% 
      as.matrix()  # dim: [4000 x K]
    
    
    
    X_beta = data.frame(model.matrix(as.formula(X_betas[1]), data = df_fastSS %>% 
                                       select(distance_size,quality,task) %>% distinct() %>% 
                                       mutate(rating_dif = list(c(0,1))) %>% unnest()))
    
    
    
    b = as.matrix(X_beta) %*% t(beta_draws) 
    
    b_long <- as.data.frame(b) %>%
      mutate(obs_id = row_number()) %>%
      pivot_longer(
        cols = -obs_id,
        names_to = "draw",
        values_to = "prediction"
      )
    
    
    b_long$draw <- as.integer(gsub("V", "", b_long$draw))
    
    
    
    conditions <- df_fastSS %>%
      select(distance_size, quality, task) %>%
      distinct() %>%
      mutate(rating_dif = list(c(0,1))) %>% unnest() %>% 
      mutate(obs_id = row_number())
    
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") %>% 
    mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
    
    
    
    LI_LI_threeway_bias = b_combined %>% filter(rating_dif == 0 & (distance_size == 1 | distance_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
      pivot_wider(names_from = c(quality,distance_size,task), values_from = prediction) %>% 
      mutate(three_way = ((warm_2_FastLI - warm_1_FastLI) - (zcold_2_FastLI - zcold_1_FastLI)) - ((warm_2_SlowLI - warm_1_SlowLI) - (zcold_2_SlowLI - zcold_1_SlowLI))) %>% 
      summarize(mean = mean(three_way), 
                q5 = HDInterval::hdi(three_way)[1],
                q95 = HDInterval::hdi(three_way)[2],
                pseudo_p_val(three_way))
    
    
    LI_LI_toway_bias = b_combined %>% filter(rating_dif == 0 & (distance_size == 1 | distance_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
      pivot_wider(names_from = c(quality,distance_size), values_from = prediction) %>% 
      group_by(task) %>%
      summarize(mean = mean((warm_2 - warm_1) - (zcold_2 - zcold_1)), 
                q5 = HDInterval::hdi((warm_2 - warm_1) - (zcold_2 - zcold_1))[1],
                q95 = HDInterval::hdi((warm_2 - warm_1) - (zcold_2 - zcold_1))[2],
                pseudo_p_val((warm_2 - warm_1) - (zcold_2 - zcold_1)))
    
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") 
    
    
    LI_main_directional_bias = b_combined %>% filter(rating_dif == 0 & (distance_size == 1 | distance_size == 2)) %>% select(-c(obs_id,rating_dif)) %>% 
      pivot_wider(names_from = c(distance_size), values_from = prediction) %>% 
      group_by(task,quality) %>%
      summarize(mean = mean(`2` - `1`), 
                q5 = HDInterval::hdi(`2` - `1`)[1],
                q95 = HDInterval::hdi(`2` - `1`)[2],
                pseudo_p_val(`2` - `1`))
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") %>% 
      mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
    
    
    LI_task_quality_toway_bias = b_combined %>% filter(rating_dif == 0 & (distance_size == 1))%>% select(-c(obs_id,rating_dif,distance_size)) %>% 
      pivot_wider(names_from = c(task,quality), values_from = prediction) %>% 
      mutate(prediction = (FastLI_warm - SlowLI_warm) - (FastLI_zcold - SlowLI_zcold)) %>% 
      summarize(mean = mean(prediction), 
                q5 = HDInterval::hdi(prediction)[1],
                q95 = HDInterval::hdi(prediction)[2],
                pseudo_p_val(prediction))
    
    

    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") 
      # mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
    
    LI_quality_main_bias = b_combined %>% filter(rating_dif == 0 & (distance_size == 1))%>% select(-c(obs_id,rating_dif,distance_size)) %>% 
      group_by(task,quality) %>%
      mutate(prediction = (prediction)) %>% 
      summarize(mean = mean(prediction), 
                q5 = HDInterval::hdi(prediction)[1],
                q95 = HDInterval::hdi(prediction)[2],
                pseudo_p_val(prediction))
    
    
    
    ## ratings:
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id") %>% 
      mutate(prediction = ifelse(quality == "warm",-prediction,prediction))
    
  
    
    
    LI_rat_threeway_bias = b_combined %>% filter(distance_size == 1) %>% select(-c(obs_id, distance_size)) %>% 
      pivot_wider(names_from = c(quality,rating_dif,task), values_from = prediction) %>% 
      mutate(three_way = ((warm_1_FastLI - warm_0_FastLI) - (zcold_1_FastLI - zcold_0_FastLI)) - ((warm_1_SlowLI - warm_0_SlowLI) - (zcold_1_SlowLI - zcold_0_SlowLI))) %>% 
      summarize(mean = mean(three_way), 
                q5 = HDInterval::hdi(three_way)[1],
                q95 = HDInterval::hdi(three_way)[2],
                pseudo_p_val(three_way))
    
    
    LI_rat_twoway_bias =  b_combined %>% filter(distance_size == 1) %>% select(-c(obs_id, distance_size)) %>% 
      pivot_wider(names_from = c(quality,rating_dif), values_from = prediction) %>% 
      group_by(task) %>%
      summarize(mean = mean((warm_1 - warm_0) - (zcold_1 - zcold_0)), 
                q5 = HDInterval::hdi((warm_1 - warm_0) - (zcold_1 - zcold_0))[1],
                q95 = HDInterval::hdi((warm_1 - warm_0) - (zcold_1 - zcold_0))[2],
                pseudo_p_val((warm_1- warm_0) - (zcold_1 - zcold_0)))
    
    
    b_combined <- b_long %>%
      left_join(conditions, by = "obs_id")
    
    LI_rat_main_bias_directional = b_combined %>% filter(distance_size == 1) %>% select(-c(obs_id, distance_size)) %>% 
      pivot_wider(names_from = c(rating_dif), values_from = prediction) %>% 
      group_by(task,quality) %>%
      summarize(mean = mean(`1` - `0`), 
                q5 = HDInterval::hdi(`1` - `0`)[1],
                q95 = HDInterval::hdi(`1` - `0`)[2],
                pseudo_p_val(`1` - `0`))
    
    
    
    
    Lateral_inhibition_bias = list(Spatial_summation_SS = list(LI_LI_threeway_bias = LI_LI_threeway_bias,
                                                              LI_LI_toway_bias = LI_LI_toway_bias,
                                                              LI_main_directional_bias = LI_main_directional_bias),
                                  
                                  Spatial_summation_task_qual = list(LI_task_quality_toway_bias = LI_task_quality_toway_bias,
                                                                     LI_quality_main_bias = LI_quality_main_bias),
                                  
                                  Spatial_summation_rat = list(LI_rat_threeway_bias,LI_rat_threeway_bias,
                                                               LI_rat_twoway_bias = LI_rat_twoway_bias,
                                                               LI_rat_main_bias_directional = LI_rat_main_bias_directional)
    )
    
  
    LI = list(Lateral_inhibition_drift = Lateral_inhibition_drift,
              Lateral_inhibition_boundary = Lateral_inhibition_boundary,
              Lateral_inhibition_ndt = Lateral_inhibition_ndt,
              Lateral_inhibition_bias = Lateral_inhibition_bias)
    
    

    
  return(list(SS,LI))
  
}

