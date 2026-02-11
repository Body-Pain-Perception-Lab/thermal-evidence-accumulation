
## script to run and estimate both STAN-models (workflow follows extracting data, setting up design matrices and then fitting the model)

set.seed(1993)

if (!requireNamespace("pacman", quietly = TRUE)) {
install.packages("pacman")
}

library(pacman)

# Define path and data
p_load(here, readr, tidyverse, dplyr, cmdstanr,loo)


# Load required functions
source(here::here("Analysis","functions","utils.R"))


# get data
df = read.csv(here::here("Data","exp1.csv")) %>%
  mutate(trial_index = row_number())

df_fastSS <- df %>% 
  filter(outlier == 0, ID != 1) %>% 
  group_by(ID, quality) %>%
  filter(corrected_RT > 0.2) %>%
  ungroup() %>% 
  dplyr::select(ID,config,area_size,corrected_RT,quality,response_cw,rating,task,trial_index) %>% 
  arrange(trial_index) %>% 
  select(-trial_index) %>% 
  drop_na()

df_fastSS = df_fastSS %>% group_by(ID,quality,area_size,task) %>% mutate(rating_dif = rating-mean(rating)) %>% ungroup()


df_fastSS = df_fastSS %>% mutate(quality = ifelse(quality == "cold","zcold","warm"))



# boundary separation
X_alpha = model.matrix(~1+task, data = df_fastSS)


# non-decision time
X_tau = model.matrix(~0+quality:task, data = df_fastSS)

# drift rate
X_deltas = c("~ 1 + area_size * quality * task + quality * rating_dif * task")

# bias
X_betas = c("~ 1 + area_size * quality * task + quality * rating_dif * task")


# compile the model and run it

# DDM = cmdstanr::cmdstan_model(here::here("STAN","stanmodels", "stanmodel_task_v2.stan"), force_recompile = T)
DDM = cmdstanr::cmdstan_model(here::here("STAN","stanmodels", "stanmodel_task.stan"), force_recompile = T)

X_delta <- model.matrix(as.formula(X_deltas[1]), data = df_fastSS)
X_beta <- model.matrix(as.formula(X_betas[1]), data = df_fastSS)

data_stan = list(
  trials = nrow(df_fastSS), # How many trials in total 
  S = length(unique(df_fastSS$ID)), # How many subjects in total 
  S_id = as.numeric(as.factor(as.numeric(df_fastSS$ID))), # Index of each ID 
  minRT_cold = as.matrix(df_fastSS %>% filter(quality == "zcold") %>% group_by(ID,task) %>% summarize(minrt = min(corrected_RT))  %>% ungroup() %>% select(task,minrt) %>% pivot_wider(names_from = "task",values_from = "minrt") %>% unnest()),
  minRT_warm = as.matrix(df_fastSS %>% filter(quality == "warm") %>% group_by(ID,task) %>% summarize(minrt = min(corrected_RT))  %>% ungroup() %>% select(task,minrt) %>% pivot_wider(names_from = "task",values_from = "minrt") %>% unnest()),
  RT = df_fastSS$corrected_RT,
  quality = df_fastSS %>% mutate(quality = ifelse(quality == "warm", 0, 1)) %>% .$quality,
  task = df_fastSS %>% mutate(task = ifelse(task == "FastSS", 0, 1)) %>% .$task,
  X_alpha = X_alpha,
  X_delta = X_delta,
  X_beta = X_beta,
  X_tau = X_tau,
  N_alpha = ncol(X_alpha), # number of cols = number of parameters 
  N_beta = ncol(X_beta),
  N_delta = ncol(X_delta),
  N_tau = ncol(X_tau),
  resp = df_fastSS %>% mutate(resp = ifelse(response_cw ==  "warm",0, 1)) %>% .$resp
)

# Run the model
fit_models_simp <- DDM$sample(
  data = data_stan,
  chains = 4,
  iter_sampling = 1000, 
  iter_warmup = 1000, 
  parallel_chains = 4,
  adapt_delta = 0.9,
  refresh = 50,
  init = 0,
  max_treedepth = 12
)

# save the fit

fit_models_simp$save_object(here("STAN","taskmodels","fit_SS.rds"))


#####################################################################
####################### For Lateral inhibition! #####################
#####################################################################


# Load required functions
source(here::here("Analysis","functions","utils.R"))

df = read.csv(here::here("Data","exp1.csv")) %>%
mutate(trial_index = row_number())

df_fastSS <- df %>% 
filter(ID != 1) %>% 
group_by(ID, quality) %>%
filter(corrected_RT > 0.2) %>%
ungroup() %>% 
dplyr::select(ID,config,distance_size,corrected_RT,quality,response_cw,rating,task,trial_index) %>% 
arrange(trial_index) %>% 
select(-trial_index) %>% 
drop_na()

df_fastSS = df_fastSS %>% group_by(ID,quality,distance_size) %>% mutate(rating_dif = rating-mean(rating)) %>% ungroup()


df_fastSS = df_fastSS %>% mutate(quality = ifelse(quality == "cold","zcold","warm"))

# boundary separation
X_alpha = model.matrix(~1+task, data = df_fastSS)


# non-decision time
X_tau = model.matrix(~0+quality*task, data = df_fastSS)


X_deltas = c("~ 1 + distance_size * quality * task + quality * rating_dif * task")


X_betas = c("~ 1 + distance_size * quality * task + quality * rating_dif * task")



# DDM = cmdstanr::cmdstan_model(here::here("STAN","stanmodels", "stanmodel_task_v2.stan"), force_recompile = T)
DDM = cmdstanr::cmdstan_model(here::here("STAN","stanmodels", "stanmodel_task.stan"), force_recompile = T)



X_delta <- model.matrix(as.formula(X_deltas[1]), data = df_fastSS)
X_beta <- model.matrix(as.formula(X_betas[1]), data = df_fastSS)

data_stan = list(
trials = nrow(df_fastSS), # How many trials in total 
S = length(unique(df_fastSS$ID)), # How many subjects in total 
S_id = as.numeric(as.factor(as.numeric(df_fastSS$ID))), # Index of each ID 
minRT_cold = as.matrix(df_fastSS %>% filter(quality == "zcold") %>% group_by(ID,task) %>% summarize(minrt = min(corrected_RT))  %>% ungroup() %>% select(task,minrt) %>% pivot_wider(names_from = "task",values_from = "minrt") %>% unnest()),
minRT_warm = as.matrix(df_fastSS %>% filter(quality == "warm") %>% group_by(ID,task) %>% summarize(minrt = min(corrected_RT))  %>% ungroup() %>% select(task,minrt) %>% pivot_wider(names_from = "task",values_from = "minrt") %>% unnest()),
RT = df_fastSS$corrected_RT,
quality = df_fastSS %>% mutate(quality = ifelse(quality == "warm", 0, 1)) %>% .$quality,
task = df_fastSS %>% mutate(task = ifelse(task == "FastLI", 0, 1)) %>% .$task,
X_alpha = X_alpha,
X_delta = X_delta,
X_beta = X_beta,
X_tau = X_tau,
N_alpha = ncol(X_alpha), # number of cols = number of parameters 
N_beta = ncol(X_beta),
N_delta = ncol(X_delta),
N_tau = ncol(X_tau),
resp = df_fastSS %>% mutate(resp = ifelse(response_cw ==  "warm",0, 1)) %>% .$resp
)



# Run the model and store in list with model name
fit_models_simp <- DDM$sample(
  data = data_stan,
  chains = 4,
  iter_sampling = 1000, # should be min 250
  iter_warmup = 1000, # should be min 250
  parallel_chains = 4,
  adapt_delta = 0.90,
  refresh = 50,
  init = 0,
  max_treedepth = 12
)

fit_models_simp$save_object(here("STAN","taskmodels","fit_LI.rds"))





