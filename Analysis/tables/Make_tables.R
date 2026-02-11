

# seed
set.seed(1993)

# Load required packages
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library("pacman")

# Define path and data
p_load(here, readr, tidyverse, dplyr)

# Visualization
p_load(ggplot2, ggpubr, gghalves, ggrain, cowplot, scales, flextable)

# Statistics  
p_load(lmerTest, MuMIn, gamlss, devtools, glmmTMB, DHARMa, ggeffects)


# Load required functions
source(here::here("Analysis","functions","effectsize_cohend.R"))
source(here::here("Analysis","functions","utils_Jesper.R"))
source(here::here("Analysis","functions","load_data.R"))


source(here::here("Analysis","tables","make_supplementary_table.R"))
## R Markdown


## Load data
df <- read_csv(here("Data","exp1.csv"), show_col_types = FALSE)

# Excluded participant 
excluded_participants = list(1)

# Remove excluded trials 
original_rows <- nrow(df)
origional_IDs <- length(unique(df$ID))
df <- df %>%
  filter(excluded == 0, !ID %in% excluded_participants) %>% 
  mutate(gender = ifelse(gender == FALSE, "F", gender))

deleted_rows <- original_rows - nrow(df)
deleted_IDs <- origional_IDs - length(unique(df$ID))

cat("Number of trials deleted:", deleted_rows, "\n")
cat("Number of participants deleted:", deleted_IDs, "\n")

# Scaling rating col from between 0 and 100 to fit between 0 and 1 
df$rating_scaled = df$rating/100


# spatial summation (rating) (AIC

# testing the hypotheses:

ss_fast_rating <- df %>% filter(task == "FastSS" | task == "SlowSS") %>% 
  dplyr::select(ID, size, area_size, quality, trial, trial_rep, rating_scaled, corrected_RT, task) %>% 
  drop_na() %>% 
  mutate(ID = as.factor(ID),  
         quality = as.factor(quality), 
         size = as.numeric(size), 
         area_size = as.numeric(area_size)) 

ss_fast_rating$quality <- relevel(ss_fast_rating$quality, ref = "warm")
ss_fast_rating$task <- relevel(as.factor(ss_fast_rating$task), ref = "SlowSS")

ss_rating_3way_log <- gamlss(rating_scaled ~ task*quality * log(area_size) +task* corrected_RT + trial  +  random(ID),
                             nu.formula = ~ task*quality * log(area_size) + corrected_RT + trial  +  random(ID),
                             sigma.formula = ~ task*quality * log(area_size) + corrected_RT + trial+  random(ID),
                             data = ss_fast_rating,
                             family = BEINF0(mu.link = "logit", sigma.link = "logit", nu.link = "logit"),
                             control = gamlss.control(n.cyc = 100, trace = T))


table_vas_ss = get_main_tables(ss_rating_3way_log)



ss_fast_rt <- df %>% filter(task == "FastSS" | task == "SlowSS") %>% 
  dplyr::select(ID, size, area_size, quality, trial, trial_rep, rating_scaled, corrected_RT, task) %>% 
  drop_na() %>% 
  mutate(ID = as.factor(ID),  
         quality = as.factor(quality), 
         size = as.numeric(size), 
         area_size = as.numeric(area_size)) 


ss_fast_rt$quality <- relevel(ss_fast_rt$quality, ref = "warm")
ss_fast_rt$task <- relevel(as.factor(ss_fast_rt$task), ref = "SlowSS")

ss_rt_3way_log <- gamlss(corrected_RT ~ log(area_size) * quality * task + trial + random(ID),
                         nu.formula = ~ log(area_size) * quality * task + trial + random(ID),
                         sigma.formula = ~ log(area_size)* quality * task + trial + random(ID),
                         data = ss_fast_rt,
                         family = GA(mu.link = "log", sigma.link = "log"),
                         control = gamlss.control(n.cyc = 100, trace = T))


table_RT_ss = get_main_tables(ss_rt_3way_log)






####### lateral inhibition
ss_fast_rating <- df %>% filter(task == "FastLI" | task == "SlowLI") %>% 
  dplyr::select(ID, size, distance_size, quality, trial, trial_rep, rating_scaled, corrected_RT, task) %>% 
  drop_na() %>% 
  mutate(ID = as.factor(ID),  
         quality = as.factor(quality), 
         size = as.numeric(size), 
         distance_size = as.numeric(distance_size)) %>% filter(rating_scaled != 1)

ss_fast_rating$quality <- relevel(ss_fast_rating$quality, ref = "warm")
ss_fast_rating$task <- relevel(as.factor(ss_fast_rating$task), ref = "SlowLI")

ss_rating_3way_li <- gamlss(rating_scaled ~ task*quality * distance_size + task* corrected_RT + trial  +  random(ID),
                            nu.formula = ~ task*quality * distance_size + task* corrected_RT + trial  +  random(ID),
                            sigma.formula = ~ task*quality * distance_size + task* corrected_RT + trial+  random(ID),
                            data = ss_fast_rating,
                            family = BEINF0(mu.link = "logit", sigma.link = "logit", nu.link = "logit"),
                            control = gamlss.control(n.cyc = 100, trace = T))



table_vas_li = get_main_tables(ss_rating_3way_li)



ss_fast_rt <- df %>% filter(task == "FastLI" | task == "SlowLI") %>% 
  dplyr::select(ID, size, distance_size, quality, trial, trial_rep, rating_scaled, corrected_RT, task) %>% 
  drop_na() %>% 
  mutate(ID = as.factor(ID),  
         quality = as.factor(quality), 
         size = as.numeric(size), 
         distance_size = as.numeric(distance_size))

ss_fast_rt$quality <- relevel(ss_fast_rt$quality, ref = "warm")
ss_fast_rt$task <- relevel(as.factor(ss_fast_rt$task), ref = "SlowLI")


ss_rt_3way_li <- gamlss(corrected_RT ~ distance_size * quality * task + trial + random(ID),
                        nu.formula = ~ distance_size * quality * task + trial + random(ID),
                        sigma.formula = ~ distance_size * quality * task + trial + random(ID),
                        data = ss_fast_rt,
                        family = GA(mu.link = "log", sigma.link = "log"),
                        control = gamlss.control(n.cyc = 100, trace = T))

table_rt_li = get_main_tables(ss_rt_3way_li)

# show the,

table_rt_li
table_vas_li
table_RT_ss
table_vas_ss

#save them

table_vars <- ls(pattern = "table_")
tables_list <- mget(table_vars, envir = .GlobalEnv)

tables_list = tables_list[-3]

saveRDS(tables_list, file = here::here("Analysis","tables","tables.RData"))
