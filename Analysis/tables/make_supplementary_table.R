#tables:


#function take takes the models used in the study and makes a table of all regression coefficients. (These are saved in xlsx files)
get_main_tables <- function(model, round = 2) {
  
  #retrieve the fixed effects of the model
  
  if("Component" %in% names(data.frame(parameters::model_parameters(model)))){
    
    fixedeffecs <- parameters::model_parameters(model) %>%
      mutate(CI = NULL, CI_low = NULL, CI_high = NULL, df_error = NULL) %>%
      dplyr::rename(parameter = Component) %>%
      dplyr::select(parameter, everything()) %>% 
      mutate(parameter = ifelse(str_detect(parameter, "conditional"), "μ", ifelse(str_detect(parameter, "sigma"), "σ", ifelse(str_detect(parameter, "tau"), "τ", "ν"))))
    #renaming
    names(fixedeffecs) <- c("parameter","contrast", "\u03B2", "SE", "t", "p")
    #formular for the model (i.e. the math)
    formular <- as.character(formula(model))
  }else{
    
    fixedeffecs <- parameters::model_parameters(model) %>%
      mutate(CI = NULL, CI_low = NULL, CI_high = NULL, df_error = NULL, Component = "conditional") %>%
      dplyr::rename(parameter = Component) %>%
      dplyr::select(parameter, everything()) %>% 
      mutate(parameter = ifelse(str_detect(parameter, "conditional"), "μ", ifelse(str_detect(parameter, "sigma"), "σ", ifelse(str_detect(parameter, "tau"), "τ", "ν"))))
    #renaming
    names(fixedeffecs) <- c("parameter","contrast", "\u03B2", "SE", "t", "p")
    #formular for the model (i.e. the math)
    formular <- as.character(formula(model))
    
    
  }
  #get family
  family = family(model)[2]
  link = model$mu.link
  
  if(family == "Beta Inflated"){
    family = "ZOIB"
  }else if(family == "Beta Inflated zero"){
    family = "ZIB"
  }
  
  #big renaming of columns: 
  ###################################### spatial summation model
  if(sum(grepl("corrected_RT", fixedeffecs$contrast)) != 0 & sum(grepl("area_size", fixedeffecs$contrast)) != 0) {
    fixedeffecs <- fixedeffecs %>% mutate(contrast = case_when(
      contrast == "taskSlowSS" ~ "Task(Target)",
      contrast == "qualitywarm" ~ "Quality(Warm)",
      contrast == "log(area_size)" ~ "log(Area)",
      contrast == "corrected_RT" ~ "Response time",
      contrast == "trial" ~ "Trial_n",
      contrast == "taskSlowSS:qualitywarm" ~ "Task(Target):Quality(Warm)",
      contrast == "taskSlowSS:log(area_size)" ~ "Task(Target):log(Area)",
      contrast == "qualitywarm:log(area_size)" ~ "Quality(Warm):log(Area)",

      contrast == "taskSlowSS:corrected_RT" ~ "Task(Target):Response time",
      contrast == "taskSlowSS:qualitywarm:log(area_size)" ~ "Task(Target):Quality(Warm):\nlog(Area)",
      TRUE ~ contrast
    ))
    references <- c(
      "Quality" = "Warm",
      "Task" = "Fixed"
    )
  }
 
  # response time:
  # spatial summation models:
  if(sum(grepl("corrected_RT", fixedeffecs$contrast)) == 0 & sum(grepl("area_size", fixedeffecs$contrast)) != 0) {
    fixedeffecs <- fixedeffecs %>% mutate(contrast = case_when(
      contrast == "taskSlowSS" ~ "Task(Target)",
      contrast == "qualitywarm" ~ "Quality(Warm)",
      contrast == "log(area_size)" ~ "log(Area)",
      contrast == "trial" ~ "Trial_n",
      contrast == "taskSlowSS:qualitywarm" ~ "Task(Target):Quality(Warm)",
      contrast == "taskSlowSS:log(area_size)" ~ "Task(Target):log(Area)",
      contrast == "qualitywarm:log(area_size)" ~ "Quality(Warm):log(Area)",
      contrast == "taskSlowSS:qualitywarm:log(area_size)" ~ "Task(Target):Quality(Warm):\nlog(Area)",
      TRUE ~ contrast
    ))
    references <- c(
      "Quality" = "Warm",
      "Task" = "Fixed"
    )
  }
  
  
  
  ###################################### spatial summation model
  if(sum(grepl("corrected_RT", fixedeffecs$contrast)) != 0) {
    fixedeffecs <- fixedeffecs %>% mutate(contrast = case_when(
      contrast == "taskSlowLI" ~ "Task(Target)",
      contrast == "qualitywarm" ~ "Quality(Warm)",
      contrast == "distance_size" ~ "Distance",
      contrast == "corrected_RT" ~ "Response time",
      contrast == "trial" ~ "Trial_n",
      contrast == "taskSlowLI:qualitywarm" ~ "Task(Target):Quality(Warm)",
      contrast == "taskSlowLI:distance_size" ~ "Task(Target):Distance",
      contrast == "qualitywarm:distance_size" ~ "Quality(Warm):Distance",
      contrast == "taskSlowLI:corrected_RT" ~ "Task(Target):Response time",
      contrast == "taskSlowLI:qualitywarm:distance_size" ~ "Task(Target):Quality(Warm):\nDistance",
      TRUE ~ contrast
    ))
    references <- c(
      "Quality" = "Warm",
      "Task" = "Target"
    )
  }
  
  # response time:
  # lateral inhibition
  if(sum(grepl("corrected_RT", fixedeffecs$contrast)) == 0 & sum(grepl("distance_size", fixedeffecs$contrast)) != 0) {
    fixedeffecs <- fixedeffecs %>% mutate(contrast = case_when(
      contrast == "taskSlowLI" ~ "Task(Target)",
      contrast == "qualitywarm" ~ "Quality(Warm)",
      contrast == "distance_size" ~ "Distance",
      contrast == "trial" ~ "Trial_n",
      contrast == "qualitywarm:taskSlowLI" ~ "Quality(Warm):Task(Target)",
      contrast == "distance_size:taskSlowLI" ~ "Distance:Task(Target)",
      contrast == "distance_size:qualitywarm" ~ "Distance:Quality(Warm)",
      contrast == "distance_size:qualitywarm:taskSlowLI" ~ "Distance:Quality(Warm):\nTask(Target)",
      TRUE ~ contrast
    ))
    references <- c(
      "Quality" = "Warm",
      "Task" = "Fixed"
    )
  }
  
  
  

  if(formular[2] == "rating_scaled" ){
    formular[2] = "Rating"
  }
  
  if(formular[2] == "corrected_RT"){
    formular[2] = "Response time"
  }
  
  
  #formating and rounding the numeric values:
  fixedeffecs[, 3:5] <- apply(fixedeffecs[, 3:5], 2, function(x) {
    formatC(x, format = "f", digits = 2)
  })
  
  library(stringr)
  
  fixedeffecs[, 6] <- ifelse(fixedeffecs[, 6] < 0.001, "< 0.001", 
                             ifelse(fixedeffecs[, 6] < 0.01, "< 0.01", 
                                    ifelse(fixedeffecs[, 6] < 0.05, "< 0.05", 
                                           formatC(fixedeffecs[, 6], format = "f", digits = 3))))
  
  # Add leading spaces for alignment
  fixedeffecs[, 6] <- str_pad(fixedeffecs[, 6], width = 8, side = "left")
  
  ft <- flextable(fixedeffecs) %>%
    add_header_row(values = paste0(formular[2], formular[1], formular[3], ", ", family, "(link = ",link,")"), colwidths = c(ncol(fixedeffecs))) %>%
    #add_header_lines(values = title) %>%
    width(j = c(1, 3:ncol(fixedeffecs)), width = 1) %>%
    width(j = 2, width = 1.8) %>%
    fontsize(size = 10, part = "all") %>%
    theme_vanilla() %>%
    align(i = 1:2, j = NULL, align = "center", part = "header")
  
  
  return(ft)
}
