# function to retrieve the data from OSF. Options for getting the data is either to get the preprocessed data i.e. rerun = FALSE
# alternatively to get the raw data and preprocess it choose rerun = TRUE. The returned list contains 3 elements, the preprocessed TPL dataframe, the TPL threshold dataframe and
# lastly a dataframe of participants that were removed from the analysis given the extreme amount of missing trials.
get_data <- function(osf_token, rerun = FALSE) {
  # get the preanalyzed data from osf:
  if (!dir.exists(here::here("Analysis"))) {
    osfr::osf_auth(token = osf_token)

    TPL <- osfr::osf_retrieve_node("https://osf.io/pw956/")

    TPL %>%
      osfr::osf_ls_files(pattern = "Analysis") %>%
      osfr::osf_download(path = here::here(), recurse = TRUE, conflicts = "overwrite", progress = TRUE)
  }


  if (rerun == FALSE) {
    # load the files from the Analysis folder
    return(list(
      data = read.csv(here::here("Analysis", "Cleaned-data", "TPL_df.csv")),
      thresholds = read.csv(here::here("Analysis", "Cleaned-data", "Thresholds.csv")),
      removers = read.csv(here::here("Analysis", "Cleaned-data", "removed.csv"))
    ))
  } else {
    osfr::osf_auth(token = osf_token)

    TPL <- osfr::osf_retrieve_node("https://osf.io/pw956/")

    if (!dir.exists(here::here("matlab", "csv_files"))) {
      dir.create(here::here("matlab", "csv_files"), showWarnings = FALSE)
    }
    # get the raw data from osf and preprocess it:
    TPL %>%
      osfr::osf_ls_files(pattern = "Data") %>%
      osfr::osf_download(path = here::here("matlab", "csv_files"), recurse = TRUE, conflicts = "overwrite", progress = TRUE)



    # read the cue structure for the two contigency spaces
    cues_even <- read.csv(here::here("matlab", "csv_files", "Data", "study-level-data", "contingency_even.csv"))
    cues_uneven <- read.csv(here::here("matlab", "csv_files", "Data", "study-level-data", "contingency_uneven.csv"))

    cues_even$sequence <- NULL
    cues_uneven$sequence <- NULL

    # get the location of the thermode for each trial for each sequence.
    loc1 <- read.csv(here::here("matlab", "csv_files", "Data", "study-level-data", "location_uneven.csv"), header = F)
    loc4 <- read.csv(here::here("matlab", "csv_files", "Data", "study-level-data", "location_even.csv"), header = F)

    # This is the directory to where the raw behavioral data is stored both the experimental TPL and the first QST thresholding

    datapath <- here::here("matlab", "csv_files", "Data")
    datafiles <- list.files(path = datapath, pattern = "*.tsv$", recursive = T)

    # initializing dataframes
    qst <- data.frame(NULL)
    data <- data.frame(NULL)
    # reading in the dataframes:
    for (i in 1:length(datafiles)) {
      if (str_detect(datafiles[i], "qst")) {
        qst1 <- read.delim(here(datapath, datafiles[i]))
        qst <- rbind(qst, qst1)
      } else if (str_detect(datafiles[i], "tpl")) {
        painlearning1 <- read.delim(here(datapath, datafiles[i]))
        if (painlearning1[1, 4] == 1) {
          painlearning1 <- cbind(painlearning1, cues_uneven)
          painlearning1$location <- loc1$V1
        } else if (painlearning1[1, 4] == 4) {
          painlearning1 <- cbind(painlearning1, cues_even)
          painlearning1$location <- loc4$V1
        } else {
          print("Error")
        }
        data <- rbind(data, painlearning1)
      }
    }
    # converting factors to factors:
    columns <- c("cohort", "stim", "predResp", "predAcc", "block_type", "sequence", "desired_prob", "blockNumber", "cue", "location", "id")
    data[, columns] <- lapply(data[, columns], as.factor)


    # function to determine participants that should be excluded. Where the argument prob is the proportion of missed trials in %
    make_remover <- function(data, prob) {
      # count the number of times a participant missed a response:
      g <- data %>%
        filter(as.character(predResp) == "NaN") %>%
        dplyr::group_by(id) %>%
        dplyr::summarize(n())

      # find all the times a participant had at least one missing values in their VAS-ratings
      q <- data %>% filter(is.nan(vasResp_1) | is.nan(vasResp_2) | is.nan(vasResp_3))

      # now counting the times where all VAS-responses were missing
      q$sort <- ifelse(is.nan(q$vasResp_1) & is.nan(q$vasResp_2) & is.nan(q$vasResp_3), 0, 1)
      q <- q %>% 
        filter(sort == 1) %>%
        dplyr::group_by(id) %>%
        dplyr::summarize(n())

      # now finding the participants that have missing trials above the threshold
      g <- data.frame(g) %>% filter(n.. > 306 * prob)
      q <- data.frame(q) %>% filter(n.. > 306 * prob)


      removers <- rbind(g, q)
    }
    # find the participants above 10% missing trials (either VAS or predictions):
    removers <- make_remover(data, 0.1)
    data <- data %>% filter(id %in% removers$id == F)
    thresh <- qst %>% filter(id %in% removers$id == F)

    # renaming factors
    levels(data$cue) <- c("low-tone", "high-tone")
    levels(data$stim) <- c("cold", "warm", "TGI", "NaN")
    levels(data$predResp) <- c("cold", "warm", "NaN")
    levels(data$predAcc) <- c("Wrong", "Right", "TGI", "Wrong")


    # making a column that gives the expectedness of the stimuli, that is given the underlying probability and the cue was the outcome to be expected (E), unexpected(UE) or neutral(N).
    # all stimuli in 0.5 blocks are neutral and all TGI trials are coded as TGI. Here its important to note that the probabilities in desired_prob reflect the probability of a
    # high-tone to give a cold stimulus
    data$expected <- ifelse(data$desired_prob == "0.5" & data$stim != "TGI", a <- "N",
      ifelse(data$desired_prob == "0.18" & data$cue == "low-tone" & data$stim == "cold", a <- "E",
        ifelse(data$desired_prob == "0.18" & data$cue == "high-tone" & data$stim == "warm", a <- "E",
          ifelse(data$desired_prob == "0.82" & data$cue == "low-tone" & data$stim == "warm", a <- "E",
            ifelse(data$desired_prob == "0.82" & data$cue == "high-tone" & data$stim == "cold", a <- "E", a <- "UE")
          )
        )
      )
    )
    
    data$expected <- ifelse(data$stim == "TGI", a <- "TGI", b <- data$expected)

    data$expected <- as.factor(data$expected)

    write.csv(removers %>% mutate(id = as.numeric(as.character(id)), n.. = as.numeric(n..)), here::here("matlab", "removers.csv"))

    return(list(data = data, thresholds = thresh, removers = removers))
  }
}
# function to check simulated residudals of a fitted model:
residual_check <- function(model) {
  model_sim <- DHARMa::simulateResiduals(model)
  return(plot(model_sim))
}

# convinient function to make p-values to desired format
make_pvalue <- function(p_value) {
  if (p_value > 0.05) {
    p_value <- round(p_value, 2)
    p <- paste("p = ", p_value)
  }
  if (p_value < 0.05) {
    p <- "p < .05 "
  }
  if (p_value < 0.01) {
    p <- "p < .01"
  }
  if (p_value < 0.001) {
    p <- "p < .001"
  }
  if (p_value < 0.0001) {
    p <- "p < .0001"
  }
  if(p_value == 0.05){
    p = "p = 0.05"
  }
  return(p)
}

# get summary statistics from a generalized linear mixed effects model.
# Arguments: model is the model to get statistics on. Coefficients is the index of the coefficients when using summary(model) (that is excluding the intercept and until the coefficients arguemnt),
# round is the rounding of the statistics.
summary_stats <- function(model, coefficients, round, intercept = FALSE) {
  if (intercept == TRUE) {
    a <- summary(model)
    coef <- a$coefficients$cond[1, 1]
    std <- a$coefficients$cond[1, 2]
    stat <- a$coefficients$cond[1, 3]
    p <- a$coefficients$cond[1, 4]

    return(list(beta = round(coef, round), std = round(std, round), stat = round(stat, round), p = p))
  }
  a <- summary(model)
  coef <- array(NA, coefficients)
  std <- array(NA, coefficients)
  stat <- array(NA, coefficients)
  p <- array(NA, coefficients)

  for (i in 1:coefficients) {
    coef[i] <- a$coefficients$cond[1 + i, 1]
    std[i] <- a$coefficients$cond[1 + i, 2]
    stat[i] <- a$coefficients$cond[1 + i, 3]
    p[i] <- a$coefficients$cond[1 + i, 4]
  }
  return(list(beta = as.numeric(format(round(coef, round), nsmall = round)), std = as.numeric(format(round(std, round), nsmall = round)), stat = as.numeric(format(round(stat, round), nsmall = round)), p = p))
}


summary_stats_new <- function(model, coefficients,modelterm, round, intercept = FALSE) {
  
  updating = as.formula(paste(". ~ . - ", names(model$frame)[modelterm]))
  
  partialmodel <- update(model, updating)
  
  full = MuMIn::r.squaredGLMM(model)
  partial = MuMIn::r.squaredGLMM(partialmodel)
  
  
  cohens_f = (full[1,1]-partial[1,1])/(1-full[1,1])
  
  if (intercept == TRUE) {
    a <- summary(model)
    coef <- a$coefficients$cond[1, 1]
    std <- a$coefficients$cond[1, 2]
    stat <- a$coefficients$cond[1, 3]
    p <- a$coefficients$cond[1, 4]
    
    return(list(beta = round(coef, round), std = round(std, round), stat = round(stat, round), p = p, cohens_f))
  }
  a <- summary(model)
  coef <- array(NA, coefficients)
  std <- array(NA, coefficients)
  stat <- array(NA, coefficients)
  p <- array(NA, coefficients)
  
  for (i in 1:coefficients) {
    coef[i] <- a$coefficients$cond[1 + i, 1]
    std[i] <- a$coefficients$cond[1 + i, 2]
    stat[i] <- a$coefficients$cond[1 + i, 3]
    p[i] <- a$coefficients$cond[1 + i, 4]
  }
  return(list(beta = as.numeric(format(round(coef, round), nsmall = round)),
              std = as.numeric(format(round(std, round), nsmall = round)),
              stat = as.numeric(format(round(stat, round), nsmall = round)),
              p = p,
              cohens_f = as.numeric(format(round(cohens_f, round), nsmall = round))))
}


# same function as above but for the zero, one inflated beta regressions fit. Here the part argument is whether the coefficients should
# be from the mean of the beta or the probability parameters for the one or zero inflated part.
summary_stats_zoib <- function(model, coefficients, round, part, intercept = FALSE) {
  if (intercept == TRUE) {
    a <- summary(model)
    coef <- a[1, 1]
    std <- a[1, 2]
    stat <- a[1, 3]
    p <- a[1, 4]

    return(list(beta = round(coef, round), std = round(std, round), stat = round(stat, round), p = p))
  }

  if (part == "mu") {
    index <- 1
  } else if (part == "nu") {
    index <- 1 + 2 + coefficients
  } else if (part == "tau") {
    index <- 2 + 2 + 2 * coefficients
  } else {
    print("Give correct part")
  }

  a <- summary(model)
  coef <- array(NA, coefficients)
  std <- array(NA, coefficients)
  stat <- array(NA, coefficients)
  p <- array(NA, coefficients)

  for (i in 1:coefficients) {
    coef[i] <- a[index + i, 1]
    std[i] <- a[index + i, 2]
    stat[i] <- a[index + i, 3]
    p[i] <- a[index + i, 4]
  }

  return(list(beta = round(coef, round), std = round(std, round), stat = round(stat, round), p = p))
}


summary_stats_zoib_new <- function(model, coefficients, round, part, intercept = FALSE) {
  
  dfdataframe = parameters::model_parameters(model)
  if(part == "mu"){
    part = "conditional"
  }
  
  dd = data.frame(dfdataframe) %>% filter(Component == part)
  if(intercept){
    dd = dd[1:(coefficients),]
  }else{
    dd = dd[2:(coefficients),]
  }

  beta = round(dd$Coefficient,round)
  se = round(dd$SE, round)
  CI_low = round(dd$CI_low, round)
  CI_high = round(dd$CI_high, round)
  stat = round(dd$t, round)
  p = round(dd$p, round)
  df = round(dd$df_error, round)
  
  
    
  return(list(beta = beta, se = se, stat = stat,CI_low = CI_low,CI_high = CI_high,df = df, p = p))
  
}


pseudo_p_val = function(p){
  2*(min(sum(p>0)/length(p),sum(p<0)/length(p)))
}


#####################################################################################
#                                                                                   #
#                   Plot raw data and model predictions                             #
#                                                                                   #
#####################################################################################
plot_raw_predictions <- function(data, task, y, x, predict_data, x_title, y_title, colors) {
  data %>%
    ggplot() +
    # Warm boxplot
    geom_boxplot(
      aes(x = as.factor(.data[[x]]), y = .data[[y]], fill = quality, col = quality),
      width = 0.3,
      position = position_nudge(x = 0),
      outlier.shape = NA,
      alpha = 0.6
    ) +
    # Cold points
    geom_point(
      aes(x = as.factor(.data[[x]]), y = .data[[y]], fill = quality, col = quality),
      position = position_nudge(x = 0),
      alpha = 0.3
    ) +
    # Lines connecting points
    geom_line(
      aes(x = as.factor(.data[[x]]), y = .data[[y]], fill = quality, col = quality, group = ID),
      alpha = 0.2
    ) +
    # Predicted lines and ribbons
    geom_line(
      data = predict_data,
      aes(x = .data[[x]], y = predicted, group = interaction(quality, task), col = quality),
      linewidth = 1
    ) +
    geom_ribbon(
      data = predict_data,
      aes(x = .data[[x]], ymin = conf.low, ymax = conf.high, group = interaction(quality, task), fill = quality),
      alpha = 0.75
    ) +
    # Custom colors
    scale_colour_manual(values = colors) +
    scale_fill_manual(values = colors) +
    # Facet vertically
    facet_grid(~quality) +
    # Titles and theme
    labs(x = x_title, y = y_title) +
    SSLI_theme +
    theme(strip.text = element_blank()) # Remove facet titles
  
  
  #####################################################################################
  #                                                                                   #
  #                                  Report statistics                                #
  #                                                                                   #
  #####################################################################################
  make_new_reporting = function(stats, number, Z){
    
    if(Z){
      z_stat = stats$stat[number]
      p_value = stats$p[number]
      ci = list(low = as.numeric(stats$beta[number])-2*as.numeric(stats$std[number]), high = as.numeric(stats$beta[number])+2*as.numeric(stats$std[number]))
      beta_stat = stats$beta[number]
      
      text = paste0("$\\", "beta", "$"," = ",beta_stat,", 95% CI = [", ci$low, ";", ci$high,"]",", Z = ", z_stat, ", ", make_pvalue(p_value))
      
      
      return(text)
      
    }else if(!Z){
      beta_stat = stats$beta[number]
      t_stat = stats$stat[number]
      p_value = stats$p[number]
      ci = list(low = stats$CI_low[number], high = stats$CI_high[number])
      df = stats$df[number]
      
      text = paste0("$\\", "beta", "$"," = ",beta_stat, ", 95% CI = [", ci$low, ";", ci$high,"]",", t(",df,") = ",t_stat,", ",make_pvalue(p_value))
      return(text)
      
    }
    
    return("error")
  }
  
  
  summary_stats_zoib_new <- function(model, coefficients, round, part, intercept = FALSE) {
    
    dfdataframe = parameters::model_parameters(model)
    if(part == "mu"){
      part = "conditional"
    }
    
    dd = data.frame(dfdataframe) %>% filter(Component == part)
    if(intercept){
      dd = dd[1:(coefficients),]
    }else{
      dd = dd[2:(coefficients),]
    }
    
    beta = round(dd$Coefficient,round)
    se = round(dd$SE, round)
    CI_low = round(dd$CI_low, round)
    CI_high = round(dd$CI_high, round)
    stat = round(dd$t, round)
    p = round(dd$p, round)
    df = round(dd$df_error, round)
    
    
    
    return(list(beta = beta, se = se, stat = stat,CI_low = CI_low,CI_high = CI_high,df = df, p = p))
    
  }
}


format_model_results <- function(model) {
  formatted_results <- mapply(function(beta, ci_low, ci_high, p) {
    # Format p-value
    if (p < 0.001) {
      p_str <- "p < .001"
    } else {
      p_str <- sprintf("p = %.3f", p)
    }
    
    # Format beta and confidence interval
    sprintf("(β = %.2f, 95%% CI = [%.2f; %.2f], %s)", beta, ci_low, ci_high, p_str)
  }, model$beta, model$CI_low, model$CI_high, model$p, SIMPLIFY = TRUE)
  
  return(formatted_results)
}
