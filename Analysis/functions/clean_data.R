########################################################################################
#                                                                                      #
#                             Function to mark outliers                                #
#                                                                                      #
########################################################################################
# Define function to remove outliers in a given data set 
# This function takes a data frame and removes outliers in the RT columns based on z scores that are higher than 3
# Input: 
# data: a data frame with data that will be summarized
# area_condition: The name of the condition which the data will be grouped by

find_outliers <- function(data, task) {
  data <- data %>% 
    # Create a new column to indicate if a button press error occurred
    mutate(button_error = if_else(accuracy == 0 & rating == 0, 1, 0)) %>%  
    dplyr::rename(corrected_RT = correctet_RT) %>% 

    # Create a column to check if the RT is an outlier 
    {if(task == 0) 
      group_by(., ID, size, quality) 
     if(task == 1) 
      group_by(., ID, area_size, distance_size, quality)
     if(task == 2)
       group_by(., id, size, intensity, quality)
      } %>%
    
    mutate(z_score = scale(corrected_RT)) %>% 
    ungroup %>%  
    mutate(outlier = if_else(z_score < 3, 0, 1)) %>%  
    
    # If the tiral has a button error or is an outlier, it should be marked as excluded 
    mutate(excluded = if_else(outlier == 1 | button_error == 1, 1, 0))
  return(data)
} 
