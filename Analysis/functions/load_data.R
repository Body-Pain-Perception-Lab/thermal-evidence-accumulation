########################################################################################
#                                                                                      #
#          Function to load BIDS tsv data into data frames based on task name          #
#                                                                                      #
########################################################################################
# Define function to load data 
# Function which will find all data sets from each participant in the BIDS format, which corresponds to a specific task name. 
# Inputs: 
# WD: The working directory pointing to the BIDS folder containing the data files 
# task: The name of the task for which data files will be loaded
load_data <- function(WD, task) {
  # Define the path to the task specific behavioral data files
  data_path <- file.path(WD, sprintf("sub-*/ses-01/beh/*-%s_beh", task))

  # Find all files that match the naming convention in BIDS
  files <- Sys.glob(data_path)
  
  # Load the first file as a data frame
  df <- read_tsv(files[1])
  
  # Bind all other files from the list to the first data frame 
  for (f in files[-1]) df <- rbind(df, read_tsv(f)) 

    # Return the combined data frame
  return(df)
} 
