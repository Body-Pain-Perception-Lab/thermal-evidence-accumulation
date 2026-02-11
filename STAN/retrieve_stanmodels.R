
# function that takes an osftoken and downloads the two STAN-models from OSF
get_models <- function(osftoken){
  # get the preanalyzed models from osf:
  if (dir.exists(here::here("results"))) {
    osfr::osf_auth(token = osftoken)
    
    temp <- osfr::osf_retrieve_node("https://osf.io/ev2mk/")
    
    downloaded = temp %>%
      osfr::osf_ls_files(pattern = "taskmodels") %>%
      osfr::osf_download(path = here::here("STAN"), recurse = TRUE, conflicts = "overwrite", progress = TRUE)
    
    
  }
  
}

#insert the path to your person OSF token here 

get_models(read_lines(here::here("osf","osf.txt"))[1])