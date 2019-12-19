# libraries
library(dplyr)
library(tidyr)
library(magrittr)
library(mice)

# make complete data sets using MICE

# import data
for (i in c("10", "20", "30", "40", "50", "60")) {
  data <- read.csv(paste0("MCAR/data_clean_missing",i,".csv"), stringsAsFactors = FALSE)
  
  # clean data
  # replace all spaces with NA
  data[data==""] = NA
  
  # remove blank rows after 2750
  data = data[1:2750,]
  
  remove gender_text, politics_text, article, deleted_responses, consent, consent2
  data %<>% select(-gender_text,
                      -politics_text,
                      -article,
                      -deleted_responses,
                      -consent,
                      -consent2)

  # condition, post_type, post_trust, post_trustworthy, post_accurate, like_post, share_post, have_facebook_q are not covariates
  covariates <- data %>% select(-condition,
                                   -post_type,
                                   -post_trust,
                                   -post_trustworthy,
                                   -post_accurate,
                                   -like_post,
                                   -share_post,
                                   -have_facebook_q)
  
  experimental_variables <- data %>% select(condition,
                                            post_type,
                                            post_trust,
                                            like_post,
                                            share_post)
  
  covariates[sapply(covariates, is.character)] <- lapply(covariates[sapply(covariates, is.character)], as.factor)
  
  # run mice
  covariates_filled <- complete(mice(data=covariates))
  
  # combine covariates with experimental_variables
  complete_df <- bind_cols(covariates_filled, experimental_variables)
  write.csv(complete_df, paste0("MCAR/data_clean_mice",i,".csv"), row.names=FALSE)
}

