# libraries
library(dplyr)
library(tidyr)
library(magrittr)
library(mice)

# make complete data sets using MICE

type <- "MAR" 
# import data
for (i in c("10", "20", "30", "40", "50", "60")) {
  data <- read.csv(paste0(type, "/data_clean_",type,"_missing_",i,".csv"), stringsAsFactors = FALSE)
  
  # clean data
  # replace all spaces with NA
  data[data==""] = NA
  
  # remove blank rows after 2750
  data = data[1:2750,]
  
  # remove gender_text, politics_text, article, deleted_responses, consent, consent2
  # data %<>% select(-gender_text, 
  #                     -politics_text, 
  #                     -article, 
  #                     -deleted_responses, 
  #                     -consent,
  #                     -consent2)
  # 
  # condition, post_type, post_trust, post_trustworthy, post_accurate, like_post, share_post, have_facebook_q are not covariates
  covariates <- data %>% select(-condition,
                                   -post_type,
                                   -post_trust,
                                   -share_post)
  
  experimental_variables <- data %>% select(condition,
                                            post_type,
                                            post_trust,
                                            share_post)
  
  covariates[sapply(covariates, is.character)] <- lapply(covariates[sapply(covariates, is.character)], as.factor)
  
  # run mice
  covariates_filled <- complete(mice(data=covariates))
  
  # combine covariates with experimental_variables
  complete_df <- bind_cols(covariates_filled, experimental_variables)
  write.csv(complete_df, paste0(type,"/data_clean_mice",i,".csv"), row.names=FALSE)
}

# reduced
for (i in c("10", "20", "30", "40", "50", "60")) {
  data <- read.csv(paste0(type, "/data_clean_",type,"_missing_",i,".csv"), stringsAsFactors = FALSE)
  
  # clean data
  # replace all spaces with NA
  data[data==""] = NA
  
  # remove blank rows after 2750
  data = data[1:2750,]
  
  # remove gender_text, politics_text, article, deleted_responses, consent, consent2
  # data %<>% select(-gender_text, 
  #                     -politics_text, 
  #                     -article, 
  #                     -deleted_responses, 
  #                     -consent,
  #                     -consent2)
  # 
  # condition, post_type, post_trust, post_trustworthy, post_accurate, like_post, share_post, have_facebook_q are not covariates
  covariates <- data %>% select(gender,
                                education,
                                conservative,
                                today_concentration_already_occured,
                                burning_oil_co2)
  
  experimental_variables <- data %>% select(condition,
                                            post_type,
                                            post_trust,
                                            share_post)
  
  covariates[sapply(covariates, is.character)] <- lapply(covariates[sapply(covariates, is.character)], as.factor)
  
  # run mice
  covariates_filled <- complete(mice(data=covariates))
  
  # combine covariates with experimental_variables
  complete_df <- bind_cols(covariates_filled, experimental_variables)
  write.csv(complete_df, paste0(type,"/data_clean_mice_reduced",i,".csv"), row.names=FALSE)
}
