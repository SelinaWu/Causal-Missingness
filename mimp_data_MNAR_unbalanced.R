# libraries
library(dplyr)
library(tidyr)
library(magrittr)
library(mice)
library(stringdist)

# make complete data sets using MIMP

for (percent in c("10", "20", "30", "40", "50", "60")) {
  data <- read.csv(paste0("MNAR_unbalanced/data_clean_MNAR_unbalanced_missing_", percent,".csv"), stringsAsFactors = FALSE)
  
  # clean data
  # replace all spaces with NA
  data[data==""] = NA
  
  # remove blank rows after 2750
  data = data[1:2750,]
  
  # remove gender_text, politics_text, article, deleted_responses, consent, consent2
  # data %<>% select(-gender_text, 
  #                  -politics_text, 
  #                  -article, 
  #                  -deleted_responses, 
  #                  -consent,
  #                  -consent2)
  
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
  
  # create individual S missing pattern columns for each covariate
  s_cols <- covariates %>% mutate_all(list(~ ifelse(is.na(.), 1, 0)))
  colnames(s_cols) <- paste0("S_", colnames(s_cols))
  s_cols_original <- s_cols
  
  # create omnibus S
  s_cols$S_omnibus <- do.call(paste0, s_cols)
  
  
  #======================================================================
  
  # cluster each pattern to make sure each count is above n=100
  s_groups <- s_cols %>% group_by(S_omnibus) %>% 
    summarise( count = n()) %>% 
    ungroup() %>% 
    arrange(desc(count)) %>%
    mutate(label = row_number())
  
  # hierarchically cluster groups to be > 100
  
  s_groups$new_label <- NA
  s_groups_num <- seq(1, nrow(s_groups))
  
  # Merge big groups with small groups
  for (i in s_groups_num) {
    label_obs <- s_groups %>% slice(i) %>% pull(count)
    merged <- s_groups %>% slice(i) %>% pull(new_label)
    pattern <- s_groups %>% slice(i) %>% pull(S_omnibus)
    current_label <- s_groups %>% slice(i) %>% pull(label)
  
    if (label_obs < 100 & is.na(merged)) {
      # filter for not merged, only fewer observations than label_ob
      candidates <- s_groups %>% filter(is.na(new_label) & count <= label_obs & label != current_label)
      # calculate edit distances
      candidates <- candidates %>% mutate(dist = stringdist(pattern, S_omnibus)) %>%
        arrange(dist, count)
      while (nrow(candidates) > 0 & label_obs < 100) {
        # take first
        label_to_add <- candidates %>% slice(1) %>% pull(label)
        # add new label
        s_groups[s_groups$label == label_to_add, "new_label"] <- current_label
        # add to label obs
        add_obs <- s_groups %>% filter(label == label_to_add) %>% pull(count)
        label_obs <- label_obs + add_obs
        candidates %<>% slice(-1)
      }
    }
  }
  s_groups %<>% mutate(new_label = ifelse(is.na(new_label), label, new_label))
  
  # Merge remaining small groups <100 with biggest groups
  remaining_groups <- s_groups %>% 
    group_by(new_label) %>% 
    summarise(count = sum(count)) %>%
    ungroup() %>%
    filter(count < 100)
  
  if (nrow(remaining_groups) > 0) {
    remaining_groups_num <- seq(1, nrow(remaining_groups))
    
    for (i in remaining_groups_num) {
      label_obs <- remaining_groups %>% slice(i) %>% pull(count)
      pattern <- s_groups %>% filter(label == remaining_groups$new_label[i]) %>% pull(S_omnibus)
      current_label <- remaining_groups %>% slice(i) %>% pull(new_label)
  
      # filter for not merged, only fewer observations than label_ob
      candidates <- s_groups %>% filter(count >= label_obs & label != current_label)
      # calculate edit distances
      candidates <- candidates %>% mutate(dist = stringdist(pattern, S_omnibus)) %>%
        arrange(dist, desc(count))
      while (nrow(candidates) > 0 & label_obs < 100) {
        # take first
        label_to_add <- candidates %>% slice(1) %>% pull(label)
        # add new label to all currently assigned to remaining group
        s_groups[s_groups$new_label == current_label, "new_label"] <- label_to_add
        # add to label obs
        add_obs <- s_groups %>% filter(new_label == label_to_add) %>% summarise(count = sum(count)) %>% pull(count)
        label_obs <- add_obs
        # remove candidate
        candidates %<>% slice(-1)
      }
    }
  }
  
  # Check groups with count < 100
  check_new_groups <- s_groups %>% 
    group_by(new_label) %>% 
    summarise(count = sum(count)) %>%
    ungroup()
  
  # create missingness pattern column
  s_cols %<>% left_join(s_groups, by = c("S_omnibus" = "S_omnibus"))
  
  # create missigness pattern
  covariates$missingness_pattern <- as.factor(s_cols$new_label)
  
  # run mice
  covariates_filled <- complete(mice(data=covariates))
  
  # combine covariates with experimental_variables
  complete_df <- bind_cols(covariates_filled, experimental_variables)
  complete_df <- bind_cols(complete_df, s_cols_original)
  write.csv(complete_df, paste0("MNAR_unbalanced/data_clean_mimp", percent,".csv"), row.names=FALSE)
}

# reduced #############################
for (percent in c("10", "20", "30", "40", "50", "60")) {
  data <- read.csv(paste0("MNAR_unbalanced/data_clean_MNAR_unbalanced_missing_", percent,".csv"), stringsAsFactors = FALSE)
  
  # clean data
  # replace all spaces with NA
  data[data==""] = NA
  
  # remove blank rows after 2750
  data = data[1:2750,]
  
  # remove gender_text, politics_text, article, deleted_responses, consent, consent2
  # data %<>% select(-gender_text, 
  #                  -politics_text, 
  #                  -article, 
  #                  -deleted_responses, 
  #                  -consent,
  #                  -consent2)
  
  # condition, post_type, post_trust, post_trustworthy, post_accurate, like_post, share_post, have_facebook_q are not covariates
  covariates <- data %>% select(gender, conservative,
                                today_concentration_already_occured,
                                burning_oil_co2,
                                education)
  
  experimental_variables <- data %>% select(condition,
                                            post_type,
                                            post_trust,
                                            share_post)
  
  covariates[sapply(covariates, is.character)] <- lapply(covariates[sapply(covariates, is.character)], as.factor)
  
  # create individual S missing pattern columns for each covariate
  s_cols <- covariates %>% mutate_all(list(~ ifelse(is.na(.), 1, 0)))
  colnames(s_cols) <- paste0("S_", colnames(s_cols))
  s_cols_original <- s_cols
  
  # create omnibus S
  s_cols$S_omnibus <- do.call(paste0, s_cols)
  
  
  #======================================================================
  
  # cluster each pattern to make sure each count is above n=100
  s_groups <- s_cols %>% group_by(S_omnibus) %>% 
    summarise( count = n()) %>% 
    ungroup() %>% 
    arrange(desc(count)) %>%
    mutate(label = row_number())
  
  # hierarchically cluster groups to be > 100
  
  s_groups$new_label <- NA
  s_groups_num <- seq(1, nrow(s_groups))
  
  # Merge big groups with small groups
  for (i in s_groups_num) {
    label_obs <- s_groups %>% slice(i) %>% pull(count)
    merged <- s_groups %>% slice(i) %>% pull(new_label)
    pattern <- s_groups %>% slice(i) %>% pull(S_omnibus)
    current_label <- s_groups %>% slice(i) %>% pull(label)
    
    if (label_obs < 100 & is.na(merged)) {
      # filter for not merged, only fewer observations than label_ob
      candidates <- s_groups %>% filter(is.na(new_label) & count <= label_obs & label != current_label)
      # calculate edit distances
      candidates <- candidates %>% mutate(dist = stringdist(pattern, S_omnibus)) %>%
        arrange(dist, count)
      while (nrow(candidates) > 0 & label_obs < 100) {
        # take first
        label_to_add <- candidates %>% slice(1) %>% pull(label)
        # add new label
        s_groups[s_groups$label == label_to_add, "new_label"] <- current_label
        # add to label obs
        add_obs <- s_groups %>% filter(label == label_to_add) %>% pull(count)
        label_obs <- label_obs + add_obs
        candidates %<>% slice(-1)
      }
    }
  }
  s_groups %<>% mutate(new_label = ifelse(is.na(new_label), label, new_label))
  
  # Merge remaining small groups <100 with biggest groups
  remaining_groups <- s_groups %>% 
    group_by(new_label) %>% 
    summarise(count = sum(count)) %>%
    ungroup() %>%
    filter(count < 100)
  
  if (nrow(remaining_groups) > 0) {
    remaining_groups_num <- seq(1, nrow(remaining_groups))
    
    for (i in remaining_groups_num) {
      label_obs <- remaining_groups %>% slice(i) %>% pull(count)
      pattern <- s_groups %>% filter(label == remaining_groups$new_label[i]) %>% pull(S_omnibus)
      current_label <- remaining_groups %>% slice(i) %>% pull(new_label)
      
      # filter for not merged, only fewer observations than label_ob
      candidates <- s_groups %>% filter(count >= label_obs & label != current_label)
      # calculate edit distances
      candidates <- candidates %>% mutate(dist = stringdist(pattern, S_omnibus)) %>%
        arrange(dist, desc(count))
      while (nrow(candidates) > 0 & label_obs < 100) {
        # take first
        label_to_add <- candidates %>% slice(1) %>% pull(label)
        # add new label to all currently assigned to remaining group
        s_groups[s_groups$new_label == current_label, "new_label"] <- label_to_add
        # add to label obs
        add_obs <- s_groups %>% filter(new_label == label_to_add) %>% summarise(count = sum(count)) %>% pull(count)
        label_obs <- add_obs
        # remove candidate
        candidates %<>% slice(-1)
      }
    }
  }
  
  # Check groups with count < 100
  check_new_groups <- s_groups %>% 
    group_by(new_label) %>% 
    summarise(count = sum(count)) %>%
    ungroup()
  
  # create missingness pattern column
  s_cols %<>% left_join(s_groups, by = c("S_omnibus" = "S_omnibus"))
  
  # create missigness pattern
  covariates$missingness_pattern <- as.factor(s_cols$new_label)
  
  # run mice
  covariates_filled <- complete(mice(data=covariates))
  
  # combine covariates with experimental_variables
  complete_df <- bind_cols(covariates_filled, experimental_variables)
  complete_df <- bind_cols(complete_df, s_cols_original)
  write.csv(complete_df, paste0("MNAR_unbalanced/data_clean_mimp_reduced", percent,".csv"), row.names=FALSE)
}
