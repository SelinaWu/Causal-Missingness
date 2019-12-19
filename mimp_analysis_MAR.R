# libraries
library(dplyr)
library(tidyr)
library(magrittr)
library(MatchIt)
library(boot)
library(ggplot2)

missing_level_df <- data.frame(missing_level=c(),
              mean_difference=c(),
              p_value=c(),
              trust_share=c()) 

type <- "MAR"

for (i in c(10, 20, 30, 40, 50, 60)) {
  # read in data
  df <- read.csv(paste0(type, "/data_clean_mimp", i, ".csv"))
  
  # remove Guidelines observations
  df %<>% filter(condition != "Guidelines" & post_type == "Fake")
  
  # create dummy variable where Enhanced Guidelines = 1 and Control = 0
  df %<>% mutate(enhanced_guidelines_yn = ifelse(condition == "Enhanced Guidelines", 1, 0))
  
  # use matchit to match
  match_out <- MatchIt::matchit(enhanced_guidelines_yn ~ 
                                  believe_climate_change_q +
                                  believer_doubter +
                                  burning_oil_co2 +
                                  nuclear_power_co2 +
                                  co2_methane +
                                  global_co2_increase +
                                  climate_change_human_activities +
                                  last_century_global_increase +
                                  today_concentration_already_occured +
                                  scientists_expect_warming_increase_melting +
                                  scientists_expect_increase_extreme_events +
                                  scientists_expect_increase_water_evap +
                                  scientists_expect_climate_change_evenly +
                                  feelings_toward_fb +
                                  gender +
                                  age +
                                  education +
                                  politics +
                                  conservative
                                  , data=df,
                       method = "optimal",
                       distance = "logit",
                       ratio = 1,
                       replace = FALSE)
  
  matched_df <- MatchIt::match.data(match_out, group = 'all',
                                     weights = "weights")
  
  # check
  summary(match_out)$nn
  
  # bootstrap matches to calculate:
  # 1. p-value
  # 2. ATT
  # 3. sd of the effect
  # 4. 95% confidence interval
  
  # first put matches in wide form
  treat <- matched_df %>% filter(enhanced_guidelines_yn == 1) %>%
    select(post_trust, share_post, subclass) %>%
    rename(`post_trust_t` = post_trust, `share_post_t` = share_post)
  
  control <- matched_df %>% filter(enhanced_guidelines_yn == 0) %>%
    select(post_trust, share_post, subclass) %>%
    rename(`post_trust_c` = "post_trust", `share_post_c` = "share_post")
  
  matched_df_wide <- left_join(treat, control, by = c("subclass"="subclass")) %>%
    mutate(trust_diff = post_trust_t - post_trust_c,
           share_diff = share_post_t - share_post_c)
  #row check
  nrow(matched_df_wide)
  
  calc_treatment_trust <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$post_trust_t))
  }
  
  calc_control_trust <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$post_trust_c))
  }
  
  calc_treatment_share <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$share_post_t))
  }
  
  calc_control_share <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$share_post_c))
  }
  
  calc_diff_trust <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$trust_diff))
  }
  
  calc_diff_share <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$share_diff))
  }
  
  
  boot_stats <- function(data, fn, label) {
    nc_boot <- boot(data, fn, 5000)
    print(paste("===========", label, "==============="))
    conf_int <- boot.ci(nc_boot, 0.95, type = "perc")
    point_estimate <- round(median(nc_boot$t), 3)
    conf_int_lower <- round(conf_int$percent[4], 3)
    conf_int_upper <- round(conf_int$percent[5], 3)
    standard_error <- round(sd(nc_boot$t), 3)   
    print(paste0("The point estimate is: ", point_estimate))
    print(paste0("The 95% conf interval (two way) is: [", conf_int_lower, ",", conf_int_upper, "]"))
    print(paste0("The standard error is: ", standard_error))
    return(list(point_estimate, standard_error, conf_int_lower, conf_int_upper, nc_boot$t))
  }
  
  treated_trust <- boot_stats(matched_df_wide, calc_treatment_trust, "Trust treatment")
  treated_share <- boot_stats(matched_df_wide, calc_treatment_share, "Share treatment")
  control_trust <- boot_stats(matched_df_wide, calc_control_trust, "Trust treatment")
  control_share <- boot_stats(matched_df_wide, calc_control_share, "Share treatment")
  diff_trust <- boot_stats(matched_df_wide, calc_diff_trust, "Trust ATT")
  diff_share <- boot_stats(matched_df_wide, calc_diff_share, "Share ATT")
  
  # Find p-value
  
  share_p_value <- sum(control_share[[5]] < treated_share[[1]])/length(control_share[[5]])
  print(paste("Share ATT p-value:", share_p_value))
  trust_p_value <- sum(control_trust[[5]] < treated_trust[[1]])/length(control_trust[[5]])
  print(paste("Trust ATT p-value:", trust_p_value))
  
  # TODO Create histogram plots for share and trust, save them
  
  # add mean difference values, p-values to data frames
  trust_entry <- data.frame(missing_level=i/100,
                                 mean_difference=diff_trust[[1]],
                                 se = diff_trust[[2]],
                                 ci_lower = diff_trust[[3]],
                                 ci_upper = diff_trust[[4]],
                                 p_value=trust_p_value,
                                 trust_share="trust")
  missing_level_df <- bind_rows(missing_level_df, trust_entry)
  share_entry <- data.frame(missing_level=i/100,
                            mean_difference=diff_share[[1]],
                            se = diff_trust[[2]],
                            ci_lower = diff_trust[[3]],
                            ci_upper = diff_trust[[4]],
                            p_value=share_p_value,
                            trust_share="share")
  missing_level_df <- bind_rows(missing_level_df, share_entry)
}

# graph missing level df
missing_level_df %<>% mutate(significant_yn = ifelse(p_value < 0.05, 1, 0))
ggplot(missing_level_df, aes(x=missing_level, y=mean_difference, group=trust_share)) +
  geom_line(aes(linetype=trust_share))+
  geom_point(aes(color=significant_yn))+
  geom_hline(yintercept=-0.24702,)+ #trust RCT ATT
  geom_hline(yintercept=-0.731479)+ # share RCT ATT
  geom_errorbar(aes(ymin=mean_difference-se, ymax=mean_difference+se), color="darkseagreen4") +
  theme(legend.position="bottom")

# write out table
write.csv(missing_level_df, paste0(type,"/mimp.csv"), row.names=FALSE)

####################################

# reduced ########################


missing_level_df <- data.frame(missing_level=c(),
                               mean_difference=c(),
                               p_value=c(),
                               trust_share=c()) 

for (i in c(10, 20, 30, 40, 50, 60)) {
  # read in data
  df <- read.csv(paste0(type, "/data_clean_mimp_reduced", i, ".csv"))
  
  # remove Guidelines observations
  df %<>% filter(condition != "Guidelines" & post_type == "Fake")
  
  # create dummy variable where Enhanced Guidelines = 1 and Control = 0
  df %<>% mutate(enhanced_guidelines_yn = ifelse(condition == "Enhanced Guidelines", 1, 0))
  
  # use matchit to match
  match_out <- MatchIt::matchit(enhanced_guidelines_yn ~ 
                                  conservative +
                                  education +
                                  today_concentration_already_occured +
                                  burning_oil_co2 +
                                  gender
                                , data=df,
                                method = "optimal",
                                distance = "logit",
                                ratio = 1,
                                replace = FALSE)
  
  matched_df <- MatchIt::match.data(match_out, group = 'all',
                                    weights = "weights")
  
  # check
  summary(match_out)$nn
  
  # bootstrap matches to calculate:
  # 1. p-value
  # 2. ATT
  # 3. sd of the effect
  # 4. 95% confidence interval
  
  # first put matches in wide form
  treat <- matched_df %>% filter(enhanced_guidelines_yn == 1) %>%
    select(post_trust, share_post, subclass) %>%
    rename(`post_trust_t` = post_trust, `share_post_t` = share_post)
  
  control <- matched_df %>% filter(enhanced_guidelines_yn == 0) %>%
    select(post_trust, share_post, subclass) %>%
    rename(`post_trust_c` = "post_trust", `share_post_c` = "share_post")
  
  matched_df_wide <- left_join(treat, control, by = c("subclass"="subclass")) %>%
    mutate(trust_diff = post_trust_t - post_trust_c,
           share_diff = share_post_t - share_post_c)
  #row check
  nrow(matched_df_wide)
  
  calc_treatment_trust <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$post_trust_t))
  }
  
  calc_control_trust <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$post_trust_c))
  }
  
  calc_treatment_share <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$share_post_t))
  }
  
  calc_control_share <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$share_post_c))
  }
  
  calc_diff_trust <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$trust_diff))
  }
  
  calc_diff_share <- function(data, boot_idx) {
    data_boot <- data[boot_idx,]
    return(mean(data_boot$share_diff))
  }
  
  
  boot_stats <- function(data, fn, label) {
    nc_boot <- boot(data, fn, 5000)
    print(paste("===========", label, "==============="))
    conf_int <- boot.ci(nc_boot, 0.95, type = "perc")
    point_estimate <- round(median(nc_boot$t), 3)
    conf_int_lower <- round(conf_int$percent[4], 3)
    conf_int_upper <- round(conf_int$percent[5], 3)
    standard_error <- round(sd(nc_boot$t), 3)   
    print(paste0("The point estimate is: ", point_estimate))
    print(paste0("The 95% conf interval (two way) is: [", conf_int_lower, ",", conf_int_upper, "]"))
    print(paste0("The standard error is: ", standard_error))
    return(list(point_estimate, standard_error, conf_int_lower, conf_int_upper, nc_boot$t))
  }
  
  treated_trust <- boot_stats(matched_df_wide, calc_treatment_trust, "Trust treatment")
  treated_share <- boot_stats(matched_df_wide, calc_treatment_share, "Share treatment")
  control_trust <- boot_stats(matched_df_wide, calc_control_trust, "Trust treatment")
  control_share <- boot_stats(matched_df_wide, calc_control_share, "Share treatment")
  diff_trust <- boot_stats(matched_df_wide, calc_diff_trust, "Trust ATT")
  diff_share <- boot_stats(matched_df_wide, calc_diff_share, "Share ATT")
  
  # Find p-value
  
  share_p_value <- sum(control_share[[5]] < treated_share[[1]])/length(control_share[[5]])
  print(paste("Share ATT p-value:", share_p_value))
  trust_p_value <- sum(control_trust[[5]] < treated_trust[[1]])/length(control_trust[[5]])
  print(paste("Trust ATT p-value:", trust_p_value))
  
  # TODO Create histogram plots for share and trust, save them
  
  # add mean difference values, p-values to data frames
  trust_entry <- data.frame(missing_level=i/100,
                            mean_difference=diff_trust[[1]],
                            se = diff_trust[[2]],
                            ci_lower = diff_trust[[3]],
                            ci_upper = diff_trust[[4]],
                            p_value=trust_p_value,
                            trust_share="trust")
  missing_level_df <- bind_rows(missing_level_df, trust_entry)
  share_entry <- data.frame(missing_level=i/100,
                            mean_difference=diff_share[[1]],
                            se = diff_trust[[2]],
                            ci_lower = diff_trust[[3]],
                            ci_upper = diff_trust[[4]],
                            p_value=share_p_value,
                            trust_share="share")
  missing_level_df <- bind_rows(missing_level_df, share_entry)
}

# graph missing level df
missing_level_df %<>% mutate(significant_yn = ifelse(p_value < 0.05, 1, 0))
ggplot(missing_level_df, aes(x=missing_level, y=mean_difference, group=trust_share)) +
  geom_line(aes(linetype=trust_share))+
  geom_point(aes(color=significant_yn))+
  geom_hline(yintercept=-0.24702,)+ #trust RCT ATT
  geom_hline(yintercept=-0.731479)+ # share RCT ATT
  geom_errorbar(aes(ymin=mean_difference-se, ymax=mean_difference+se), color="darkseagreen4") +
  theme(legend.position="bottom")

# write out table
write.csv(missing_level_df, paste0(type,"/mimp_reduced.csv"), row.names=FALSE)

