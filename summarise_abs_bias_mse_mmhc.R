#Setup
packages <- c("tidyverse","future.apply","fst",
              "AIPW","tmle","here","progressr","broom","RColorBrewer",
              "lme4",'lmerTest',"forcats"
)
#remotes::install_github("yqzhong7/AIPW")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

#############-----------MMHC Default----------############
#load fst
est_df_param <- read.fst(here("Analysis","est_df_param_cd.fst"))

est_df_param_unadj <- read.fst(here("Analysis","est_df_param_wo_cd_all.fst")) %>%
  filter(CovSet == 'Unadjusted') %>%
  mutate(listNo = unlist(lapply(as.character(1:100), rep, 8640*3))) %>%
  select(listNo, simNo, ATE, est) %>%
  rename(unadj_est = est)

#load param
load(here("Analysis","param_df.RData"))
load(here("Analysis",'param_list_df.RData'))

#summarise over 2000MC per DAG
rct_sum_mmhc_recode <- function(df = est_df_param){
  rct_sum <- df %>%
    filter(ATE=='RD' & !(CovSet %in% c("Parent","Swapped Parent"))) %>%
    select(listNo:Method, starts_with("OR"), m_bias, ATE, est, true_est) %>%
    full_join(est_df_param_unadj, by = c('listNo','simNo','ATE')) %>%
    mutate(unadj = ifelse(is.na(est),1,0),
           est = ifelse(unadj==1, unadj_est,est),
           #report ancestor set only
           CovSet = fct_recode(CovSet, Ancestor="Empty")
    ) 
  #replace with unadjusted estimates when adjustment set is empty
  rct_sum_empty_replacement <- rct_sum %>%
    filter(unadj==1) %>%
    bind_rows(.,.,.,.) %>% #dont ask why im doing it..
    mutate(Method = case_when(row_number() %in% 1:(nrow(.)/4) ~ "gComp",
                              row_number() %in% (1+nrow(.)/4):(nrow(.)/2) ~ "IPW",
                              row_number() %in% (1+nrow(.)/2):(3*nrow(.)/4) ~ "TMLE",
                              row_number() %in% (1+3*nrow(.)/4):nrow(.) ~ "AIPW"
    ))
  
  rct_sum_replaced <- rct_sum %>%
    filter(unadj != 1) %>%
    bind_rows(.,rct_sum_empty_replacement) %>%
    group_by_at(vars(CovSet, Method, ATE, starts_with("OR"),m_bias)) %>%
    summarise(count = n(),
              mean_est = mean(est,na.rm = T),
              true_est = mean(true_est),
              Bias = mean(est,na.rm = T) - true_est,
              abs_Bias = abs(mean(est,na.rm = T) - true_est),
              Bias_se = sd(est,na.rm = T) / sqrt(count),
              MSE = mean((est - true_est)^2,na.rm = T),
              unadj = sum(unadj)) %>%
    left_join(param_df, by = colnames(param_df)[1:7]) %>%
    mutate_at(vars(starts_with('OR')), function(x) case_when(x<1 ~ 'neg',
                                                             x==1 ~'null',
                                                             x>1 ~ 'pos')) %>%
    mutate_at(vars(CovSet:m_bias), factor)
  
  return(rct_sum_replaced)
}

rct_sum_factor_mmhc <- rct_sum_mmhc_recode()
save(rct_sum_factor_mmhc, file = here("Analysis","summarised_output","MMHC_default.RData"))
# rm(rct_sum_factor_mmhc)
# gc()

#############-----------MMHC Tuned----------############

#load fst
est_df_param_cd <- read.fst(here("Analysis","est_df_param_cd_tuned.fst"))

# use the function above
rct_sum_factor_mmhc_tuned <- rct_sum_mmhc_recode(df = est_df_param_cd)
save(rct_sum_factor_mmhc_tuned, file = here("Analysis","summarised_output","MMHC_tuned.RData"))
