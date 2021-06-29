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

sim_file_names <- grep(".RData",list.files(here("Simulation","2021-03-27T125752EDT")), value = T)
load(here("Simulation","2021-03-27T125752EDT",sim_file_names[1]))

# #--------load parameters---------####
param_list <- lapply(simList, function(x) attributes(x)$param)
param_list_df <- matrix(unlist(param_list), length(param_list), byrow=T) %>%
  data.frame() %>%
  rename_all(~names(param_list[[1]])) %>%
  mutate(simNo = 1:nrow(.))

param_list_df %>% unique()
save(param_list_df, file = here("Analysis",'param_list_df.RData'))

param_df <- param_list_df %>% select(starts_with("OR")) %>% unique() %>% mutate(DGM_param=1:nrow(.))
save(param_df, file = here("Analysis",'param_df.RData'))


rm(param_list)
rm(simList)

#-------load true values--------####
load(here("Analysis",'param_df.RData'))
load(here("Simulation","sim_rct_true_est2021-03-27T125752EDT.RData"))
true_est_df <- bind_rows(sim_true_est_list)
true_est_df_long <- true_est_df %>%
  rename_at(vars(starts_with('true')), ~c('RD','logRR','logOR')) %>%
  pivot_longer(cols = RD:logOR, names_to = 'ATE', values_to = 'true_est')

#--------load output------------####
est_file_names <- grep(".RData",list.files(here("Estimation","CD_MMHC_tuned"), full.names = T), value = T)
est_df_list <- list()
mmhc_tuned_out <- list()
(s_time <- Sys.time())
for (i in 1:length(est_file_names)){
  print(i)
  si_time <- Sys.time()
  load(est_file_names[i])
  #extract mmhc output
  mmhc_tuned_out[[i]] <- lapply(estList, function(x) as.data.frame(attributes(x)$edge_list))
  mmhc_tuned_out[[i]] <- bind_rows(mmhc_tuned_out[[i]], .id='simNo')
  #extract estimates
  est_df_list[[i]] <- bind_rows(estList, .id= "simNo")

  print(Sys.time()-si_time)
  }

#for estimates
est_df <- bind_rows(est_df_list, .id = 'listNo')

est_df_param_cd <- est_df %>%
  mutate(simNo = as.numeric(simNo)) %>%
  left_join(param_list_df, by = "simNo") %>%
  mutate_at(vars(RD:logOR), as.numeric) %>%
  pivot_longer(cols = RD:logOR, names_to = 'ATE', values_to = "est") %>%
  left_join(true_est_df_long, by = c(colnames(param_list_df)[1:7],"ATE")) %>%
  mutate(m_bias = as.factor(as.numeric(ORyc1 == 1 & ORxc1 ==1)))

#for mmhc
mmhc_tuned_df <- bind_rows(mmhc_tuned_out, .id = 'listNo') %>%
  mutate(simNo = as.numeric(simNo)) %>%
  left_join(param_list_df, by = "simNo")
(e_time <- Sys.time() - s_time)

#combine and save into fst
write.fst(mmhc_tuned_df,  here("Analysis","mmhc_tuned_df.fst"))
write.fst(est_df_param_cd, here("Analysis","est_df_param_cd_tuned.fst"))


#load fst
est_df_param_cd <- read.fst(here("Analysis","est_df_param_cd_tuned.fst"))

est_df_param_unadj <- read.fst(here("Analysis","est_df_param_wo_cd_all.fst")) %>%
  filter(CovSet == 'Unadjusted') %>%
  mutate(listNo = unlist(lapply(as.character(1:100), rep, 8640*3))) %>%
  select(listNo, simNo, ATE, est) %>%
  rename(unadj_est = est)


#load param
load(here("Analysis","param_df.RData"))
load(here("Analysis",'param_list_df.RData'))
#load mmhc
mmhc_tuned_df <- read.fst(here("Analysis","mmhc_tuned_df.fst"))

#--------------Summary_CD----------####
### block backdoor == 1
# 1. m-bias==1:  C1 not adjusted OR (C1 + (loss_num and/or hs_CRP_level))
# 2. m-bias==0 & ORxc1 == 1: (C1 and hs_CRP_level) OR loss_num
# 3. m-bias==0 & ORyc1 == 1: (C1 and loss_num) OR hs_CRP_level
# 4. else: C1 and (loss_num and/or hs_CRP_level)
# helper function for recoding
mmhc_recode <- function(df=mmhc_df, child = c('X',"Y")){
  mmhc_df_xy <- df %>%
    filter(to %in% c(child) & !(from %in% c("treatment","X","Y"))) %>%
    dplyr::select(-to,-N,-marginalProb) 
  
  if (length(child) == 2){
    mmhc_df_xy <- mmhc_df_xy %>%distinct() 
  } 
  
  
  mmhc_df_xy_wide <- mmhc_df_xy %>%
    pivot_wider(id_cols = c(listNo, simNo, ORyx:ORc1c3),names_from = c(from), values_from = from) %>%
    mutate_at(vars(C1:last_col()),function(x) as.numeric(!is.na(x))) %>%
    mutate(DGM_type = case_when(ORxc1==1 & ORyc1== 1 ~ 'M-bias',
                                ORxc1==1 & ORyc1!= 1 ~ 'Right-tri',
                                ORxc1!=1 & ORyc1== 1 ~ 'Left-tri',
                                ORxc1!=1 & ORyc1!= 1 ~ 'Butterfly') %>% factor())
  if (is.null(mmhc_df_xy_wide$hsCRP_level)){
    mmhc_df_xy_wide$hsCRP_level <- 0
  }
  mmhc_df_xy_wide <-  mmhc_df_xy_wide %>% mutate (block_backdoor = ifelse( 
    (DGM_type == 'M-bias' & (C1==0 | (C1==1 & (loss_num==1 | hsCRP_level==1)) )) | 
      (DGM_type == 'Right-tri' & (loss_num==1  | (C1==1 & hsCRP_level==1)  )) | 
      (DGM_type == 'Left-tri' & (hsCRP_level==1   | (C1==1 & loss_num==1) )) |
      (DGM_type == 'Butterfly' & (C1==1 & (loss_num==1 | hsCRP_level==1)))
    , 1,0))
  
  
  if (length(child) == 1){
    names(mmhc_df_xy_wide)[ncol(mmhc_df_xy_wide)]<- paste0('block_backdoor.',child)
  } 
  
  return(mmhc_df_xy_wide)
}

mmhc_df_x <- mmhc_recode(df = mmhc_tuned_df, child = 'X') 
mmhc_df_y <- mmhc_recode(df = mmhc_tuned_df, child = 'Y') 
mmhc_df_xy <- mmhc_recode(df = mmhc_tuned_df, child = c('X','Y')) 


#create a NA df for unadjusted scenarios
mmhc_df_unadj <- crossing(listNo = as.character(1:100), simNo = 1:8640 ) %>%
  left_join(param_list_df, by = 'simNo') %>%
  select(-N, -marginalProb) %>%
  mutate(DGM_type = case_when(ORxc1==1 & ORyc1== 1 ~ 'M-bias',
                              ORxc1==1 & ORyc1!= 1 ~ 'Right-tri',
                              ORxc1!=1 & ORyc1== 1 ~ 'Left-tri',
                              ORxc1!=1 & ORyc1!= 1 ~ 'Butterfly') %>% factor()
  ) 

mmhc_df_xy_sum <-  mmhc_df_xy %>%
  left_join(mmhc_df_x , by =colnames(.[1:16])) %>%
  left_join(mmhc_df_y, by =colnames(.[1:16])) %>%
  full_join(mmhc_df_unadj, by = c(colnames(.[1:9]),'DGM_type')) %>%
  mutate(across(c(C1:employed,starts_with('block_')), ~replace_na(.x, 0)),
         unadj = if_else((C1==0 & loss_num==0 & hsCRP_level==0 & eligibility==0 & smoke==0 & employed == 0), 1, 0)
         ) %>%
  mutate_at(vars(starts_with('block_')), function(x) ifelse(.$DGM_type=='M-bias', x+.$unadj, x))

mmhc_df_xy_sum_tuned <-mmhc_df_xy_sum
save(mmhc_df_xy_sum_tuned, file = here("Analysis","summarised_output","accuracy_MMHC_tuned.RData"))

mean(mmhc_df_xy_sum$block_backdoor)

