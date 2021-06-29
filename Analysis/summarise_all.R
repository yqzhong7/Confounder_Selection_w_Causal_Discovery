#Setup
packages <- c("tidyverse","future.apply","fst",
              "AIPW","tmle","here","progressr","broom","RColorBrewer",
              "lme4",'lmerTest','quantreg', 'knitr','kableExtra','scales', 'tableone',
              "gt", "gtsummary","forcats","sjlabelled"
)
#remotes::install_github("yqzhong7/AIPW")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}


#####---------------TABLE 1-----------------#####
eager_raw <- read_csv(here("Data","eager_base_imputed.csv")) 
eager_bl <- eager_raw %>%
  mutate(meanAP = (BPS + 2*BPD)/3,
         live_birth = as.numeric(outcome=="live birth"),
         time_try_pregnant=if_else(as.numeric(time_try_pregnant)>10,10,as.numeric(time_try_pregnant)),
         eligibility = as.numeric(eligibility=="new"),
         #needs to be changed?
         # doi: 10.1111/ppe.12409 [0,2), [2,10), [10,inf)
         # hsCRP_level = cut(hsCRP, breaks=c(-Inf, 10, 20, Inf), labels=c("low/medium","high","extreme high"))
         hsCRP_level = cut(hsCRP, breaks=c(-Inf, 2, 10, Inf), labels=c("low","medium","high"), right = F, include.lowest = T)
  ) %>%
  dplyr::select(loss_num,hsCRP_level,treatment,eligibility,alcohol,smoke,employed,age,time_try_pregnant,BMI,meanAP) 

covariates <- colnames(eager_bl[3:ncol(eager_bl)])
eager_bl
summary(eager_bl)
table(eager_bl$loss_num)
table(eager_bl$hsCRP_level)
tab1 <- tableone::CreateTableOne(data=eager_bl,
                         factorVars = colnames(eager_bl %>% select(-age,-time_try_pregnant,-BMI,-meanAP))) %>%
  print(showAllLevels = T)

kable(tab1,booktab=T, 'latex',
      caption = "Covariates from the EAGeR trial used for simulations") %>%
  kable_styling() 

## load summarised data
# 
load( file = here("Analysis","summarised_output","wo_CD.RData"))
load( file = here("Analysis","summarised_output","MMHC_default.RData"))
load( file = here("Analysis","summarised_output","MMHC_tuned.RData"))
#
# # create table for medians of abs bias and MSE
rct_sum_all <- bind_rows(rct_sum_factor,rct_sum_factor_mmhc,rct_sum_factor_mmhc_tuned, .id = c("CovSel")) %>%
  ungroup() %>%
  mutate(DGM_type = case_when(ORxc1=='null' & ORyc1== 'null' ~ 'M',
                              ORxc1=='null' & ORyc1!= 'null' ~ 'Right-triangle',
                              ORxc1!='null' & ORyc1== 'null' ~ 'Left-triangle',
                              ORxc1!='null' & ORyc1!= 'null' ~ 'Butterfly') %>% 
           factor(levels = c("Butterfly","M","Left-triangle","Right-triangle")) ,
         CovSel = factor(CovSel, levels = 1:3, labels = c("Manual Selection","MMHC default","MMHC tuned")),
         CovSet = factor(CovSet, levels = unique(CovSet),
                         labels = c("Direct Causes","Swapped Direct Causes",
                                    "All Causes", "All Covariates", "Collider Only",
                                    "Empty Set / Unadjusted"))
         )
# save(rct_sum_all, file = here("Analysis","summarised_output","rct_sum_all.RData"))

# load(file = here("Analysis","summarised_output","rct_sum_all.RData"))

#overall magnitude
rct_sum_all %>% 
  filter(ATE == 'RD') %>%
  summarise(across(c(mean_est,true_est, abs_Bias, MSE), list(mean = mean, min=min,max=max),.names = "{.col}.{.fn}"))


median_table <- rct_sum_all %>% 
  filter(ATE == 'RD') %>%
  group_by( CovSet, CovSel) %>% 
  summarise(#abs_Bias = mean(mean_est-true_est),
            across(.cols = c(abs_Bias,MSE), list(~median(.x), ~quantile(.x,.25),  ~quantile(.x,.75)))) %>%
  # pivot_wider(id_cols = c(CovSet, CovSel), names_from = ATE, values_from = c(abs_Bias_1:MSE_3)) %>%
  # select(CovSel, CovSet, ends_with('RD')#, ends_with('logRR'), ends_with('logOR')
  #        ) %>%
  arrange(CovSel)%>%
  ungroup() %>%
  mutate(across(matches("abs_Bias"), ~.x/sum((row_number() == 3)*.x), .names = "rel_{.col}")) %>%
  mutate(across(matches("^abs_Bias"), ~label_number(scale = 1000, accuracy = 0.01)(.x)),
         across(matches("MSE"), ~label_number(scale = 1000, accuracy = 0.01)(.x)),
         across(matches("rel_abs_Bias"),  ~label_number(accuracy = 0.01)(.x)),
         abs_Bias = paste0(abs_Bias_1, " (",abs_Bias_2, ", ", abs_Bias_3,")"),
         rel_abs_Bias = paste0(rel_abs_Bias_1, " (",rel_abs_Bias_2, ", ", rel_abs_Bias_3,")"),
         MSE = paste0(MSE_1, " (",MSE_2, ", ", MSE_3,")")
  ) %>%
  select( CovSet, abs_Bias, rel_abs_Bias, MSE)  %>%
  filter(!(CovSet %in% c("Direct Causes", "Swapped Direct Causes")))

# rct_sum_all %>% 
#   group_by(ATE, CovSet, CovSel,DGM_type) %>% 
#   summarise(across(.cols = c(abs_Bias,MSE), 
#                    list(~median(.x), ~quantile(.x,.25),  ~quantile(.x,.75)))) %>%
#   filter(ATE=='RD', CovSet == 'C1 Only')
# 
# rct_sum_all %>% 
#   group_by(ATE, CovSet, CovSel) %>% 
#   summarise(across(.cols = c(abs_Bias,MSE), 
#                    list(~median(.x), ~quantile(.x,.25),  ~quantile(.x,.75)))) %>%
#   filter(ATE=='RD', CovSet == 'C1 Only')

median_table %>%  kable('latex', booktabs = T,
        caption = "Absolute bias and MSE of average treatment effect from different covariate adjustment sets selected by different methods",
        col.names = c("CovSet", rep("Median (Q1, Q3)",3))
  ) %>%
  add_header_above(c(" " = 1, "Abs. Bias"= 1, "Relative Abs. Bias" = 1,  "MSE"= 1#, "Abs. Bias"= 1,  "MSE"= 1, "Abs. Bias"= 1,  "MSE"= 1
                     )) %>%
  # add_header_above(c(" " = 1, "RD"=2#, "logRR"=2,"logOR"=2
  #                    )) %>%
  pack_rows("Manual Covariate Selection", 1, 4) %>%
  pack_rows("Automated Covariate Selection: MMHC default", 5, 5) %>%
  pack_rows("Automated Covariate Selection: MMHC tuned", 6, 5) 


######-------Plotting abs bias and MSE by scenario method--------######
rct_sum_plot_df <- rct_sum_all %>%
  filter(!(CovSet %in% c("Direct Causes", "Swapped Direct Causes")) & ATE == "RD") %>%
  mutate(abs_Bias = set_label(abs_Bias, label = "Absolute Bias"),
         Method = fct_relevel(Method,"AIPW",after = 2),
         MSE = set_label(MSE, label = "Mean Squared Error"),
         CovSet = set_label(CovSet, label = "Adjustment Scenario")
         ) %>%
  select(CovSel:Method, DGM_type, abs_Bias,MSE) %>%
  rename(Estimator = Method)


my.cols <- brewer.pal(8, "Paired")[c(7,8,3,4,6)]
violin_plot_wrap <- function(df, y){
  g <- ggplot(data=df, aes_string(y=y, x='CovSet', fill='Estimator')) +
    geom_violin(aes(alpha = 0.9), position = position_dodge(width = 0.8), show.legend = F) +
    geom_boxplot( width=0.08, outlier.colour=NA, position = position_dodge(width = 0.8)) +
    facet_grid(DGM_type ~ CovSel, scales = "free", space='free_x') +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    theme_bw() +
    theme(legend.position = 'bottom',
          # strip.text.x = element_text(size = 8),
          # axis.text.x = element_text(angle = 35, hjust = 1),
          text = element_text(size=12)
    ) +
    labs(x=attr(df %>% select(CovSet) %>% pull(),"label"),
         y=attr(df %>% select(all_of(y)) %>% pull(),"label")) +
    scale_fill_manual(values=my.cols) 
  return(g)
}

pdf(here("Analysis","Medium_violin_plot_abs_bias_mse_all.pdf"), width = 10, height = 8)
print(violin_plot_wrap(df = rct_sum_plot_df, y= "abs_Bias"))
print(violin_plot_wrap(df = rct_sum_plot_df %>%filter(!(CovSet %in% "Empty Set / Unadjusted")), y= "abs_Bias"))
print(violin_plot_wrap(df = rct_sum_plot_df, y= "MSE"))
print(violin_plot_wrap(df = rct_sum_plot_df %>%filter(!(CovSet %in% "Empty Set / Unadjusted")), y= "MSE"))
dev.off()

########output

pdf(here("Analysis","Medium_violin_plot_sum_cd_tuned.pdf"), width = 22, height = 10)
#measures
perf_measure <- c("abs_Bias","Bias","MSE")
for (i in perf_measure){
  print(i)
  print(violin_plot_wrap(df = rct_sum, y= i))
  print(violin_plot_wrap(df = rct_sum%>% filter(Method != 'Empty'), y= i))
}
dev.off()





######-------quantile regressions--------######
# ATE_scale = c("RD", "logRR", "logOR")
ATE_scale = c("RD")
# form1 <- as.formula(paste0(outcome, "* 1000 ~ CovSet + Method  + DGM_type" ))

extract_rq <- function(rq_obj, num_tau) {
  rq_summ <- summary(rq_obj)
  if (num_tau==1){
    rq_df <- rq_summ$coefficients %>% 
      data.frame() %>% 
      mutate(v = rownames(.), 
             across(-v, round, 2),
             est_se = paste0(Value, " (", Std..Error, ") ")) %>%
      rename(p_val = Pr...t..) %>%
      select(v, est_se, p_val)
  } else{
    rq_df <- bind_cols(lapply(rq_summ, function(x) x$coefficients %>% 
                       data.frame() %>% 
                       mutate(v = rownames(.), 
                              across(-v, round, 2),
                              est_se = paste0(Value, " (", Std..Error, ") ")) %>%
                         rename(p_val = Pr...t..) %>%
                       select(v, est_se, p_val)
                       )) 
  }
  
  return(rq_df)
}


tidy_rq <- function(df = rct_sum_all %>% filter(CovSel==1), form = "abs_Bias*1000 ~  CovSet + Method  + DGM_type", tau_val = c(0.5, 0.75)){
  rq_obj_list <- lapply(ATE_scale, 
                        function(x) extract_rq(
                          rq(data = df %>% filter(ATE == x)), form, tau = tau_val), 
                          num_tau = length(tau_val) %>%
                          data.frame() %>%
                          rename_with( ~paste0(x, "_",.x))
                        )
  #create output table
  rq_out <- bind_cols(rq_obj_list) %>% 
    mutate(coeff = RD_v...1) %>% 
    select(-matches('_v[[:punct:]]')) %>% 
    relocate(coeff,.before=RD_est_se...2)
  return(rq_out)
}

###tbl_rq
# tbl_rq <- function(df = rct_sum_all %>% filter(CovSel==1), tau_val = c(0.5, 0.75)){
#   #abs bias
#   rq_bias_list <- lapply(tau_val, function(tau){
#                           rq(data = df %>% filter(ATE == "RD" & !(CovSet %in% c('Unadjusted', 'Empty'))), 
#                              "abs_Bias*1000 ~  CovSet + Method  + DGM_type", tau = tau) %>%
#                             tbl_regression(intercept = T) %>%
#                             bold_p() 
#                         }) %>% 
#     tbl_merge(tab_spanner = c("**50%**", "**75%**"))
#   #mse
#   rq_mse_list <- lapply(tau_val, function(tau){
#     rq(data = df %>% filter(ATE == "RD" & !(CovSet %in% c('Unadjusted', 'Empty'))), 
#        "MSE*1000 ~  CovSet + Method  + DGM_type", tau = tau) %>%
#       tbl_regression(intercept = T) %>%
#       bold_p() 
#   })%>% 
#     tbl_merge(tab_spanner = c("**50%**", "**75%**"))
#   #merge outpt 
#   rq_out <- tbl_merge(list(rq_bias_list, rq_mse_list), tab_spanner = c("**Abs. Bias**", "**MSE**"))
#   return(rq_out)
# }
# rq_performance_rd <- lapply(1:3, function(i) tbl_rq(df = rct_sum_all %>% filter(CovSel==i))) %>% 
#   tbl_stack(group_header = c("Manual", "MMHC Default", "MMHC Tuned")) 
# 
# rq_performance_rd %>% as_gt() 

#absolute bias
rct_sum_rq_df <- rct_sum_all %>%
  filter(!(CovSet %in% c("Direct Causes", "Swapped Direct Causes"#, "Empty Set / Unadjusted"
                         )) &
           # Method != "Empty" & 
           ATE == "RD") %>%
  # mutate(across(c(Method,DGM_type,CovSet), as.character),
  #        Method = recode(Method, Empty = "Unadjusted")) %>%
  select(CovSel:Method, DGM_type, abs_Bias,MSE) 

median_table_detail <- rct_sum_rq_df %>% 
  group_by(CovSel, Method, DGM_type, CovSet) %>%
  summarise(across(c(abs_Bias,MSE), list(median = median, Q3=~quantile(.x,probs = c(0.75)), max = max ))) %>%
  mutate(across(abs_Bias_median:MSE_max, label_number(scale = 1000, accuracy = 0.01)))

rq_bias_all <- bind_rows(tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==1), form = "abs_Bias*1000 ~   Method  + DGM_type * CovSet"),
                         tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==2), form = "abs_Bias*1000 ~   Method  + DGM_type"),
                         tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==3), form = "abs_Bias*1000 ~   Method  + DGM_type"), .id='method')

rq_bias_all_interaction <- bind_rows(tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==1), form = "abs_Bias*1000 ~  CovSet + Method  + DGM_type * CovSet"),
                                     tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==2), form = "abs_Bias*1000 ~  CovSet + Method  + DGM_type * CovSet"),
                                     tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==3), form = "abs_Bias*1000 ~  CovSet + Method  + DGM_type * CovSet"), 
                                     .id='method')

#MSE
rq_mse_all <- bind_rows(tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==1), form = "MSE*1000 ~  CovSet + Method  + DGM_type"),
                        tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==2), form = "MSE*1000 ~  CovSet + Method  + DGM_type"),
                        tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==3), form = "MSE*1000 ~  CovSet + Method  + DGM_type"),
                        .id='method')
rq_mse_all_interaction <- bind_rows(tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==1), form = "MSE*1000 ~  CovSet + Method  + DGM_type * CovSet"),
                                    tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==2), form = "MSE*1000 ~  CovSet + Method  + DGM_type * CovSet"),
                                    tidy_rq(df = rct_sum_rq_df %>% filter(CovSel==3), form = "MSE*1000 ~  CovSet + Method  + DGM_type * CovSet"), 
                                    .id='method')

#performance in RD
rq_performance_rd <- rq_bias_all %>%
  left_join(rq_mse_all ,by = c('method','coeff')) 

rq_performance_rd %>% 
  select(coeff, starts_with("RD_")) %>%
  kable('latex', booktabs = T,
        caption = "Quantile regressions on absolute bias and MSE (1e-4) of average treatment effect from different covariate adjustment sets selected by different methods",
        col.names = c("Coefficients", rep(c("Est. (SE)", 'p'),4))
  ) %>%
  add_header_above(c(" " = 1, "50%"= 2,  "75%"= 2, "50%"= 2,  "75%"= 2)) %>%
  add_header_above(c(" " = 1, "Abs. Bias ($10^{-4}$)"=4, "MSE ($10^{-4}$)"=4)) %>%
  pack_rows("Manual Covariate Selection", 1, 11) %>%
  pack_rows("Causal Discovery Algorithm: MMHC default", 12, 20) %>%
  pack_rows("Causal Discovery Algorithm: MMHC tuned", 21, 29)








#####-------------Accuracy of MMHC-----------#####
load( file = here("Analysis","summarised_output","accuracy_MMHC_default.RData"))
load( file = here("Analysis","summarised_output","accuracy_MMHC_tuned.RData"))
mean(mmhc_df_xy_sum$block_backdoor)
mean(mmhc_df_xy_sum_tuned$block_backdoor)

mmhc_all_sum <- mmhc_df_xy_sum %>%
  bind_rows(mmhc_df_xy_sum_tuned, .id='tuned') %>%
  mutate(tuned = factor(tuned, levels = 1:2, c("Default","Tuned")),
         C1only = ifelse(C1==1 & hsCRP_level==0 & loss_num ==0 & eligibility==0 & smoke ==0 & employed == 0,1,0),
         C2only = ifelse(C1==0 & hsCRP_level==0 & loss_num ==1 & eligibility==0 & smoke ==0 & employed == 0,1,0),
         C3only = ifelse(C1==0 & hsCRP_level==1 & loss_num ==0 & eligibility==0 & smoke ==0 & employed == 0,1,0),
         C1C2 = ifelse(C1==1 & hsCRP_level==0 & loss_num ==1 & eligibility==0 & smoke ==0 & employed == 0,1,0),
         C1C3 = ifelse(C1==1 & hsCRP_level==1 & loss_num ==0 & eligibility==0 & smoke ==0 & employed == 0,1,0),
         C2C3 = ifelse(C1==0 & hsCRP_level==1 & loss_num ==1 & eligibility==0 & smoke ==0 & employed == 0,1,0),
         C1C2C3 = ifelse(C1==1 & hsCRP_level==1 & loss_num ==1 & eligibility==0 & smoke ==0 & employed == 0,1,0),
         Empty = ifelse(C1==0 & hsCRP_level==0 & loss_num ==0 & eligibility==0 & smoke ==0 & employed == 0,1,0)) %>%
  ungroup() 

mmhc_var_selected_all <- mmhc_all_sum %>% 
  group_by(tuned) %>% 
  summarise(across(c(C1:employed), list(Percent = ~mean(.x), se = ~sd(.x)/sqrt(n())), .names = "{col}.{fn}")) %>%
  ungroup() %>%
  mutate(DGM_type = "Overall")

mmhc_var_selected <- mmhc_all_sum %>% 
  group_by(tuned, DGM_type) %>% 
  summarise(across(c(C1:employed), list(Percent = ~mean(.x), se = ~sd(.x)/sqrt(n())), .names = "{col}.{fn}")) %>%
  ungroup() %>%
  bind_rows(mmhc_var_selected_all) %>%
  mutate(across(-tuned:-DGM_type, ~round(.x*100,2))) %>%
  pivot_longer(cols = -tuned:-DGM_type) %>%
  separate(name, into = c("Covariate","stat"), sep = "\\.") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(Covariate = factor(Covariate,levels = c("C1", "loss_num", "hsCRP_level", "eligibility", "smoke", "employed","Empty"), 
                      labels = c("C1", "C2" ,"C3", "B (Eligibility)", "B (Smoke)", "B (Employed)","Empty")) %>%
           forcats::fct_rev(),
         DGM_type = gsub("-tri","-triangle",DGM_type),
         DGM_type = gsub("-bias","",DGM_type),
         DGM_type = factor(DGM_type, levels = c("Overall","Butterfly","M","Left-triangle","Right-triangle"),
                           labels = c("Overall (n=432)","Butterfly DAG (n=192)","M-structure DAG (n=48)","Left-triangle DAG (n=96)","Right-triangle DAG (n=96)")),
         tuned = fct_rev(tuned)) %>%
  rename(MMHC = tuned,
         DAG_type = DGM_type)


pdf(here("Analysis","prop_cov_selected.pdf"), width = 10,height = 10)
ggplot(data= mmhc_var_selected %>% filter(MMHC=='Tuned'), aes(x = Covariate , y = Percent)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_text(stat='identity', aes(label= paste0(round(Percent,1), "(", round(se,1), ")")),
            position = position_dodge(width = 1.1), size=3 , hjust = 0, color = "black") +
  facet_wrap(~DAG_type,ncol = 1, strip.position = "right") +
  scale_fill_grey(start = 0.4, end = .7)+
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0,100)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'bottom',
        text = element_text(size=12)#,
        # axis.text.x = element_text(size=10,angle = 25, hjust = 1)
        ) +
  labs(y = "Probability of being selected")

# by MMHC tuned
ggplot(data= mmhc_var_selected, aes(x = Covariate , y = Percent, group = MMHC, fill = MMHC)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_text(stat='identity', aes(label= paste0(round(Percent,1), "(", round(se,1), ")")),
            position = position_dodge(width = 1.1), size=3 , hjust = 0, color = "black") +
  facet_wrap(~DAG_type,ncol = 1, strip.position = "right") +
  scale_fill_grey(start = 0.4, end = .7)+
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0,100)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'bottom',
        text = element_text(size=12)#,
        # axis.text.x = element_text(size=10,angle = 25, hjust = 1)
  ) +
  labs(y = "Probability of being selected")
dev.off()

#------overall------
mmhc_overall <- mmhc_all_sum %>% 
  group_by(tuned) %>%
  summarise(n = n()/2000, 
            across(c(block_backdoor:block_backdoor.Y, unadj:C1C2C3), 
                   list(mean = ~mean(.x), se = ~sd(.x)/sqrt(n())), .names = "{col}.{fn}"),
            DGM_type='Overall') %>%
  relocate(DGM_type, .after = n)
#by DGM 
mmhc_dgm <- mmhc_all_sum %>% 
  group_by(tuned,DGM_type) %>%
  summarise(n = n()/2000, 
            across(c(block_backdoor:block_backdoor.Y, unadj:C1C2C3), 
                   list(mean = ~mean(.x), se = ~sd(.x)/sqrt(n())), .names = "{col}.{fn}"))

mmhc_acc_tab <- mmhc_overall %>%
  bind_rows(mmhc_dgm) %>%
  select(tuned:block_backdoor.se, unadj.mean:C1C2C3.se)%>% 
  relocate(c(unadj.mean:unadj.se), .after = C1C2C3.se) %>%
  pivot_longer(cols = -tuned:-DGM_type) %>%
  pivot_wider(names_from = tuned, values_from = value) %>%
  separate(name, into = c("adjset","stat"), sep = '\\.') %>%
  mutate(adjset = factor(adjset, levels = unique(adjset),
                         labels = c("Accuracy",
                                    "C1 Only",'C2 Only', 'C3 Only',
                                    '{C1,C2}','{C1,C3}','{C2,C3}',
                                    '{C1,C2,C3}','Unadjusted')))

mmhc_acc_tab_wide <- mmhc_acc_tab %>%
  pivot_wider(names_from = stat, values_from = Default:Tuned) %>%
  mutate(across(Default_mean:Tuned_se, ~label_percent(accuracy = 0.01)(.x)),
         Default = paste0(Default_mean, " (", Default_se,")"),
         Tuned = paste0(Tuned_mean, " (", Tuned_se,")")) %>%
  select(-Default_mean:-Tuned_se) 

mmhc_acc_tab_wide %>% 
  filter(DGM_type %in% c("Overall","Butterfly","M-bias")) %>%
    kable("latex",align = "c", booktabs = T, digits = 1,
          caption = "Accuracy of MMHC algorithm in selecting correct covariate adjustment set",
          col.names = c("Num. DAGs", "DAG type", "Accuracy / Adjustment Set", rep(c("Prop (SE)"),2))) %>%
  kable_styling() %>%
  add_header_above(c(" " = 3, "Default" = 1, "Tuned" = 1))%>%
  collapse_rows(columns = 1:2) #%>%
  # pack_rows("Selected Adjustment Set",2,9) %>%
  # pack_rows("Selected Adjustment Set",11,18) %>%
  # pack_rows("Selected Adjustment Set",20,27)
  
mmhc_acc_tab_wide %>% 
  filter(!(DGM_type %in% c("Overall","Butterfly","M-bias"))) %>%
  kable("latex",align = "c", booktabs = T, digits = 1,
        caption = "Accuracy of MMHC algorithm in selecting correct covariate adjustment set",
        col.names = c("Num. DAGs", "DAG type", "Adjustment Set", rep(c("Prop (SE)"),2))) %>%
  kable_styling() %>%
  add_header_above(c(" " = 3, "Default" = 1, "Tuned" = 1))%>%
  collapse_rows(columns = 1:2)

#plot
mmhc_acc_plot <- mmhc_acc_tab %>%
  pivot_longer(Default:Tuned, names_to = "MMHC") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(adjset = fct_rev(adjset),
         DGM_type = gsub("-tri","-triangle",DGM_type),
         DGM_type = gsub("-bias","",DGM_type),
         DGM_type = factor(DGM_type, levels = c("Overall","Butterfly","M","Left-triangle","Right-triangle"),
                           labels = c("Overall (n=432)","Butterfly DAG (n=192)","M-structure DAG (n=48)","Left-triangle DAG (n=96)","Right-triangle DAG (n=96)")),
         MMHC = fct_rev(MMHC)
         ) %>%
  rename(DAG_type = DGM_type,
         Percent = mean) 

pdf(here("Analysis","AccMMHC_fig.pdf"), width = 10,height = 12)
covset_fonts <- c("bold.italic", rep("plain",3),rep("bold",2),"plain","bold","plain") %>% rev()
ggplot(data= mmhc_acc_plot %>% filter(MMHC=='Tuned'), aes(x = adjset , y = Percent)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_text(stat='identity', aes(label= paste0(label_number(accuracy = 0.01,scale = 100)(Percent), 
                                               "(", label_number(accuracy = 0.01,scale = 100)(se), ")")),
            position = position_dodge(width = 1.1), size=3 , hjust = 0, color = "black") +
  facet_wrap(~DAG_type, labeller = 'label_value',ncol = 1, strip.position = "right") +
  scale_fill_grey(start = 0.4, end = .7)+
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), 
                     limits = c(0,1.04), n.breaks = 9) +
  theme_bw() +
  theme(legend.position = 'bottom',
        text = element_text(size=12),
        axis.text.y = element_text(face=covset_fonts)
  ) +
  labs(x = "Accuracy / Confounder Adjustment Set",
       y = "Probability")

## by MMHC Tuned
ggplot(data= mmhc_acc_plot, aes(x = adjset , y = Percent, group = MMHC, fill = MMHC)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_text(stat='identity', aes(label= paste0(label_number(accuracy = 0.01,scale = 100)(Percent), 
                                               "(", label_number(accuracy = 0.01,scale = 100)(se), ")")),
            position = position_dodge(width = 1.1), size=3 , hjust = 0, color = "black") +
  facet_wrap(~DAG_type, labeller = 'label_value',ncol = 1, strip.position = "right") +
  scale_fill_grey(start = 0.4, end = .7)+
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), 
                     limits = c(0,1.04), n.breaks = 9) +
  theme_bw() +
  theme(legend.position = 'bottom',
        text = element_text(size=12),
        axis.text.y = element_text(face=covset_fonts)
        ) +
  labs(x = "Accuracy / Confounder Adjustment Set",
       y = "Probability")
dev.off()
