#Setup
packages <- c("tidyverse","future.apply","fst",
              "AIPW","tmle","here","progressr","broom","RColorBrewer",
              "lme4",'lmerTest','quantreg', 'knitr','kableExtra'
)
#remotes::install_github("yqzhong7/AIPW")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

sim_file_names <- grep(".RData",list.files(here("Simulation","2021-03-27T125752EDT")), value = T)
load(here("Simulation","2021-03-27T125752EDT",sim_file_names[1]))

# --------load parameters---------####
param_list <- lapply(simList, function(x) attributes(x)$param)
param_list_df <- matrix(unlist(param_list), length(param_list), byrow=T) %>%
  data.frame() %>%
  rename_all(~names(param_list[[1]])) %>%
  mutate(simNo = 1:nrow(.))
param_list_df %>% unique()

param_df <- param_list_df %>% select(starts_with("OR")) %>% unique() %>% mutate(DGM_param=1:nrow(.))
save(param_df, file = here("Analysis",'param_df.RData'))


rm(param_list)
rm(simList)

# #-------load true values--------####
load(here("Simulation","sim_rct_true_est2021-03-27T125752EDT.RData"))
true_est_df <- bind_rows(sim_true_est_list)
true_est_df_long <- true_est_df %>%
  rename_at(vars(starts_with('true')), ~c('RD','logRR','logOR')) %>%
  pivot_longer(cols = RD:logOR, names_to = 'ATE', values_to = 'true_est')

# # #--------load output------------####
est_file_names <- grep(".RData",list.files(here("Estimation","Full_Knowledge"), full.names = T), value = T)
est_df_list <- list()
(s_time <- Sys.time())
for (i in 1:length(est_file_names)){
  print(i)
  si_time <- Sys.time()
  load(est_file_names[i])
  est_df_list[[i]] <- bind_rows(estList, .id= "simNo")
  print(Sys.time()-si_time)
  }

est_df <- bind_rows(est_df_list)

est_df_param <- est_df %>%
  mutate(simNo = as.numeric(simNo)) %>%
  left_join(param_list_df, by = "simNo") %>%
  mutate_at(vars(RD:logOR), as.numeric) %>%
  pivot_longer(cols = RD:logOR, names_to = 'ATE', values_to = "est") %>%
  left_join(true_est_df_long, by = c(colnames(param_list_df)[1:7],"ATE")) %>%
  mutate(m_bias = as.factor(as.numeric(ORyc1 == 1 & ORxc1 ==1)),
         DGM_type = case_when(ORxc1==1 & ORyc1== 1 ~ 'M-bias',
                              ORxc1==1 & ORyc1!= 1 ~ 'Right-tri',
                              ORxc1!=1 & ORyc1== 1 ~ 'Left-tri',
                              ORxc1!=1 & ORyc1!= 1 ~ 'Butterfly') %>% factor())
(e_time <- Sys.time() - s_time)

# combine and save into fst
write.fst(est_df_param, here("Analysis","est_df_param_wo_cd_all.fst"))


#load fst
est_df_param <- read.fst(here("Analysis","est_df_param_wo_cd_all.fst"))

#load param
load(here("Analysis","param_df.RData"))

  
#--------------Summary----------####
my.cols <- brewer.pal(8,"Paired")[5:8]
violin_plot_wrap <- function(df, y){
  df_mean <- df %>% group_by(CovSet,Method,ATE,DGM_type) %>%summarise(across(.cols = all_of(y), mean))
  g <- ggplot(data=df, aes_string(y=y, x='CovSet', fill='DGM_type')) +
      geom_violin(aes(alpha=0.9), position = position_dodge(width = 0.4)) +
      geom_boxplot(width=.1, outlier.colour=NA, position = position_dodge(width = 0.4)) +
      # geom_point(data = df_mean, aes(x = y), shape=20, size=5, color="red", fill="red") +
      facet_grid(ATE~Method, scales = 'free') +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme_bw() +
      theme(legend.position = 'bottom',
            # strip.text.x = element_text(size = 8),
            axis.text.x = element_text(angle = 35, hjust = 1),
            text = element_text(size=18)
            ) +
    scale_fill_manual(values=my.cols) +
    # scale_fill_brewer('OrRd') +
    ggtitle(y)
  return(g)
}


## summarizing over 2000 nSim
rct_sum <- est_df_param %>%
  select(simNo:Method, starts_with("OR"), DGM_type, ATE, est, true_est) %>%
  group_by_at(vars(CovSet, Method, ATE, starts_with("OR"),DGM_type)) %>%
  summarise(count = n(),
            mean_est = mean(est),
            true_est = mean(true_est),
            Bias = mean(est) - true_est,
            abs_Bias = abs(mean(est) - true_est),
            Bias_se = sd(est) / sqrt(count),
            MSE = mean((est - true_est)^2)) %>%
  left_join(param_df, by = colnames(param_df)[1:7])

sum_min_max <- rct_sum %>%
  group_by(ATE) %>% 
  summarise_at(vars(Bias:MSE), funs(min(.),max(.), IQR(.))) %>%
  select(ATE, sort(colnames(.)))
rct_sum %>%
  group_by(Method,ATE) %>% 
  summarise_at(vars(Bias:MSE), funs(min(.),max(.), IQR(.))) %>%
  select(Method, ATE,sort(colnames(.))) %>%
  arrange(ATE)


###from summarised output
pdf(here("Analysis","Medium_violin_plot_sum.pdf"), width = 22, height = 10)
#measures
perf_measure <- c("abs_Bias","MSE")
for (i in perf_measure){
  print(i)
  print(violin_plot_wrap(df = rct_sum, y= i))
  print(violin_plot_wrap(df = rct_sum%>% filter(Method != 'Unadjusted'), y= i))
}
dev.off()



###--------- find poor performance cases--------- ####
rct_sum_factor <- rct_sum %>%
  # filter(Method != 'Unadjusted') %>%
  mutate_at(vars(starts_with('OR')), function(x) case_when(x<1 ~ 'neg',
                                                            x==1 ~'null',
                                                            x>1 ~ 'pos')) %>%
  mutate_at(vars(CovSet:DGM_type), factor)

save(rct_sum_factor, file = here("Analysis","summarised_output","wo_CD.RData"))



