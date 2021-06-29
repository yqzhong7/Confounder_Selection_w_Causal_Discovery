packages <- c("tidyverse","future.apply",
              "AIPW","tmle","here","progressr"
)
#remotes::install_github("yqzhong7/AIPW")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}


#-------------Setup--------------------####
source(here("est_helper.R"))


#---------------Helper functions-------------####
estimator_wrapper <- function(data, node_list, scale=scale){
  g <- g_formula.logistic(data,node_list,scale)
  ipw <- ipw_fun(data,node_list,scale)
  tmle_aipw <- tmle_aipw_wrap(data,node_list, scale)
  
  output <- data.frame(Method = c("gComp","IPW","TMLE","AIPW"),
                  rbind(g,ipw,tmle_aipw))
  
  colnames(output) <- c("Method","RD","logRR","logOR")
  invisible(output)
}



#---------------Estimation function--------------####
estimation <- function(df_list){
  iter <- 1:length(df_list)
  pb <- progressr::progressor(along=iter)
  estList <- future.apply::future_lapply(
    iter, 
    future.seed = T,
    future.packages = packages,
    future.label = T, 
    function(i,...) {
      df = df_list[[i]]
      
      param = attributes(df)$param
      #define node list
      node_list <- list(A = "X",
                        Y = "Y",
                        W.g = c("loss_num"),
                        W.Q = c("hsCRP_level"))
      #parents
      parent_node_list <- node_list
      if (param$ORyc1 != 1 & param$ORxc1 != 1) {
        parent_node_list$W.Q <- append(node_list$W.Q, "C1")
        parent_node_list$W.g <- append(node_list$W.g, "C1")
      } else if (param$ORyc1 == 1) {
        parent_node_list$W.g <- append(node_list$W.g, "C1")
      } else if (param$ORxc1 == 1) {
        parent_node_list$W.Q <- append(node_list$W.Q, "C1")
      } 
      
      parent_est <- estimator_wrapper(data=df, node_list = parent_node_list, 'All')
      parent_est <- cbind(CovSet='Parent',parent_est)
      
      #swapped parents
      s_paraent_node_list <- node_list
      s_paraent_node_list$W.Q <- parent_node_list$W.g
      s_paraent_node_list$W.g <- parent_node_list$W.Q
      
      s_parent_est <- estimator_wrapper(data=df, node_list = s_paraent_node_list, 'All')
      s_parent_est <- cbind(CovSet='Swapped Parent',s_parent_est)
      
      #ancestors (All causes with collider (not in M))
      ancestor_node_list <- node_list
      ancestor_node_list$W.Q <- ancestor_node_list$W.g <- unique(c(parent_node_list$W.Q,parent_node_list$W.g))
      
      ancestor_est <- estimator_wrapper(data=df, node_list = ancestor_node_list, 'All')
      ancestor_est <- cbind(CovSet='Ancestor',ancestor_est)
      
      #all
      all_node_list <- node_list
      all_node_list$W.Q <- all_node_list$W.g <- c("eligibility", "loss_num", "age","time_try_pregnant" ,"hsCRP_level", "BMI",
                                                  "meanAP","alcohol",   "smoke", "employed", "C1" )
      
      all_est <- estimator_wrapper(data=df, node_list = all_node_list, 'All')
      all_est <- cbind(CovSet='All',all_est)
      
      #c1 only
      c1_node_list <- node_list
      c1_node_list$W.Q <- c1_node_list$W.g <- c("C1")
      
      c1_est <- estimator_wrapper(data=df, node_list = c1_node_list, 'All')
      c1_est <- cbind(CovSet='C1 Only',c1_est)
      
      #unadjusted
      unadj_est <- matrix(c("Unadjusted","Unadjusted", get_all_est(prop.table(table(df$X,df$Y), margin = 1))), nrow = 1, byrow = T)
      colnames(unadj_est) <- colnames(parent_est)
      
      #combine
      output <- rbind(parent_est, s_parent_est, ancestor_est, all_est, c1_est, unadj_est)
      rownames(output) <- c()
      # output[3:5] <- apply(output[3:5], 2, as.numeric)
      
      pb(message=sprintf("No.%g", i))
      return(output)
    })
}
# #function testing
# suppressWarnings({
# with_progress(out <- estimation(df = sliced_list[1:10]))
# })
# tibble(out)



#-------------------Estimation w/ future----------########
#setup
sim_file_names <- grep(".RData",list.files(here("Simulation","2020-11-16T193628EST_sliced")), value = T)
set.seed(123)
handlers("progress")
#timestamp
(date_time <- format(Sys.time(),format = "%Y-%m-%dT%H%M%S%Z"))

#multiprocess
availableCores()
plan(multisession,workers=10,gc=T)
# plan(sequential)

#increase memory
options(future.globals.maxSize= 1.2*1024^3)

#create a folder
dir.create(here("Estimation","Full_Knowledge",date_time))

#simulation
s_time<- Sys.time()
for (i in sim_file_names[1]){
  (si_time <- Sys.time())
  #load data
  load(here("Simulation","2020-11-16T193628EST_sliced",i))
  trunck_num <- gsub("\\.","", str_extract(i, "_[[:digit:]]{1,2}_[[:digit:]]\\."))
  print(trunck_num)

  progressr::with_progress(estList <- estimation(sliced_list))  

  save(estList, file = here("Estimation","Full_Knowledge",date_time,paste0("est_wo_cd",trunck_num,".RData")))
  (ei_time <- Sys.time() - si_time)
  print(ei_time)
  gc()
}
(e_time <-Sys.time() - s_time)

