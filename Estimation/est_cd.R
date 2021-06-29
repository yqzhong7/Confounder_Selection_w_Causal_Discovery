packages <- c("tidyverse","future.apply",
              "AIPW","tmle","here","progressr","bnlearn"
)
#remotes::install_github("yqzhong7/AIPW")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}


#-------------Setup--------------------####
setwd(here())
source("./est_helper.R")
sessionInfo()


#---------------Helper functions-------------####
estimator_wrapper <- function(data, node_list, scale=scale){
  if (!is.na(node_list$W.Q) ){
    g <- g_formula.logistic(data,node_list,scale)
  } else{
    g <- rep(NA, 3)
  }
  
  if (!is.na(node_list$W.g)){
    ipw <- ipw_fun(data,node_list,scale)
  } else{
    ipw <- rep(NA, 3)
  }
  
  if (!is.na(node_list$W.Q) & !is.na(node_list$W.g)){
    tmle_aipw <- tmle_aipw_wrap(data,node_list, scale)
  } else{
    tmle_aipw <- matrix(rep(NA, 6), ncol = 3)
  }
  
  output <- data.frame(Method = c("gComp","IPW","TMLE","AIPW"),
                  rbind(g,ipw,tmle_aipw))
  
  colnames(output) <- c("Method","RD","logRR","logOR")
  invisible(output)
}
  

#---------------CD setup--------------####
covar = c( "eligibility","loss_num","age","time_try_pregnant", "hsCRP_level","BMI","meanAP","alcohol","smoke","employed","C1")
nCov = length(covar)
cd_constrain = list(
  wl = matrix(c('treatment','X','X','Y'), ncol = 2, byrow = T),
  bl = rbind(matrix(c(covar,rep("treatment",nCov*2),covar), ncol = 2, byrow = F),
             matrix(c(rep("X",nCov),covar), ncol = 2, byrow = F),
             matrix(c(rep("Y",nCov),covar), ncol = 2, byrow = F))
)

#---------------Estimation function--------------####
estimation_cd <- function(df_list, cd_constrain){
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
      
      #---mmhc---#
      #discretization
      df[,c("age", "time_try_pregnant", "BMI", "meanAP")] <- bnlearn::discretize(df[,c("age", "time_try_pregnant", "BMI", "meanAP")],
                                                                                 breaks = 5, method = 'interval')
      df <- df %>% mutate_all(as.factor) %>% as.data.frame()
      
      #cd
      bn <- mmhc(df, whitelist = cd_constrain$wl, blacklist = cd_constrain$bl)
      edge_list <- bn$arcs %>% 
        as.data.frame() %>% 
        filter(to %in% c("X","Y") & !(from %in% c('treatment','X','Y'))) %>%
        mutate_all(as.character)
      
      #--estimation--#
      #remove discretization
      df <- df_list[[i]]
      
      #define node list
      node_list <- list(A = "X",
                        Y = "Y")
     
      #parents
      parent_node_list <- node_list
      if (nrow(edge_list)==0){
        empty_mat <- data.frame(matrix(c("Empty","Empty",NA,NA,NA), nrow=1))
        empty_mat[,3:5] <- apply(empty_mat[,3:5], 1, as.numeric)
        attr(empty_mat,'edge_list') <- bn$arcs
        names(empty_mat) <- c("CovSet","Method","RD","logRR","logOR")
        return(empty_mat)
      } else{
        parent_node_list$W.Q <- ifelse(nrow(edge_list[edge_list$to=='Y',])>0,edge_list[edge_list$to=='Y',]$from, NA)
        parent_node_list$W.g <- ifelse(nrow(edge_list[edge_list$to=='X',])>0,edge_list[edge_list$to=='X',]$from, NA)
      }
      
      parent_est <- estimator_wrapper(data=df, node_list = parent_node_list, 'All')
      parent_est <- cbind(CovSet='Parent',parent_est)
      
      #swapped parents
      s_paraent_node_list <- node_list
      s_paraent_node_list$W.Q <- parent_node_list$W.g
      s_paraent_node_list$W.g <- parent_node_list$W.Q

      s_parent_est <- estimator_wrapper(data=df, node_list = s_paraent_node_list, 'All')
      s_parent_est <- cbind(CovSet='Swapped Parent',s_parent_est)
      
      #ancestors
      ancestor_node_list <- node_list
      ancestor_node_list$W.Q <- ancestor_node_list$W.g <- unique(edge_list$from)
      
      ancestor_est <- estimator_wrapper(data=df, node_list = ancestor_node_list, 'All')
      ancestor_est <- cbind(CovSet='Ancestor',ancestor_est)
      
      #combine
      output <- rbind(parent_est, s_parent_est, ancestor_est)
      rownames(output) <- c()
      attr(output, 'edge_list') <- bn$arcs
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
sim_file_names <- grep(".RData",list.files("./Simulation/2020-11-17T180306EST"), value = T)
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
file_path_dt <- paste0("./Estimation/CD/",date_time)
dir.create(file_path_dt)

#simulation
s_time<- Sys.time()
for (i in sim_file_names[2]){
  (si_time <- Sys.time())
  #load data
  load(paste0("Simulation/2020-11-17T180306EST/",i))
  trunck_num <- gsub("\\.","", str_extract(i, "_[[:digit:]]{1,3}\\."))
  print(trunck_num)
  
  try({
    progressr::with_progress(estList <- estimation_cd(simList, cd_constrain = cd_constrain))  
    save(estList, file = paste0("Estimation/CD/",date_time,paste0("/est_cd",trunck_num,".RData")))
  }, silent = T)
  (ei_time <- Sys.time() - si_time)
  print(ei_time)
  gc()
}
(e_time <-Sys.time() - s_time)

#####output sessioninfo#####
fileConn<-file(paste0("Estimation/CD/",date_time,paste0("/sessionInfo_",date_time,".txt")))
writeLines(capture.output(sessionInfo()), fileConn)
close(fileConn)

