####Setup####
packages <- c("purrr","future","future.apply","here","speedglm","tidyverse")

for (package in packages) {
  library(package, character.only=T)
}

set.seed(123)
eager_raw <- read_csv(here("Data","eager_base_imputed.csv")) 
eager_bl <- eager_raw %>%
  mutate(meanAP = (BPS + 2*BPD)/3,
         live_birth = as.numeric(outcome=="live birth"),
         time_try_pregnant=if_else(as.numeric(time_try_pregnant)>10,10,as.numeric(time_try_pregnant)),
         eligibility = as.numeric(eligibility=="new"),
         hsCRP_level = cut(hsCRP, breaks=c(-Inf, 2, 10, Inf), labels=c("low","medium","high"), right = F, include.lowest = T)
         ) %>%
  dplyr::select(live_birth,treatment,eligibility,loss_num,age,time_try_pregnant,hsCRP_level,BMI,meanAP,alcohol,smoke,employed)
covariates <- colnames(eager_bl[3:ncol(eager_bl)])
eager_bl
summary(eager_bl)
table(eager_bl$loss_num)
table(eager_bl$hsCRP_level)


#--------------------helper functions-----------------#
source(here("est_helper.R"))
expit <- function(x) exp(x)/(exp(x)+1)
inv_expit <- function(y) log(y/(1-y))
#' Title Generate probablity for binomial distributiokn
#'
#' @param N sample size
#' @param intercept intercept
#' @param varname data.frame, parent variables
#' @param OR vector, Correspoinding odds ratio of the variables to the child node
#' @param marginalProb predefined marginal probability 
genProb <- function(df, N, intercept=NULL,varname, OR, marginalProb = NULL) {
  vars = df %>% dplyr::select(all_of(varname)) %>% mutate_all(as.numeric) %>% as.matrix()
  coeff <- vars%*%log(OR)
  if (is.null(intercept) & !is.null(marginalProb)){
    intercept = inv_expit(marginalProb) - mean(coeff)
  }
  return(expit(rep(intercept,N) + coeff))
}


start_time <- Sys.time()
####-------------------RCT simulation-------------------####
sim_rct_func <- function(N=1228, marginalProb = 0.5, ORyx = 1.75, 
                         ORxc1 = 2, ORxc2 = 5, ORc1c2 = 10,
                         ORyc1 = 2, ORyc3 = 5, ORc1c3 = 5) {
  #resample from eager 
  eager_boot_id <- sample(1:nrow(eager_bl),size = N, replace = T)
  eager_boot <- eager_bl[eager_boot_id, 2:ncol(eager_bl)]
  
  # #C1
  pC1 <- genProb(df = eager_boot, N, varname = c("loss_num","hsCRP_level"), OR=c(ORc1c2,ORc1c3), marginalProb = 0.5, #intercept = -2 
                 )
  eager_boot$C1 <- rbinom(N,1,pC1)
  
  #X
  pX <- genProb(df = eager_boot, N, varname = c("treatment","C1","loss_num"),  OR=c(0.8, ORxc1, ORxc2), marginalProb = 0.5)
  eager_boot$X <- ifelse(eager_boot$treatment, rbinom(N,1,pX), 0)

  #Y
  pY <- genProb(df = eager_boot, N, varname = c("X","C1","hsCRP_level"), OR=c(ORyx, OR=c(ORyc1, ORyc3)), marginalProb = 0.5)
  eager_boot$Y <- rbinom(N,1,pY)
  
  #simulated dataset
  attr(eager_boot,"avg_prob") <- list(avg_pX=mean(pX),avg_pY=mean(pY),avg_C1=mean(pC1))
  return(eager_boot)
}


true<- g_formula.logistic(data= sim_rct_func(N=1e+6,
                                             ORxc1=1, ORyc1=1, ORxc2 = 2, ORc1c2 = 2,ORyc3 = 2, ORc1c3 = 2),
                   node_list = list(A = "X",Y = "Y", W.g = c("loss_num"),W.Q = c("hsCRP_level")),
                   scale = 'All')
est <- lapply(1:1000, function(x) g_formula.logistic(data= sim_rct_func(
                                                                        ORxc1=1, ORyc1=1, ORxc2 = 2, ORc1c2 = 2,ORyc3 = 2, ORc1c3 = 2),
                          node_list = list(A = "X",Y = "Y", W.g = c("C1"),W.Q = c("C1")),
                          scale = 'All'))
est_df <- do.call(rbind,est)
mean(abs(est_df[,1]-true[1]))
mean(abs(est_df[,2]-true[2]))
mean(abs(est_df[,3]-true[3]))

#--------------------True Est.-------------------#
#parameters for the cartesian product
param_cross <- purrr::cross(list(
  ORyx = c(1/1.75, 1, 1.75),
  ORxc1 = c(0.5, 1, 2), ORxc2 = c(0.5, 2), ORc1c2 = c(0.5, 2),
  ORyc1 = c(0.5, 1, 2), ORyc3 = c(0.5, 2), ORc1c3 = c(0.5,2)
))
#param that does not vary
param_true <- lapply(param_cross, append,
                list(N = 1e+6,  marginalProb=c(0.5)) 
)
#true estimate function
sim_rct_true_est_func <- function(param){
  #simulate population
  true_df <- do.call(sim_rct_func, param)
  #define node list
  node_list <- list(A = "X",
                    Y = "Y",
                    W.g = c("loss_num"),
                    W.Q = c("hsCRP_level"))
  if (param$ORyc1 != 1 & param$ORxc1 != 1) {
    node_list$W.Q <- append(node_list$W.Q, "C1")
    node_list$W.g <- append(node_list$W.g, "C1")
  } else if (param$ORyc1 == 1) {
    node_list$W.g <- append(node_list$W.g, "C1")
  } else if (param$ORxc1 == 1) {
    node_list$W.Q <- append(node_list$W.Q, "C1")
  } 
  #marginal estimation
  rct_true_est <- g_formula.logistic(data=true_df, node_list = node_list, 'All')
  output <- append(param, list(
    true_rd = rct_true_est[1], true_logrr = rct_true_est[2], true_logor = rct_true_est[3], 
    marginalProb_X = mean(true_df$X), marginalProb_Y = mean(true_df$Y), marginalProb_C1= mean(true_df$C1))
  )
  return(output)
}
#use future est
plan(multisession,workers=10,gc=T)
s_time_true <- Sys.time()
sim_true_est_list <- future.apply::future_lapply(
  param_true, future.seed = T,
  function(x) do.call(sim_rct_true_est_func,list(x))
  )
(e_time_true <- Sys.time() - s_time_true)

#timestamp
date_time <- format(Sys.time(),format = "%Y-%m-%dT%H%M%S%Z")
save(sim_true_est_list,file=here("Simulation",paste0("sim_rct_true_est",date_time,".RData")))



#---------------------Simulation-------------------#
set.seed(123)
#param that does not vary
param <- lapply(param_cross, append,
                     list(N = 1228,  marginalProb=c(0.5)) )
#number of simulations per param set
nSim = 2000
nSlice = 100
#timestamp
(date_time <- format(Sys.time(),format = "%Y-%m-%dT%H%M%S%Z"))
#multiprocess
availableCores()
plan(multisession,workers=10,gc=T)
#create a folder
dir.create(here("Simulation",date_time))

#simulation
s_time<- Sys.time()
for (i in 1:nSlice){
  simList <- future.apply::future_sapply(
    param, future.seed = T, 
    function(x) lapply(1:(nSim/nSlice), function(y) {
      df <- do.call(sim_rct_func,x)
      attr(df,"param") <- x
      return(df)
      }))
  save(simList, file = here("Simulation",date_time,paste0("simList_rct_",date_time,"_",i,".RData")))
  gc()
}
(e_time <-Sys.time() - s_time)


#####output sessioninfo#####
fileConn<-file(here("Simulation",date_time,paste0("sessionInfo_",date_time,".txt")))
writeLines(capture.output(sessionInfo()), fileConn)
close(fileConn)
