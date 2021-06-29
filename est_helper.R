#' Title get Odds ratio
#'
#' @param x a 2*2 table 
#'
#' @return Odds ratio
get_OR <- function(x) {
  mat <- as.matrix(x)
  OR <- (mat[2,2]*mat[1,1]) / (mat[2,1]*mat[1,2])
  return(OR)
}

get_RD <- function(x) {
  mat <- as.matrix(x)
  RD <- (mat[2,2]/sum(mat[2,]) )-( mat[1,2]/sum(mat[1,]))
  return(RD)
}

get_RR <- function(x) {
  mat <- as.matrix(x)
  RR <- (mat[2,2]/sum(mat[2,])) / (mat[1,2]/sum(mat[1,]))
  return(RR)
}

get_all_est <- function(x){
  return(c(get_RD(x), log(get_RR(x)), log(get_OR(x))))
}

################Estimation function wrappers################
#' Title Estimating OR_bn from a mutilated graph with bnlearn function 
#'
#' @param dag a pre specify DAG using bnlearn
#' @param data dataframe
#' @param node_list a list denotes Y (outcome), A(exposure) and W(covariates)
#'
#' @return \psi_bn
bn_mut_jt_fun <- function(dag, data, node_list, scale = c("RD","OR","RR")){
  require(bnlearn)
  data = dplyr::mutate_all(data, as.factor)
  #######fit an bn with mle (pre-intervention (P(Y,C,X)))#######
  sim_bn_mle.fit <- bnlearn::bn.fit(x=dag, data=data, method="mle")
  
  #######do operator (or graph mutilation)#######
  # do operator to set treatment to x
  operation_list <- list(0)
  names(operation_list) <- node_list$A
  mle_mut <- bnlearn::mutilated(bn.net(sim_bn_mle.fit),operation_list)
  # refit bn with the mutilated graph (post-intervention:P(Y,C, do(X)))
  mle_mut <- bnlearn::bn.fit(mle_mut,data)
  
  ##########junction tree#########
  # exact inference with junction tree algo with gRain on the mutilated graph
  sim_mle.grain = bnlearn::as.grain(mle_mut)
  # joint prob query (marginal prob of X and Y: P(Y|X) = \sum_C P(Y|X,C)P(C))
  joint_prob_mle <- gRain::querygrain(sim_mle.grain, nodes= c(node_list$A,node_list$Y), type="joint")
  
  # calculate marginal RD / OR 
  if (scale=="RD"){
    bn_mle_rd = get_RD(joint_prob_mle)
    return(bn_mle_rd)
  }else if(scale == "OR"){
    bn_mle_OR <-get_OR(joint_prob_mle)
    return(log(bn_mle_OR))
  } else{
    bn_mle_rr = get_RR(joint_prob_mle)
    return(log(bn_mle_rr))
  }

}


#' Title Estimating \hat{\psi}
#'
#' @param data dataframe
#' @param node_list a list denotes Y (outcome), A(exposure) and W(covariates)
#' @param RD logical, return RD
#'
#' @return \hat{\psi}
#' parametric g computation
g_formula.logistic <- function(data,node_list, scale = c("RD","OR","RR","All")){
  formula <- as.formula(paste0(node_list$Y,"~",paste(c(node_list$A,node_list$W.Q), collapse= "+")))
  X <- node_list$A
  #g(X,C)
  lm <- glm(data=data, formula, family = "binomial")
  
  lm_tx <- mean(predict(lm, newdata = transform(data,X = 1), type="response"))
  lm_untx<- mean(predict(lm, newdata = transform(data,X = 0), type="response"))
  
  rd <- lm_tx - lm_untx
  logrr <- log(lm_tx / lm_untx)
  logor <- log((lm_tx/(1-lm_tx)) / (lm_untx/(1-lm_untx)))
  
  if (scale == "RD"){
    return(rd)
  } else if(scale == "RR") {
    return(logrr)
  }  else if(scale == "OR") {
    return(logor)
  } else {
    return(c(rd, logrr, logor))
  }
  
}
#' IP-Weighting with stablized weights
ipw_fun <- function(data, node_list, scale = c("RD","OR","RR","All")){
  exposure_formula <- as.formula(paste0(node_list$A,"~",paste(node_list$W.g, collapse= "+")))
  outcome_formula <- as.formula(paste0(node_list$Y,"~",paste(node_list$A, collapse= "+")))
  
  p_score <- glm(exposure_formula, data=data,family=binomial("logit"))$fitted.values
  X <- get(node_list$A, data)
  data$ate_sw <- as.numeric( (mean(X))*(X/p_score) + (mean(1-X))*((1-X)/(1-p_score)) )
  if (scale == "RD"){
    ipw <- glm(data=data, outcome_formula, weights=ate_sw,family=gaussian(link = "identity"))
    ipw_est <- as.numeric(ipw$coefficients[2])
  } else if (scale == "RR"){
    ipw <- glm(data=data, outcome_formula, weights=ate_sw,family=poisson(link = "log"))
    ipw_est <- as.numeric(ipw$coefficients[2])
  }else if (scale == "OR"){
    ipw <- glm(data=data, outcome_formula, weights=ate_sw,family=binomial(link = "logit"))
    ipw_est <- as.numeric(ipw$coefficients[2])
  }  else{
    ipw1 <- glm(data=data, outcome_formula, weights=ate_sw,family=gaussian(link = "identity"))
    ipw2 <- glm(data=data, outcome_formula, weights=ate_sw,family=poisson(link = "log"))
    ipw3 <- glm(data=data, outcome_formula, weights=ate_sw,family=binomial(link = "logit"))
    ipw_est <- c(as.numeric(ipw1$coefficients[2]),as.numeric(ipw2$coefficients[2]),as.numeric(ipw3$coefficients[2]))
  }
  return(ipw_est)
}

#' AIPW and TMLE using tmle3 function (sample splitting k = 10 and bound = 0.025)
tmle_aipw_wrap <- function(data,node_list, scale = c("RD","OR","RR","All")){
  # require(tmle)
  #https://yqzhong7.github.io/AIPW/
  #remotes::install_github("yqzhong7/AIPW")
  # require(SuperLearner)
  # require(AIPW)
  data <- data %>% 
    select(node_list$Y, node_list$A,node_list$W.g,node_list$W.Q) %>%
    rename(Y = node_list$Y,
           A = node_list$A)
  
  exposure_formula <- paste0("A ~ ",paste(node_list$W.g, collapse= "+"))
  outcome_formula <- paste0("Y ~ A+",paste(node_list$W.Q, collapse= "+"))

  W = data %>% select(node_list$W.g,node_list$W.Q) %>% as.data.frame()
  
  #manual input
  Q.form = glm(data = data , outcome_formula, "binomial")
  Q1 <- predict(Q.form, newdata = transform(data, A=1), type = "response")
  Q0 <- predict(Q.form, newdata = transform(data, A=0), type = "response")
  Q <- cbind(Q0,Q1)
  g1W <- glm(data = data, exposure_formula, "binomial")$fitted.values
  
  tmle_fit <- tmle(Y=data$Y,
                   A=data$A,
                   W=W,
                   Q = Q,
                   g1W = g1W,
                   # Q.SL.library = "SL.glm",
                   # g.SL.library = "SL.glm",
                   family = "binomial",
                   cvQinit = F,
                   gbound = 0.025,
                   V = 10,
                   automate = FALSE)
  
  
  #aipw 
  suppressMessages({
    # aipw <- AIPW::AIPW$new(A=A, 
    #                       Y=Y, 
    #                       W=data %>% select(node_list$W) %>% as.matrix(),
    #                       Q.SL.library = "SL.glm",
    #                       g.SL.library = "SL.glm",
    #                       k_split = 1,
    #                       verbose=FALSE)$fit()$summary()
    aipw <- AIPW::AIPW_tmle$new(A=data$A, 
                           Y=data$Y, 
                           tmle_fit = tmle_fit,
                           verbose=FALSE)$summary()
  })
  
  #return log OR
  if (scale == "RD"){
    tmle_RD <- tmle_fit$estimates$ATE$psi
    aipw_RD <- aipw$estimates$RD[1]
    return(c(tmle_RD,aipw_RD))
  }else if (scale=="RR"){
    tmle_logRR <- tmle_fit$estimates$RR$log.psi
    aipw_logRR <- log(aipw$estimates$RR[1])
    return(c(tmle_logRR,aipw_logRR))
  }else if (scale=="OR"){
    tmle_logOR <- tmle_fit$estimates$OR$log.psi
    aipw_logOR <- log(aipw$estimates$OR[1])
    return(c(tmle_logOR,aipw_logOR))
  } else {
    return(rbind(c(tmle_fit$estimates$ATE$psi,tmle_fit$estimates$RR$log.psi,tmle_fit$estimates$OR$log.psi),
                 c(aipw$estimates$RD[1],log(aipw$estimates$RR[1]),log(aipw$estimates$OR[1]))))
  }
  
}

