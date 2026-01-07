rm(list = ls())
library(MASS)
library(survival)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(DPCri)
library(mvtnorm)
library(Landmarking)
library(nlme)
library(JMbayes2)
library(refund)

#setwd("/Users/user/Desktop/Rcodes/PP_TM/New_stan/11th")#
setwd("C:\\Users\\p80744tb\\Desktop\\Desktop 2026\\intensive")
source("iauc_sd.R")
NN=100
Esti=CovR=matrix(0,NN,3)
Results=list()
for(kkk in 36:NN){ 
  set.seed(kkk)
  
  # ----------------------------
  # Simulation parameters
  # ----------------------------
  nsample <- 1000
  alpha   <- -0.5
  Beta    <- c(-0.5, -0.5, 0.5, 0.5)
  sigma   <- 1
  gamma_w <- c(0.2, -0.2)
  followup <- 2
  t <- seq(from = 0, to = followup, by = 0.005)
  
  # Covariates and random effects
  x1 <- x2 <- w1 <- w2 <- rep(NA_real_, nsample)
  death <- rep(0L, nsample)
  Y1 <- matrix(NA_real_, nrow = nsample, ncol = length(t))
  time <- rep(NA_real_, nsample)
  id   <- seq_len(nsample)
  ct   <- rep(NA_integer_, nsample)
  evtim <- td <- td_censor <- rep(NA_real_, nsample)
  D <- matrix(c(1,0.5,0.5,1),2,2)
  u <- rmvnorm(nsample, mean = c(0,0), sigma = D)
  
  for (i in 1:nsample) {
    x1[i] <- rnorm(1)
    x2[i] <- rbinom(1,1,0.5)
    w1[i] <- rnorm(1)
    w2[i] <- rbinom(1,1,0.5)
  }
  
  haz_fun_i <- function(s, x1i, x2i, w1i, w2i, ui) {
    eta <- gamma_w[1]*w1i + gamma_w[2]*w2i +
      alpha*(Beta[1] + Beta[2]*x1i + Beta[3]*x2i + Beta[4]*s + ui[1] + ui[2]*s)
    exp(eta)
  }
  
  cumhaz_i <- function(T, x1i, x2i, w1i, w2i, ui) {
    if (T <= 0) return(0)
    integrate(function(s) haz_fun_i(s,x1i,x2i,w1i,w2i,ui),
              lower = 0, upper = T, rel.tol = 1e-6, abs.tol = 1e-8)$value
  }
  
  surv_i <- function(T, x1i, x2i, w1i, w2i, ui) exp(-cumhaz_i(T,x1i,x2i,w1i,w2i,ui))
  
  sim_event_time_i <- function(x1i,x2i,w1i,w2i,ui,Tmax=followup) {
    U <- runif(1)
    S_Tmax <- surv_i(Tmax,x1i,x2i,w1i,w2i,ui)
    if (S_Tmax > U) return(Inf)
    f <- function(T) surv_i(T,x1i,x2i,w1i,w2i,ui)-U
    res <- tryCatch(uniroot(f, lower=.Machine$double.eps, upper=Tmax)$root,
                    error=function(e) Tmax)
    res
  }
  
  # ----------------------------
  # Generate longitudinal and survival data
  # ----------------------------
  for (i in 1:nsample) {
    evt <- sim_event_time_i(x1[i], x2[i], w1[i], w2[i], u[i,])
    td[i] <- rexp(1,1)
    evtim[i] <- evt
    
    obs_time <- min(evtim[i], td[i], followup)
    time[i] <- obs_time
    death[i] <- as.integer(is.finite(evtim[i]) && evtim[i] <= td[i] && evtim[i] <= followup)
    
    for (j in seq_along(t)) {
      if (t[j] <= obs_time) {
        Y1[i,j] <- Beta[1]+Beta[2]*x1[i]+Beta[3]*x2[i]+Beta[4]*t[j] +
          u[i,1]+u[i,2]*t[j]+rnorm(1,0,sigma)
      } else Y1[i,j] <- NA_real_
    }
    ct[i] <- sum(!is.na(Y1[i,]))
  }
  
  # ----------------------------
  # Pack data into long format
  # ----------------------------
  mat <- matrix(NA_real_, ncol = 7, nrow = nsample*length(t))
  colnames(mat) <- c("obstime","id","x1","x2","Event","survtime","Y1")
  row_idx <- 1L
  for (i in 1:nsample) {
    if (ct[i] > 0) {
      idxs <- which(!is.na(Y1[i,]))
      k <- length(idxs)
      rows <- row_idx:(row_idx+k-1L)
      mat[rows,1] <- t[idxs]
      mat[rows,2] <- i
      mat[rows,3] <- x1[i]
      mat[rows,4] <- x2[i]
      mat[rows,5] <- death[i]
      mat[rows,6] <- time[i]
      mat[rows,7] <- Y1[i,idxs]
      row_idx <- row_idx+k
    }
  }
  mat <- mat[1:(row_idx-1L), , drop=FALSE]
  long.data <- as.data.frame(mat)
  surv.data <- data.frame(id=id, w1=w1, w2=w2, survtime=time, Event=death)
  
  # ----------------------------
  # Train/validation split
  # ----------------------------
  train_ids <- 1:500
  valid_ids <- setdiff(id, train_ids)
  long_train <- subset(long.data, id %in% train_ids)
  long_valid <- subset(long.data, id %in% valid_ids)
  surv_train <- subset(surv.data, id %in% train_ids)
  surv_valid <- subset(surv.data, id %in% valid_ids)
  
  
  
  # Quick checks
  
  mean(table(long.data$id))
  median(table(long.data$id))
  min(table(long.data$id))
  max(table(long.data$id))
  
  mean(surv.data$survtime)
  mean(surv.data$Event)
  
  
  ################################################################################
  ################################################################################
  ################################################################################
  
   
  #===============================================================================
  #============================ TS  ==============================================
  #===============================================================================
  start1<- Sys.time()
  
  
  library(lme4)
  library(survival)
  library(dplyr)
  
  # ---- Stage 1: Fit linear mixed model on training data ----
  lme_fit <- lmer(Y1 ~ obstime + x1 + x2 + (1 + obstime | id), data = long_train)
  
  # Extract fixed and random effects
  fixef_vals <- fixef(lme_fit)
  ranef_df <- ranef(lme_fit)$id
  colnames(ranef_df) <- c("b0", "b1")
  ranef_df$id <- as.numeric(rownames(ranef_df))
  
  # ---- Stage 2: Construct predicted longitudinal values ----
  long_train_td <- long_train %>%
    left_join(ranef_df, by = "id") %>%
    mutate(
      Y1_pred = b0 + b1 * obstime + fixef_vals["x1"] * x1 + fixef_vals["x2"] * x2
    )
  
  # ---- Prepare survival data for tmerge ----
  surv_train_td <- surv_train %>%
    rename(Id = id, deathTimes = survtime, Event = Event)
  
  # ---- Baseline tmerge object ----
  td_base <- tmerge(
    data1 = surv_train_td, data2 = surv_train_td,
    id = Id, endpt = event(deathTimes, Event)
  )
  
  # ---- Add longitudinal time-dependent covariate ----
  # Make sure your long_train_td has: Id, time (obstime), Y1_pred
  long_train_td <- long_train_td %>% rename(Id = id, time = obstime)
  
  long_data_td <- tmerge(
    td_base, long_train_td,
    id = Id,
    lY1 = tdc(time, Y1_pred)
  )
  
  # ---- Fit extended Cox model ----
  cox_mts <- coxph(
    Surv(tstart, tstop, endpt) ~ lY1 + w1 + w2,
    data = long_data_td,
    id = Id,
    cluster = Id
  )
  
  Esti[kkk,]=summary(cox_mts)$coefficients[, "coef"]
  
  # 2. Extract standard errors
  se <- summary(cox_mts)$coefficients[, "se(coef)"]
  Real=c(alpha,gamma_w)
  # 3. 95% confidence intervals
  ci <- confint(cox_mts)  # default is 95%
  ci
  for(ii in 1:length(Real)){ 
    if((Real[ii]>ci[ii,1]) & (Real[ii]<ci[ii,2]))(CovR[kkk,ii]=1)
  }
  
  
  # ---- Dynamic predictions for validation using BLUPs ----
  dp_results <- list()
  Sigma_b <- as.matrix(VarCorr(lme_fit)$id)[1:2, 1:2]  # 2x2 covariance of random effects
  sigma2  <- sigma(lme_fit)^2
  beta_hat <- fixef(lme_fit)
  
  S <- c(0.25, 0.5,.75, 1)   # landmark times
  t_pred <- 0.25         # horizon
  est4=sd4=matrix(0,length(S),2)
  
  r=0
  for (s in S) {
    r=r+1
    # Subjects at risk at landmark s
    valid_at_risk <- subset(surv_valid, survtime >= s)
    valid_ids <- valid_at_risk$id
    Y_pred_valid <- numeric(length(valid_ids))
    
    for (i in seq_along(valid_ids)) {
      idv <- valid_ids[i]
      subj_long <- subset(long_valid, id == idv & obstime <= s)
      
      # baseline covariates
      x1i <- if (nrow(subj_long) > 0) subj_long$x1[1] else valid_at_risk$x1[valid_at_risk$id == idv]
      x2i <- if (nrow(subj_long) > 0) subj_long$x2[1] else valid_at_risk$x2[valid_at_risk$id == idv]
      
      # design for fixed effects at landmark
      X_landmark <- model.matrix(~ obstime + x1 + x2, 
                                 data = data.frame(obstime = s, x1 = x1i, x2 = x2i))
      
      if (nrow(subj_long) == 0) {
        # No history -> fixed effects only
        Y_pred_valid[i] <- as.numeric(X_landmark %*% beta_hat)
      } else {
        # With history -> compute BLUP
        X_s <- model.matrix(~ obstime + x1 + x2, data = subj_long)
        Z_s <- cbind(1, subj_long$obstime)
        Y_s <- subj_long$Y1
        
        V <- Z_s %*% Sigma_b %*% t(Z_s) + sigma2 * diag(nrow(Z_s))
        b_hat <- Sigma_b %*% t(Z_s) %*% solve(V, (Y_s - X_s %*% beta_hat))
        
        Z_landmark <- c(1, s)
        Y_pred_valid[i] <- as.numeric(X_landmark %*% beta_hat + Z_landmark %*% b_hat)
      }
    }
    
    # merge predictions with covariates for Cox prediction
    valid_landmark <- merge(valid_at_risk, 
                            data.frame(id = valid_ids, lY1 = Y_pred_valid),
                            by = "id", sort = FALSE)
    
    # baseline cumulative hazard
    H0_df <- basehaz(cox_mts, centered = FALSE)
    H0_tpred <- approx(H0_df$time, H0_df$hazard, xout = s + t_pred, rule = 2)$y
    
    # linear predictor and survival
    lp_valid <- predict(cox_mts, newdata = valid_landmark, type = "lp")
    surv_pred <- exp(-H0_tpred * exp(lp_valid))
    
    # store results
    dp_results[[paste0("S_", s)]] <- data.frame(
      id = valid_ids,
      surv_pred = surv_pred,
      landmark = s,
      horizon = s + t_pred
    )
    
    
    
    
    
    
    Crit <- Criteria(
      s = s,
      t = t_pred,
      Survt = valid_at_risk$survtime,
      CR = valid_at_risk$Event,
      P = 1-dp_results[[paste0("S_", s)]][,2],
      cause = 1
    )
    
    est4[r,] <- Crit$Cri[, 1]
    sd4[r,]  <- Crit$Cri[, 2]
    
  }
  
  
  
  res_inc4 <- compute_iAUC_iBS(
    s = S,
    auc_mat=cbind(est4[,1],sd4[,1]), bs_mat=cbind(est4[,2],sd4[,2]), 
    survtime=surv_valid$survtime, death=surv_valid$Event, t_pred=t_pred,
    type = "incident",     # or "cumulative"
    estimate_sd = TRUE, 
    B = 200
  )
  
  
  c(res_inc4$iAUC,res_inc4$sd_iAUC)
  c(res_inc4$iBS,res_inc4$sd_iBS)
  est4
  
  end1 <- Sys.time()
  time_TS=difftime(end1,start1,units ="mins")
  
  Results[[kkk]]=list(time_TS=time_TS,est4=est4,AUCi_ts=res_inc4$iAUC,BSi_ts=res_inc4$iBS,
                      Esti=Esti,
                      CovR=CovR)
  
  print(rep(kkk,10))
  
}













Time1=matrix(0,NN,2)
Cr1=Cr2=matrix(0,NN,10)
Est=CR=matrix(0,NN,2)
for(kkk in c(1:NN)[-35]){ 
  
  Time1[kkk,1]=Results[[kkk]]$time_TS

  Cr1[kkk,]=c(as.numeric(Results[[kkk]]$est4),c(Results[[kkk]]$AUCi_ts,Results[[kkk]]$BSi_ts))

  
}

apply(Results[[kkk]]$Esti[-35,],2,mean)
apply(Results[[kkk]]$CovR[-35,],2,mean)



apply(Cr1,2,mean)
apply(Cr2,2,mean)


apply(Time1,2,mean)



cbind(apply(Time1,2,mean),apply(Time1,2,sd))


cbind(apply(Cr1,2,mean),apply(Cr1,2,sd))
cbind(apply(Cr2,2,mean),apply(Cr2,2,sd))

 

MSE <- matrix(0, NN, length(Real))
for(i in c(1:NN)[-35]){ 
  MSE[i, ] <- (Results[[kkk]]$Esti[i,] - Real)^2
}

# Compute summary metrics
Res <- cbind(
  Real,
  apply(Results[[kkk]]$Esti[-35,],2,mean),
  apply(Results[[kkk]]$Esti[-35,],2,sd),
  (apply(Results[[kkk]]$Esti[-35,],2,mean) - Real) / Real, # Relative bias
  sqrt(apply(MSE, 2, mean)),     # RMSE
  apply(Results[[kkk]]$CovR[-35,],2,mean)
)
colnames(Res) <- c("Real", "Est", "sd", "Rbias", "RMSE", "CR")

xtable::xtable(Res,digits=3)

 
