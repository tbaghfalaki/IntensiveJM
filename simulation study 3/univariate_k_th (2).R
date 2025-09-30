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
library(dplyr)

library(JMbayes2)
library(refund)

setwd("/Users/user/Desktop/Rcodes/PP_TM/New_stan/11th")#
#setwd("C:\\Users\\p80744tb\\Desktop\\PP_TM\\PP_TM\\New_stan\\11th")
source("iauc_sd.R")
# ----------------------------
# Simulation parameters
# ----------------------------

NN=100
Estijm=CovRjm=matrix(0,NN,3)
Results1=list()
for(kkk in 1:NN){ 
  set.seed(kkk)
  
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
      alpha*(Beta[1] + Beta[2]*x1i + Beta[3]*x2i + sin(pi*s)+cos(pi*s) + ui[1] + ui[2]*s)
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
        Y1[i,j] <- Beta[1]+Beta[2]*x1[i]+Beta[3]*x2[i]+ sin(pi*t[j])+cos(pi*t[j])   +
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
  
  
  library(dplyr)
  
  
  downsample_long <- function(data, kk = 3) {
    data %>%
      group_by(id) %>%
      arrange(obstime, .by_group = TRUE) %>%
      mutate(idx = row_number()) %>%
      filter((idx - 1) %% kk == 0) %>%  # pick 0, kk, 2kk, ...
      select(-idx) %>%
      ungroup()
  }
  
  # Example usage
  long.data.down3 <- downsample_long(long.data, kk = 3)
  long.data.down5 <- downsample_long(long.data, kk = 5)
  
  
  long.data=long.data.down5
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
  start1<- Sys.time()
  
  # ----------------------------
  # Fit longitudinal mixed model
  # ----------------------------
  # ----------------------------
  fm1 <- lme(fixed = Y1 ~ x1 + x2 + ns(obstime, df=4),
             random = ~ obstime | id,
             data = long_train)
  
  
  # ----------------------------
  # Fit Cox model
  # ----------------------------
  fCox1 <- coxph(Surv(survtime, Event) ~ w1 + w2,
                 data = surv_train)
  
  # ----------------------------
  # Fit joint model
  # ----------------------------
  joint_model_fit_1 <- jm(fCox1, fm1, time_var="obstime",
                          n_chains=1L, n_iter=2000L, n_burnin=1000L)
  
  
  
  Real=c(alpha,gamma_w)
  
  # Extract MCMC samples
  mcmc_samples <- joint_model_fit_1$mcmc
  
  ci=rbind(summary(mcmc_samples$alphas)[2]$quantiles[c(1,5)],
           summary(mcmc_samples$gammas)[2]$quantiles[,c(1,5)])
  
  
  
  Estijm[kkk,]=c(summary(mcmc_samples$alphas)[1]$statistics[1],summary(mcmc_samples$gammas)[1]$statistics[,1])
  
  
  for(ii in 1:length(Real)){ 
    if((Real[ii]>ci[ii,1]) & (Real[ii]<ci[ii,2]))(CovRjm[kkk,ii]=1)
  }
  
  
  # ----------------------------
  # Dynamic predictions for id
  # ----------------------------
  
  S <- c(0.25, 0.5,.75, 1)   # landmark times
  t_pred <- 0.25         # horizon
  
  # container for predictions
  all_preds <- list()
  
  for (id_i in valid_ids) {
    for (s in S) {
      # longitudinal up to landmark
      this_long <- subset(long_valid, id == id_i & obstime <= s)
      if (nrow(this_long) == 0) next
      
      # add baseline covariates from surv_valid
      this_long$w1 <- surv_valid$w1[surv_valid$id == id_i]
      this_long$w2 <- surv_valid$w2[surv_valid$id == id_i]
      this_long$survtime <- s
      this_long$Event <- 0
      
      pred_times <- s + t_pred
      
      pr <- tryCatch(
        predict(joint_model_fit_1,
                newdata = this_long,
                times = pred_times,
                process = "event",
                type = "subject_specific",
                return_newdata = TRUE),
        error = function(e) NULL
      )
      
      if (!is.null(pr)) {
        all_preds[[paste0("id_", id_i, "_s_", s)]] <- list(
          id = id_i,
          landmark = s,
          pred_obj = pr
        )
      }
    }
  }
  
  # ----------------------------
  # Convert all to tidy tibble
  # ----------------------------
  tidy_preds_all <- bind_rows(lapply(all_preds, function(el) {
    pr <- el$pred_obj
    if (is.null(pr)) return(NULL)
    
    pr %>%
      transmute(
        id = el$id,
        landmark = el$landmark,
        time = survtime,
        horizon = survtime - el$landmark,
        pred = pred_CIF,
        lower = low_CIF,
        upper = upp_CIF
      )
  }))
  
  # final tidy predictions
  #print(tidy_preds_all)
  
  # --- keep only horizon predictions (exclude time == landmark)
  tidy_preds_all <- bind_rows(lapply(all_preds, function(el) {
    pr <- el$pred_obj
    if (is.null(pr)) return(NULL)
    
    pr %>%
      filter(survtime > el$landmark) %>%   # remove landmark == time
      transmute(
        id = el$id,
        landmark = el$landmark,
        time = survtime,
        horizon = survtime - el$landmark,
        pred = pred_CIF
      )
  }))
  
  # --- split by landmark time
  preds_by_landmark <- split(tidy_preds_all, tidy_preds_all$landmark)
  
  # Example: predictions for each landmark
  preds_025 <- preds_by_landmark[["0.25"]] %>% select(id, pred)
  preds_050 <- preds_by_landmark[["0.5"]]  %>% select(id, pred)
  preds_075 <- preds_by_landmark[["0.75"]]  %>% select(id, pred)
  preds_100 <- preds_by_landmark[["1"]]    %>% select(id, pred)
  
  # check outputs
  head(preds_025)
  head(preds_050)
  head(preds_075)
  head(preds_100)
  
  S=c(0.25,0.5,0.75,1)
  t_pred=0.25
  est2=sd2=matrix(0,length(S),2)
  
  
  Crit <- Criteria(
    s = .25,
    t = t_pred,
    Survt = surv_valid$survtime,
    CR = surv_valid$Event,
    P = as.numeric(unlist(preds_025[,2])),
    cause = 1
  )
  
  est2[1,] <- Crit$Cri[, 1]
  sd2[1,]  <- Crit$Cri[, 2]
  
  
  Crit <- Criteria(
    s = .5,
    t = t_pred,
    Survt = surv_valid$survtime,
    CR = surv_valid$Event,
    P = as.numeric(unlist(preds_050[,2])),
    cause = 1
  )
  
  est2[2, ] <- Crit$Cri[, 1]
  sd2[2,]  <- Crit$Cri[, 2]
  
  
  Crit <- Criteria(
    s = .75,
    t = t_pred,
    Survt = surv_valid$survtime,
    CR = surv_valid$Event,
    P = as.numeric(unlist(preds_075[,2])),
    cause = 1
  )
  
  est2[3,] <- Crit$Cri[, 1]
  sd2[3, ]  <- Crit$Cri[, 2]
  
  Crit <- Criteria(
    s = 1,
    t = t_pred,
    Survt = surv_valid$survtime,
    CR = surv_valid$Event,
    P = as.numeric(unlist(preds_075[,2])),
    cause = 1
  )
  
  est2[4,] <- Crit$Cri[, 1]
  sd2[4, ]  <- Crit$Cri[, 2]
  
  
  
  res_inc2 <- compute_iAUC_iBS(
    s = S,
    auc_mat=cbind(est2[,1],sd2[,1]), bs_mat=cbind(est2[,2],sd2[,2]), 
    survtime=surv_valid$survtime, death=surv_valid$Event, t_pred=t_pred,
    type = "incident",     # or "cumulative"
    estimate_sd = TRUE, 
    B = 200
  )
  
  
  c(res_inc2$iAUC,res_inc2$sd_iAUC)
  c(res_inc2$iBS,res_inc2$sd_iBS)
  est2
  
  end1 <- Sys.time()
  TimeMulti=difftime(end1,start1,units ="mins")
  
  
  
  
  Results1[[kkk]]=list(iAUC_jm=res_inc2$iAUC,
                       iBS_jm=res_inc2$iBS,
                       est2=est2,
                       TimeMulti=TimeMulti,Estijm=Estijm,CovRjm=CovRjm)
  
  
  print(rep(kkk,10))
  
  
}




Time1=c()
Cr1=Cr2=matrix(0,NN,10)
Est=CR=matrix(0,NN,2)
for(kkk in c(1:NN)){ 
  
  Time1[kkk]=Results1[[kkk]]$TimeMulti
  
  Cr1[kkk,]=c(as.numeric(Results1[[kkk]]$est2),c(Results1[[kkk]]$iAUC_jm,Results1[[kkk]]$iBS_jm))
  
  
}

apply(Results1[[kkk]]$Esti,2,mean)
apply(Results1[[kkk]]$CovR,2,mean)
AA=cbind(apply(Cr1,2,mean),apply(Cr1,2,sd))
mean(Time1[-59]);sd(Time1[-59]);
xtable::xtable(AA,digits=3)

MSE <- matrix(0, NN, length(Real))
for (i in 1:NN) {
  MSE[i, ] <- (Results1[[kkk]]$Esti[i,] - Real)^2
}

# Compute summary metrics
Res <- cbind(
  Real,
  apply(Results1[[kkk]]$Esti,2,mean),
  apply(Results1[[kkk]]$Esti,2,sd),
  (apply(Results1[[kkk]]$Esti,2,mean) - Real) / Real, # Relative bias
  sqrt(apply(MSE, 2, mean)),     # RMSE
  apply(Results1[[kkk]]$CovR,2,mean)
)
colnames(Res) <- c("Real", "Est", "sd", "Rbias", "RMSE", "CR")

xtable::xtable(Res,digits=3)



mean(Time1[-59]);sd(Time1[-59]);
[1] 6.57772
[1] 1.983472
> xtable::xtable(AA,digits=3)
% latex table generated in R 4.5.1 by xtable 1.8-4 package
% Mon Sep 29 10:42:13 2025
\begin{table}[ht]
\centering
\begin{tabular}{rrr}
\hline
& 1 & 2 \\ 
\hline
1 & 0.694 & 0.044 \\ 
2 & 0.709 & 0.043 \\ 
3 & 0.726 & 0.051 \\ 
4 & 0.722 & 0.061 \\ 
5 & 0.112 & 0.016 \\ 
6 & 0.140 & 0.018 \\ 
7 & 0.170 & 0.022 \\ 
8 & 0.195 & 0.034 \\ 
9 & 0.708 & 0.025 \\ 
10 & 0.142 & 0.010 \\ 
\hline
\end{tabular}
\end{table}
> xtable::xtable(Res,digits=3)
% latex table generated in R 4.5.1 by xtable 1.8-4 package
% Mon Sep 29 10:42:16 2025
\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrr}
\hline
& Real & Est & sd & Rbias & RMSE & CR \\ 
\hline
1 & -0.500 & -0.498 & 0.058 & -0.005 & 0.057 & 0.910 \\ 
2 & 0.200 & 0.204 & 0.072 & 0.019 & 0.071 & 0.980 \\ 
3 & -0.200 & -0.241 & 0.142 & 0.204 & 0.147 & 0.940 \\ 
\hline
\end{tabular}
\end{table}
