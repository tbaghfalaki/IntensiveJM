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
  
  
  
  #===============================================================================
  ################################################################################
  ################################ real values of DP ############################# 
  ################################################################################
  S <- c(0.25, 0.5, 0.75, 1)
  t_pred <- 0.25  # horizon
  
  estr=sdr=matrix(0,length(S),2)
  
  dp_results <- list()
  r=0
  for (s in S) {
    r=r+1
    # Subjects still at risk at landmark
    valid_at_risk <- subset(surv_valid, survtime >= s)
    valid_ids <- valid_at_risk$id
    
    Y_pred_valid <- numeric(length(valid_ids))
    
    for (i in seq_along(valid_ids)) {
      idv <- valid_ids[i]
      
      # Last observed longitudinal value before landmark
      subj_long <- subset(long_valid, id == idv & obstime <= s)
      
      if (nrow(subj_long) == 0) {
        # No history: use fixed effects only
        X_landmark <- c(1, x1[idv], x2[idv], s)  # intercept, time, x1, x2
        Y_pred_valid[i] <- sum(Beta * X_landmark)  # Beta[1]+Beta[2]*x1+...
      } else {
        # Use true random effects (u) + fixed effects
        X_landmark <- c(1, x1[idv], x2[idv], s)
        Z_landmark <- c(1, s)
        Y_pred_valid[i] <- sum(Beta * X_landmark) + sum(Z_landmark * u[idv,])
      }
    }
    
    # Compute “true” survival probability at horizon s+t_pred
    surv_pred <- numeric(length(valid_ids))
    
    for (i in seq_along(valid_ids)) {
      idv <- valid_ids[i]
      
      # Integrate hazard from s to s + t_pred
      haz_fun <- function(t) {
        gamma_w[1]*w1[idv] + gamma_w[2]*w2[idv] + alpha*(
          Beta[1]  + Beta[2]*x1[idv] + Beta[2]*x2[idv]+ Beta[4]*t  + 
            u[idv,1] + u[idv,2]*t
        )
      }
      
      cumhaz <- integrate(function(t) exp(haz_fun(t)), lower = s, upper = s + t_pred)$value
      surv_pred[i] <- exp(-cumhaz)
    }
    
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
    
    estr[r,] <- Crit$Cri[, 1]
    sdr[r,]  <- Crit$Cri[, 2]
  }
  
  
  estr
  
  # Compute incremental iAUC and iBS
  res_incr <- compute_iAUC_iBS(
    s = S,
    auc_mat = cbind(estr[, 1], sdr[, 1]),
    bs_mat  = cbind(estr[, 2], sdr[, 2]), 
    survtime = surv_valid$survtime,
    death    = surv_valid$Event,
    t_pred   = t_pred,
    type     = "incident",
    estimate_sd = TRUE, 
    B = 200
  )
  
  # Results
  c(res_incr$iAUC, res_incr$sd_iAUC)
  c(res_incr$iBS,  res_incr$sd_iBS)
  estr
  
  
  #===============================================================================
  
  start1<- Sys.time()
  
  ############# laandmarking locf
  S=c(0.25,0.5,0.75,1)
  t_pred=0.25
  est=sd=matrix(0,length(S),2)
  
  # ---------------- Landmarking (LOCF only) predictions on validation (CORRECTED) ----------------
  dp_results <- vector("list", length(S))
  names(dp_results) <- paste0("LM_", S)
  for(j in seq_along(S)){
    s <- S[j]
    
    ## subjects at risk at the landmark (train / valid)
    train_at_risk <- subset(surv_train, survtime >= s)
    valid_at_risk <- subset(surv_valid, survtime >= s)
    
    ## longitudinal observations up to the landmark
    long_train_s <- subset(long_train, id %in% train_at_risk$id & obstime <= s)
    long_valid_s <- subset(long_valid, id %in% valid_at_risk$id & obstime <= s)
    
    ## LOCF summary for training (one row per id: last obs before or at s)
    if(nrow(long_train_s) > 0){
      tmp_train <- lapply(split(long_train_s, long_train_s$id),
                          function(df) df[which.max(df$obstime), c("id","Y1")])
      long_locf_train <- do.call(rbind, tmp_train)
      long_locf_train <- as.data.frame(long_locf_train, stringsAsFactors = FALSE)
      names(long_locf_train) <- c("id","Y1")
      long_locf_train$id <- as.integer(long_locf_train$id)
      long_locf_train$Y1 <- as.numeric(long_locf_train$Y1)
    } else {
      long_locf_train <- data.frame(id=integer(0), Y1=numeric(0))
    }
    
    ## LOCF summary for validation
    if(nrow(long_valid_s) > 0){
      tmp_valid <- lapply(split(long_valid_s, long_valid_s$id),
                          function(df) df[which.max(df$obstime), c("id","Y1")])
      long_locf_valid <- do.call(rbind, tmp_valid)
      long_locf_valid <- as.data.frame(long_locf_valid, stringsAsFactors = FALSE)
      names(long_locf_valid) <- c("id","Y1")
      long_locf_valid$id <- as.integer(long_locf_valid$id)
      long_locf_valid$Y1 <- as.numeric(long_locf_valid$Y1)
    } else {
      long_locf_valid <- data.frame(id=integer(0), Y1=numeric(0))
    }
    
    ## Merge LOCF summaries with survival (keep order of train_at_risk / valid_at_risk)
    train_landmark <- merge(train_at_risk, long_locf_train, by = "id", sort = FALSE)
    valid_landmark <- merge(valid_at_risk, long_locf_valid, by = "id", sort = FALSE)
    
    ## If not enough data to fit or no valid subjects at risk, return NA preds
    if(nrow(train_landmark) < 1 || nrow(valid_at_risk) == 0){
      out_all <- data.frame(id = valid_at_risk$id,
                            surv_pred = NA_real_,
                            landmark = s,
                            horizon = s + t_pred)
      dp_results[[j]] <- out_all
      next
    }
    if(nrow(valid_landmark) == 0){
      # none of the validation-at-risk subjects had a prior measurement -> all NA
      out_all <- data.frame(id = valid_at_risk$id,
                            surv_pred = NA_real_,
                            landmark = s,
                            horizon = s + t_pred)
      dp_results[[j]] <- out_all
      next
    }
    
    ## Fit Cox on training landmark data (time reset at s)
    cox_fit <- coxph(Surv(survtime - s, Event) ~ Y1 + w1 + w2, data = train_landmark, ties = "breslow")
    
    ## Baseline cumulative hazard H0(t) from the fitted Cox model
    H0_df <- basehaz(cox_fit, centered = FALSE)   # columns: time, hazard
    # get H0 at t_pred (time since landmark)
    H0_tpred <- approx(x = H0_df$time, y = H0_df$hazard, xout = t_pred, rule = 2)$y
    # rule = 2 -> use last value if t_pred is beyond last event time (mimics extend=TRUE)
    
    ## Subject-specific linear predictors for validation subjects (must match model vars)
    lp_valid <- predict(cox_fit, newdata = valid_landmark, type = "lp")
    
    ## predicted survival at t_pred for each validation subject: S_i(t) = exp(-H0(t) * exp(lp_i))
    preds <- as.numeric(exp(- H0_tpred * exp(lp_valid)))
    
    ## Put predictions back into the full set of validation-at-risk subjects (NA for those without LOCF)
    out_all <- data.frame(id = valid_at_risk$id, surv_pred = NA_real_)
    m <- match(valid_landmark$id, out_all$id)
    out_all$surv_pred[m] <- preds
    
    ## annotate and store
    out_all$landmark <- s
    out_all$horizon  <- s + t_pred
    
    dp_results[[j]] <- out_all
  }
  
  # Combined results for all landmarks (validation)
  dp_valid_all <- do.call(rbind, dp_results)
  
  # Quick check
  head(dp_valid_all, 20)
  
  
  for(kk in 1:length(S)){ 
    s=S[kk]
    # take predictions for s = 0.25
    tmp <- subset(dp_valid_all, landmark == s)
    head(tmp)
    # keep only id and event probability
    tmp <- tmp[, c("id", "surv_pred")]
    
    # merge with survival data to align
    merged <- merge(surv_valid, tmp, by = "id", all.x = TRUE)
    head(merged)
    
    # subjects who failed before s=0.25 have no prediction → set NA or 0
    merged$surv_pred[is.na(merged$surv_pred)] <- 0  
    
    # now compute Criteria
    Crit <- Criteria(
      s = s,
      t = t_pred,
      Survt = merged$survtime,
      CR = merged$Event,
      P = 1-merged$surv_pred,
      cause = 1
    )
    
    Crit
    est[kk,]=Crit$Cri[,1]
    sd[kk,]=Crit$Cri[,2]
    
    # print(rep(kk,9))
  }
  
  
  colnames(est)=c("est-AUC","est-BS")
  colnames(sd)=c("sd-AUC","sd-BS")
  
  rownames(est)=rownames(sd)=S
  est
  sd
  
  
  
  
  res_inc <- compute_iAUC_iBS(
    s = S,
    auc_mat=cbind(est[,1],sd[,1]), bs_mat=cbind(est[,2],sd[,2]), 
    survtime=surv_valid$survtime, death=surv_valid$Event, t_pred=t_pred,
    type = "incident",     # or "cumulative"
    estimate_sd = TRUE, 
    B = 200
  )
  
  c(res_inc$iAUC,res_inc$sd_iAUC)
  c(res_inc$iBS,res_inc$sd_iBS)
  
  end1 <- Sys.time()
  landmarking_locf=difftime(end1,start1,units ="mins")
  
  
  ############ ############ ############ ############ ############ ############ ############ 
  ############ ############ ############  landmarking based on lmm ############ ############ 
  ############ ############ ############ ############ ############ ############ ############ 
  library(lme4)
  library(survival)
  start1<- Sys.time()
  
  S <- c(0.25, 0.5, 0.75, 1)
  t_pred <- 0.25
  dp_results_mixed <- vector("list", length(S))
  names(dp_results_mixed) <- paste0("LM_Mixed_", S)
  
  est1 <- sd1 <- matrix(NA_real_, length(S), 2)
  
  for (j in seq_along(S)) {
    s <- S[j]
    
    train_at_risk <- subset(surv_train, survtime >= s)
    valid_at_risk <- subset(surv_valid, survtime >= s)
    
    long_train_s <- subset(long_train, id %in% train_at_risk$id & obstime <= s)
    long_valid_s <- subset(long_valid, id %in% valid_at_risk$id & obstime <= s)
    
    if (nrow(long_train_s) < 1 || nrow(valid_at_risk) == 0) {
      dp_results_mixed[[j]] <- data.frame(
        id = valid_at_risk$id,
        surv_pred = NA_real_,
        landmark = s,
        horizon = s + t_pred
      )
      next
    }
    
    # Fit LMM on training longitudinal up to s
    lmm_fit <- lmer(Y1 ~ obstime + x1 + x2 + (1 + obstime | id),
                    data = long_train_s, REML = FALSE)
    summary(lmm_fit)
    # --- Predictions for training subjects (BLUP or fixed if none) ---
    train_ids <- unique(train_at_risk$id)
    
    Y_pred_train=c()
    for(i in 1:dim(train_at_risk)[1]){ 
      x1i <- long_train_s$x1[long_train_s$id==train_at_risk$id[i]][1]
      x2i <-  long_train_s$x2[long_train_s$id==train_at_risk$id[i]][1]
      idv <- train_at_risk$id[i]
      newdat <- data.frame(obstime = s, x1 = x1i, x2 = x2i, id = idv)
      xnewdat <- c(1,obstime = s, x1 = x1i, x2 = x2i)
      b <- lme4::ranef(lmm_fit)$id
      
      Y_pred_train1=lmm_fit@beta%*%c(1,obstime = s, x1 = x1i, x2 = x2i)+
        b[i,1]+b[i,2]*s
      Y_pred_train <- append(Y_pred_train,Y_pred_train1)
    }
    
    train_at_risk$Ypred=Y_pred_train
    names(train_at_risk)
    # Remove any rows with NA Y1_pred
    
    # Fit Cox model
    cox_fit <- coxph(Surv(survtime - s, Event) ~ Ypred + w1 + w2,
                     data = train_at_risk, ties = "breslow")
    summary(cox_fit)
    
    
    # --- Predictions for validation (for subjects at risk) ---
    valid_ids_atrisk <- valid_at_risk$id
    
    
    #############
    Y_pred_valid <- c()
    
    Sigma_b <- as.matrix(VarCorr(lmm_fit)$id)[1:2, 1:2]  # 2x2 G matrix
    sigma2  <- sigma(lmm_fit)^2
    beta_hat <- lmm_fit@beta
    
    for (i in 1:nrow(valid_at_risk)) { 
      idv <- valid_at_risk$id[i]
      subj_long <- subset(long_valid_s, id == idv & obstime <= s)
      
      # baseline covariates (constant per subject)
      x1i <- if (nrow(subj_long) > 0) subj_long$x1[1] else valid_at_risk$x1[i]
      x2i <- if (nrow(subj_long) > 0) subj_long$x2[1] else valid_at_risk$x2[i]
      
      # Design for prediction at time s
      X_landmark <- model.matrix(~ obstime + x1 + x2,
                                 data = data.frame(obstime = s, x1 = x1i, x2 = x2i))
      
      if (nrow(subj_long) == 0) {
        ## No history -> fixed effects only
        Y_pred_valid1 <- as.numeric(X_landmark %*% beta_hat)
      } else {
        ## With history -> compute BLUP
        X_s <- model.matrix(~ obstime + x1 + x2, data = subj_long)
        Z_s <- cbind(1, subj_long$obstime)
        Y_s <- subj_long$Y1
        
        V <- Z_s %*% Sigma_b %*% t(Z_s) + sigma2 * diag(nrow(Z_s))
        b_hat <- Sigma_b %*% t(Z_s) %*% solve(V, (Y_s - X_s %*% beta_hat))
        
        Z_landmark <- c(1, s)
        Y_pred_valid1 <- as.numeric(X_landmark %*% beta_hat + Z_landmark %*% b_hat)
      }
      
      Y_pred_valid <- c(Y_pred_valid, Y_pred_valid1)
    }
    
    
    
    
    valid_landmark <- data.frame(id = valid_ids_atrisk, Ypred = Y_pred_valid)
    # merge other covariates from valid_at_risk (w1,w2,survtime,Event)
    valid_landmark <- merge(valid_at_risk, valid_landmark, by = "id", sort = FALSE)
    
    # Baseline cumulative hazard
    H0_df <- basehaz(cox_fit, centered = FALSE)
    H0_tpred <- approx(H0_df$time, H0_df$hazard, xout = t_pred, rule = 2)$y
    
    # Linear predictor for validation
    lp_valid <- predict(cox_fit, newdata = valid_landmark, type = "lp")
    surv_pred <- exp(-H0_tpred * exp(lp_valid))   # survival prob to horizon (from landmark)
    
    # Build output aligned to valid_at_risk order
    out_all <- data.frame(id = valid_at_risk$id, surv_pred = NA_real_)
    m <- match(valid_landmark$id, out_all$id)
    out_all$surv_pred[m] <- surv_pred
    out_all$landmark <- s
    out_all$horizon <- s + t_pred
    dp_results_mixed[[j]] <- out_all
    
    # --- Evaluate Criteria: must evaluate on same set & order
    # Merge predictions with surv_valid restricted to at-risk subjects, in same order
    eval_df <- merge(valid_at_risk, out_all, by = "id", sort = FALSE)
    # event probability (P) = 1 - survival probability
    P_event <- 1 - eval_df$surv_pred
    
    Crit <- Criteria(
      s = s,
      t = t_pred,
      Survt = eval_df$survtime,
      CR = eval_df$Event,
      P = P_event,
      cause = 1
    )
    
    est1[j, ] <- Crit$Cri[, 1]
    sd1[j, ]  <- Crit$Cri[, 2]
  }
  
  est1
  sd1
  
  
  
  res_inc1 <- compute_iAUC_iBS(
    s = S,
    auc_mat=cbind(est1[,1],sd1[,1]), bs_mat=cbind(est1[,2],sd1[,2]), 
    survtime=surv_valid$survtime, death=surv_valid$Event, t_pred=t_pred,
    type = "incident",     # or "cumulative"
    estimate_sd = TRUE, 
    B = 200
  )
  
  
  
  end1 <- Sys.time()
  landmarking_lmm=difftime(end1,start1,units ="mins")
  
  
  Results1[[kkk]]=list(iAUC_lmm=res_inc1$iAUC,
                       iBS_lmm=res_inc1$iBS,
                       est1=est1,
                       landmarking_lmm=landmarking_lmm,
                       
                       iAUC_locf=res_inc$iAUC,
                       iBS_locf=res_inc$iBS,
                       est=est,
                       landmarking_locf=landmarking_locf,
                       
                       Real=estr,
                       iAUC_real=res_incr$iAUC, 
                       iBS_real=res_incr$iBS)
  
  
  print(rep(kkk,10))
  
  
}









Time1=matrix(0,NN,2)
Cr0=Cr1=Cr2=matrix(0,NN,10)
for(kkk in 1:NN){ 
  
  Time1[kkk,1]=Results1[[kkk]]$landmarking_lmm
  Time1[kkk,2]=Results1[[kkk]]$landmarking_locf
  
  Cr0[kkk,]=c(as.numeric(Results1[[kkk]]$Real),c(Results1[[kkk]]$iAUC_real,Results1[[kkk]]$iBS_real))
  Cr1[kkk,]=c(as.numeric(Results1[[kkk]]$est1),c(Results1[[kkk]]$iAUC_lmm,Results1[[kkk]]$iBS_lmm))
  Cr2[kkk,]=c(as.numeric(Results1[[kkk]]$est),c(Results1[[kkk]]$iAUC_locf,Results1[[kkk]]$iBS_locf))
}

cbind(apply(Time1,2,mean),apply(Time1,2,sd))


cbind(apply(Cr0,2,mean),apply(Cr0,2,sd))
cbind(apply(Cr1,2,mean),apply(Cr1,2,sd))
cbind(apply(Cr2,2,mean),apply(Cr2,2,sd))






