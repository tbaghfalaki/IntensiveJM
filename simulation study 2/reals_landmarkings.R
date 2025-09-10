# ----------------------------
# Libraries
# ----------------------------
rm(list = ls())
library(MASS)
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvtnorm)
library(JMbayes2)
library(DPCri)



source("C:/Users/p80744tb/Desktop/FDA/Sep_2025/iauc_sd.R")
# ----------------------------
# Simulation parameters
# ----------------------------
NN=100

Estijm=CovRjm=matrix(0,NN,5)
Results1=list()


for(kkk in 1:NN){ 
  set.seed(kkk)
  nsample <- 1000
  followup <- 2
  t <- seq(0, followup, by = 0.005)
  
  # Association and fixed effect parameters
  alpha <- c(-0.5, 0.3, -0.2)
  Beta1 <- c(-0.5, -0.5, 0.5, 0.5)
  Beta2 <- c(0.2,  0.3, -0.3, 0.4)
  Beta3 <- c(1.0, -0.2,  0.1, 0.6)
  sigma <- 1
  gamma_w <- c(0.2, -0.2)
  
  # ----------------------------
  # Covariates and random effects
  # ----------------------------
  id <- 1:nsample
  x1 <- rnorm(nsample)
  x2 <- rbinom(nsample, 1, 0.5)
  w1 <- rnorm(nsample)
  w2 <- rbinom(nsample, 1, 0.5)
  death <- rep(0L, nsample)
  time <- rep(NA_real_, nsample)
  
  # Random effects: 3 markers, intercept + slope each
  D_block <- matrix(c(1,0.6,0.6,1), 2,2)
  D <- matrix(0.4,6,6)
  D[1:2,1:2] <- D_block
  D[3:4,3:4] <- D_block
  D[5:6,5:6] <- D_block
  u <- rmvnorm(nsample, mean = rep(0,6), sigma = D)
  
  # ----------------------------
  # Hazard functions
  # ----------------------------
  haz_fun_i <- function(s, x1i, x2i, w1i, w2i, ui) {
    mu1 <- Beta1[1] + Beta1[2]*x1i + Beta1[3]*x2i + Beta1[4]*s + ui[1] + ui[2]*s
    mu2 <- Beta2[1] + Beta2[2]*x1i + Beta2[3]*x2i + Beta2[4]*s + ui[3] + ui[4]*s
    mu3 <- Beta3[1] + Beta3[2]*x1i + Beta3[3]*x2i + Beta3[4]*s + ui[5] + ui[6]*s
    eta <- gamma_w[1]*w1i + gamma_w[2]*w2i + alpha[1]*mu1 + alpha[2]*mu2 + alpha[3]*mu3
    exp(eta)
  }
  
  cumhaz_i <- function(T, x1i, x2i, w1i, w2i, ui) {
    if(T<=0) return(0)
    integrate(function(s) haz_fun_i(s,x1i,x2i,w1i,w2i,ui), lower=0, upper=T, rel.tol=1e-6)$value
  }
  
  surv_i <- function(T, x1i, x2i, w1i, w2i, ui) exp(-cumhaz_i(T,x1i,x2i,w1i,w2i,ui))
  
  sim_event_time_i <- function(x1i, x2i, w1i, w2i, ui, Tmax=followup){
    U <- runif(1)
    S_Tmax <- surv_i(Tmax,x1i,x2i,w1i,w2i,ui)
    if(S_Tmax>U) return(Inf)
    f <- function(T) surv_i(T,x1i,x2i,w1i,w2i,ui)-U
    res <- tryCatch(uniroot(f, lower=.Machine$double.eps, upper=Tmax)$root, error=function(e) Tmax)
    res
  }
  
  # ----------------------------
  # Generate longitudinal and survival data
  # ----------------------------
  Y1 <- Y2 <- Y3 <- matrix(NA_real_, nsample, length(t))
  for(i in 1:nsample){
    evt <- sim_event_time_i(x1[i],x2[i],w1[i],w2[i], u[i,])
    td <- rexp(1,1)
    obs_time <- min(evt, td, followup)
    time[i] <- obs_time
    death[i] <- as.integer(is.finite(evt) && (evt<=td) && (evt<=followup))
    
    for(j in seq_along(t)){
      if(t[j]<=obs_time){
        Y1[i,j] <- Beta1[1] + Beta1[2]*x1[i] + Beta1[3]*x2[i] + Beta1[4]*t[j] + u[i,1] + u[i,2]*t[j] + rnorm(1,0,sigma)
        Y2[i,j] <- Beta2[1] + Beta2[2]*x1[i] + Beta2[3]*x2[i] + Beta2[4]*t[j] + u[i,3] + u[i,4]*t[j] + rnorm(1,0,sigma)
        Y3[i,j] <- Beta3[1] + Beta3[2]*x1[i] + Beta3[3]*x2[i] + Beta3[4]*t[j] + u[i,5] + u[i,6]*t[j] + rnorm(1,0,sigma)
      }
    }
  }
  
  # ----------------------------
  # Pack into long format with Y1,Y2,Y3 columns
  # ----------------------------
  long.list <- vector("list", nsample)
  for(i in 1:nsample){
    idxs <- which(!is.na(Y1[i,]))
    if(length(idxs)==0) next
    long.list[[i]] <- data.frame(
      id = i,
      obstime = t[idxs],
      x1 = x1[i],
      x2 = x2[i],
      w1 = w1[i],
      w2 = w2[i],
      Event = death[i],
      survtime = time[i],
      Y1 = Y1[i, idxs],
      Y2 = Y2[i, idxs],
      Y3 = Y3[i, idxs]
    )
  }
  
  long.data <- do.call(rbind, long.list)
  
  surv.data <- data.frame(id=id, x1=x1, x2=x2, w1=w1, w2=w2, survtime=time, Event=death)
  
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
  ############################### Real value #################################################
  # ----------------------------
  # Dynamic Prediction for 3 longitudinal markers
  # ----------------------------
  
  S <- c(0.25, 0.5, 0.75, 1)  # landmark times
  t_pred <- 0.25             # prediction horizon
  
  estr <- sdr <- matrix(0, length(S), 2)
  dp_results <- list()
  
  for (r in seq_along(S)) {
    s <- S[r]
    
    # Subjects at risk at landmark s
    valid_at_risk <- subset(surv_valid, survtime >= s)
    valid_ids <- valid_at_risk$id
    
    surv_pred <- numeric(length(valid_ids))
    
    for (i in seq_along(valid_ids)) {
      idv <- valid_ids[i]
      
      # Subset longitudinal data for this subject up to landmark
      subj_long <- subset(long_valid, id == idv & obstime <= s)
      
      # Initialize marker predictions at landmark
      mu_landmark <- numeric(3)
      
      for (m in 1:3) {  # loop over 3 markers
        # Pick Beta and random effects for marker m
        Beta_m <- switch(m, Beta1, Beta2, Beta3)
        u_m <- u[idv, (2*m-1):(2*m)]  # intercept + slope for marker m
        
        if (nrow(subj_long) == 0) {
          # No history: fixed effects only
          X_landmark <- c(1, x1[idv], x2[idv], s)
          mu_landmark[m] <- sum(Beta_m * X_landmark)
        } else {
          # Include random effects
          X_landmark <- c(1, x1[idv], x2[idv], s)
          Z_landmark <- c(1, s)
          mu_landmark[m] <- sum(Beta_m * X_landmark) + sum(Z_landmark * u_m)
        }
      }
      
      # Hazard integration for survival prediction
      haz_fun <- function(t) {
        # contributions from the three markers
        marker_contrib <- 0
        for (m in 1:3) {
          Beta_m <- switch(m, Beta1, Beta2, Beta3)
          u_m <- u[idv, (2*m-1):(2*m)]
          mu_t <- Beta_m[1] + Beta_m[2]*x1[idv] + Beta_m[3]*x2[idv] + Beta_m[4]*t +
            u_m[1] + u_m[2]*t
          marker_contrib <- marker_contrib + alpha[m] * mu_t
        }
        
        # full linear predictor
        eta <- gamma_w[1]*w1[idv] + gamma_w[2]*w2[idv] + marker_contrib
        exp(eta)
      }
      
      cumhaz <- integrate(haz_fun, lower = s, upper = s + t_pred, rel.tol = 1e-6)$value
      surv_pred[i] <- exp(-cumhaz)
    }
    
    dp_results[[paste0("S_", s)]] <- data.frame(
      id = valid_ids,
      surv_pred = surv_pred,
      landmark = s,
      horizon = s + t_pred
    )
    
    # Evaluate criteria (iAUC/iBS)
    Crit <- Criteria(
      s = s,
      t = t_pred,
      Survt = valid_at_risk$survtime,
      CR = valid_at_risk$Event,
      P = 1 - dp_results[[paste0("S_", s)]][, "surv_pred"],
      cause = 1
    )
    
    estr[r, ] <- Crit$Cri[, 1]
    sdr[r, ]  <- Crit$Cri[, 2]
  }
  
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
  list(
    iAUC = c(res_incr$iAUC, res_incr$sd_iAUC),
    iBS  = c(res_incr$iBS,  res_incr$sd_iBS),
    estr  = estr
  )
  
  #===============================================================================
  
  library(lme4)
  library(survival)
  
  # Start time
  start1 <- Sys.time()
  
  # Landmarks and prediction horizon
  S <- c(0.25, 0.5, 0.75, 1)
  t_pred <- 0.25
  
  # Store results
  dp_results_mixed <- vector("list", length(S))
  names(dp_results_mixed) <- paste0("LM_Mixed_", S)
  est_lmm <- sd_lmm <- matrix(NA_real_, length(S), 2)
  
  # Loop over landmarks
  for (j in seq_along(S)) {
    s <- S[j]
    
    # Subjects at risk at the landmark
    train_at_risk <- subset(surv_train, survtime >= s)
    valid_at_risk <- subset(surv_valid, survtime >= s)
    
    # Longitudinal data up to the landmark
    long_train_s <- subset(long_train, id %in% train_at_risk$id & obstime <= s)
    long_valid_s <- subset(long_valid, id %in% valid_at_risk$id & obstime <= s)
    
    if (nrow(train_at_risk) < 1 || nrow(valid_at_risk) < 1) {
      dp_results_mixed[[j]] <- data.frame(
        id = valid_at_risk$id,
        surv_pred = NA_real_,
        landmark = s,
        horizon = s + t_pred
      )
      next
    }
    
    # Fit LMM for longitudinal marker
    lmm_fit <- lmer(Y1 ~ obstime + x1 + x2 + (1 + obstime | id),
                    data = long_train_s, REML = FALSE)
    
    beta_hat <- lmm_fit@beta
    Sigma_b <- as.matrix(VarCorr(lmm_fit)$id)[1:2, 1:2]
    sigma2 <- sigma(lmm_fit)^2
    
    # --- Predictions for training subjects ---
    Ypred_train <- numeric(nrow(train_at_risk))
    for (i in seq_along(train_at_risk$id)) {
      idv <- train_at_risk$id[i]
      subj_long <- subset(long_train_s, id == idv & obstime <= s)
      
      X_landmark <- model.matrix(~ obstime + x1 + x2,
                                 data = data.frame(obstime = s,
                                                   x1 = train_at_risk$x1[i],
                                                   x2 = train_at_risk$x2[i]))
      
      if (nrow(subj_long) == 0) {
        # No prior measurement -> fixed effects only
        Ypred_train[i] <- as.numeric(X_landmark %*% beta_hat)
      } else {
        X_s <- model.matrix(~ obstime + x1 + x2, data = subj_long)
        Z_s <- cbind(1, subj_long$obstime)
        Y_s <- subj_long$Y1
        V <- Z_s %*% Sigma_b %*% t(Z_s) + sigma2 * diag(nrow(Z_s))
        b_hat <- Sigma_b %*% t(Z_s) %*% solve(V, (Y_s - X_s %*% beta_hat))
        Z_landmark <- c(1, s)
        Ypred_train[i] <- as.numeric(X_landmark %*% beta_hat + Z_landmark %*% b_hat)
      }
    }
    train_at_risk$Ypred <- Ypred_train
    
    # Fit Cox model
    cox_fit <- coxph(Surv(survtime - s, Event) ~ Ypred + w1 + w2,
                     data = train_at_risk, ties = "breslow")
    
    # --- Predictions for validation subjects ---
    Ypred_valid <- numeric(nrow(valid_at_risk))
    for (i in seq_along(valid_at_risk$id)) {
      idv <- valid_at_risk$id[i]
      subj_long <- subset(long_valid_s, id == idv & obstime <= s)
      
      X_landmark <- model.matrix(~ obstime + x1 + x2,
                                 data = data.frame(obstime = s,
                                                   x1 = valid_at_risk$x1[i],
                                                   x2 = valid_at_risk$x2[i]))
      if (nrow(subj_long) == 0) {
        Ypred_valid[i] <- as.numeric(X_landmark %*% beta_hat)
      } else {
        X_s <- model.matrix(~ obstime + x1 + x2, data = subj_long)
        Z_s <- cbind(1, subj_long$obstime)
        Y_s <- subj_long$Y1
        V <- Z_s %*% Sigma_b %*% t(Z_s) + sigma2 * diag(nrow(Z_s))
        b_hat <- Sigma_b %*% t(Z_s) %*% solve(V, (Y_s - X_s %*% beta_hat))
        Z_landmark <- c(1, s)
        Ypred_valid[i] <- as.numeric(X_landmark %*% beta_hat + Z_landmark %*% b_hat)
      }
    }
    
    valid_landmark <- data.frame(id = valid_at_risk$id, Ypred = Ypred_valid)
    valid_landmark <- merge(valid_at_risk, valid_landmark, by = "id", sort = FALSE)
    
    # Baseline cumulative hazard
    H0_df <- basehaz(cox_fit, centered = FALSE)
    H0_tpred <- approx(H0_df$time, H0_df$hazard, xout = t_pred, rule = 2)$y
    
    # Linear predictor for validation
    lp_valid <- predict(cox_fit, newdata = valid_landmark, type = "lp")
    surv_pred <- exp(-H0_tpred * exp(lp_valid))
    
    # Build output
    out_all <- data.frame(id = valid_at_risk$id, surv_pred = surv_pred,
                          landmark = s, horizon = s + t_pred)
    dp_results_mixed[[j]] <- out_all
    
    # --- Evaluate Criteria ---
    eval_df <- merge(valid_at_risk, out_all, by = "id", sort = FALSE)
    P_event <- 1 - eval_df$surv_pred
    
    Crit <- Criteria(
      s = s,
      t = t_pred,
      Survt = eval_df$survtime,
      CR = eval_df$Event,
      P = P_event,
      cause = 1
    )
    
    est_lmm[j, ] <- as.numeric(Crit$Cri[, 1])
    sd_lmm[j, ]  <- as.numeric(Crit$Cri[, 2])
  }
  
  colnames(est_lmm) <- c("est-AUC", "est-BS")
  colnames(sd_lmm)  <- c("sd-AUC", "sd-BS")
  rownames(est_lmm) <- rownames(sd_lmm) <- S
  
  # Compute iAUC and iBS
  res_inc_lmm <- compute_iAUC_iBS(
    s = S,
    auc_mat = cbind(est_lmm[,1], sd_lmm[,1]),
    bs_mat  = cbind(est_lmm[,2], sd_lmm[,2]),
    survtime = surv_valid$survtime,
    death = surv_valid$Event,
    t_pred = t_pred,
    type = "incident",
    estimate_sd = TRUE,
    B = 200
  )
  
  end1 <- Sys.time()
  landmarking_lmm <- difftime(end1, start1, units = "mins")
  
  # Outputs
  est_lmm
  sd_lmm
  res_inc_lmm
  landmarking_lmm
  
  
  ###############################
  
  
  S <- c(0.25, 0.5, 0.75, 1)   # Landmarks
  t_pred <- 0.25                # Prediction horizon
  
  ###########################
  # 2. LOCF Landmarking
  ###########################
  start1 <- Sys.time()
  est_locf <- sd_locf <- matrix(0, length(S), 2)
  dp_results_locf <- vector("list", length(S))
  names(dp_results_locf) <- paste0("LM_", S)
  
  for(j in seq_along(S)) {
    s <- S[j]
    
    # Subjects at risk at landmark
    train_at_risk <- subset(surv_train, survtime >= s)
    valid_at_risk <- subset(surv_valid, survtime >= s)
    
    # Longitudinal data up to landmark
    long_train_s <- subset(long_train, id %in% train_at_risk$id & obstime <= s)
    long_valid_s <- subset(long_valid, id %in% valid_at_risk$id & obstime <= s)
    
    # LOCF summaries
    locf_summary <- function(df, varname) {
      if(nrow(df) == 0) return(data.frame(id=integer(0), Y= numeric(0)))
      tmp <- lapply(split(df, df$id), function(d) d[which.max(d$obstime), c("id", varname)])
      res <- do.call(rbind, tmp)
      res <- as.data.frame(res)
      names(res) <- c("id", "Y")
      res$id <- as.integer(res$id)
      res$Y <- as.numeric(res$Y)
      res
    }
    
    long_locf_train <- locf_summary(long_train_s, "Y1")
    long_locf_valid <- locf_summary(long_valid_s, "Y1")
    
    # Merge with survival
    train_landmark <- merge(train_at_risk, long_locf_train, by = "id", sort = FALSE)
    valid_landmark <- merge(valid_at_risk, long_locf_valid, by = "id", sort = FALSE)
    
    # Skip if not enough data
    if(nrow(train_landmark) < 1 || nrow(valid_at_risk) == 0) {
      out_all <- data.frame(id = valid_at_risk$id,
                            surv_pred = NA_real_,
                            landmark = s,
                            horizon = s + t_pred)
      dp_results_locf[[j]] <- out_all
      next
    }
    
    # Fit Cox model
    cox_fit <- coxph(Surv(survtime - s, Event) ~ Y + w1 + w2, data = train_landmark, ties = "breslow")
    H0_df <- basehaz(cox_fit, centered = FALSE)
    H0_tpred <- approx(x = H0_df$time, y = H0_df$hazard, xout = t_pred, rule = 2)$y
    lp_valid <- predict(cox_fit, newdata = valid_landmark, type = "lp")
    surv_pred <- exp(- H0_tpred * exp(lp_valid))
    
    # Align predictions
    out_all <- data.frame(id = valid_at_risk$id, surv_pred = NA_real_)
    m <- match(valid_landmark$id, out_all$id)
    out_all$surv_pred[m] <- surv_pred
    out_all$landmark <- s
    out_all$horizon <- s + t_pred
    dp_results_locf[[j]] <- out_all
    
    # Compute criteria
    merged <- merge(valid_at_risk, out_all, by = "id", all.x = TRUE)
    merged$surv_pred[is.na(merged$surv_pred)] <- 0
    Crit <- Criteria(s = s, t = t_pred, Survt = merged$survtime, CR = merged$Event,
                     P = 1 - merged$surv_pred, cause = 1)
    est_locf[j, ] <- as.numeric(Crit$Cri[,1])
    sd_locf[j, ]  <- as.numeric(Crit$Cri[,2])
  }
  
  res_inc_locf <- compute_iAUC_iBS(
    s = S,
    auc_mat = cbind(est_locf[,1], sd_locf[,1]),
    bs_mat  = cbind(est_locf[,2], sd_locf[,2]),
    survtime = surv_valid$survtime,
    death = surv_valid$Event,
    t_pred = t_pred,
    type = "incident",
    estimate_sd = TRUE,
    B = 200
  )
  end1 <- Sys.time()
  landmarking_locf <- difftime(end1, start1, units = "mins")
  
  
  
 
  
  
  Results1[[kkk]]=list(iAUC_lmm=res_inc_locf$iAUC,
                       iBS_lmm=res_inc_locf$iBS,
                       est1=est_lmm,
                       landmarking_lmm=landmarking_lmm,
                       
                       iAUC_locf=res_inc_locf$iAUC,
                       iBS_locf=res_inc_locf$iBS,
                       est=est_locf,
                       landmarking_locf=landmarking_locf,
                       
                       
                       
                       iAUC_real = res_incr$iAUC, 
                       iBS_real  = res_incr$iBS,  
                       Real  = estr )
  
  
  
  
  
  
  
  # Outputs
  
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

apply(Time1,2,mean)


apply(Cr0,2,mean)
apply(Cr1,2,mean)
apply(Cr2,2,mean)





