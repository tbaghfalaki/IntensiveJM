

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

#setwd("/Users/user/Desktop/Rcodes/PP_TM/New_stan/11th")#
setwd("C:\\Users\\p80744tb\\Desktop\\PP_TM\\PP_TM\\New_stan\\11th")

source("iauc_sd.R")
# ----------------------------
# Simulation parameters
# ----------------------------

NN=100
Esti=CovR=matrix(0,NN,3)
Results=list()


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
        Y1[i,j] <- Beta[1]+Beta[2]*x1[i]+Beta[3]*x2[i]+Beta[4]*sin(pi*t[j])+cos(pi*t[j])  +
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
  #============================ MFPCCOX ===========================================
  #===============================================================================

#===============================================================================
#============================ MFPCCOX (Corrected) ===========================================
#===============================================================================
start1 <- Sys.time()
source("C:/Users/p80744tb/Desktop/intensive/functions.R")

################################################################################
# --- FPCA on training data ---
################################################################################
Ytrain_mat <- matrix(NA, nrow = length(train_ids), ncol = length(t))
for (i in seq_along(train_ids)) {
  tmp <- subset(long_train, id == train_ids[i])
  idx <- sapply(tmp$obstime, function(tt) which.min(abs(t - tt))) # safer matching
  Ytrain_mat[i, idx] <- tmp$Y1
}

ufpca <- fpca.sc(Ytrain_mat, argvals = t, pve = 0.99, var = TRUE)
K <- ncol(ufpca$scores)

sc_train <- as.data.frame(ufpca$scores)
colnames(sc_train) <- paste0("rho", 1:K)
surv_train_model <- surv_train %>% select(-id)  # remove id from Cox model
surv_train_model <- cbind(surv_train_model, sc_train)

CoxFit <- coxph(Surv(survtime, Event) ~ ., data = surv_train_model, x = TRUE)

################################################################################
# --- Prediction for validation set ---
################################################################################
S <- landmarks <- c(0.25, 0.5, 0.75, 1)
t_pred <- horizons <- 0.25
delta_t <- t[2] - t[1]  # time step for projection

project_scores <- function(fpca, newY) {
  # newY: vector of length t, may contain NAs
  idx <- !is.na(newY)
  if (sum(idx) == 0) return(rep(0, ncol(fpca$efunctions)))
  
  Y_centered <- newY[idx] - fpca$mu[idx]
  scores <- t(Y_centered) %*% fpca$efunctions[idx, , drop=FALSE] * delta_t
  as.numeric(scores)
}

results <- list()

for (t0 in landmarks) {
  for (dt in horizons) {
    pred_list <- lapply(valid_ids, function(iid) {
      this_long <- subset(long_valid, id == iid & obstime <= t0)
      this_surv <- subset(surv_valid, id == iid)
      if (nrow(this_surv) == 0) return(NULL)
      
      Ytest <- rep(NA, length(t))
      idx <- sapply(this_long$obstime, function(tt) which.min(abs(t - tt)))
      Ytest[idx] <- this_long$Y1
      
      sc <- project_scores(ufpca, Ytest)
      sc_full <- rep(0, K)
      sc_full[1:length(sc)] <- sc
      sc_df <- as.data.frame(as.list(sc_full))
      colnames(sc_df) <- paste0("rho", 1:K)
      
      newdata <- cbind(this_surv %>% select(-id), sc_df)
      sf <- survfit(CoxFit, newdata = newdata)
      p <- tryCatch(summary(sf, times = t0 + dt)$surv, error = function(e) NA)
      
      data.frame(id = iid, landmark = t0, horizon = dt, pred = p)
    })
    results[[paste(t0, dt)]] <- do.call(rbind, pred_list)
  }
}

pred_all <- do.call(rbind, results)

# Optional: split by landmark with correct order
pred_by_landmark <- split(pred_all, factor(pred_all$landmark, levels = S))

# iAUC / iBS computation
est3 <- sd3 <- matrix(0, length(S), 2)

for (i in seq_along(S)) {
  Crit <- Criteria(
    s = S[i],
    t = t_pred,
    Survt = surv_valid$survtime,
    CR = surv_valid$Event,
    P = 1 - pred_by_landmark[[i]]$pred,
    cause = 1
  )
  est3[i, ] <- Crit$Cri[, 1]
  sd3[i, ]  <- Crit$Cri[, 2]
}

res_inc3 <- compute_iAUC_iBS(
  s = S,
  auc_mat = cbind(est3[, 1], sd3[, 1]),
  bs_mat = cbind(est3[, 2], sd3[, 2]),
  survtime = surv_valid$survtime,
  death = surv_valid$Event,
  t_pred = t_pred,
  type = "incident",
  estimate_sd = TRUE,
  B = 200
)

c(res_inc3$iAUC, res_inc3$sd_iAUC)
c(res_inc3$iBS, res_inc3$sd_iBS)
est3

end1 <- Sys.time()
time_MFPCCOX <- difftime(end1, start1, units = "mins")
time_MFPCCOX



res_inc3 <- compute_iAUC_iBS(
  s = S,
  auc_mat=cbind(est3[,1],sd3[,1]), bs_mat=cbind(est3[,2],sd3[,2]), 
  survtime=surv_valid$survtime, death=surv_valid$Event, t_pred=t_pred,
  type = "incident",     # or "cumulative"
  estimate_sd = TRUE, 
  B = 200
)


c(res_inc3$iAUC,res_inc3$sd_iAUC)
c(res_inc3$iBS,res_inc3$sd_iBS)
est3
end1 <- Sys.time()
time_MFPCCOX=difftime(end1,start1,units ="mins")



end1 <- Sys.time()
time_TS=difftime(end1,start1,units ="mins")

Results[[kkk]]=list(
  
                    est3=est3,AUCi_mfpcox=res_inc3$iAUC,BSi_mfpcox=res_inc3$iBS,
                    time_MFPCCOX=time_MFPCCOX)

print(rep(kkk,10))

}





Time1=matrix(0,NN,2)
Cr1=Cr2=matrix(0,NN,10)
Est=CR=matrix(0,NN,2)
for(kkk in c(1:NN)){ 
  
  Time1[kkk,2]=Results[[kkk]]$time_MFPCCOX
  
  Cr2[kkk,]=c(as.numeric(Results[[kkk]]$est3),c(Results[[kkk]]$AUCi_mfpcox,Results[[kkk]]$BSi_mfpcox))
  
  
}



cbind(apply(Time1,2,mean),apply(Time1,2,sd))#time
cbind(apply(Cr2,2,mean),apply(Cr2,2,sd))#mfpccox

