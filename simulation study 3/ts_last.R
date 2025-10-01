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
      Y1[i,j] <- Beta[1]+Beta[2]*x1[i]+Beta[3]*x2[i]+sin(pi*t[j])+cos(pi*t[j])  +
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
#============================ TS  ==============================================
#===============================================================================
start1 <- Sys.time()

############ Fast Two-Stage GAMM Landmarking ############

library(mgcv)
library(survival)
library(dplyr)
library(purrr)

start_time <- Sys.time()

# ---- user settings ----
S <- c(0.25, 0.5, 0.75, 1)    # landmark times
t_pred <- 0.25                # prediction horizon
min_obs_for_blup <- 2         # minimum history for BLUP

# ---- containers ----
dp_results_ts <- vector("list", length(S))
names(dp_results_ts) <- paste0("TS_Mixed_", S)
est_ts <- sd_ts <- matrix(NA_real_, length(S), 2)

# ---- Stage 1: fit GAMM on training longitudinal data ----
gamm_fit <- gamm(Y1 ~ s(obstime, k = 5) + x1 + x2,
                 random = list(id = ~ 1 + obstime),
                 data = long_train)

Sigma_b <- as.matrix(gamm_fit$lme$modelStruct$reStruct$id,2,2)  # 'id' grouping
sigma2  <- gamm_fit$lme$sigma^2                     # residual variance
beta_hat <- coef(gamm_fit$gam)[c("(Intercept)", "x1", "x2")]

# ---- Precompute lpmatrix for all longitudinal rows ----
B_lpmatrix_train <- predict(gamm_fit$gam, newdata = long_train, type = "lpmatrix")
B_lpmatrix_valid <- predict(gamm_fit$gam, newdata = long_valid, type = "lpmatrix")
smooth_idx <- grep("s\\(obstime\\)", colnames(B_lpmatrix_train))

# ---- Helper: compute Ypred using BLUP and smooth ----
compute_Ypred <- function(subj_hist, t_val, x1_val, x2_val, B_lpmatrix_subj, smooth_idx) {
  # smooth at t_val
  B_row <- predict(gamm_fit$gam, newdata = data.frame(obstime = t_val, x1 = x1_val, x2 = x2_val), type = "lpmatrix")
  f_t <- sum(B_row[, smooth_idx] * coef(gamm_fit$gam)[smooth_idx])
  
  fixed_part <- sum(c(1, x1_val, x2_val) * beta_hat)
  
  if (nrow(subj_hist) < min_obs_for_blup) return(f_t + fixed_part)
  
  # BLUP computation
  Z_i <- cbind(1, subj_hist$obstime)
  X_i <- cbind(1, subj_hist$x1, subj_hist$x2)
  Y_i <- subj_hist$Y1
  f_i <- rowSums(B_lpmatrix_subj[, smooth_idx, drop=FALSE] * coef(gamm_fit$gam)[smooth_idx])
  resid <- Y_i - X_i %*% beta_hat - f_i
  V_i <- Z_i %*% Sigma_b %*% t(Z_i) + sigma2 * diag(nrow(Z_i))
  b_hat <- Sigma_b %*% t(Z_i) %*% solve(V_i, resid)
  
  Z_s <- c(1, t_val)
  f_t + fixed_part + sum(Z_s * b_hat)
}

# ---- Landmark loop ----
for (j in seq_along(S)) {
  s <- S[j]
  
  # Subjects at risk
  train_at_risk <- subset(surv_train, survtime >= s)
  valid_at_risk <- subset(surv_valid, survtime >= s)
  
  long_train_s <- subset(long_train, id %in% train_at_risk$id & obstime <= s)
  long_valid_s <- subset(long_valid, id %in% valid_at_risk$id & obstime <= s)
  
  if (nrow(long_train_s) < 1 || nrow(valid_at_risk) == 0) {
    dp_results_ts[[j]] <- data.frame(
      id = valid_at_risk$id, surv_pred = NA_real_,
      landmark = s, horizon = s + t_pred, Ypred = NA_real_
    )
    next
  }
  
  # --- TRAIN: compute Ypred at s using BLUP ---
  Ypred_train <- sapply(train_at_risk$id, function(idv) {
    subj_hist <- subset(long_train_s, id == idv)
    x1v <- if(nrow(subj_hist) > 0) subj_hist$x1[1] else train_at_risk$x1[train_at_risk$id==idv]
    x2v <- if(nrow(subj_hist) > 0) subj_hist$x2[1] else train_at_risk$x2[train_at_risk$id==idv]
    B_subj <- B_lpmatrix_train[long_train$id %in% subj_hist$id & long_train$obstime <= s, , drop=FALSE]
    compute_Ypred(subj_hist, s, x1v, x2v, B_subj, smooth_idx)
  })
  train_at_risk$Ypred <- Ypred_train
  
  # --- Fit Cox at landmark ---
  cox_fit <- coxph(Surv(survtime - s, Event) ~ Ypred + w1 + w2,
                   data = train_at_risk, ties = "breslow")
  
  
  
  Esti[kkk,]=summary(cox_fit)$coefficients[, "coef"]
  
  # 2. Extract standard errors
  se <- summary(cox_fit)$coefficients[, "se(coef)"]
  Real=c(alpha,gamma_w)
  # 3. 95% confidence intervals
  ci <- confint(cox_fit)  # default is 95%
  ci
  for(ii in 1:length(Real)){ 
    if((Real[ii]>ci[ii,1]) & (Real[ii]<ci[ii,2]))(CovR[kkk,ii]=1)
  }
  
  
  
  
  
  # --- VALIDATION: compute Ypred at s ---
  Ypred_valid <- sapply(valid_at_risk$id, function(idv) {
    subj_hist <- subset(long_valid_s, id == idv)
    x1v <- if(nrow(subj_hist) > 0) subj_hist$x1[1] else valid_at_risk$x1[valid_at_risk$id==idv]
    x2v <- if(nrow(subj_hist) > 0) subj_hist$x2[1] else valid_at_risk$x2[valid_at_risk$id==idv]
    B_subj <- B_lpmatrix_valid[long_valid$id %in% subj_hist$id & long_valid$obstime <= s, , drop=FALSE]
    compute_Ypred(subj_hist, s, x1v, x2v, B_subj, smooth_idx)
  })
  valid_at_risk$Ypred <- Ypred_valid
  
  # --- Predicted survival over horizon t_pred ---
  bh <- basehaz(cox_fit, centered = FALSE)
  H0_tpred <- approx(bh$time, bh$hazard, xout = t_pred, rule = 2)$y
  lp_valid <- predict(cox_fit, newdata = valid_at_risk, type = "lp")
  surv_pred <- exp(-H0_tpred * exp(lp_valid))
  P_event <- 1 - surv_pred
  
  # --- Criteria ---
  Crit <- Criteria(s = s, t = t_pred,
                   Survt = valid_at_risk$survtime,
                   CR = valid_at_risk$Event,
                   P = P_event,
                   cause = 1)
  est_ts[j, ] <- Crit$Cri[,1]
  sd_ts[j, ]  <- Crit$Cri[,2]
  
  # store results
  dp_results_ts[[j]] <- data.frame(
    id = valid_at_risk$id,
    surv_pred = P_event,
    landmark = s,
    horizon = s + t_pred,
    Ypred = Ypred_valid
  )
  
  message(sprintf("landmark=%.2f: n_train=%d, n_valid=%d, mean(Ypred_valid)=%.4g, mean(P_event)=%.4g",
                  s, nrow(train_at_risk), nrow(valid_at_risk),
                  mean(Ypred_valid, na.rm=TRUE), mean(P_event, na.rm=TRUE)))
}

end_time <- Sys.time()
message("Elapsed: ", round(difftime(end_time, start_time, units="secs"),2))







res_inc4 <- compute_iAUC_iBS(
  s = S,
  auc_mat=cbind(est_ts[,1],sd_ts[,1]), bs_mat=cbind(est_ts[,2],sd_ts[,2]), 
  survtime=surv_valid$survtime, death=surv_valid$Event, t_pred=t_pred,
  type = "incident",     # or "cumulative"
  estimate_sd = TRUE, 
  B = 200
)

 

end1 <- Sys.time()
time_TS=difftime(end1,start1,units ="mins")

Results[[kkk]]=list(time_TS=time_TS,est4=est_ts,AUCi_ts=res_inc4$iAUC,BSi_ts=res_inc4$iBS,
                    Esti=Esti,
                    CovR=CovR)

print(rep(kkk,10))

}













Time1=matrix(0,NN,2)
Cr1=Cr2=matrix(0,NN,10)
Est=CR=matrix(0,NN,2)
for(kkk in c(1:NN)){ 
  
  Time1[kkk,1]=Results[[kkk]]$time_TS
  Time1[kkk,2]=Results[[kkk]]$time_MFPCCOX
  
  Cr1[kkk,]=c(as.numeric(Results[[kkk]]$est4),c(Results[[kkk]]$AUCi_ts,Results[[kkk]]$BSi_ts))
  Cr2[kkk,]=c(as.numeric(Results[[kkk]]$est3),c(Results[[kkk]]$AUCi_mfpcox,Results[[kkk]]$BSi_mfpcox))
  
  
}

apply(Results[[kkk]]$Esti,2,mean)
apply(Results[[kkk]]$CovR,2,mean)



apply(Cr1,2,mean)
apply(Cr2,2,mean)


apply(Time1,2,mean)



cbind(apply(Time1,2,mean),apply(Time1,2,sd))#time
cbind(apply(Cr1,2,mean),apply(Cr1,2,sd))#ts
cbind(apply(Cr2,2,mean),apply(Cr2,2,sd))#mfpccox






MSE <- matrix(0, NN, length(Real))
for(i in c(1:NN)){ 
  MSE[i, ] <- (Results[[kkk]]$Esti[i,] - Real)^2
}

# Compute summary metrics
Res <- cbind(
  Real,
  apply(Results[[kkk]]$Esti,2,mean),
  apply(Results[[kkk]]$Esti,2,sd),
  (apply(Results[[kkk]]$Esti,2,mean) - Real) / Real, # Relative bias
  sqrt(apply(MSE, 2, mean)),     # RMSE
  apply(Results[[kkk]]$CovR,2,mean)
)
colnames(Res) <- c("Real", "Est", "sd", "Rbias", "RMSE", "CR")

xtable::xtable(Res,digits=3)

1 & -0.500 & -0.568 & 0.140 & 0.136 & 0.155 & 0.990 \\ 
2 & 0.200 & 0.203 & 0.202 & 0.013 & 0.201 & 1.000 \\ 
3 & -0.200 & -0.205 & 0.375 & 0.023 & 0.373 & 1.000 \\ 

cbind(apply(Time1,2,mean),apply(Time1,2,sd))#time
[,1]       [,2]
[1,] 1.599588 0.13694957
[2,] 0.824840 0.08284549
> cbind(apply(Cr1,2,mean),apply(Cr1,2,sd))#ts
[,1]       [,2]
[1,] 0.6991214 0.04629667
[2,] 0.7103508 0.04336872
[3,] 0.7236069 0.04919513
[4,] 0.7198950 0.06586897
[5,] 0.1114756 0.01563058
[6,] 0.1413372 0.01728922
[7,] 0.1847241 0.02344440
[8,] 0.2335963 0.04454819
[9,] 0.7098783 0.02526949
[10,] 0.1487989 0.01053831
> cbind(apply(Cr2,2,mean),apply(Cr2,2,sd))#mfpccox
[,1]       [,2]
[1,] 0.5640021 0.11224492
[2,] 0.5607214 0.14879507
[3,] 0.5679936 0.17770629
[4,] 0.5582636 0.19205070
[5,] 0.1403033 0.02470205
[6,] 0.2306634 0.05926054
[7,] 0.3635401 0.10220864
[8,] 0.4716685 0.12122130
[9,] 0.5635038 0.13768960
[10,] 0.2478829 0.05539092
