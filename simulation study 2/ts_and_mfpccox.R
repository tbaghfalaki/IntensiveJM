
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
NN=3

 Results=list()

Esti=CovR=matrix(0,NN,5)

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


# ----------------------------
# MFPCCox for three markers
# ----------------------------
start1 <- Sys.time()
source("C:/Users/p80744tb/Desktop/intensive/functions.R")

t <- sort(unique(long_train$obstime))
train_ids <- unique(long_train$id)
valid_ids <- unique(long_valid$id)

# Initialize matrices for FPCA
Ytrain_mat1 <- Ytrain_mat2 <- Ytrain_mat3 <- matrix(NA, nrow = length(train_ids), ncol = length(t))

for (i in seq_along(train_ids)) {
  tmp <- subset(long_train, id == train_ids[i] & obstime %in% t)
  if(nrow(tmp) > 0) {
    tmp <- tmp %>% group_by(obstime) %>% summarize(Y1 = mean(Y1, na.rm=TRUE),
                                                   Y2 = mean(Y2, na.rm=TRUE),
                                                   Y3 = mean(Y3, na.rm=TRUE))
    idx <- match(tmp$obstime, t)
    Ytrain_mat1[i, idx] <- tmp$Y1
    Ytrain_mat2[i, idx] <- tmp$Y2
    Ytrain_mat3[i, idx] <- tmp$Y3
  }
}

# FPCA for each marker
fpca1 <- fpca.sc(Ytrain_mat1, argvals = t, pve = 0.99)
fpca2 <- fpca.sc(Ytrain_mat2, argvals = t, pve = 0.99)
fpca3 <- fpca.sc(Ytrain_mat3, argvals = t, pve = 0.99)

K1 <- ncol(fpca1$scores)
K2 <- ncol(fpca2$scores)
K3 <- ncol(fpca3$scores)

# Combine FPCA scores into survival training set
surv_train <- cbind(surv_train,
                    setNames(as.data.frame(fpca1$scores), paste0("rho1_", 1:K1)),
                    setNames(as.data.frame(fpca2$scores), paste0("rho2_", 1:K2)),
                    setNames(as.data.frame(fpca3$scores), paste0("rho3_", 1:K3)))

# Fit Cox model
cox_model <- coxph(Surv(survtime, Event) ~ w1 + w2 + ., data = surv_train, x = TRUE)
summary(cox_model)

# Project FPCA scores for validation
project_scores <- function(fpca, newY) {
  Y_centered <- newY - fpca$mu
  scores <- colSums(outer(Y_centered, fpca$efunctions) * (1/length(newY)), na.rm = TRUE)
  scores[1:ncol(fpca$efunctions)]
}

S <- landmarks <- c(0.25, 0.5, 0.75, 1)
t_pred <- horizons <- 0.25
results <- list()

for (t0 in landmarks) {
  for (dt in horizons) {
    pred_list <- lapply(valid_ids, function(iid) {
      this_long <- subset(long_valid, id == iid & obstime <= t0)
      this_surv <- subset(surv_valid, id == iid)
      if (nrow(this_surv) == 0) return(NULL)
      
      Y1test <- Y2test <- Y3test <- rep(NA, length(t))
      
      if (nrow(this_long) > 0) {
        tmp <- this_long %>% group_by(obstime) %>%
          summarize(Y1 = mean(Y1, na.rm=TRUE),
                    Y2 = mean(Y2, na.rm=TRUE),
                    Y3 = mean(Y3, na.rm=TRUE))
        idx <- match(tmp$obstime, t)
        Y1test[idx] <- tmp$Y1
        Y2test[idx] <- tmp$Y2
        Y3test[idx] <- tmp$Y3
      }
      
      sc1 <- project_scores(fpca1, Y1test)
      sc2 <- project_scores(fpca2, Y2test)
      sc3 <- project_scores(fpca3, Y3test)
      
      sc_full <- c(sc1, sc2, sc3)
      sc_df <- as.data.frame(as.list(sc_full))
      colnames(sc_df) <- c(paste0("rho1_", 1:length(sc1)),
                           paste0("rho2_", 1:length(sc2)),
                           paste0("rho3_", 1:length(sc3)))
      
      newdata <- cbind(this_surv, sc_df)
      sf <- survfit(cox_model, newdata = newdata)
      p <- tryCatch(summary(sf, times = t0 + dt)$surv, error = function(e) NA)
      
      data.frame(id = iid, landmark = t0, horizon = dt, pred = p)
    })
    results[[paste(t0, dt)]] <- do.call(rbind, pred_list)
  }
}

pred_all <- do.call(rbind, results)
pred_by_landmark <- split(pred_all, pred_all$landmark)

# Compute iAUC and iBS
est3 <- sd3 <- matrix(0, length(S), 2)

for (i in seq_along(S)) {
  Crit <- Criteria(
    s = S[i],
    t = t_pred,
    Survt = surv_valid$survtime,
    CR = surv_valid$Event,
    P =  1-pred_by_landmark[[i]][, "pred"],
    cause = 1
  )
  
  est3[i, ] <- Crit$Cri[, 1]
  sd3[i, ]  <- Crit$Cri[, 2]
}

res_inc3 <- compute_iAUC_iBS(
  s = S,
  auc_mat = cbind(est3[,1], sd3[,1]),
  bs_mat  = cbind(est3[,2], sd3[,2]),
  survtime = surv_valid$survtime,
  death    = surv_valid$Event,
  t_pred   = t_pred,
  type     = "incident",
  estimate_sd = TRUE,
  B = 200
)

c(res_inc3$iAUC, res_inc3$sd_iAUC)
c(res_inc3$iBS, res_inc3$sd_iBS)
est3

end1 <- Sys.time()
time_MFPCCOX <- difftime(end1, start1, units ="mins")
time_MFPCCOX


#===============================================================================
#============================ Two-Stage Joint Model (3 markers) ==============
#===============================================================================

library(lme4)
library(survival)
library(dplyr)

start1 <- Sys.time()

#================= Stage 1: Fit LMMs for each marker ==========================
# Marker 1
lme1 <- lmer(Y1 ~ obstime + x1 + x2 + (1 + obstime | id), data = long_train)
fixef1 <- fixef(lme1)
ranef1 <- ranef(lme1)$id
colnames(ranef1) <- c("b0_1", "b1_1")
ranef1$id <- as.numeric(rownames(ranef1))

# Marker 2
lme2 <- lmer(Y2 ~ obstime + x1 + x2 + (1 + obstime | id), data = long_train)
fixef2 <- fixef(lme2)
ranef2 <- ranef(lme2)$id
colnames(ranef2) <- c("b0_2", "b1_2")
ranef2$id <- as.numeric(rownames(ranef2))

# Marker 3
lme3 <- lmer(Y3 ~ obstime + x1 + x2 + (1 + obstime | id), data = long_train)
fixef3 <- fixef(lme3)
ranef3 <- ranef(lme3)$id
colnames(ranef3) <- c("b0_3", "b1_3")
ranef3$id <- as.numeric(rownames(ranef3))

#================= Stage 2: Construct predicted longitudinal values ===========

long_train_td <- long_train %>%
  left_join(ranef1, by = "id") %>%
  mutate(Y1_pred = b0_1 + b1_1 * obstime + fixef1["x1"] * x1 + fixef1["x2"] * x2) %>%
  left_join(ranef2, by = "id") %>%
  mutate(Y2_pred = b0_2 + b1_2 * obstime + fixef2["x1"] * x1 + fixef2["x2"] * x2) %>%
  left_join(ranef3, by = "id") %>%
  mutate(Y3_pred = b0_3 + b1_3 * obstime + fixef3["x1"] * x1 + fixef3["x2"] * x2)

#================= Prepare survival data for Cox model ========================

surv_train_td <- surv_train %>%
  rename(Id = id, deathTimes = survtime, Event = Event)

# Baseline tmerge object
td_base <- tmerge(
  data1 = surv_train_td, data2 = surv_train_td,
  id = Id, endpt = event(deathTimes, Event)
)

# Rename for time-dependent covariates
long_train_td <- long_train_td %>% rename(Id = id, time = obstime)

long_data_td <- tmerge(
  td_base, long_train_td,
  id = Id,
  lY1 = tdc(time, Y1_pred),
  lY2 = tdc(time, Y2_pred),
  lY3 = tdc(time, Y3_pred)
)

#================= Fit extended Cox model ====================================

cox_mts <- coxph(
  Surv(tstart, tstop, endpt) ~ lY1 + lY2 + lY3 + w1 + w2,
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

#================= Dynamic predictions at landmarks ==========================

S <- c(0.25, 0.5, 0.75, 1)  # landmark times
t_pred <- 0.25              # prediction horizon
dp_results <- list()
est4 <- sd4 <- matrix(0, length(S), 2)

# Covariance matrices and residual variance
Sigma_b1 <- as.matrix(VarCorr(lme1)$id)[1:2, 1:2]
Sigma_b2 <- as.matrix(VarCorr(lme2)$id)[1:2, 1:2]
Sigma_b3 <- as.matrix(VarCorr(lme3)$id)[1:2, 1:2]
sigma2_1 <- sigma(lme1)^2
sigma2_2 <- sigma(lme2)^2
sigma2_3 <- sigma(lme3)^2

r <- 0
for (s in S) {
  r <- r + 1
  valid_at_risk <- subset(surv_valid, survtime >= s)
  valid_ids <- valid_at_risk$id
  
  Y1_pred_valid <- Y2_pred_valid <- Y3_pred_valid <- numeric(length(valid_ids))
  
  for (i in seq_along(valid_ids)) {
    idv <- valid_ids[i]
    
    # Subject's longitudinal history
    subj_long <- subset(long_valid, id == idv & obstime <= s)
    
    # Fixed covariates
    x1i <- if (nrow(subj_long) > 0) subj_long$x1[1] else valid_at_risk$x1[valid_at_risk$id == idv]
    x2i <- if (nrow(subj_long) > 0) subj_long$x2[1] else valid_at_risk$x2[valid_at_risk$id == idv]
    
    # Design matrices
    X_landmark <- model.matrix(~ obstime + x1 + x2, data = data.frame(obstime = s, x1 = x1i, x2 = x2i))
    Z_landmark <- c(1, s)
    
    # Marker 1
    if (nrow(subj_long) == 0) {
      Y1_pred_valid[i] <- as.numeric(X_landmark %*% fixef1)
    } else {
      X_s <- model.matrix(~ obstime + x1 + x2, data = subj_long)
      Z_s <- cbind(1, subj_long$obstime)
      Y_s <- subj_long$Y1
      V <- Z_s %*% Sigma_b1 %*% t(Z_s) + sigma2_1 * diag(nrow(Z_s))
      b_hat <- Sigma_b1 %*% t(Z_s) %*% solve(V, Y_s - X_s %*% fixef1)
      Y1_pred_valid[i] <- as.numeric(X_landmark %*% fixef1 + Z_landmark %*% b_hat)
    }
    
    # Marker 2
    if (nrow(subj_long) == 0) {
      Y2_pred_valid[i] <- as.numeric(X_landmark %*% fixef2)
    } else {
      X_s <- model.matrix(~ obstime + x1 + x2, data = subj_long)
      Z_s <- cbind(1, subj_long$obstime)
      Y_s <- subj_long$Y2
      V <- Z_s %*% Sigma_b2 %*% t(Z_s) + sigma2_2 * diag(nrow(Z_s))
      b_hat <- Sigma_b2 %*% t(Z_s) %*% solve(V, Y_s - X_s %*% fixef2)
      Y2_pred_valid[i] <- as.numeric(X_landmark %*% fixef2 + Z_landmark %*% b_hat)
    }
    
    # Marker 3
    if (nrow(subj_long) == 0) {
      Y3_pred_valid[i] <- as.numeric(X_landmark %*% fixef3)
    } else {
      X_s <- model.matrix(~ obstime + x1 + x2, data = subj_long)
      Z_s <- cbind(1, subj_long$obstime)
      Y_s <- subj_long$Y3
      V <- Z_s %*% Sigma_b3 %*% t(Z_s) + sigma2_3 * diag(nrow(Z_s))
      b_hat <- Sigma_b3 %*% t(Z_s) %*% solve(V, Y_s - X_s %*% fixef3)
      Y3_pred_valid[i] <- as.numeric(X_landmark %*% fixef3 + Z_landmark %*% b_hat)
    }
  }
  
  # Merge predictions with covariates
  valid_landmark <- merge(
    valid_at_risk,
    data.frame(id = valid_ids, lY1 = Y1_pred_valid, lY2 = Y2_pred_valid, lY3 = Y3_pred_valid),
    by = "id", sort = FALSE
  )
  
  # Baseline hazard
  H0_df <- basehaz(cox_mts, centered = FALSE)
  H0_tpred <- approx(H0_df$time, H0_df$hazard, xout = s + t_pred, rule = 2)$y
  
  # Linear predictor and survival
  lp_valid <- predict(cox_mts, newdata = valid_landmark, type = "lp")
  surv_pred <- exp(-H0_tpred * exp(lp_valid))
  
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
    P = 1 - dp_results[[paste0("S_", s)]][,2],
    cause = 1
  )
  
  est4[r, ] <- as.numeric(Crit$Cri[, 1])
  sd4[r, ]  <- as.numeric(Crit$Cri[, 2])
}

#================= Compute iAUC and iBS ======================================

res_inc4 <- compute_iAUC_iBS(
  s = S,
  auc_mat = cbind(est4[,1], sd4[,1]),
  bs_mat  = cbind(est4[,2], sd4[,2]), 
  survtime = surv_valid$survtime,
  death    = surv_valid$Event,
  t_pred   = t_pred,
  type     = "incident",
  estimate_sd = TRUE, 
  B = 200
)

c(res_inc4$iAUC, res_inc4$sd_iAUC)
c(res_inc4$iBS, res_inc4$sd_iBS)
est4

end1 <- Sys.time()
time_TS <- difftime(end1, start1, units = "mins")
time_TS




Results[[kkk]]=list(time_TS=time_TS,est4=est4,AUCi_ts=res_inc4$iAUC,BSi_ts=res_inc4$iBS,
                    Esti=Esti,
                    CovR=CovR,
                    est3=est3,AUCi_mfpcox=res_inc3$iAUC,BSi_mfpcox=res_inc3$iBS,
                    time_MFPCCOX=time_MFPCCOX)

print(rep(kkk,10))

}





Time1=matrix(0,NN,2)
Cr1=Cr2=matrix(0,NN,10)
Est=CR=matrix(0,NN,2)
for(kkk in c(1:NN)[-35]){ 
  
  Time1[kkk,1]=Results[[kkk]]$time_TS
  Time1[kkk,2]=Results[[kkk]]$time_MFPCCOX
  
  Cr1[kkk,]=c(as.numeric(Results[[kkk]]$est4),c(Results[[kkk]]$AUCi_ts,Results[[kkk]]$BSi_ts))
  Cr2[kkk,]=c(as.numeric(Results[[kkk]]$est3),c(Results[[kkk]]$AUCi_mfpcox,Results[[kkk]]$BSi_mfpcox))
  
  
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
