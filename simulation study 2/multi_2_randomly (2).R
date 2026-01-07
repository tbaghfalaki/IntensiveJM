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



setwd("/Users/user/Desktop/Rcodes/PP_TM/New_stan/11th")#
source("iauc_sd.R")

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
D_block <- matrix(c(1,0.5,0.5,1), 2,2)
D <- matrix(0.5,6,6)
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



library(dplyr)

downsample_long_fraction <- function(data, frac = 0.2) {
  data %>%
    group_by(id) %>%
    arrange(obstime, .by_group = TRUE) %>%
    mutate(idx = row_number()) %>%
    {
      # Keep first observation
      first_obs <- filter(., idx == 1)
      # Sample 30% of remaining
      remaining <- filter(., idx > 1)
      sampled <- remaining %>% slice_sample(prop = frac)
      bind_rows(first_obs, sampled)
    } %>%
    arrange(id, obstime) %>%  # preserve order
    select(-idx) %>%
    ungroup()
}

# Example usage
long.data.down30 <- downsample_long_fraction(long.data, frac = 0.2)





long.data=long.data.down30



# ----------------------------
# Train/validation split
# ----------------------------
train_ids <- 1:500
valid_ids <- setdiff(id, train_ids)

long_train <- subset(long.data, id %in% train_ids)
long_valid <- subset(long.data, id %in% valid_ids)

surv_train <- subset(surv.data, id %in% train_ids)
surv_valid <- subset(surv.data, id %in% valid_ids)
########========================================================================
# outcomes and mixed_model() for categorical outcomes.
start1<- Sys.time()

fm1 <- lme(fixed = Y1 ~ x1 + x2 + obstime,
           random = ~ obstime | id,
           data = long_train)

fm2 <- lme(fixed = Y2 ~ x1 + x2 + obstime,
           random = ~ obstime | id,
           data = long_train)


fm3 <- lme(fixed = Y3 ~ x1 + x2 + obstime,
           random = ~ obstime | id,
           data = long_train)



# [2] Save all the fitted mixed-effects models in a list.
Mixed <- list(fm1, fm2, fm3)

# [3] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(survtime, Event) ~ w1 + w2,
               data = surv_train)
# [4] The joint model is fitted using a call to jm() i.e.,
joint_model_fit_1 <- jm(fCox1, Mixed, time_var = "obstime",
                        n_chains = 1L, n_iter = 1000L, n_burnin = 500L)


# ----------------------------
# Dynamic predictions for id=501
# ----------------------------
library(dplyr)

# ----------------------------
# Dynamic predictions for id=501
# ----------------------------
library(dplyr)

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


library(dplyr)

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
  P = as.numeric(unlist(preds_100[,2])),
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


Real=c(alpha,gamma_w)

# Extract MCMC samples
mcmc_samples <- joint_model_fit_1$mcmc

ci=rbind(summary(mcmc_samples$alphas)[2]$quantiles[,c(1,5)],
         summary(mcmc_samples$gammas)[2]$quantiles[,c(1,5)])



Estijm[kkk,]=c(summary(mcmc_samples$alphas)[1]$statistics[,1],summary(mcmc_samples$gammas)[1]$statistics[,1])


for(ii in 1:length(Real)){ 
  if((Real[ii]>ci[ii,1]) & (Real[ii]<ci[ii,2]))(CovRjm[kkk,ii]=1)
}







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
apply(Cr1,2,mean)
mean(Time1);sd(Time1);

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
