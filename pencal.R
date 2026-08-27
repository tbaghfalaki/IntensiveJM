rm(list = ls())
library(MASS)
library(mvtnorm)
library(lme4)
library(survival)
library(DPCri)
library(pencal)

NN <- 100
nsample <- 1000
followup <- 2
S <- c(0.25, 0.50, 0.75, 1.00)
t_pred <- 0.25
# =============================================================================

AUC_pen  <- matrix(NA_real_, NN, length(S))
BS_pen   <- matrix(NA_real_, NN, length(S))

# Computational time
TIME_pen  <- numeric(NN)

# Optional: store all predictions
Results1 <- vector("list", NN)

 
for (kkk in seq_len(NN)) {
  
  set.seed(kkk)
  
  start_total <- Sys.time()
  
  # ===========================================================================
  # Simulation parameters
  # ===========================================================================
  
  t <- seq(0, followup, by = 0.005)
  
  # ---------------------------------------------------------------------------
  # Association parameters
  # ---------------------------------------------------------------------------
  
  alpha <- c(-0.5, 0.3, -0.2)
  
  # ---------------------------------------------------------------------------
  # Fixed effects
  # ---------------------------------------------------------------------------
  
  Beta1 <- c(-0.5, -0.5, 0.5, 0.5)
  Beta2 <- c( 0.2,  0.3, -0.3, 0.4)
  Beta3 <- c( 1.0, -0.2,  0.1, 0.6)
  
  sigma <- 1
  
  gamma_w <- c(0.2, -0.2)
  
  # ===========================================================================
  # Covariates
  # ===========================================================================
  
  id <- 1:nsample
  
  x1 <- rnorm(nsample)
  x2 <- rbinom(nsample, 1, 0.5)
  
  w1 <- rnorm(nsample)
  w2 <- rbinom(nsample, 1, 0.5)
  
  death <- rep(0L, nsample)
  time <- rep(NA_real_, nsample)
  
  # ===========================================================================
  # Random effects
  # ===========================================================================
  
  D_block <- matrix(
    c(
      1,   0.6,
      0.6, 1
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
  
  D <- matrix(0.5, nrow = 6, ncol = 6)
  
  D[1:2, 1:2] <- D_block
  D[3:4, 3:4] <- D_block
  D[5:6, 5:6] <- D_block
  
  u <- rmvnorm(
    n = nsample,
    mean = rep(0, 6),
    sigma = D
  )
  
  # ===========================================================================
  # Hazard function
  # ===========================================================================
  
  haz_fun_i <- function(s, x1i, x2i, w1i, w2i, ui) {
    
    mu1 <-
      Beta1[1] +
      Beta1[2] * x1i +
      Beta1[3] * x2i +
      Beta1[4] * s +
      ui[1] +
      ui[2] * s
    
    mu2 <-
      Beta2[1] +
      Beta2[2] * x1i +
      Beta2[3] * x2i +
      Beta2[4] * s +
      ui[3] +
      ui[4] * s
    
    mu3 <-
      Beta3[1] +
      Beta3[2] * x1i +
      Beta3[3] * x2i +
      Beta3[4] * s +
      ui[5] +
      ui[6] * s
    
    eta <-
      gamma_w[1] * w1i +
      gamma_w[2] * w2i +
      alpha[1] * mu1 +
      alpha[2] * mu2 +
      alpha[3] * mu3
    
    exp(eta)
  }
  
  # ===========================================================================
  # Cumulative hazard
  # ===========================================================================
  
  cumhaz_i <- function(T, x1i, x2i, w1i, w2i, ui) {
    
    if (T <= 0) {
      return(0)
    }
    
    integrate(
      function(s) {
        haz_fun_i(
          s = s,
          x1i = x1i,
          x2i = x2i,
          w1i = w1i,
          w2i = w2i,
          ui = ui
        )
      },
      lower = 0,
      upper = T,
      rel.tol = 1e-6
    )$value
  }
  
  # ===========================================================================
  # Survival function
  # ===========================================================================
  
  surv_i <- function(T, x1i, x2i, w1i, w2i, ui) {
    
    exp(
      -cumhaz_i(
        T = T,
        x1i = x1i,
        x2i = x2i,
        w1i = w1i,
        w2i = w2i,
        ui = ui
      )
    )
  }
  
  # ===========================================================================
  # Simulate event time
  # ===========================================================================
  
  sim_event_time_i <- function(
    x1i,
    x2i,
    w1i,
    w2i,
    ui,
    Tmax = followup
  ) {
    
    U <- runif(1)
    
    S_Tmax <- surv_i(
      T = Tmax,
      x1i = x1i,
      x2i = x2i,
      w1i = w1i,
      w2i = w2i,
      ui = ui
    )
    
    if (S_Tmax > U) {
      return(Inf)
    }
    
    f <- function(T) {
      
      surv_i(
        T = T,
        x1i = x1i,
        x2i = x2i,
        w1i = w1i,
        w2i = w2i,
        ui = ui
      ) - U
    }
    
    tryCatch(
      uniroot(
        f,
        lower = .Machine$double.eps,
        upper = Tmax
      )$root,
      error = function(e) Tmax
    )
  }
  
  # ===========================================================================
  # Generate longitudinal and survival data
  # ===========================================================================
  
  Y1 <- matrix(NA_real_, nsample, length(t))
  Y2 <- matrix(NA_real_, nsample, length(t))
  Y3 <- matrix(NA_real_, nsample, length(t))
  
  for (i in seq_len(nsample)) {
    
    evt <- sim_event_time_i(
      x1i = x1[i],
      x2i = x2[i],
      w1i = w1[i],
      w2i = w2[i],
      ui = u[i, ]
    )
    
    # Independent censoring
    td <- rexp(1, rate = 1)
    
    obs_time <- min(evt, td, followup)
    
    time[i] <- obs_time
    
    death[i] <- as.integer(
      is.finite(evt) &&
        evt <= td &&
        evt <= followup
    )
    
    # -------------------------------------------------------------------------
    # Longitudinal measurements
    # -------------------------------------------------------------------------
    
    for (j in seq_along(t)) {
      
      if (t[j] <= obs_time) {
        
        # Marker 1
        Y1[i, j] <-
          Beta1[1] +
          Beta1[2] * x1[i] +
          Beta1[3] * x2[i] +
          Beta1[4] * t[j] +
          u[i, 1] +
          u[i, 2] * t[j] +
          rnorm(1, 0, sigma)
        
        # Marker 2
        Y2[i, j] <-
          Beta2[1] +
          Beta2[2] * x1[i] +
          Beta2[3] * x2[i] +
          Beta2[4] * t[j] +
          u[i, 3] +
          u[i, 4] * t[j] +
          rnorm(1, 0, sigma)
        
        # Marker 3
        Y3[i, j] <-
          Beta3[1] +
          Beta3[2] * x1[i] +
          Beta3[3] * x2[i] +
          Beta3[4] * t[j] +
          u[i, 5] +
          u[i, 6] * t[j] +
          rnorm(1, 0, sigma)
      }
    }
  }
  
  # ===========================================================================
  # Long format
  # ===========================================================================
  
  long.list <- vector("list", nsample)
  
  for (i in seq_len(nsample)) {
    
    idxs <- which(!is.na(Y1[i, ]))
    
    if (length(idxs) == 0) {
      next
    }
    
    long.list[[i]] <- data.frame(
      id       = i,
      obstime  = t[idxs],
      x1       = x1[i],
      x2       = x2[i],
      w1       = w1[i],
      w2       = w2[i],
      Event    = death[i],
      survtime = time[i],
      Y1       = Y1[i, idxs],
      Y2       = Y2[i, idxs],
      Y3       = Y3[i, idxs]
    )
  }
  
  long.data <- do.call(rbind, long.list)
  
  # ===========================================================================
  # Survival data
  # ===========================================================================
  
  surv.data <- data.frame(
    id       = id,
    x1       = x1,
    x2       = x2,
    w1       = w1,
    w2       = w2,
    survtime = time,
    Event    = death
  )
  
  # ===========================================================================
  # Train / validation split
  # ===========================================================================
  
  train_ids <- 1:500
  valid_ids <- setdiff(id, train_ids)
  
  long_train <- subset(
    long.data,
    id %in% train_ids
  )
  
  long_valid <- subset(
    long.data,
    id %in% valid_ids
  )
  
  surv_train <- subset(
    surv.data,
    id %in% train_ids
  )
  
  surv_valid <- subset(
    surv.data,
    id %in% valid_ids
  )
  # ===========================================================================
  # ===========================================================================
  # ===========================================================================
  
  start_total=Sys.time()
  dp_results_pen <- vector("list", length(S))
  
  
  # =============================================================================
  # Loop over landmark times
  # =============================================================================
  AUC_pen=BS_pen=c()
  for (j in seq_along(S)) {
    
    s <- S[j]
    Delta_t <- t_pred
    
    
    start_landmark <- Sys.time()
    
    
    train_lm_ids <- surv_train$id[
      surv_train$survtime >= s
    ]
    
    valid_lm_ids <- surv_valid$id[
      surv_valid$survtime >= s
    ]
    
    cat(
      "Training subjects at landmark:",
      length(train_lm_ids),
      "\n"
    )
    
    cat(
      "Validation subjects at landmark:",
      length(valid_lm_ids),
      "\n"
    )
    
    
    # ===========================================================================
    # Check sample size
    # ===========================================================================
    
    if (length(train_lm_ids) < 20 ||
        length(valid_lm_ids) < 10) {
      
      warning(
        paste(
          "Too few subjects at landmark",
          s,
          "- prediction skipped."
        )
      )
      
      AUC_pen[kkk, j] <- NA_real_
      BS_pen[kkk, j] <- NA_real_
      
      dp_results_pen[[j]] <- NULL
      
      next
    }
    
    
    
    long_train_lm <- subset(
      long_train,
      id %in% train_lm_ids &
        obstime <= s
    )
    
    long_train_lm <- as.data.frame(long_train_lm)
    
    
    
    surv_train_lm <- subset(
      surv_train,
      id %in% train_lm_ids
    )
    
    surv_train_lm <- data.frame(
      id = surv_train_lm$id,
      time = surv_train_lm$survtime,
      event = surv_train_lm$Event,
      x1 = surv_train_lm$x1,
      x2 = surv_train_lm$x2,
      w1 = surv_train_lm$w1,
      w2 = surv_train_lm$w2
    )
    
    
    
    common_train_ids <- intersect(
      unique(long_train_lm$id),
      unique(surv_train_lm$id)
    )
    
    long_train_lm <- long_train_lm[
      long_train_lm$id %in% common_train_ids,
      ,
      drop = FALSE
    ]
    
    surv_train_lm <- surv_train_lm[
      surv_train_lm$id %in% common_train_ids,
      ,
      drop = FALSE
    ]
    
    
    # ===========================================================================
    # 5. Sort training data
    # ===========================================================================
    
    long_train_lm <- long_train_lm[
      order(
        long_train_lm$id,
        long_train_lm$obstime
      ),
      ,
      drop = FALSE
    ]
    
    surv_train_lm <- surv_train_lm[
      order(surv_train_lm$id),
      ,
      drop = FALSE
    ]
    
    
    step1_pen <- fit_lmms(
      y.names = c(
        "Y1",
        "Y2",
        "Y3"
      ),
      
      fixefs = ~ x1 + x2 + obstime,
      
      ranefs = ~ obstime | id,
      
      t.from.base = obstime,
      
      long.data = long_train_lm,
      
      surv.data = surv_train_lm,
      
      n.boots = 0,
      
      n.cores = 1,
      
      verbose = FALSE
    )
    
    
    # ===========================================================================
    # 7. PRC-LMM: STEP 2
    # ===========================================================================
    #
    # Calculate subject-specific predicted random effects.
    # ===========================================================================
    
    step2_pen <- summarize_lmms(
      object = step1_pen,
      
      n.cores = 1,
      
      verbose = FALSE
    )
    
    
    
    step3_pen <- fit_prclmm(
      object = step2_pen,
      
      surv.data = surv_train_lm,
      
      baseline.covs = ~ x1 + x2 + w1 + w2,
      
      penalty = "ridge",
      
      standardize = TRUE,
      
      pfac.base.covs = 0,
      
      n.cores = 1,
      
      verbose = FALSE
    )
    
    
    # ===========================================================================
    # 9. VALIDATION LONGITUDINAL DATA
    # ===========================================================================
    #
    # Only measurements available up to landmark s.
    # ===========================================================================
    
    long_valid_lm <- subset(
      long_valid,
      id %in% valid_lm_ids &
        obstime <= s
    )
    
    long_valid_lm <- as.data.frame(long_valid_lm)
    
    
    
    surv_valid_lm <- subset(
      surv_valid,
      id %in% valid_lm_ids
    )
    
    surv_valid_lm <- data.frame(
      id = surv_valid_lm$id,
      x1 = surv_valid_lm$x1,
      x2 = surv_valid_lm$x2,
      w1 = surv_valid_lm$w1,
      w2 = surv_valid_lm$w2
    )
    
    
    # ===========================================================================
    # 11. Ensure common validation subjects
    # ===========================================================================
    
    common_valid_ids <- intersect(
      unique(long_valid_lm$id),
      unique(surv_valid_lm$id)
    )
    
    long_valid_lm <- long_valid_lm[
      long_valid_lm$id %in% common_valid_ids,
      ,
      drop = FALSE
    ]
    
    surv_valid_lm <- surv_valid_lm[
      surv_valid_lm$id %in% common_valid_ids,
      ,
      drop = FALSE
    ]
    
    
    # ===========================================================================
    # 12. Sort validation data
    # ===========================================================================
    
    long_valid_lm <- long_valid_lm[
      order(
        long_valid_lm$id,
        long_valid_lm$obstime
      ),
      ,
      drop = FALSE
    ]
    
    surv_valid_lm <- surv_valid_lm[
      order(surv_valid_lm$id),
      ,
      drop = FALSE
    ]
    
    
    # ===========================================================================
    # 13. Dynamic prediction
    # ===========================================================================
    #
    # pencal directly returns:
    #
    # P(T > t | T > s, longitudinal history up to s)
    #
    # Therefore we specify:
    #
    # t = s + Delta_t
    #
    # and DO NOT divide by another survival probability.
    # ===========================================================================
    
    prediction_time <- s + Delta_t
    
    pred_pen <- survpred_prclmm(
      step1 = step1_pen,
      
      step2 = step2_pen,
      
      step3 = step3_pen,
      
      times = prediction_time,
      
      new.longdata = long_valid_lm,
      
      new.basecovs = surv_valid_lm
    )
    
    
    # ===========================================================================
    # 14. Extract prediction matrix
    # ===========================================================================
    
    pred_df <- as.data.frame(
      pred_pen$predicted_survival
    )
    
    cat("\nPrediction output:\n")
    print(head(pred_df))
    
    
    # ===========================================================================
    # 15. Identify prediction column
    # ===========================================================================
    
    prediction_column <- setdiff(
      names(pred_df),
      "id"
    )
    
    if (length(prediction_column) == 0) {
      
      stop(
        "Could not identify predicted survival column."
      )
      
    }
    
    pred_survival <- pred_df[[prediction_column[1]]]
    
    
    # ===========================================================================
    # IMPORTANT:
    #
    # The above can be affected by console formatting.
    # Use [[ ]] explicitly:
    # ===========================================================================
    
    pred_survival <- pred_df[[prediction_column[1]]]
    
    
    # ===========================================================================
    # 16. Prediction data
    # ===========================================================================
    
    prediction_data <- data.frame(
      id = pred_df$id,
      
      pred_survival = pred_survival
    )
    
    
    # ===========================================================================
    # 17. Add observed validation survival information
    # ===========================================================================
    #
    # Original survival time and event are needed for evaluation.
    # ===========================================================================
    
    observed_valid <- subset(
      surv_valid,
      id %in% prediction_data$id
    )
    
    observed_valid <- observed_valid[
      ,
      c(
        "id",
        "survtime",
        "Event"
      ),
      drop = FALSE
    ]
    
    prediction_data <- merge(
      prediction_data,
      observed_valid,
      by = "id",
      sort = FALSE
    )
    
    
    # ===========================================================================
    # 18. Define prediction horizon
    # ===========================================================================
    
    prediction_data$prediction_horizon <- prediction_time
    
    
    
    prediction_data$event_horizon <- as.integer(
      prediction_data$Event == 1 &
        prediction_data$survtime > s &
        prediction_data$survtime <= prediction_time
    )
    
    
    
    prediction_data$pred_risk <-
      1 - prediction_data$pred_survival
    
    
    # ===========================================================================
    # 21. Store predictions
    # ===========================================================================
    
    dp_results_pen[[j]] <- prediction_data
    
    Crit_pen <- Criteria(
      s = s,
      t = t_pred,
      Survt = prediction_data$survtime,
      CR = prediction_data$Event,
      P = prediction_data$pred_risk,
      cause = 1
    )
    
    AUC_pen[j] <-Crit_pen$Cri[1, 1]
    
    
    BS_pen[j] <-Crit_pen$Cri[2, 1]
    
    print(j)
  }
  
  cbind(AUC_pen,BS_pen)
  # =============================================================================
  # TOTAL COMPUTATIONAL TIME
  # =============================================================================
  
  TIME_pen[kkk] <- as.numeric(
    difftime(
      Sys.time(),
      start_total,
      units = "mins"
    )
  )
  
  
  # =============================================================================
  # STORE RESULTS
  # =============================================================================
  
  Results1[[kkk]] <- list(
    
    AUC_pen =  AUC_pen ,
    
    BS_pen = BS_pen,
    
    TIME_pen = TIME_pen[kkk] 
    
  )
  
  
  print(kkk)
}

auc=bs=matrix(0,NN,4)
ct=c()
for(kkk in 1:NN){  

auc[kkk,]=Results1[[kkk]]$AUC_pen
bs[kkk,]=Results1[[kkk]]$BS_pen
ct[kkk]=Results1[[kkk]]$TIME_pen
}

cbind(apply(auc,2,mean),apply(auc,2,sd))
cbind(apply(bs,2,mean),apply(bs,2,sd))
mean(ct);sd(ct)

 