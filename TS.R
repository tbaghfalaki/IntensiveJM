
# =============================================================================
# Libraries
# =============================================================================

rm(list = ls())

library(mvtnorm)
library(DPCri)
library(MASS)
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvtnorm)
library(JMbayes2)
library(DPCri)
library(lme4)

# =============================================================================
# Simulation parameters
# =============================================================================

NN <- 100

Results1 <- vector("list", NN)

# =============================================================================
# Simulation loop
# =============================================================================

for (kkk in 74:NN) {
  
  set.seed(kkk)
  
  nsample <- 1000
  followup <- 2
  t <- seq(0, followup, by = 0.005)
  
  # ---------------------------------------------------------------------------
  # Association and fixed-effect parameters
  # ---------------------------------------------------------------------------
  
  alpha <- c(-0.5, 0.3, -0.2)
  
  Beta1 <- c(-0.5, -0.5, 0.5, 0.5)
  Beta2 <- c(0.2, 0.3, -0.3, 0.4)
  Beta3 <- c(1.0, -0.2, 0.1, 0.6)
  
  sigma <- 1
  
  gamma_w <- c(0.2, -0.2)
  
  # ---------------------------------------------------------------------------
  # Covariates and survival variables
  # ---------------------------------------------------------------------------
  
  id <- 1:nsample
  
  x1 <- rnorm(nsample)
  x2 <- rbinom(nsample, 1, 0.5)
  
  w1 <- rnorm(nsample)
  w2 <- rbinom(nsample, 1, 0.5)
  
  death <- rep(0L, nsample)
  time <- rep(NA_real_, nsample)
  
  # ---------------------------------------------------------------------------
  # Random effects
  #
  # Three longitudinal markers.
  #
  # Each marker has:
  #   random intercept
  #   random slope
  #
  # Intercept/slope correlation = 0.6 within each marker.
  #
  # IMPORTANT:
  # The following covariance construction is kept exactly as provided.
  # ---------------------------------------------------------------------------
  
  D_block <- matrix(
    c(
      1, 0.6,
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
  
  # ---------------------------------------------------------------------------
  # Individual hazard function
  # ---------------------------------------------------------------------------
  
  haz_fun_i <- function(
    s,
    x1i,
    x2i,
    w1i,
    w2i,
    ui
  ) {
    
    # Marker 1
    mu1 <-
      Beta1[1] +
      Beta1[2] * x1i +
      Beta1[3] * x2i +
      Beta1[4] * s +
      ui[1] +
      ui[2] * s
    
    # Marker 2
    mu2 <-
      Beta2[1] +
      Beta2[2] * x1i +
      Beta2[3] * x2i +
      Beta2[4] * s +
      ui[3] +
      ui[4] * s
    
    # Marker 3
    mu3 <-
      Beta3[1] +
      Beta3[2] * x1i +
      Beta3[3] * x2i +
      Beta3[4] * s +
      ui[5] +
      ui[6] * s
    
    # Survival linear predictor
    eta <-
      gamma_w[1] * w1i +
      gamma_w[2] * w2i +
      alpha[1] * mu1 +
      alpha[2] * mu2 +
      alpha[3] * mu3
    
    exp(eta)
  }
  
  # ---------------------------------------------------------------------------
  # Cumulative hazard
  # ---------------------------------------------------------------------------
  
  cumhaz_i <- function(
    T,
    x1i,
    x2i,
    w1i,
    w2i,
    ui
  ) {
    
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
  
  # ---------------------------------------------------------------------------
  # Survival function
  # ---------------------------------------------------------------------------
  
  surv_i <- function(
    T,
    x1i,
    x2i,
    w1i,
    w2i,
    ui
  ) {
    
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
  
  # ---------------------------------------------------------------------------
  # Simulate event time
  # ---------------------------------------------------------------------------
  
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
    
    res <- tryCatch(
      
      uniroot(
        f,
        lower = .Machine$double.eps,
        upper = Tmax
      )$root,
      
      error = function(e) Tmax
    )
    
    res
  }
  
  # =============================================================================
  # Generate longitudinal and survival data
  # =============================================================================
  
  Y1 <- matrix(
    NA_real_,
    nsample,
    length(t)
  )
  
  Y2 <- matrix(
    NA_real_,
    nsample,
    length(t)
  )
  
  Y3 <- matrix(
    NA_real_,
    nsample,
    length(t)
  )
  
  for (i in 1:nsample) {
    
    # True event time
    evt <- sim_event_time_i(
      x1i = x1[i],
      x2i = x2[i],
      w1i = w1[i],
      w2i = w2[i],
      ui = u[i, ]
    )
    
    # Independent censoring time
    td <- rexp(
      1,
      rate = 1
    )
    
    # Observed time
    obs_time <- min(
      evt,
      td,
      followup
    )
    
    time[i] <- obs_time
    
    # Event indicator
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
          rnorm(
            1,
            0,
            sigma
          )
        
        # Marker 2
        Y2[i, j] <-
          Beta2[1] +
          Beta2[2] * x1[i] +
          Beta2[3] * x2[i] +
          Beta2[4] * t[j] +
          u[i, 3] +
          u[i, 4] * t[j] +
          rnorm(
            1,
            0,
            sigma
          )
        
        # Marker 3
        Y3[i, j] <-
          Beta3[1] +
          Beta3[2] * x1[i] +
          Beta3[3] * x2[i] +
          Beta3[4] * t[j] +
          u[i, 5] +
          u[i, 6] * t[j] +
          rnorm(
            1,
            0,
            sigma
          )
      }
    }
  }
  
  # =============================================================================
  # Long format
  # =============================================================================
  
  long.list <- vector(
    "list",
    nsample
  )
  
  for (i in 1:nsample) {
    
    idxs <- which(
      !is.na(Y1[i, ])
    )
    
    if (length(idxs) == 0) {
      next
    }
    
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
  
  long.data <- do.call(
    rbind,
    long.list
  )
  
  # =============================================================================
  # Survival data
  # =============================================================================
  
  surv.data <- data.frame(
    id = id,
    x1 = x1,
    x2 = x2,
    w1 = w1,
    w2 = w2,
    survtime = time,
    Event = death
  )
  
  # =============================================================================
  # Train / validation split
  # =============================================================================
  
  train_ids <- 1:500
  
  valid_ids <- setdiff(
    id,
    train_ids
  )
  
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
  
  # =============================================================================
  # Quick checks
  # =============================================================================
  
  if (kkk == 1) {
    
    cat("\n====================================================\n")
    cat("DATA GENERATION CHECKS\n")
    cat("====================================================\n")
    
    cat(
      "Mean number of longitudinal observations:",
      mean(table(long.data$id)),
      "\n"
    )
    
    cat(
      "Median number of longitudinal observations:",
      median(table(long.data$id)),
      "\n"
    )
    
    cat(
      "Minimum number of longitudinal observations:",
      min(table(long.data$id)),
      "\n"
    )
    
    cat(
      "Maximum number of longitudinal observations:",
      max(table(long.data$id)),
      "\n"
    )
    
    cat(
      "Mean survival time:",
      mean(surv.data$survtime),
      "\n"
    )
    
    cat(
      "Event rate:",
      mean(surv.data$Event),
      "\n"
    )
  }
  
  # =============================================================================
  # TWO-STAGE APPROACH
  #
  # Stage 1:
  # Separate LMMs for Y1, Y2 and Y3.
  #
  # Stage 2:
  # Extended Cox model using subject-specific LMM trajectories.
  #
  # Dynamic prediction:
  # Landmark times = 0.25, 0.50, 0.75, 1.00
  # Prediction horizon = 0.25
  # =============================================================================
  
  start_TS <- Sys.time()
  
  # =============================================================================
  # Stage 1: LMMs
  # =============================================================================
  
  lmm1 <- lmer(
    Y1 ~ obstime + x1 + x2 +
      (1 + obstime | id),
    data = long_train,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 2e5
      )
    )
  )
  
  lmm2 <- lmer(
    Y2 ~ obstime + x1 + x2 +
      (1 + obstime | id),
    data = long_train,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 2e5
      )
    )
  )
  
  lmm3 <- lmer(
    Y3 ~ obstime + x1 + x2 +
      (1 + obstime | id),
    data = long_train,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 2e5
      )
    )
  )
  
  # =============================================================================
  # Helper function:
  # Extract random-effects covariance matrix from LMM
  # =============================================================================
  
  get_D_lmm <- function(lmm_fit) {
    
    Dhat <- as.matrix(
      VarCorr(lmm_fit)$id
    )
    
    Dhat[1:2, 1:2]
  }
  
  # =============================================================================
  # Helper function:
  #
  # Predict subject-specific longitudinal value.
  #
  # Only observations <= landmark are used to calculate the BLUP.
  # =============================================================================
  
  predict_lmm_value <- function(
    idv,
    pred_time,
    landmark,
    marker,
    lmm_fit,
    history_data
  ) {
    
    # -------------------------------------------------------------------------
    # Subject's observed history up to landmark
    # -------------------------------------------------------------------------
    
    subj <- history_data[
      history_data$id == idv &
        history_data$obstime <= landmark,
    ]
    
    # -------------------------------------------------------------------------
    # Fixed effects
    # -------------------------------------------------------------------------
    
    beta_hat <- fixef(lmm_fit)
    
    # -------------------------------------------------------------------------
    # Random-effects covariance
    # -------------------------------------------------------------------------
    
    Sigma_b <- get_D_lmm(
      lmm_fit
    )
    
    # -------------------------------------------------------------------------
    # Residual variance
    # -------------------------------------------------------------------------
    
    sigma2 <- sigma(lmm_fit)^2
    
    # -------------------------------------------------------------------------
    # Subject covariates
    # -------------------------------------------------------------------------
    
    x1i <- history_data$x1[
      history_data$id == idv
    ][1]
    
    x2i <- history_data$x2[
      history_data$id == idv
    ][1]
    
    # -------------------------------------------------------------------------
    # Fixed-effect prediction
    # -------------------------------------------------------------------------
    
    newdata_pred <- data.frame(
      obstime = pred_time,
      x1 = x1i,
      x2 = x2i
    )
    
    X_pred <- model.matrix(
      ~ obstime + x1 + x2,
      data = newdata_pred
    )
    
    fixed_pred <- as.numeric(
      X_pred %*% beta_hat
    )
    
    # -------------------------------------------------------------------------
    # If no history is available
    # -------------------------------------------------------------------------
    
    if (nrow(subj) == 0) {
      return(fixed_pred)
    }
    
    # -------------------------------------------------------------------------
    # Observed marker values
    # -------------------------------------------------------------------------
    
    Y <- subj[[marker]]
    
    # -------------------------------------------------------------------------
    # Fixed-effect design matrix
    # -------------------------------------------------------------------------
    
    X <- model.matrix(
      ~ obstime + x1 + x2,
      data = subj
    )
    
    # -------------------------------------------------------------------------
    # Random-effect design matrix
    # -------------------------------------------------------------------------
    
    Z <- cbind(
      1,
      subj$obstime
    )
    
    # -------------------------------------------------------------------------
    # Marginal covariance
    # -------------------------------------------------------------------------
    
    V <- Z %*%
      Sigma_b %*%
      t(Z) +
      sigma2 * diag(
        nrow(Z)
      )
    
    # -------------------------------------------------------------------------
    # BLUP of random effects
    # -------------------------------------------------------------------------
    
    b_hat <- Sigma_b %*%
      t(Z) %*%
      solve(
        V,
        Y - X %*% beta_hat
      )
    
    b_hat <- as.numeric(
      b_hat
    )
    
    # -------------------------------------------------------------------------
    # Subject-specific prediction
    # -------------------------------------------------------------------------
    
    Z_pred <- c(
      1,
      pred_time
    )
    
    pred <-
      fixed_pred +
      as.numeric(
        Z_pred %*% b_hat
      )
    
    pred
  }
  
  # =============================================================================
  # Stage 2:
  # Construct predicted longitudinal trajectories in training data
  #
  # Each predicted value at time t uses only observations available up to t.
  # =============================================================================
  
  long_td <- long_train %>%
    dplyr::select(
      id,
      obstime
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(
      id,
      obstime
    )
  
  # =============================================================================
  # Training predicted trajectory: Y1
  # =============================================================================
  
  long_td$Y1_pred <- mapply(
    
    function(idv, tt) {
      
      predict_lmm_value(
        idv = idv,
        pred_time = tt,
        landmark = tt,
        marker = "Y1",
        lmm_fit = lmm1,
        history_data = long_train
      )
      
    },
    
    long_td$id,
    long_td$obstime
  )
  
  # =============================================================================
  # Training predicted trajectory: Y2
  # =============================================================================
  
  long_td$Y2_pred <- mapply(
    
    function(idv, tt) {
      
      predict_lmm_value(
        idv = idv,
        pred_time = tt,
        landmark = tt,
        marker = "Y2",
        lmm_fit = lmm2,
        history_data = long_train
      )
      
    },
    
    long_td$id,
    long_td$obstime
  )
  
  # =============================================================================
  # Training predicted trajectory: Y3
  # =============================================================================
  
  long_td$Y3_pred <- mapply(
    
    function(idv, tt) {
      
      predict_lmm_value(
        idv = idv,
        pred_time = tt,
        landmark = tt,
        marker = "Y3",
        lmm_fit = lmm3,
        history_data = long_train
      )
      
    },
    
    long_td$id,
    long_td$obstime
  )
  
  # =============================================================================
  # Stage 2: Create extended survival data
  # =============================================================================
  
  td_base <- tmerge(
    data1 = surv_train,
    data2 = surv_train,
    id = id,
    death = event(
      survtime,
      Event
    )
  )
  
  td_data <- tmerge(
    td_base,
    long_td,
    id = id,
    
    Y1_pred_td = tdc(
      obstime,
      Y1_pred
    ),
    
    Y2_pred_td = tdc(
      obstime,
      Y2_pred
    ),
    
    Y3_pred_td = tdc(
      obstime,
      Y3_pred
    )
  )
  
  # =============================================================================
  # Extended Cox model
  # =============================================================================
  
  cox_ext <- coxph(
    
    Surv(
      tstart,
      tstop,
      death
    ) ~
      Y1_pred_td +
      Y2_pred_td +
      Y3_pred_td +
      w1 +
      w2,
    
    data = td_data,
    
    id = id,
    
    cluster = id,
    
    ties = "breslow"
  )
  
  # =============================================================================
  # Dynamic prediction settings
  # =============================================================================
  
  S <- c(
    0.25,
    0.50,
    0.75,
    1.00
  )
  
  t_pred <- 0.25
  
  # =============================================================================
  # Helper function:
  # Baseline cumulative hazard
  # =============================================================================
  
  get_H0 <- function(
    cox_fit,
    tt
  ) {
    
    bh <- basehaz(
      cox_fit,
      centered = FALSE
    )
    
    if (nrow(bh) == 0) {
      return(0)
    }
    
    if (tt <= 0) {
      return(0)
    }
    
    if (tt < min(bh$time)) {
      return(0)
    }
    
    H0 <- approx(
      x = bh$time,
      y = bh$hazard,
      xout = tt,
      method = "constant",
      f = 1,
      rule = 2
    )$y
    
    as.numeric(H0)
  }
  
  # =============================================================================
  # Helper function:
  # Dynamic prediction
  #
  # P(T <= s+t | T>s, history up to s)
  # =============================================================================
  
  predict_dynamic_TS <- function(
    cox_fit,
    valid_at_risk,
    long_valid,
    lmm1,
    lmm2,
    lmm3,
    s,
    t_pred
  ) {
    
    ids <- valid_at_risk$id
    
    # -------------------------------------------------------------------------
    # Prediction grid
    # -------------------------------------------------------------------------
    
    pred_grid <- seq(
      s,
      s + t_pred,
      by = 0.005
    )
    
    if (
      tail(pred_grid, 1) <
      s + t_pred
    ) {
      
      pred_grid <- c(
        pred_grid,
        s + t_pred
      )
    }
    
    # -------------------------------------------------------------------------
    # Result vector
    # -------------------------------------------------------------------------
    
    DP <- numeric(
      length(ids)
    )
    
    # -------------------------------------------------------------------------
    # Loop over validation subjects
    # -------------------------------------------------------------------------
    
    for (ii in seq_along(ids)) {
      
      idv <- ids[ii]
      
      # -----------------------------------------------------------------------
      # Subject-specific covariates
      # -----------------------------------------------------------------------
      
      subj_info <- valid_at_risk[
        valid_at_risk$id == idv,
      ]
      
      w1i <- subj_info$w1[1]
      w2i <- subj_info$w2[1]
      
      # -----------------------------------------------------------------------
      # Predicted Y1 trajectory
      #
      # Only history <= s is used.
      # -----------------------------------------------------------------------
      
      Y1_grid <- sapply(
        
        pred_grid,
        
        function(tt) {
          
          predict_lmm_value(
            idv = idv,
            pred_time = tt,
            landmark = s,
            marker = "Y1",
            lmm_fit = lmm1,
            history_data = long_valid
          )
          
        }
      )
      
      # -----------------------------------------------------------------------
      # Predicted Y2 trajectory
      # -----------------------------------------------------------------------
      
      Y2_grid <- sapply(
        
        pred_grid,
        
        function(tt) {
          
          predict_lmm_value(
            idv = idv,
            pred_time = tt,
            landmark = s,
            marker = "Y2",
            lmm_fit = lmm2,
            history_data = long_valid
          )
          
        }
      )
      
      # -----------------------------------------------------------------------
      # Predicted Y3 trajectory
      # -----------------------------------------------------------------------
      
      Y3_grid <- sapply(
        
        pred_grid,
        
        function(tt) {
          
          predict_lmm_value(
            idv = idv,
            pred_time = tt,
            landmark = s,
            marker = "Y3",
            lmm_fit = lmm3,
            history_data = long_valid
          )
          
        }
      )
      
      # -----------------------------------------------------------------------
      # Cox prediction data
      # -----------------------------------------------------------------------
      
      newdata_grid <- data.frame(
        
        Y1_pred_td = Y1_grid,
        
        Y2_pred_td = Y2_grid,
        
        Y3_pred_td = Y3_grid,
        
        w1 = w1i,
        
        w2 = w2i
      )
      
      # -----------------------------------------------------------------------
      # Cox linear predictor
      # -----------------------------------------------------------------------
      
      lp_grid <- predict(
        cox_fit,
        newdata = newdata_grid,
        type = "lp"
      )
      
      # -----------------------------------------------------------------------
      # Baseline cumulative hazard
      # -----------------------------------------------------------------------
      
      H_grid <- sapply(
        
        pred_grid,
        
        function(tt) {
          
          get_H0(
            cox_fit = cox_fit,
            tt = tt
          )
          
        }
      )
      
      # -----------------------------------------------------------------------
      # Incremental baseline cumulative hazard
      # -----------------------------------------------------------------------
      
      H_start <- get_H0(
        cox_fit,
        s
      )
      
      dH <- diff(
        c(
          H_start,
          H_grid
        )
      )
      
      # -----------------------------------------------------------------------
      # Integrated cumulative hazard
      # -----------------------------------------------------------------------
      
      cumulative_hazard <- sum(
        dH * exp(lp_grid)
      )
      
      # -----------------------------------------------------------------------
      # Dynamic event probability
      # -----------------------------------------------------------------------
      
      DP[ii] <- 1 -
        exp(
          -cumulative_hazard
        )
    }
    
    DP
  }
  
  # =============================================================================
  # Dynamic prediction and evaluation
  # =============================================================================
  
  AUC_TS <- rep(
    NA_real_,
    length(S)
  )
  
  BS_TS <- rep(
    NA_real_,
    length(S)
  )
  
  dp_results <- vector(
    "list",
    length(S)
  )
  
  names(dp_results) <- paste0(
    "TS_",
    S
  )
  
  # =============================================================================
  # Loop over landmarks
  # =============================================================================
  
  for (j in seq_along(S)) {
    
    s <- S[j]
    
    cat(
      "Simulation:",
      kkk,
      "| Landmark:",
      s,
      "\n"
    )
    
    # -------------------------------------------------------------------------
    # Subjects event-free at landmark
    # -------------------------------------------------------------------------
    
    valid_at_risk <- subset(
      surv_valid,
      survtime >= s
    )
    
    # -------------------------------------------------------------------------
    # Dynamic prediction
    # -------------------------------------------------------------------------
    
    DP <- predict_dynamic_TS(
      
      cox_fit = cox_ext,
      
      valid_at_risk = valid_at_risk,
      
      long_valid = long_valid,
      
      lmm1 = lmm1,
      lmm2 = lmm2,
      lmm3 = lmm3,
      
      s = s,
      
      t_pred = t_pred
    )
    
    # -------------------------------------------------------------------------
    # Store predictions
    # -------------------------------------------------------------------------
    
    dp_results[[j]] <- data.frame(
      
      id = valid_at_risk$id,
      
      landmark = s,
      
      horizon = t_pred,
      
      DP = DP,
      
      survtime = valid_at_risk$survtime,
      
      Event = valid_at_risk$Event
    )
    
    # -------------------------------------------------------------------------
    # AUC and Brier score
    # -------------------------------------------------------------------------
    
    Crit <- tryCatch(
      
      Criteria(
        
        s = s,
        
        t = t_pred,
        
        Survt = valid_at_risk$survtime,
        
        CR = valid_at_risk$Event,
        
        P = DP,
        
        cause = 1
        
      ),
      
      error = function(e) {
        
        message(
          "Criteria failed at landmark ",
          s,
          ": ",
          e$message
        )
        
        NULL
      }
    )
    
    # -------------------------------------------------------------------------
    # Save AUC and Brier score
    # -------------------------------------------------------------------------
    
    if (!is.null(Crit)) {
      
      AUC_TS[j] <- as.numeric(
        Crit$Cri[1, 1]
      )
      
      BS_TS[j] <- as.numeric(
        Crit$Cri[2, 1]
      )
    }
  }
  
  # =============================================================================
  # Computational time
  # =============================================================================
  
  end_TS <- Sys.time()
  
  time_TS <- as.numeric(
    difftime(
      end_TS,
      start_TS,
      units = "mins"
    )
  )
  
  # =============================================================================
  # Store simulation results
  # =============================================================================
  
  Results1[[kkk]] <- list(
    
    AUC = AUC_TS,
    
    BS = BS_TS,
    
    time = time_TS,
    
    predictions = dp_results,
    
    lmm1 = lmm1,
    
    lmm2 = lmm2,
    
    lmm3 = lmm3,
    
    cox = cox_ext
  )
  
  cat(
    "Simulation",
    kkk,
    "completed in",
    round(time_TS, 2),
    "minutes.\n"
  )
}

# =============================================================================
# SUMMARY ACROSS 100 SIMULATIONS
# =============================================================================

S <- c(
  0.25,
  0.50,
  0.75,
  1.00
)

# =============================================================================
# Extract AUC
# =============================================================================

AUC_matrix <- do.call(
  rbind,
  lapply(
    Results1,
    function(x) x$AUC
  )
)

# =============================================================================
# Extract Brier Score
# =============================================================================

BS_matrix <- do.call(
  rbind,
  lapply(
    Results1,
    function(x) x$BS
  )
)

# =============================================================================
# Extract computational time
# =============================================================================

Time_vector <- sapply(
  Results1,
  function(x) x$time
)

# =============================================================================
# Mean and SD AUC
# =============================================================================

mean_AUC_TS <- colMeans(
  AUC_matrix,
  na.rm = TRUE
)

sd_AUC_TS <- apply(
  AUC_matrix,
  2,
  sd,
  na.rm = TRUE
)

# =============================================================================
# Mean and SD Brier Score
# =============================================================================

mean_BS_TS <- colMeans(
  BS_matrix,
  na.rm = TRUE
)

sd_BS_TS <- apply(
  BS_matrix,
  2,
  sd,
  na.rm = TRUE
)

# =============================================================================
# Mean and SD computational time
# =============================================================================

mean_time_TS <- mean(
  Time_vector,
  na.rm = TRUE
)

sd_time_TS <- sd(
  Time_vector,
  na.rm = TRUE
)

# =============================================================================
# Summary table
# =============================================================================

Summary_TS <- data.frame(
  
  Landmark = S,
  
  Mean_AUC = mean_AUC_TS,
  
  SD_AUC = sd_AUC_TS,
  
  Mean_BS = mean_BS_TS,
  
  SD_BS = sd_BS_TS
)

# =============================================================================
# Computational time summary
# =============================================================================

Time_summary_TS <- data.frame(
  
  Method = "Two-stage LMM",
  
  Mean_Time = mean_time_TS,
  
  SD_Time = sd_time_TS
)

# =============================================================================
# Print results
# =============================================================================

cat("\n====================================================\n")
cat("TWO-STAGE LMM DYNAMIC PREDICTION\n")
cat("====================================================\n\n")

print(
  Summary_TS
)

cat("\n====================================================\n")
cat("COMPUTATIONAL TIME (MINUTES)\n")
cat("====================================================\n\n")

print(
  Time_summary_TS
)

# =============================================================================
# Save results
# =============================================================================

save(
  Results1,
  Summary_TS,
  Time_summary_TS,
  file = "TS_LMM_results.RData"
)

 