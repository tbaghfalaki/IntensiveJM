rm(list = ls())

# =============================================================================
# Libraries
# =============================================================================

library(MASS)
library(mvtnorm)
library(lme4)
library(survival)
library(DPCri)

# =============================================================================
# Simulation settings
# =============================================================================

NN <- 100

nsample <- 1000
followup <- 2

S <- c(0.25, 0.50, 0.75, 1.00)
t_pred <- 0.25

# =============================================================================
# Containers for simulation results
# =============================================================================

# Each row = one simulation
# Columns = four landmark times

AUC_LMM  <- matrix(NA_real_, NN, length(S))
BS_LMM   <- matrix(NA_real_, NN, length(S))

AUC_LOCF <- matrix(NA_real_, NN, length(S))
BS_LOCF  <- matrix(NA_real_, NN, length(S))

# Computational time
TIME_LMM  <- numeric(NN)
TIME_LOCF <- numeric(NN)

# Optional: store all predictions
Results1 <- vector("list", NN)

# =============================================================================
# Simulation loop
# =============================================================================

for (kkk in seq_len(NN)) {
  
  cat("\n============================================================\n")
  cat("Simulation:", kkk, "of", NN, "\n")
  cat("============================================================\n")
  
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
  # Helper: BLUP at landmark
  # ===========================================================================
  
  get_blup_at_landmark <- function(
    subj_data,
    s,
    lmm_fit
  ) {
    
    beta_hat <- fixef(lmm_fit)
    
    Sigma_b <- as.matrix(
      VarCorr(lmm_fit)$id
    )
    
    sigma2 <- sigma(lmm_fit)^2
    
    # -------------------------------------------------------------------------
    # Design matrix at landmark
    # -------------------------------------------------------------------------
    
    newdata_s <- data.frame(
      obstime = s,
      x1 = subj_data$x1[1],
      x2 = subj_data$x2[1]
    )
    
    X_landmark <- model.matrix(
      ~ obstime + x1 + x2,
      data = newdata_s
    )
    
    # -------------------------------------------------------------------------
    # No observations
    # -------------------------------------------------------------------------
    
    if (nrow(subj_data) == 0) {
      
      return(
        as.numeric(
          X_landmark %*% beta_hat
        )
      )
    }
    
    # -------------------------------------------------------------------------
    # Subject-specific history
    # -------------------------------------------------------------------------
    
    X <- model.matrix(
      ~ obstime + x1 + x2,
      data = subj_data
    )
    
    Z <- cbind(
      1,
      subj_data$obstime
    )
    
    Y <- subj_data$Y
    
    V <-
      Z %*%
      Sigma_b %*%
      t(Z) +
      sigma2 * diag(nrow(Z))
    
    # BLUP
    b_hat <-
      Sigma_b %*%
      t(Z) %*%
      solve(
        V,
        Y - X %*% beta_hat
      )
    
    # Random effects evaluated at landmark
    Z_landmark <- c(1, s)
    
    as.numeric(
      X_landmark %*% beta_hat +
        Z_landmark %*% b_hat
    )
  }
  
  # ===========================================================================
  # Helper: baseline cumulative hazard
  # ===========================================================================
  
  get_H0 <- function(cox_fit, t_pred) {
    
    H0 <- basehaz(
      cox_fit,
      centered = FALSE
    )
    
    if (nrow(H0) == 0) {
      return(0)
    }
    
    # Cumulative hazard at t_pred
    if (t_pred < min(H0$time)) {
      return(0)
    }
    
    H0_t <- approx(
      x = H0$time,
      y = H0$hazard,
      xout = t_pred,
      method = "constant",
      f = 1,
      rule = 2
    )$y
    
    as.numeric(H0_t)
  }
  
  # ===========================================================================
  # Helper: LOCF
  # ===========================================================================
  
  get_locf <- function(df, marker) {
    
    if (nrow(df) == 0) {
      
      return(
        data.frame(
          id = integer(0),
          value = numeric(0)
        )
      )
    }
    
    split_data <- split(
      df,
      df$id
    )
    
    out <- lapply(
      split_data,
      function(d) {
        
        d <- d[order(d$obstime), ]
        
        last <- d[nrow(d), ]
        
        data.frame(
          id = last$id,
          value = last[[marker]]
        )
      }
    )
    
    out <- do.call(rbind, out)
    
    rownames(out) <- NULL
    
    out
  }
  
  # ===========================================================================
  # Containers for this simulation
  # ===========================================================================
  
  dp_results_lmm <- vector("list", length(S))
  dp_results_locf <- vector("list", length(S))
  
  # ===========================================================================
  # Start LMM timing
  # ===========================================================================
  
  start_lmm <- Sys.time()
  
  # ===========================================================================
  # 1. LMM LANDMARKING
  # ===========================================================================
  
  for (j in seq_along(S)) {
    
    s <- S[j]
    
    cat("  LMM landmark =", s, "\n")
    
    # -------------------------------------------------------------------------
    # Subjects at risk
    # -------------------------------------------------------------------------
    
    train_at_risk <- subset(
      surv_train,
      survtime >= s
    )
    
    valid_at_risk <- subset(
      surv_valid,
      survtime >= s
    )
    
    if (
      nrow(train_at_risk) < 10 ||
      nrow(valid_at_risk) == 0
    ) {
      next
    }
    
    # -------------------------------------------------------------------------
    # Longitudinal history up to landmark
    # -------------------------------------------------------------------------
    
    long_train_s <- subset(
      long_train,
      id %in% train_at_risk$id &
        obstime <= s
    )
    
    long_valid_s <- subset(
      long_valid,
      id %in% valid_at_risk$id &
        obstime <= s
    )
    
    # ========================================================================
    # Fit three separate LMMs
    # ========================================================================
    
    lmm1 <- lmer(
      Y1 ~ obstime + x1 + x2 +
        (1 + obstime | id),
      data = long_train_s,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    )
    
    lmm2 <- lmer(
      Y2 ~ obstime + x1 + x2 +
        (1 + obstime | id),
      data = long_train_s,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    )
    
    lmm3 <- lmer(
      Y3 ~ obstime + x1 + x2 +
        (1 + obstime | id),
      data = long_train_s,
      REML = FALSE,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    )
    
    # ========================================================================
    # Training predictions
    # ========================================================================
    
    Ypred_train <- matrix(
      NA_real_,
      nrow = nrow(train_at_risk),
      ncol = 3
    )
    
    for (i in seq_len(nrow(train_at_risk))) {
      
      idv <- train_at_risk$id[i]
      
      subj <- subset(
        long_train_s,
        id == idv &
          obstime <= s
      )
      
      # Marker 1
      subj1 <- subj[, c(
        "id", "Y1", "obstime", "x1", "x2"
      )]
      
      names(subj1)[names(subj1) == "Y1"] <- "Y"
      
      Ypred_train[i, 1] <-
        get_blup_at_landmark(
          subj1,
          s,
          lmm1
        )
      
      # Marker 2
      subj2 <- subj[, c(
        "id", "Y2", "obstime", "x1", "x2"
      )]
      
      names(subj2)[names(subj2) == "Y2"] <- "Y"
      
      Ypred_train[i, 2] <-
        get_blup_at_landmark(
          subj2,
          s,
          lmm2
        )
      
      # Marker 3
      subj3 <- subj[, c(
        "id", "Y3", "obstime", "x1", "x2"
      )]
      
      names(subj3)[names(subj3) == "Y3"] <- "Y"
      
      Ypred_train[i, 3] <-
        get_blup_at_landmark(
          subj3,
          s,
          lmm3
        )
    }
    
    # ========================================================================
    # Landmark Cox model
    # ========================================================================
    
    train_lm <- train_at_risk
    
    train_lm$Y1_pred <- Ypred_train[, 1]
    train_lm$Y2_pred <- Ypred_train[, 2]
    train_lm$Y3_pred <- Ypred_train[, 3]
    
    train_lm <- train_lm[
      complete.cases(
        train_lm[, c(
          "Y1_pred",
          "Y2_pred",
          "Y3_pred"
        )]
      ),
    ]
    
    if (nrow(train_lm) < 10) {
      next
    }
    
    cox_lmm <- coxph(
      Surv(
        survtime - s,
        Event
      ) ~
        Y1_pred +
        Y2_pred +
        Y3_pred +
        w1 +
        w2,
      data = train_lm,
      ties = "breslow",
      x = TRUE
    )
    
    # ========================================================================
    # Validation predictions
    # ========================================================================
    
    Ypred_valid <- matrix(
      NA_real_,
      nrow = nrow(valid_at_risk),
      ncol = 3
    )
    
    for (i in seq_len(nrow(valid_at_risk))) {
      
      idv <- valid_at_risk$id[i]
      
      subj <- subset(
        long_valid_s,
        id == idv &
          obstime <= s
      )
      
      # Marker 1
      subj1 <- subj[, c(
        "id", "Y1", "obstime", "x1", "x2"
      )]
      
      names(subj1)[names(subj1) == "Y1"] <- "Y"
      
      Ypred_valid[i, 1] <-
        get_blup_at_landmark(
          subj1,
          s,
          lmm1
        )
      
      # Marker 2
      subj2 <- subj[, c(
        "id", "Y2", "obstime", "x1", "x2"
      )]
      
      names(subj2)[names(subj2) == "Y2"] <- "Y"
      
      Ypred_valid[i, 2] <-
        get_blup_at_landmark(
          subj2,
          s,
          lmm2
        )
      
      # Marker 3
      subj3 <- subj[, c(
        "id", "Y3", "obstime", "x1", "x2"
      )]
      
      names(subj3)[names(subj3) == "Y3"] <- "Y"
      
      Ypred_valid[i, 3] <-
        get_blup_at_landmark(
          subj3,
          s,
          lmm3
        )
    }
    
    # ========================================================================
    # Dynamic prediction
    # ========================================================================
    
    valid_lm <- valid_at_risk
    
    valid_lm$Y1_pred <- Ypred_valid[, 1]
    valid_lm$Y2_pred <- Ypred_valid[, 2]
    valid_lm$Y3_pred <- Ypred_valid[, 3]
    
    complete_valid <- complete.cases(
      valid_lm[, c(
        "Y1_pred",
        "Y2_pred",
        "Y3_pred"
      )]
    )
    
    DP <- rep(
      NA_real_,
      nrow(valid_lm)
    )
    
    if (sum(complete_valid) > 0) {
      
      valid_complete <- valid_lm[
        complete_valid,
      ]
      
      H0_tpred <- get_H0(
        cox_lmm,
        t_pred
      )
      
      lp_valid <- predict(
        cox_lmm,
        newdata = valid_complete,
        type = "lp"
      )
      
      surv_pred <- exp(
        -H0_tpred *
          exp(lp_valid)
      )
      
      DP_complete <- 1 - surv_pred
      
      DP[
        match(
          valid_complete$id,
          valid_lm$id
        )
      ] <- DP_complete
    }
    
    # ========================================================================
    # Store LMM predictions
    # ========================================================================
    
    dp_results_lmm[[j]] <- data.frame(
      id = valid_lm$id,
      Y1_pred = valid_lm$Y1_pred,
      Y2_pred = valid_lm$Y2_pred,
      Y3_pred = valid_lm$Y3_pred,
      DP = DP,
      landmark = s,
      horizon = s + t_pred
    )
  }
  
  # ===========================================================================
  # End LMM timing
  # ===========================================================================
  
  end_lmm <- Sys.time()
  
  TIME_LMM[kkk] <- as.numeric(
    difftime(
      end_lmm,
      start_lmm,
      units = "mins"
    )
  )
  
  # ===========================================================================
  # 2. LOCF LANDMARKING
  # ===========================================================================
  
  start_locf <- Sys.time()
  
  for (j in seq_along(S)) {
    
    s <- S[j]
    
    cat("  LOCF landmark =", s, "\n")
    
    # -------------------------------------------------------------------------
    # Subjects at risk
    # -------------------------------------------------------------------------
    
    train_at_risk <- subset(
      surv_train,
      survtime >= s
    )
    
    valid_at_risk <- subset(
      surv_valid,
      survtime >= s
    )
    
    if (
      nrow(train_at_risk) < 10 ||
      nrow(valid_at_risk) == 0
    ) {
      next
    }
    
    # -------------------------------------------------------------------------
    # Longitudinal history
    # -------------------------------------------------------------------------
    
    long_train_s <- subset(
      long_train,
      id %in% train_at_risk$id &
        obstime <= s
    )
    
    long_valid_s <- subset(
      long_valid,
      id %in% valid_at_risk$id &
        obstime <= s
    )
    
    # ========================================================================
    # Training LOCF
    # ========================================================================
    
    locf1_train <- get_locf(
      long_train_s,
      "Y1"
    )
    
    locf2_train <- get_locf(
      long_train_s,
      "Y2"
    )
    
    locf3_train <- get_locf(
      long_train_s,
      "Y3"
    )
    
    names(locf1_train)[2] <- "Y1_locf"
    names(locf2_train)[2] <- "Y2_locf"
    names(locf3_train)[2] <- "Y3_locf"
    
    train_lm <- merge(
      train_at_risk,
      locf1_train,
      by = "id",
      all.x = TRUE
    )
    
    train_lm <- merge(
      train_lm,
      locf2_train,
      by = "id",
      all.x = TRUE
    )
    
    train_lm <- merge(
      train_lm,
      locf3_train,
      by = "id",
      all.x = TRUE
    )
    
    train_lm <- train_lm[
      complete.cases(
        train_lm[, c(
          "Y1_locf",
          "Y2_locf",
          "Y3_locf"
        )]
      ),
    ]
    
    if (nrow(train_lm) < 10) {
      next
    }
    
    # ========================================================================
    # Landmark Cox model
    # ========================================================================
    
    cox_locf <- coxph(
      Surv(
        survtime - s,
        Event
      ) ~
        Y1_locf +
        Y2_locf +
        Y3_locf +
        w1 +
        w2,
      data = train_lm,
      ties = "breslow",
      x = TRUE
    )
    
    # ========================================================================
    # Validation LOCF
    # ========================================================================
    
    locf1_valid <- get_locf(
      long_valid_s,
      "Y1"
    )
    
    locf2_valid <- get_locf(
      long_valid_s,
      "Y2"
    )
    
    locf3_valid <- get_locf(
      long_valid_s,
      "Y3"
    )
    
    names(locf1_valid)[2] <- "Y1_locf"
    names(locf2_valid)[2] <- "Y2_locf"
    names(locf3_valid)[2] <- "Y3_locf"
    
    valid_lm <- merge(
      valid_at_risk,
      locf1_valid,
      by = "id",
      all.x = TRUE
    )
    
    valid_lm <- merge(
      valid_lm,
      locf2_valid,
      by = "id",
      all.x = TRUE
    )
    
    valid_lm <- merge(
      valid_lm,
      locf3_valid,
      by = "id",
      all.x = TRUE
    )
    
    # ========================================================================
    # Dynamic prediction
    # ========================================================================
    
    DP <- rep(
      NA_real_,
      nrow(valid_lm)
    )
    
    complete_valid <- complete.cases(
      valid_lm[, c(
        "Y1_locf",
        "Y2_locf",
        "Y3_locf"
      )]
    )
    
    if (sum(complete_valid) > 0) {
      
      valid_complete <- valid_lm[
        complete_valid,
      ]
      
      H0_tpred <- get_H0(
        cox_locf,
        t_pred
      )
      
      lp_valid <- predict(
        cox_locf,
        newdata = valid_complete,
        type = "lp"
      )
      
      surv_pred <- exp(
        -H0_tpred *
          exp(lp_valid)
      )
      
      DP_complete <- 1 - surv_pred
      
      DP[
        match(
          valid_complete$id,
          valid_lm$id
        )
      ] <- DP_complete
    }
    
    # ========================================================================
    # Store LOCF predictions
    # ========================================================================
    
    dp_results_locf[[j]] <- data.frame(
      id = valid_lm$id,
      Y1_locf = valid_lm$Y1_locf,
      Y2_locf = valid_lm$Y2_locf,
      Y3_locf = valid_lm$Y3_locf,
      DP = DP,
      landmark = s,
      horizon = s + t_pred
    )
  }
  
  # ===========================================================================
  # End LOCF timing
  # ===========================================================================
  
  end_locf <- Sys.time()
  
  TIME_LOCF[kkk] <- as.numeric(
    difftime(
      end_locf,
      start_locf,
      units = "mins"
    )
  )
  
  # ===========================================================================
  # Evaluate AUC and Brier Score
  # ===========================================================================
  
  for (j in seq_along(S)) {
    
    s <- S[j]
    
    # -------------------------------------------------------------------------
    # LMM
    # -------------------------------------------------------------------------
    
    if (!is.null(dp_results_lmm[[j]])) {
      
      pred_lmm <- dp_results_lmm[[j]]
      
      eval_lmm <- merge(
        surv_valid,
        pred_lmm[, c("id", "DP")],
        by = "id",
        all.x = FALSE
      )
      
      eval_lmm <- eval_lmm[
        is.finite(eval_lmm$DP),
      ]
      
      if (nrow(eval_lmm) > 0) {
        
        Crit_lmm <- tryCatch(
          Criteria(
            s = s,
            t = t_pred,
            Survt = eval_lmm$survtime,
            CR = eval_lmm$Event,
            P = eval_lmm$DP,
            cause = 1
          ),
          error = function(e) NULL
        )
        
        if (!is.null(Crit_lmm)) {
          
          AUC_LMM[kkk, j] <-
            as.numeric(
              Crit_lmm$Cri[1, 1]
            )
          
          BS_LMM[kkk, j] <-
            as.numeric(
              Crit_lmm$Cri[2, 1]
            )
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # LOCF
    # -------------------------------------------------------------------------
    
    if (!is.null(dp_results_locf[[j]])) {
      
      pred_locf <- dp_results_locf[[j]]
      
      eval_locf <- merge(
        surv_valid,
        pred_locf[, c("id", "DP")],
        by = "id",
        all.x = FALSE
      )
      
      eval_locf <- eval_locf[
        is.finite(eval_locf$DP),
      ]
      
      if (nrow(eval_locf) > 0) {
        
        Crit_locf <- tryCatch(
          Criteria(
            s = s,
            t = t_pred,
            Survt = eval_locf$survtime,
            CR = eval_locf$Event,
            P = eval_locf$DP,
            cause = 1
          ),
          error = function(e) NULL
        )
        
        if (!is.null(Crit_locf)) {
          
          AUC_LOCF[kkk, j] <-
            as.numeric(
              Crit_locf$Cri[1, 1]
            )
          
          BS_LOCF[kkk, j] <-
            as.numeric(
              Crit_locf$Cri[2, 1]
            )
        }
      }
    }
  }
  
  # ===========================================================================
  # Save all results from this simulation
  # ===========================================================================
  
  Results1[[kkk]] <- list(
    AUC_LMM = AUC_LMM[kkk, ],
    BS_LMM = BS_LMM[kkk, ],
    AUC_LOCF = AUC_LOCF[kkk, ],
    BS_LOCF = BS_LOCF[kkk, ],
    TIME_LMM = TIME_LMM[kkk],
    TIME_LOCF = TIME_LOCF[kkk],
    DP_LMM = dp_results_lmm,
    DP_LOCF = dp_results_locf
  )
  
  cat(
    "Time LMM  :",
    round(TIME_LMM[kkk], 3),
    "minutes\n"
  )
  
  cat(
    "Time LOCF :",
    round(TIME_LOCF[kkk], 3),
    "minutes\n"
  )
}

# =============================================================================
# FINAL SIMULATION SUMMARIES
# =============================================================================

# -----------------------------------------------------------------------------
# Mean and SD of AUC
# -----------------------------------------------------------------------------

Mean_AUC_LMM <- colMeans(
  AUC_LMM,
  na.rm = TRUE
)

SD_AUC_LMM <- apply(
  AUC_LMM,
  2,
  sd,
  na.rm = TRUE
)

Mean_AUC_LOCF <- colMeans(
  AUC_LOCF,
  na.rm = TRUE
)

SD_AUC_LOCF <- apply(
  AUC_LOCF,
  2,
  sd,
  na.rm = TRUE
)

# -----------------------------------------------------------------------------
# Mean and SD of Brier Score
# -----------------------------------------------------------------------------

Mean_BS_LMM <- colMeans(
  BS_LMM,
  na.rm = TRUE
)

SD_BS_LMM <- apply(
  BS_LMM,
  2,
  sd,
  na.rm = TRUE
)

Mean_BS_LOCF <- colMeans(
  BS_LOCF,
  na.rm = TRUE
)

SD_BS_LOCF <- apply(
  BS_LOCF,
  2,
  sd,
  na.rm = TRUE
)

# -----------------------------------------------------------------------------
# Mean and SD of computational time
# -----------------------------------------------------------------------------

Mean_TIME_LMM <- mean(
  TIME_LMM,
  na.rm = TRUE
)

SD_TIME_LMM <- sd(
  TIME_LMM,
  na.rm = TRUE
)

Mean_TIME_LOCF <- mean(
  TIME_LOCF,
  na.rm = TRUE
)

SD_TIME_LOCF <- sd(
  TIME_LOCF,
  na.rm = TRUE
)

# =============================================================================
# SUMMARY TABLE: LMM
# =============================================================================

Summary_LMM <- data.frame(
  Landmark = S,
  Mean_AUC = Mean_AUC_LMM,
  SD_AUC = SD_AUC_LMM,
  Mean_BS = Mean_BS_LMM,
  SD_BS = SD_BS_LMM
)

# =============================================================================
# SUMMARY TABLE: LOCF
# =============================================================================

Summary_LOCF <- data.frame(
  Landmark = S,
  Mean_AUC = Mean_AUC_LOCF,
  SD_AUC = SD_AUC_LOCF,
  Mean_BS = Mean_BS_LOCF,
  SD_BS = SD_BS_LOCF
)

# =============================================================================
# TIME SUMMARY
# =============================================================================

Summary_Time <- data.frame(
  Method = c(
    "LMM Landmarking",
    "LOCF Landmarking"
  ),
  Mean_Time_Minutes = c(
    Mean_TIME_LMM,
    Mean_TIME_LOCF
  ),
  SD_Time_Minutes = c(
    SD_TIME_LMM,
    SD_TIME_LOCF
  )
)



Summary_LMM
Summary_LOCF
Summary_Time


xtable::xtable(Summary_LMM, digits = 3)
xtable::xtable(Summary_LOCF, digits = 3)
xtable::xtable(Summary_Time, digits = 3)


 