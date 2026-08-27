rm(list = ls())

# ==============================================================================
# Libraries
# ==============================================================================

library(MASS)
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvtnorm)
library(refund)
library(DPCri)


# ==============================================================================
# Simulation parameters
# ==============================================================================

NN <- 100

Results1 <- vector("list", NN)


# ==============================================================================
# START SIMULATION
# ==============================================================================

for (kkk in 100:NN) {
  
  set.seed(kkk)
  
  cat("\n============================================================\n")
  cat("Simulation:", kkk, "of", NN, "\n")
  cat("============================================================\n")
  
  
  # ============================================================================
  # 1. Simulation parameters
  # ============================================================================
  
  nsample <- 1000
  
  followup <- 2
  
  t <- seq(
    0,
    followup,
    by = 0.005
  )
  
  
  # ============================================================================
  # 2. Association and fixed-effect parameters
  # ============================================================================
  
  alpha <- c(
    -0.5,
    0.3,
    -0.2
  )
  
  
  Beta1 <- c(
    -0.5,
    -0.5,
    0.5,
    0.5
  )
  
  
  Beta2 <- c(
    0.2,
    0.3,
    -0.3,
    0.4
  )
  
  
  Beta3 <- c(
    1.0,
    -0.2,
    0.1,
    0.6
  )
  
  
  sigma <- 1
  
  
  gamma_w <- c(
    0.2,
    -0.2
  )
  
  
  # ============================================================================
  # 3. Baseline covariates
  # ============================================================================
  
  id <- 1:nsample
  
  x1 <- rnorm(
    nsample,
    0,
    1
  )
  
  x2 <- rbinom(
    nsample,
    1,
    0.5
  )
  
  w1 <- rnorm(
    nsample,
    0,
    1
  )
  
  w2 <- rbinom(
    nsample,
    1,
    0.5
  )
  
  
  death <- rep(
    0L,
    nsample
  )
  
  time <- rep(
    NA_real_,
    nsample
  )
  
  
  # ============================================================================
  # 4. Random effects
  #
  # Three markers:
  #
  # Marker 1: random intercept + random slope
  # Marker 2: random intercept + random slope
  # Marker 3: random intercept + random slope
  #
  # Within-marker intercept/slope correlation = 0.6
  # Different markers independent
  # ============================================================================
  
  D_block <- matrix(
    c(
      1.0, 0.6,
      0.6, 1.0
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE
  )
  
  
  D <- matrix(
    0,
    nrow = 6,
    ncol = 6
  )
  
  
  D[1:2, 1:2] <- D_block
  D[3:4, 3:4] <- D_block
  D[5:6, 5:6] <- D_block
  
  
  u <- rmvnorm(
    n = nsample,
    mean = rep(0, 6),
    sigma = D
  )
  
  
  # ============================================================================
  # 5. Individual hazard function
  # ============================================================================
  
  haz_fun_i <- function(
    s,
    x1i,
    x2i,
    w1i,
    w2i,
    ui
  ) {
    
    # --------------------------------------------------------------------------
    # Marker 1
    # --------------------------------------------------------------------------
    
    mu1 <-
      Beta1[1] +
      Beta1[2] * x1i +
      Beta1[3] * x2i +
      Beta1[4] * s +
      ui[1] +
      ui[2] * s
    
    
    # --------------------------------------------------------------------------
    # Marker 2
    # --------------------------------------------------------------------------
    
    mu2 <-
      Beta2[1] +
      Beta2[2] * x1i +
      Beta2[3] * x2i +
      Beta2[4] * s +
      ui[3] +
      ui[4] * s
    
    
    # --------------------------------------------------------------------------
    # Marker 3
    # --------------------------------------------------------------------------
    
    mu3 <-
      Beta3[1] +
      Beta3[2] * x1i +
      Beta3[3] * x2i +
      Beta3[4] * s +
      ui[5] +
      ui[6] * s
    
    
    # --------------------------------------------------------------------------
    # Survival linear predictor
    # --------------------------------------------------------------------------
    
    eta <-
      gamma_w[1] * w1i +
      gamma_w[2] * w2i +
      alpha[1] * mu1 +
      alpha[2] * mu2 +
      alpha[3] * mu3
    
    
    exp(eta)
  }
  
  
  # ============================================================================
  # 6. Cumulative hazard
  # ============================================================================
  
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
  
  
  # ============================================================================
  # 7. Survival function
  # ============================================================================
  
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
  
  
  # ============================================================================
  # 8. Simulate event time
  # ============================================================================
  
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
    
    
    # Event occurs after follow-up
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
  
  
  # ============================================================================
  # 9. Generate longitudinal data
  # ============================================================================
  
  Y1 <- matrix(
    NA_real_,
    nrow = nsample,
    ncol = length(t)
  )
  
  Y2 <- matrix(
    NA_real_,
    nrow = nsample,
    ncol = length(t)
  )
  
  Y3 <- matrix(
    NA_real_,
    nrow = nsample,
    ncol = length(t)
  )
  
  
  for (i in 1:nsample) {
    
    # --------------------------------------------------------------------------
    # True event time
    # --------------------------------------------------------------------------
    
    evt <- sim_event_time_i(
      x1i = x1[i],
      x2i = x2[i],
      w1i = w1[i],
      w2i = w2[i],
      ui = u[i, ]
    )
    
    
    # --------------------------------------------------------------------------
    # Independent censoring
    # --------------------------------------------------------------------------
    
    td <- rexp(
      1,
      rate = 1
    )
    
    
    # --------------------------------------------------------------------------
    # Observed time
    # --------------------------------------------------------------------------
    
    obs_time <- min(
      evt,
      td,
      followup
    )
    
    
    time[i] <- obs_time
    
    
    # --------------------------------------------------------------------------
    # Event indicator
    # --------------------------------------------------------------------------
    
    death[i] <- as.integer(
      is.finite(evt) &&
        evt <= td &&
        evt <= followup
    )
    
    
    # --------------------------------------------------------------------------
    # Longitudinal measurements
    # --------------------------------------------------------------------------
    
    for (j in seq_along(t)) {
      
      if (t[j] <= obs_time) {
        
        # ======================================================================
        # Marker 1
        # ======================================================================
        
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
        
        
        # ======================================================================
        # Marker 2
        # ======================================================================
        
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
        
        
        # ======================================================================
        # Marker 3
        # ======================================================================
        
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
  
  
  # ============================================================================
  # 10. Convert to long format
  # ============================================================================
  
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
  
  
  # ============================================================================
  # 11. Survival data
  # ============================================================================
  
  surv.data <- data.frame(
    
    id = id,
    
    x1 = x1,
    x2 = x2,
    
    w1 = w1,
    w2 = w2,
    
    survtime = time,
    
    Event = death
    
  )
  
  
  # ============================================================================
  # 12. Train / validation split
  # ============================================================================
  
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
  
  
  # ============================================================================
  # 13. Diagnostics
  # ============================================================================
  
  cat(
    "Number of subjects:",
    nrow(surv.data),
    "\n"
  )
  
  cat(
    "Number of longitudinal observations:",
    nrow(long.data),
    "\n"
  )
  
  cat(
    "Average observations per subject:",
    mean(table(long.data$id)),
    "\n"
  )
  
  cat(
    "Median observations per subject:",
    median(table(long.data$id)),
    "\n"
  )
  
  cat(
    "Minimum observations:",
    min(table(long.data$id)),
    "\n"
  )
  
  cat(
    "Maximum observations:",
    max(table(long.data$id)),
    "\n"
  )
  
  cat(
    "Mean survival time:",
    mean(surv.data$survtime),
    "\n"
  )
  
  cat(
    "Event proportion:",
    mean(surv.data$Event),
    "\n"
  )
  
  
  # ============================================================================
  # 14. MFPCCox
  # ============================================================================
  
  start1 <- Sys.time()
  
  
  # ============================================================================
  # 14.1 Time grid and IDs
  # ============================================================================
  
  t_grid <- seq(
    0,
    followup,
    by = 0.005
  )
  
  
  train_ids <- sort(
    unique(long_train$id)
  )
  
  
  valid_ids <- sort(
    unique(long_valid$id)
  )
  
  
  # ============================================================================
  # 14.2 Training matrices
  # ============================================================================
  
  Ytrain_mat1 <- matrix(
    NA_real_,
    nrow = length(train_ids),
    ncol = length(t_grid)
  )
  
  Ytrain_mat2 <- matrix(
    NA_real_,
    nrow = length(train_ids),
    ncol = length(t_grid)
  )
  
  Ytrain_mat3 <- matrix(
    NA_real_,
    nrow = length(train_ids),
    ncol = length(t_grid)
  )
  
  
  for (i in seq_along(train_ids)) {
    
    iid <- train_ids[i]
    
    
    tmp <- long_train %>%
      filter(
        id == iid
      ) %>%
      group_by(
        obstime
      ) %>%
      summarise(
        Y1 = mean(Y1, na.rm = TRUE),
        Y2 = mean(Y2, na.rm = TRUE),
        Y3 = mean(Y3, na.rm = TRUE),
        .groups = "drop"
      )
    
    
    idx <- match(
      tmp$obstime,
      t_grid
    )
    
    
    Ytrain_mat1[i, idx] <- tmp$Y1
    Ytrain_mat2[i, idx] <- tmp$Y2
    Ytrain_mat3[i, idx] <- tmp$Y3
  }
  
  
  # ============================================================================
  # 14.3 FPCA
  # ============================================================================
  
  fpca1 <- fpca.sc(
    Y = Ytrain_mat1,
    argvals = t_grid,
    pve = 0.99,
    var = TRUE
  )
  
  
  fpca2 <- fpca.sc(
    Y = Ytrain_mat2,
    argvals = t_grid,
    pve = 0.99,
    var = TRUE
  )
  
  
  fpca3 <- fpca.sc(
    Y = Ytrain_mat3,
    argvals = t_grid,
    pve = 0.99,
    var = TRUE
  )
  
  
  K1 <- fpca1$npc
  K2 <- fpca2$npc
  K3 <- fpca3$npc
  
  
  cat(
    "Number of FPCs for marker 1:",
    K1,
    "\n"
  )
  
  cat(
    "Number of FPCs for marker 2:",
    K2,
    "\n"
  )
  
  cat(
    "Number of FPCs for marker 3:",
    K3,
    "\n"
  )
  
  
  # ============================================================================
  # 14.4 FPC names
  # ============================================================================
  
  rho1_names <- paste0(
    "rho1_",
    seq_len(K1)
  )
  
  rho2_names <- paste0(
    "rho2_",
    seq_len(K2)
  )
  
  rho3_names <- paste0(
    "rho3_",
    seq_len(K3)
  )
  
  
  all_rho_names <- c(
    rho1_names,
    rho2_names,
    rho3_names
  )
  
  
  # ============================================================================
  # 14.5 Training FPC scores
  # ============================================================================
  
  score_train <- data.frame(
    
    id = train_ids,
    
    fpca1$scores[
      ,
      seq_len(K1),
      drop = FALSE
    ],
    
    fpca2$scores[
      ,
      seq_len(K2),
      drop = FALSE
    ],
    
    fpca3$scores[
      ,
      seq_len(K3),
      drop = FALSE
    ]
    
  )
  
  
  names(score_train)[
    2:(K1 + 1)
  ] <- rho1_names
  
  
  names(score_train)[
    (K1 + 2):(K1 + K2 + 1)
  ] <- rho2_names
  
  
  names(score_train)[
    (K1 + K2 + 2):(K1 + K2 + K3 + 1)
  ] <- rho3_names
  
  
  # ============================================================================
  # 14.6 Survival training data + FPC scores
  # ============================================================================
  
  surv_train_mfp <- surv_train %>%
    select(
      id,
      w1,
      w2,
      survtime,
      Event
    ) %>%
    left_join(
      score_train,
      by = "id"
    )
  
  
  # ============================================================================
  # 14.7 Cox model
  # ============================================================================
  
  cox_formula <- as.formula(
    paste(
      "Surv(survtime, Event) ~ w1 + w2 +",
      paste(
        all_rho_names,
        collapse = " + "
      )
    )
  )
  
  
  cox_model <- coxph(
    formula = cox_formula,
    data = surv_train_mfp,
    x = TRUE
  )
  
  
  print(
    summary(cox_model)
  )
  
  
  # ============================================================================
  # 14.8 Project partially observed trajectories onto FPCs
  # ============================================================================
  
  project_scores <- function(
    fpca,
    newY
  ) {
    
    obs <- !is.na(newY)
    
    
    if (sum(obs) == 0) {
      
      return(
        rep(
          NA_real_,
          fpca$npc
        )
      )
    }
    
    
    Zcur <- fpca$efunctions[
      obs,
      seq_len(fpca$npc),
      drop = FALSE
    ]
    
    
    Ytilde <- newY[obs] -
      fpca$mu[obs]
    
    
    D.inv <- diag(
      1 / fpca$evalues[
        seq_len(fpca$npc)
      ]
    )
    
    
    A <- crossprod(Zcur) +
      fpca$sigma2 * D.inv
    
    
    b <- crossprod(
      Zcur,
      Ytilde
    )
    
    
    scores <- solve(
      A,
      b
    )
    
    
    as.numeric(scores)
  }
  
  
  # ============================================================================
  # 14.9 Baseline cumulative hazard
  # ============================================================================
  
  bh <- basehaz(
    cox_model,
    centered = FALSE
  )
  
  
  H0_at <- function(tt) {
    
    if (tt <= 0) {
      return(0)
    }
    
    
    j <- findInterval(
      tt,
      bh$time
    )
    
    
    if (j == 0) {
      return(0)
    }
    
    
    bh$hazard[j]
  }
  
  
  # ============================================================================
  # 14.10 Landmark times and prediction horizon
  # ============================================================================
  
  S <- c(
    0.25,
    0.50,
    0.75,
    1.00
  )
  
  
  t_pred <- 0.25
  
  
  results <- list()
  
  
  # ============================================================================
  # 14.11 Dynamic prediction
  # ============================================================================
  
  for (t0 in S) {
    
    cat(
      "\nLandmark:",
      t0,
      "\n"
    )
    
    
    # --------------------------------------------------------------------------
    # Subjects at risk at landmark
    # --------------------------------------------------------------------------
    
    risk_ids <- surv_valid$id[
      surv_valid$survtime > t0
    ]
    
    
    cat(
      "Number at risk:",
      length(risk_ids),
      "\n"
    )
    
    
    if (length(risk_ids) == 0) {
      next
    }
    
    
    # --------------------------------------------------------------------------
    # Predictions
    # --------------------------------------------------------------------------
    
    pred_list <- lapply(
      risk_ids,
      function(iid) {
        
        this_long <- long_valid %>%
          filter(
            id == iid,
            obstime <= t0
          )
        
        
        this_surv <- surv_valid %>%
          filter(
            id == iid
          )
        
        
        # ======================================================================
        # Partially observed trajectories
        # ======================================================================
        
        Y1test <- rep(
          NA_real_,
          length(t_grid)
        )
        
        Y2test <- rep(
          NA_real_,
          length(t_grid)
        )
        
        Y3test <- rep(
          NA_real_,
          length(t_grid)
        )
        
        
        if (nrow(this_long) > 0) {
          
          tmp <- this_long %>%
            group_by(
              obstime
            ) %>%
            summarise(
              Y1 = mean(Y1, na.rm = TRUE),
              Y2 = mean(Y2, na.rm = TRUE),
              Y3 = mean(Y3, na.rm = TRUE),
              .groups = "drop"
            )
          
          
          idx <- match(
            tmp$obstime,
            t_grid
          )
          
          
          Y1test[idx] <- tmp$Y1
          Y2test[idx] <- tmp$Y2
          Y3test[idx] <- tmp$Y3
        }
        
        
        # ======================================================================
        # FPC scores
        # ======================================================================
        
        sc1 <- project_scores(
          fpca1,
          Y1test
        )
        
        
        sc2 <- project_scores(
          fpca2,
          Y2test
        )
        
        
        sc3 <- project_scores(
          fpca3,
          Y3test
        )
        
        
        # ======================================================================
        # New data for Cox model
        # ======================================================================
        
        sc_df <- data.frame(
          matrix(
            c(
              sc1,
              sc2,
              sc3
            ),
            nrow = 1
          )
        )
        
        
        names(sc_df) <- all_rho_names
        
        
        newdata <- cbind(
          
          this_surv %>%
            select(
              id,
              w1,
              w2,
              survtime,
              Event
            ),
          
          sc_df
          
        )
        
        
        # ======================================================================
        # Uncentered linear predictor
        # ======================================================================
        
        lp <- predict(
          cox_model,
          newdata = newdata,
          type = "lp",
          reference = "zero"
        )
        
        
        # ======================================================================
        # Conditional survival from t0 to t0 + t_pred
        # ======================================================================
        
        H0_start <- H0_at(
          t0
        )
        
        
        H0_end <- H0_at(
          t0 + t_pred
        )
        
        
        delta_H0 <- max(
          H0_end - H0_start,
          0
        )
        
        
        surv_cond <- exp(
          -delta_H0 * exp(lp)
        )
        
        
        # ======================================================================
        # Predicted event probability
        # ======================================================================
        
        pred_risk <- 1 - surv_cond
        
        
        data.frame(
          
          id = iid,
          
          landmark = t0,
          
          horizon = t_pred,
          
          pred = as.numeric(
            pred_risk
          )
          
        )
      }
    )
    
    
    results[[as.character(t0)]] <- do.call(
      rbind,
      pred_list
    )
  }
  
  
  # ============================================================================
  # 14.12 Combine predictions
  # ============================================================================
  
  pred_all <- do.call(
    rbind,
    results
  )
  
  
  rownames(pred_all) <- NULL
  
  
  # ============================================================================
  # 14.13 Prediction diagnostics
  # ============================================================================
  
  cat(
    "\nPrediction summary:\n"
  )
  
  print(
    summary(
      pred_all$pred
    )
  )
  
  
  cat(
    "\nPrediction range:\n"
  )
  
  print(
    range(
      pred_all$pred,
      na.rm = TRUE
    )
  )
  
  
  # ============================================================================
  # 14.14 AUC and Brier Score
  # ============================================================================
  
  est3 <- matrix(
    NA_real_,
    nrow = length(S),
    ncol = 2
  )
  
  
  sd3 <- matrix(
    NA_real_,
    nrow = length(S),
    ncol = 2
  )
  
  
  rownames(est3) <- S
  rownames(sd3) <- S
  
  
  for (i in seq_along(S)) {
    
    s0 <- S[i]
    
    
    # --------------------------------------------------------------------------
    # Predictions at this landmark
    # --------------------------------------------------------------------------
    
    pred_i <- pred_all %>%
      filter(
        landmark == s0,
        horizon == t_pred
      )
    
    
    # --------------------------------------------------------------------------
    # Survival data corresponding exactly to prediction subjects
    # --------------------------------------------------------------------------
    
    surv_i_df <- surv_valid %>%
      filter(
        id %in% pred_i$id
      ) %>%
      arrange(
        match(
          id,
          pred_i$id
        )
      )
    
    
    if (
      nrow(pred_i) != nrow(surv_i_df)
    ) {
      
      stop(
        paste(
          "Mismatch between prediction and survival data at landmark",
          s0
        )
      )
    }
    
    
    # --------------------------------------------------------------------------
    # Criteria
    # --------------------------------------------------------------------------
    
    Crit <- Criteria(
      s = s0,
      t = t_pred,
      Survt = surv_i_df$survtime,
      CR = surv_i_df$Event,
      P = pred_i$pred,
      cause = 1
    )
    
    
    est3[i, ] <- Crit$Cri[, 1]
    sd3[i, ] <- Crit$Cri[, 2]
  }
  
  
  # ============================================================================
  # 14.15 Runtime
  # ============================================================================
  
  end1 <- Sys.time()
  
  
  time_MFPCCOX <- difftime(
    end1,
    start1,
    units = "mins"
  )
  
  
  cat(
    "\nMFPCCox computation time:",
    time_MFPCCOX,
    "minutes\n"
  )
  
  
  # ============================================================================
  # 14.16 Print results for this replicate
  # ============================================================================
  
  cat(
    "\nEstimated criteria:\n"
  )
  
  print(
    est3
  )
  
  
  cat(
    "\nStandard errors:\n"
  )
  
  print(
    sd3
  )
  
  
  # ============================================================================
  # 14.17 Save replicate results
  # ============================================================================
  
  Results1[[kkk]] <- list(
    
    MFP = est3,
    
    ctime = as.numeric(
      time_MFPCCOX
    )
    
  )
  
  
  cat(
    "\nCompleted simulation:",
    kkk,
    "\n"
  )
}


# ==============================================================================
# FINAL SIMULATION RESULTS
# ==============================================================================

Time1 <- numeric(NN)


Cr1 <- matrix(
  0,
  nrow = NN,
  ncol = 8
)


for (kkk in c(1:NN)[-c(92,93,97,98,99,100)]) {
  
  Time1[kkk] <-
    Results1[[kkk]]$ctime
  
  
  Cr1[kkk, ] <-
    as.numeric(
      Results1[[kkk]]$MFP
    )
}


# ==============================================================================
# Mean and SD
# ==============================================================================
CC=matrix(Cr1[-c(92,93,97,98,99,100),],nrow=94)
MFPCCox_results <- cbind(
  
  Mean = apply(
    CC,
    2,
    mean,
    na.rm = TRUE
  ),
  
  SD = apply(
    CC,
    2,
    sd,
    na.rm = TRUE
  )
)


rownames(
  MFPCCox_results
) <- 1:8


cat(
  "\n============================================================\n"
)

cat(
  "MFPCCox SIMULATION RESULTS\n"
)

cat(
  "============================================================\n"
)


print(
  MFPCCox_results
)


# ==============================================================================
# LaTeX table
# ==============================================================================

cat(
  "\nLaTeX table:\n"
)


print(
  xtable::xtable(
    MFPCCox_results,
    digits = 3
  )
)


# ==============================================================================
# Runtime
# ==============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "Runtime\n"
)

cat(
  "============================================================\n"
)


cat(
  "Mean time:",
  mean(
    Time1,
    na.rm = TRUE
  ),
  "minutes\n"
)


cat(
  "SD time:",
  sd(
    Time1,
    na.rm = TRUE
  ),
  "minutes\n"
)