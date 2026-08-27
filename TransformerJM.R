# ============================================================
# DeepSurv-LM simulation
# Landmarking + LME BLUP + DeepSurv + DPCri Criteria
# ============================================================

rm(list = ls())

library(MASS)
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvtnorm)
library(nlme)
library(DPCri)
library(torch)

# First time only:
# install.packages("torch")
# torch::install_torch()


# ============================================================
# 1. BLUP current value at landmark
# ============================================================

blup_current_value <- function(
    fit,
    dat_i,
    s,
    fixed_form  = ~ x1 + x2 + obstime,
    random_form = ~ obstime
) {
  
  beta <- nlme::fixef(fit)
  
  D <- as.matrix(
    nlme::getVarCov(
      fit,
      individuals = 1,
      type = "random.effects"
    )
  )
  
  sigma <- fit$sigma
  
  Xi <- model.matrix(
    fixed_form,
    data = dat_i
  )
  
  Zi <- model.matrix(
    random_form,
    data = dat_i
  )
  
  yname <- all.vars(formula(fit))[1]
  
  yi <- dat_i[[yname]]
  
  keep <- complete.cases(
    Xi,
    Zi,
    yi
  )
  
  Xi <- Xi[keep, , drop = FALSE]
  Zi <- Zi[keep, , drop = FALSE]
  yi <- yi[keep]
  
  if (length(yi) == 0) {
    return(NA_real_)
  }
  
  Vi <- Zi %*% D %*% t(Zi) +
    sigma^2 * diag(nrow(Zi))
  
  # ----------------------------------------------------------
  # BLUP of random effects
  # ----------------------------------------------------------
  
  bi <- D %*%
    t(Zi) %*%
    solve(
      Vi,
      yi - Xi %*% beta
    )
  
  # ----------------------------------------------------------
  # Prediction at landmark s
  # ----------------------------------------------------------
  
  newrow <- dat_i[1, , drop = FALSE]
  
  newrow$obstime <- s
  
  x_s <- model.matrix(
    fixed_form,
    data = newrow
  )
  
  z_s <- model.matrix(
    random_form,
    data = newrow
  )
  
  as.numeric(
    x_s %*% beta +
      z_s %*% bi
  )
}


# ============================================================
# 2. Construct landmark features
# ============================================================

build_landmark_features <- function(
    long_dat,
    surv_dat,
    s,
    fits,
    marker_names,
    fixed_form  = ~ x1 + x2 + obstime,
    random_form = ~ obstime
) {
  
  # ----------------------------------------------------------
  # Subjects at risk at landmark
  # ----------------------------------------------------------
  
  at_risk_ids <- surv_dat$id[
    surv_dat$survtime >= s
  ]
  
  # ----------------------------------------------------------
  # Subjects having longitudinal observations by s
  # ----------------------------------------------------------
  
  ids_long <- unique(
    long_dat$id[
      long_dat$obstime <= s
    ]
  )
  
  ids <- intersect(
    ids_long,
    at_risk_ids
  )
  
  if (length(ids) == 0) {
    return(NULL)
  }
  
  # ----------------------------------------------------------
  # Calculate BLUP current values
  # ----------------------------------------------------------
  
  feats <- lapply(
    ids,
    function(idi) {
      
      di <- long_dat[
        long_dat$id == idi &
          long_dat$obstime <= s,
        ,
        drop = FALSE
      ]
      
      if (nrow(di) == 0) {
        return(NULL)
      }
      
      yhat <- vapply(
        seq_along(fits),
        function(k) {
          
          blup_current_value(
            fit = fits[[k]],
            dat_i = di,
            s = s,
            fixed_form = fixed_form,
            random_form = random_form
          )
          
        },
        numeric(1)
      )
      
      names(yhat) <- paste0(
        "Yhat_",
        marker_names
      )
      
      wi <- surv_dat[
        surv_dat$id == idi,
        c("w1", "w2"),
        drop = FALSE
      ]
      
      if (nrow(wi) != 1) {
        return(NULL)
      }
      
      data.frame(
        id = idi,
        
        w1 = wi$w1,
        w2 = wi$w2,
        
        Yhat_Y1 = yhat[1],
        Yhat_Y2 = yhat[2],
        Yhat_Y3 = yhat[3],
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  feats <- feats[
    !vapply(
      feats,
      is.null,
      logical(1)
    )
  ]
  
  if (length(feats) == 0) {
    return(NULL)
  }
  
  feats <- do.call(
    rbind,
    feats
  )
  
  rownames(feats) <- NULL
  
  feats
}


# ============================================================
# 3. Safe LME fitting
# ============================================================

fit_lme_safe <- function(
    form,
    data,
    marker_name = ""
) {
  
  if (nrow(data) == 0) {
    stop(
      "No observations available for LME: ",
      marker_name
    )
  }
  
  if (length(unique(data$id)) < 5) {
    stop(
      "Too few subjects for LME: ",
      marker_name
    )
  }
  
  # ----------------------------------------------------------
  # Keep subjects with >= 2 observations
  # ----------------------------------------------------------
  
  nobs_id <- table(data$id)
  
  keep_ids <- names(
    nobs_id[
      nobs_id >= 2
    ]
  )
  
  data <- data[
    data$id %in% keep_ids,
    ,
    drop = FALSE
  ]
  
  if (length(unique(data$id)) < 5) {
    stop(
      "Too few subjects with >=2 observations for LME: ",
      marker_name
    )
  }
  
  # ----------------------------------------------------------
  # First attempt
  # ----------------------------------------------------------
  
  fit1 <- try(
    lme(
      fixed = form,
      random = ~ obstime | id,
      data = data,
      method = "REML",
      control = lmeControl(
        maxIter = 200,
        msMaxIter = 200,
        niterEM = 50,
        tolerance = 1e-6,
        msTol = 1e-7,
        returnObject = TRUE
      )
    ),
    silent = TRUE
  )
  
  if (!inherits(fit1, "try-error")) {
    return(fit1)
  }
  
  message(
    "First LME attempt failed for ",
    marker_name,
    ". Retrying..."
  )
  
  # ----------------------------------------------------------
  # Second attempt
  # ----------------------------------------------------------
  
  fit2 <- try(
    lme(
      fixed = form,
      random = ~ obstime | id,
      data = data,
      method = "REML",
      control = lmeControl(
        maxIter = 500,
        msMaxIter = 500,
        niterEM = 100,
        tolerance = 1e-5,
        msTol = 1e-6,
        returnObject = TRUE
      )
    ),
    silent = TRUE
  )
  
  if (!inherits(fit2, "try-error")) {
    return(fit2)
  }
  
  # ----------------------------------------------------------
  # Third attempt: ML
  # ----------------------------------------------------------
  
  message(
    "Second LME attempt failed for ",
    marker_name,
    ". Retrying using ML..."
  )
  
  fit3 <- try(
    lme(
      fixed = form,
      random = ~ obstime | id,
      data = data,
      method = "ML",
      control = lmeControl(
        maxIter = 500,
        msMaxIter = 500,
        niterEM = 100,
        tolerance = 1e-5,
        msTol = 1e-6,
        returnObject = TRUE
      )
    ),
    silent = TRUE
  )
  
  if (!inherits(fit3, "try-error")) {
    return(fit3)
  }
  
  stop(
    "\nLME failed for marker = ",
    marker_name,
    "\nNumber of subjects = ",
    length(unique(data$id)),
    "\nNumber of observations = ",
    nrow(data),
    "\n"
  )
}


# ============================================================
# 4. DeepSurv
# ============================================================

fit_deepsurv <- function(
    feat,
    surv_dat,
    s,
    hidden = c(8, 8),
    epochs = 300,
    lr = 1e-3,
    weight_decay = 1e-4,
    dropout = 0.1,
    seed = 1
) {
  
  if (is.null(feat) || nrow(feat) == 0) {
    stop(
      "No landmark features available at s = ",
      s
    )
  }
  
  torch_manual_seed(seed)
  
  # ----------------------------------------------------------
  # Add survival information
  # ----------------------------------------------------------
  
  surv_info <- surv_dat[
    ,
    c("id", "survtime", "Event"),
    drop = FALSE
  ]
  
  feat <- merge(
    feat,
    surv_info,
    by = "id",
    sort = FALSE
  )
  
  feat <- feat[
    feat$survtime >= s,
    ,
    drop = FALSE
  ]
  
  if (nrow(feat) < 10) {
    stop(
      "Too few subjects for DeepSurv at landmark s = ",
      s
    )
  }
  
  # ----------------------------------------------------------
  # Features
  # ----------------------------------------------------------
  
  xcols <- c(
    "Yhat_Y1",
    "Yhat_Y2",
    "Yhat_Y3",
    "w1",
    "w2"
  )
  
  if (!all(xcols %in% names(feat))) {
    stop(
      "Missing DeepSurv feature columns."
    )
  }
  
  Xraw <- as.matrix(
    feat[
      ,
      xcols,
      drop = FALSE
    ]
  )
  
  if (any(!is.finite(Xraw))) {
    stop(
      "Non-finite landmark features at s = ",
      s
    )
  }
  
  # ----------------------------------------------------------
  # Standardisation using TRAINING landmark data
  # ----------------------------------------------------------
  
  x_center <- apply(
    Xraw,
    2,
    mean
  )
  
  x_scale <- apply(
    Xraw,
    2,
    sd
  )
  
  x_scale[
    !is.finite(x_scale) |
      x_scale == 0
  ] <- 1
  
  Xs <- sweep(
    Xraw,
    2,
    x_center,
    "-"
  )
  
  Xs <- sweep(
    Xs,
    2,
    x_scale,
    "/"
  )
  
  # ----------------------------------------------------------
  # Survival time after landmark
  # ----------------------------------------------------------
  
  time_res <- feat$survtime - s
  
  event <- as.numeric(
    feat$Event
  )
  
  # ----------------------------------------------------------
  # Sort decreasingly
  #
  # For each event at T_i, the risk set is:
  #
  # T_j >= T_i
  #
  # Therefore log-cumsum-exp on decreasing time is correct.
  # ----------------------------------------------------------
  
  ord <- order(
    time_res,
    decreasing = TRUE
  )
  
  Xs <- Xs[
    ord,
    ,
    drop = FALSE
  ]
  
  time_o <- time_res[
    ord
  ]
  
  event_o <- event[
    ord
  ]
  
  feat_o <- feat[
    ord,
    ,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # Torch tensors
  # ----------------------------------------------------------
  
  Xt <- torch_tensor(
    Xs,
    dtype = torch_float()
  )
  
  event_t <- torch_tensor(
    event_o,
    dtype = torch_float()
  )
  
  # ----------------------------------------------------------
  # Network
  # ----------------------------------------------------------
  
  layers <- list()
  
  d_in <- ncol(Xs)
  
  for (h in hidden) {
    
    layers <- c(
      layers,
      nn_linear(
        d_in,
        h
      ),
      nn_selu(),
      nn_dropout(
        p = dropout
      )
    )
    
    d_in <- h
  }
  
  layers <- c(
    layers,
    nn_linear(
      d_in,
      1
    )
  )
  
  net <- do.call(
    nn_sequential,
    layers
  )
  
  opt <- optim_adam(
    net$parameters,
    lr = lr,
    weight_decay = weight_decay
  )
  
  # ----------------------------------------------------------
  # Training
  # ----------------------------------------------------------
  
  loss_history <- numeric(
    epochs
  )
  
  for (ep in seq_len(epochs)) {
    
    net$train()
    
    opt$zero_grad()
    
    risk <- net(
      Xt
    )$squeeze(2)
    
    log_cumsum <- torch_logcumsumexp(
      risk,
      dim = 1
    )
    
    loss <- -torch_mean(
      (risk - log_cumsum) *
        event_t
    )
    
    loss$backward()
    
    opt$step()
    
    loss_history[ep] <- loss$item()
  }
  
  # ----------------------------------------------------------
  # Training risk
  # ----------------------------------------------------------
  
  net$eval()
  
  risk_train <- with_no_grad({
    
    as.numeric(
      net(
        Xt
      )$squeeze(2)
    )
    
  })
  
  # ----------------------------------------------------------
  # Diagnostic
  # ----------------------------------------------------------
  
  diagnostic_train <- data.frame(
    id = feat_o$id,
    survtime = feat_o$survtime,
    Event = feat_o$Event,
    
    Yhat_Y1 = feat_o$Yhat_Y1,
    Yhat_Y2 = feat_o$Yhat_Y2,
    Yhat_Y3 = feat_o$Yhat_Y3,
    
    w1 = feat_o$w1,
    w2 = feat_o$w2,
    
    risk = risk_train
  )
  
  list(
    
    net = net,
    
    x_center = x_center,
    x_scale = x_scale,
    xcols = xcols,
    
    time_sorted = time_o,
    event_sorted = event_o,
    
    risk_train = risk_train,
    
    diagnostic_train = diagnostic_train,
    
    loss_history = loss_history
  )
}


# ============================================================
# 5. Breslow estimator
# ============================================================

compute_breslow <- function(
    time,
    event,
    risk_scores
) {
  
  event_times <- sort(
    unique(
      time[event == 1]
    )
  )
  
  if (length(event_times) == 0) {
    
    return(
      data.frame(
        time = numeric(0),
        dLambda = numeric(0)
      )
    )
  }
  
  out <- lapply(
    event_times,
    function(tt) {
      
      d_k <- sum(
        event == 1 &
          time == tt
      )
      
      risk_set <- time >= tt
      
      denom <- sum(
        exp(
          risk_scores[risk_set]
        )
      )
      
      data.frame(
        time = tt,
        dLambda = d_k / denom
      )
    }
  )
  
  do.call(
    rbind,
    out
  )
}


# ============================================================
# 6. Dynamic event probability
# ============================================================

predict_event_prob <- function(
    risk_i,
    breslow_tab,
    h
) {
  
  if (nrow(breslow_tab) == 0) {
    return(0)
  }
  
  idx <- breslow_tab$time <= h
  
  if (!any(idx)) {
    return(0)
  }
  
  cumhaz <- exp(
    risk_i
  ) *
    sum(
      breslow_tab$dLambda[idx]
    )
  
  cumhaz <- max(
    0,
    min(
      cumhaz,
      700
    )
  )
  
  1 - exp(
    -cumhaz
  )
}


# ============================================================
# 7. Full DeepSurv-LM dynamic prediction
# ============================================================

deepsurv_lm_dynpred <- function(
    long_train,
    surv_train,
    long_new,
    surv_new,
    s,
    horizon,
    marker_names,
    fixed_form = ~ x1 + x2 + obstime,
    random_form = ~ obstime,
    hidden = c(8, 8),
    epochs = 300,
    lr = 1e-3,
    weight_decay = 1e-4,
    dropout = 0.1,
    seed = 1
) {
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "Landmark s = ",
    s,
    "\n",
    sep = ""
  )
  
  # ==========================================================
  # TRAINING DATA
  # ==========================================================
  
  train_ids_s <- surv_train$id[
    surv_train$survtime >= s
  ]
  
  dtrain_s <- long_train[
    long_train$id %in% train_ids_s &
      long_train$obstime <= s,
    ,
    drop = FALSE
  ]
  
  cat(
    "Training subjects at risk = ",
    length(
      unique(
        dtrain_s$id
      )
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Training longitudinal observations = ",
    nrow(dtrain_s),
    "\n",
    sep = ""
  )
  
  if (nrow(dtrain_s) == 0) {
    stop(
      "No training longitudinal observations at s = ",
      s
    )
  }
  
  # ----------------------------------------------------------
  # Diagnostics
  # ----------------------------------------------------------
  
  cat(
    "\nObservations per subject:\n"
  )
  
  print(
    summary(
      table(
        dtrain_s$id
      )
    )
  )
  
  cat(
    "\nMarker SD:\n"
  )
  
  print(
    sapply(
      dtrain_s[
        ,
        marker_names,
        drop = FALSE
      ],
      sd,
      na.rm = TRUE
    )
  )
  
  # ==========================================================
  # Keep subjects with >= 2 observations
  # ==========================================================
  
  nobs_id <- table(
    dtrain_s$id
  )
  
  keep_ids <- as.integer(
    names(
      nobs_id[
        nobs_id >= 2
      ]
    )
  )
  
  dtrain_s <- dtrain_s[
    dtrain_s$id %in% keep_ids,
    ,
    drop = FALSE
  ]
  
  cat(
    "Subjects retained for LME = ",
    length(
      unique(
        dtrain_s$id
      )
    ),
    "\n",
    sep = ""
  )
  
  # ==========================================================
  # LME for each marker
  # ==========================================================
  
  fits <- lapply(
    marker_names,
    function(mk) {
      
      cat(
        "\nFitting LME for ",
        mk,
        " at s = ",
        s,
        "\n",
        sep = ""
      )
      
      form <- as.formula(
        paste(
          mk,
          "~ x1 + x2 + obstime"
        )
      )
      
      fit <- fit_lme_safe(
        form = form,
        data = dtrain_s,
        marker_name = mk
      )
      
      cat(
        "SUCCESS: ",
        mk,
        "\n",
        sep = ""
      )
      
      fit
    }
  )
  
  # ==========================================================
  # TRAIN landmark features
  # ==========================================================
  
  feat_train <- build_landmark_features(
    long_dat = dtrain_s,
    surv_dat = surv_train,
    s = s,
    fits = fits,
    marker_names = marker_names,
    fixed_form = fixed_form,
    random_form = random_form
  )
  
  if (
    is.null(feat_train) ||
    nrow(feat_train) == 0
  ) {
    
    stop(
      "No landmark training features at s = ",
      s
    )
  }
  
  cat(
    "\nLandmark training subjects = ",
    nrow(feat_train),
    "\n",
    sep = ""
  )
  
  # ==========================================================
  # DeepSurv
  # ==========================================================
  
  dsfit <- fit_deepsurv(
    feat = feat_train,
    surv_dat = surv_train,
    s = s,
    hidden = hidden,
    epochs = epochs,
    lr = lr,
    weight_decay = weight_decay,
    dropout = dropout,
    seed = seed
  )
  
  # ==========================================================
  # Breslow
  # ==========================================================
  
  breslow_tab <- compute_breslow(
    time = dsfit$time_sorted,
    event = dsfit$event_sorted,
    risk_scores = dsfit$risk_train
  )
  
  # ==========================================================
  # VALIDATION DATA
  # ==========================================================
  
  new_ids_s <- surv_new$id[
    surv_new$survtime >= s
  ]
  
  dnew_s <- long_new[
    long_new$id %in% new_ids_s &
      long_new$obstime <= s,
    ,
    drop = FALSE
  ]
  
  feat_new <- build_landmark_features(
    long_dat = dnew_s,
    surv_dat = surv_new,
    s = s,
    fits = fits,
    marker_names = marker_names,
    fixed_form = fixed_form,
    random_form = random_form
  )
  
  if (
    is.null(feat_new) ||
    nrow(feat_new) == 0
  ) {
    
    warning(
      "No validation subjects available at s = ",
      s
    )
    
    return(
      list(
        
        predictions = data.frame(
          id = integer(0),
          pred = numeric(0)
        ),
        
        true_vs_pred = data.frame(),
        
        landmark_features_train = feat_train,
        
        landmark_features_new = data.frame(),
        
        deep_surv = dsfit,
        
        breslow = breslow_tab,
        
        fits = fits
      )
    )
  }
  
  # ==========================================================
  # VALIDATION FEATURES
  # ==========================================================
  
  Xnew <- as.matrix(
    feat_new[
      ,
      dsfit$xcols,
      drop = FALSE
    ]
  )
  
  if (any(!is.finite(Xnew))) {
    stop(
      "Non-finite validation features at s = ",
      s
    )
  }
  
  # ----------------------------------------------------------
  # Use TRAINING scaling
  # ----------------------------------------------------------
  
  Xnew <- sweep(
    Xnew,
    2,
    dsfit$x_center,
    "-"
  )
  
  Xnew <- sweep(
    Xnew,
    2,
    dsfit$x_scale,
    "/"
  )
  
  Xnew_t <- torch_tensor(
    Xnew,
    dtype = torch_float()
  )
  
  # ==========================================================
  # VALIDATION RISK
  # ==========================================================
  
  risk_new <- with_no_grad({
    
    as.numeric(
      dsfit$net(
        Xnew_t
      )$squeeze(2)
    )
    
  })
  
  # ==========================================================
  # Dynamic event probability
  # ==========================================================
  
  pred <- vapply(
    risk_new,
    predict_event_prob,
    numeric(1),
    breslow_tab = breslow_tab,
    h = horizon
  )
  
  pred <- pmin(
    pmax(
      pred,
      0
    ),
    1
  )
  
  # ==========================================================
  # TRUE latent current value
  #
  # Since the simulation stores Y*_k at the exact observation
  # times, obtain the latent value at landmark s.
  # ==========================================================
  
  true_values <- lapply(
    feat_new$id,
    function(ii) {
      
      rr <- long_new[
        long_new$id == ii &
          abs(long_new$obstime - s) < 1e-10,
        ,
        drop = FALSE
      ]
      
      if (nrow(rr) == 0) {
        
        # If exact s is not found, use the last available
        # latent value before s.
        
        rr <- long_new[
          long_new$id == ii &
            long_new$obstime <= s,
          ,
          drop = FALSE
        ]
        
        if (nrow(rr) == 0) {
          return(
            c(
              Y1_true = NA_real_,
              Y2_true = NA_real_,
              Y3_true = NA_real_
            )
          )
        }
        
        rr <- rr[
          which.max(rr$obstime),
          ,
          drop = FALSE
        ]
        
      } else {
        
        rr <- rr[1, , drop = FALSE]
        
      }
      
      c(
        Y1_true = rr$Y1_true,
        Y2_true = rr$Y2_true,
        Y3_true = rr$Y3_true
      )
    }
  )
  
  true_values <- do.call(
    rbind,
    true_values
  )
  
  true_values <- as.data.frame(
    true_values
  )
  
  # ==========================================================
  # TRUE vs BLUP + survival prediction
  # ==========================================================
  
  true_vs_pred <- data.frame(
    
    id = feat_new$id,
    
    Y1_true = true_values$Y1_true,
    Y2_true = true_values$Y2_true,
    Y3_true = true_values$Y3_true,
    
    Yhat_Y1 = feat_new$Yhat_Y1,
    Yhat_Y2 = feat_new$Yhat_Y2,
    Yhat_Y3 = feat_new$Yhat_Y3,
    
    w1 = feat_new$w1,
    w2 = feat_new$w2,
    
    risk = risk_new,
    
    pred = pred
  )
  
  # ==========================================================
  # Return everything
  # ==========================================================
  
  list(
    
    predictions = data.frame(
      id = feat_new$id,
      pred = pred
    ),
    
    true_vs_pred = true_vs_pred,
    
    landmark_features_train = feat_train,
    
    landmark_features_new = feat_new,
    
    deep_surv = dsfit,
    
    breslow = breslow_tab,
    
    fits = fits
  )
}


# ============================================================
# 8. Simulation parameters
# ============================================================

NN <- 100

Results1 <- vector(
  "list",
  NN
)


# ============================================================
# 9. Monte Carlo simulation
# ============================================================
for (kkk in seq_len(NN)) {
  
  set.seed(
    kkk
  )
  
  cat(
    "\n\n####################################\n"
  )
  
  cat(
    "Monte Carlo replication = ",
    kkk,
    "\n",
    sep = ""
  )
  
  cat(
    "####################################\n"
  )
  
  nsample <- 1000
  
  followup <- 2
  
  t <- seq(
    0,
    followup,
    by = 0.005
  )
  
  # ==========================================================
  # Association parameters
  # ==========================================================
  
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
  
  # ==========================================================
  # Covariates
  # ==========================================================
  
  id <- 1:nsample
  
  x1 <- rnorm(
    nsample
  )
  
  x2 <- rbinom(
    nsample,
    1,
    0.5
  )
  
  w1 <- rnorm(
    nsample
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
  
  # ==========================================================
  # Random effects
  # ==========================================================
  
  D_block <- matrix(
    c(
      1,
      0.6,
      0.6,
      1
    ),
    nrow = 2,
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
  
  # ==========================================================
  # Individual hazard
  # ==========================================================
  
  haz_fun_i <- function(
    s,
    x1i,
    x2i,
    w1i,
    w2i,
    ui
  ) {
    
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
    
    exp(
      eta
    )
  }
  
  # ==========================================================
  # Cumulative hazard
  # ==========================================================
  
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
      function(ss) {
        
        haz_fun_i(
          ss,
          x1i,
          x2i,
          w1i,
          w2i,
          ui
        )
        
      },
      lower = 0,
      upper = T,
      rel.tol = 1e-6
    )$value
  }
  
  # ==========================================================
  # Survival
  # ==========================================================
  
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
        T,
        x1i,
        x2i,
        w1i,
        w2i,
        ui
      )
    )
  }
  
  # ==========================================================
  # Event time
  # ==========================================================
  
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
      Tmax,
      x1i,
      x2i,
      w1i,
      w2i,
      ui
    )
    
    if (S_Tmax > U) {
      return(Inf)
    }
    
    f <- function(T) {
      
      surv_i(
        T,
        x1i,
        x2i,
        w1i,
        w2i,
        ui
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
  
  # ==========================================================
  # Generate longitudinal data
  # ==========================================================
  
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
  
  # True latent trajectories
  
  Y1_true <- matrix(
    NA_real_,
    nsample,
    length(t)
  )
  
  Y2_true <- matrix(
    NA_real_,
    nsample,
    length(t)
  )
  
  Y3_true <- matrix(
    NA_real_,
    nsample,
    length(t)
  )
  
  # ==========================================================
  # Generate subjects
  # ==========================================================
  
  for (i in seq_len(nsample)) {
    
    evt <- sim_event_time_i(
      x1[i],
      x2[i],
      w1[i],
      w2[i],
      u[i, ]
    )
    
    td <- rexp(
      1,
      rate = 1
    )
    
    obs_time <- min(
      evt,
      td,
      followup
    )
    
    time[i] <- obs_time
    
    death[i] <- as.integer(
      is.finite(evt) &&
        evt <= td &&
        evt <= followup
    )
    
    for (j in seq_along(t)) {
      
      if (t[j] <= obs_time) {
        
        # ----------------------------------------------------
        # TRUE latent marker 1
        # ----------------------------------------------------
        
        Y1_true[i, j] <-
          Beta1[1] +
          Beta1[2] * x1[i] +
          Beta1[3] * x2[i] +
          Beta1[4] * t[j] +
          u[i, 1] +
          u[i, 2] * t[j]
        
        # ----------------------------------------------------
        # TRUE latent marker 2
        # ----------------------------------------------------
        
        Y2_true[i, j] <-
          Beta2[1] +
          Beta2[2] * x1[i] +
          Beta2[3] * x2[i] +
          Beta2[4] * t[j] +
          u[i, 3] +
          u[i, 4] * t[j]
        
        # ----------------------------------------------------
        # TRUE latent marker 3
        # ----------------------------------------------------
        
        Y3_true[i, j] <-
          Beta3[1] +
          Beta3[2] * x1[i] +
          Beta3[3] * x2[i] +
          Beta3[4] * t[j] +
          u[i, 5] +
          u[i, 6] * t[j]
        
        # ----------------------------------------------------
        # Observed marker 1
        # ----------------------------------------------------
        
        Y1[i, j] <-
          Y1_true[i, j] +
          rnorm(
            1,
            0,
            sigma
          )
        
        # ----------------------------------------------------
        # Observed marker 2
        # ----------------------------------------------------
        
        Y2[i, j] <-
          Y2_true[i, j] +
          rnorm(
            1,
            0,
            sigma
          )
        
        # ----------------------------------------------------
        # Observed marker 3
        # ----------------------------------------------------
        
        Y3[i, j] <-
          Y3_true[i, j] +
          rnorm(
            1,
            0,
            sigma
          )
      }
    }
  }
  
  # ==========================================================
  # Long format
  # ==========================================================
  
  long.list <- vector(
    "list",
    nsample
  )
  
  for (i in seq_len(nsample)) {
    
    idxs <- which(
      !is.na(
        Y1[i, ]
      )
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
      Y3 = Y3[i, idxs],
      
      Y1_true = Y1_true[i, idxs],
      Y2_true = Y2_true[i, idxs],
      Y3_true = Y3_true[i, idxs]
    )
  }
  
  long.data <- do.call(
    rbind,
    long.list
  )
  
  surv.data <- data.frame(
    
    id = id,
    
    x1 = x1,
    x2 = x2,
    
    w1 = w1,
    w2 = w2,
    
    survtime = time,
    Event = death
  )
  
  # ==========================================================
  # Train / validation split
  # ==========================================================
  
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
  
  # ==========================================================
  # DeepSurv-LM
  # ==========================================================
  
  start1 <- Sys.time()
  
  S <- c(
    0.25,
    0.50,
    0.75,
    1.00
  )
  
  t_pred <- 0.25
  
  marker_names <- c(
    "Y1",
    "Y2",
    "Y3"
  )
  
  preds_by_S <- list()
  
  true_vs_pred_by_S <- list()
  
  # ==========================================================
  # Landmark loop
  # ==========================================================
  
  for (s in S) {
    
    fit_s <- deepsurv_lm_dynpred(
      
      long_train = long_train,
      surv_train = surv_train,
      
      long_new = long_valid,
      surv_new = surv_valid,
      
      s = s,
      
      horizon = t_pred,
      
      marker_names = marker_names,
      
      hidden = c(8, 8),
      
      epochs = 300,
      
      lr = 1e-3,
      
      weight_decay = 1e-4,
      
      dropout = 0.1,
      
      seed = kkk
    )
    
    # ========================================================
    # IMPORTANT FIX:
    #
    # DO NOT put fit_s$predictions inside another list.
    #
    # Before:
    # list(fit_s$predictions)
    #
    # Now:
    # fit_s$predictions
    # ========================================================
    
    preds_by_S[[as.character(s)]] <-
      fit_s$predictions
    
    true_vs_pred_by_S[[as.character(s)]] <-
      fit_s$true_vs_pred
  }
  
  # ==========================================================
  # Extract predictions
  # ==========================================================
  
  preds_025 <- preds_by_S[["0.25"]]
  
  preds_050 <- preds_by_S[["0.5"]]
  
  preds_075 <- preds_by_S[["0.75"]]
  
  preds_100 <- preds_by_S[["1"]]
  
  # ==========================================================
  # Check that predictions are DATA FRAMES
  # ==========================================================
  
  cat(
    "\nPrediction object classes:\n"
  )
  
  print(
    c(
      s025 = class(preds_025)[1],
      s050 = class(preds_050)[1],
      s075 = class(preds_075)[1],
      s100 = class(preds_100)[1]
    )
  )
  
  # ==========================================================
  # Criteria
  # ==========================================================
  
  est2 <- matrix(
    NA_real_,
    nrow = length(S),
    ncol = 2
  )
  
  colnames(est2) <- c(
    "AUC",
    "BS"
  )
  
  rownames(est2) <- paste0(
    "s",
    S
  )
  
  # ==========================================================
  # Criteria helper
  # ==========================================================
  
  get_criteria <- function(
    pred_dat,
    s,
    horizon
  ) {
    
    # --------------------------------------------------------
    # Safety checks
    # --------------------------------------------------------
    
    if (!is.data.frame(pred_dat)) {
      
      stop(
        "pred_dat is not a data.frame at s = ",
        s,
        ". Class = ",
        paste(
          class(pred_dat),
          collapse = ", "
        )
      )
    }
    
    if (!all(
      c("id", "pred") %in%
      names(pred_dat)
    )) {
      
      stop(
        "Prediction data must contain id and pred."
      )
    }
    
    # --------------------------------------------------------
    # Match survival information
    # --------------------------------------------------------
    
    idx <- match(
      pred_dat$id,
      surv_valid$id
    )
    
    surv_tmp <- surv_valid[
      idx,
      ,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # Complete cases
    # --------------------------------------------------------
    
    keep <- complete.cases(
      surv_tmp$survtime,
      surv_tmp$Event,
      pred_dat$pred
    )
    
    surv_tmp <- surv_tmp[
      keep,
      ,
      drop = FALSE
    ]
    
    p_tmp <- pred_dat$pred[
      keep
    ]
    
    if (nrow(surv_tmp) == 0) {
      
      warning(
        "No complete cases for Criteria at s = ",
        s
      )
      
      return(
        c(
          AUC = NA_real_,
          BS = NA_real_
        )
      )
    }
    
    # --------------------------------------------------------
    # Check prediction range
    # --------------------------------------------------------
    
    p_tmp <- pmin(
      pmax(
        as.numeric(p_tmp),
        0
      ),
      1
    )
    
    # --------------------------------------------------------
    # Criteria
    # --------------------------------------------------------
    
    Crit <- try(
      Criteria(
        s = s,
        t = horizon,
        Survt = surv_tmp$survtime,
        CR = surv_tmp$Event,
        P = p_tmp,
        cause = 1
      ),
      silent = TRUE
    )
    
    if (inherits(
      Crit,
      "try-error"
    )) {
      
      stop(
        "\nDPCri::Criteria failed at landmark s = ",
        s,
        "\n\n",
        as.character(Crit)
      )
    }
    
    # --------------------------------------------------------
    # Inspect output once
    # --------------------------------------------------------
    
    if (s == S[1]) {
      
      cat(
        "\n====================================\n"
      )
      
      cat(
        "DPCri::Criteria output\n"
      )
      
      cat(
        "====================================\n"
      )
      
      print(
        Crit
      )
      
      cat(
        "\nStructure of Criteria output:\n"
      )
      
      str(
        Crit
      )
    }
    
    # --------------------------------------------------------
    # Extract two criteria
    # --------------------------------------------------------
    
    if (
      is.null(Crit$Cri)
    ) {
      
      stop(
        "Criteria() did not return $Cri at s = ",
        s
      )
    }
    
    cri_mat <- as.matrix(
      Crit$Cri
    )
    
    # --------------------------------------------------------
    # The package should return two criteria.
    # We explicitly take the first two values.
    # --------------------------------------------------------
    
    vals <- as.numeric(
      cri_mat[, 1]
    )
    
    if (length(vals) < 2) {
      
      # Some versions/structures may have the criteria
      # as columns rather than rows.
      
      vals2 <- as.numeric(
        cri_mat[1, ]
      )
      
      if (length(vals2) >= 2) {
        vals <- vals2
      }
    }
    
    if (length(vals) < 2) {
      
      stop(
        "\nCould not extract two criteria from Criteria(). ",
        "\nPlease inspect the printed Crit object above."
      )
    }
    
    c(
      AUC = vals[1],
      BS = vals[2]
    )
  }
  
  
  # ==========================================================
  # Calculate criteria at four landmarks
  # ==========================================================
  
  est2[1, ] <- get_criteria(
    pred_dat = preds_025,
    s = 0.25,
    horizon = t_pred
  )
  
  est2[2, ] <- get_criteria(
    pred_dat = preds_050,
    s = 0.50,
    horizon = t_pred
  )
  
  est2[3, ] <- get_criteria(
    pred_dat = preds_075,
    s = 0.75,
    horizon = t_pred
  )
  
  est2[4, ] <- get_criteria(
    pred_dat = preds_100,
    s = 1.00,
    horizon = t_pred
  )
  
  # ==========================================================
  # PRINT criteria for this replication
  # ==========================================================
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "Criteria - replication ",
    kkk,
    "\n",
    sep = ""
  )
  
  cat(
    "====================================\n"
  )
  
  print(
    est2
  )
  
  # ==========================================================
  # Runtime
  # ==========================================================
  
  end1 <- Sys.time()
  
  TimeMulti <- difftime(
    end1,
    start1,
    units = "mins"
  )
  
 
  
  # ==========================================================
  # Store everything
  # ==========================================================
  
  Results1[[kkk]] <- list(
    
    est2 = est2,
    
    TimeMulti = TimeMulti,
    
    true_vs_pred =
      true_vs_pred_by_S,
    
    predictions =
      preds_by_S
  )
  
  cat(
    "\nCompleted replication ",
    kkk,
    "\n",
    sep = ""
  )
}


# ============================================================
# 10. Final results
# ============================================================

Time1 <- numeric(
  NN
)

# ------------------------------------------------------------
# 4 landmarks x 2 criteria = 8 columns
# ------------------------------------------------------------

Cr1 <- matrix(
  NA_real_,
  nrow = NN,
  ncol = 8
)

colnames(Cr1) <- c(
  "AUC_s025",
  "AUC_s050",
  "AUC_s075",
  "AUC_s100",
  
  "BS_s025",
  "BS_s050",
  "BS_s075",
  "BS_s100"
)

for (kkk in seq_len(NN)) {
  
  Time1[kkk] <-
    as.numeric(
      Results1[[kkk]]$TimeMulti
    )
  
  # ----------------------------------------------------------
  # IMPORTANT:
  #
  # as.numeric(est2) is column-major:
  #
  # AUC s=.25
  # AUC s=.50
  # AUC s=.75
  # AUC s=1
  # BS  s=.25
  # BS  s=.50
  # BS  s=.75
  # BS  s=1
  # ----------------------------------------------------------
  
  Cr1[kkk, ] <-
    as.numeric(
      Results1[[kkk]]$est2
    )
}


# ============================================================
# 11. Summary
# ============================================================

MeanCriteria <- apply(
  Cr1,
  2,
  mean,
  na.rm = TRUE
)

SDCriteria <- apply(
  Cr1,
  2,
  sd,
  na.rm = TRUE
)


# ============================================================
# 12. Final table
# ============================================================

FinalResults <- data.frame(
  
  Landmark = c(
    "0.25",
    "0.50",
    "0.75",
    "1.00"
  ),
  
  AUC_Mean = MeanCriteria[
    1:4
  ],
  
  AUC_SD = SDCriteria[
    1:4
  ],
  
  BS_Mean = MeanCriteria[
    5:8
  ],
  
  BS_SD = SDCriteria[
    5:8
  ]
)

cat(
  "\n\n====================================\n"
)

cat(
  "FINAL RESULTS\n"
)

cat(
  "====================================\n"
)

print(
  FinalResults
)


# ============================================================
# 13. Raw 8-column result
# ============================================================

cat(
  "\n\n====================================\n"
)

cat(
  "Mean and SD for all 8 criteria\n"
)

cat(
  "====================================\n"
)

FinalResults8 <- cbind(
  Mean = MeanCriteria,
  SD = SDCriteria
)

print(
  FinalResults8
)


# ============================================================
# 14. Runtime
# ============================================================

cat(
  "\n\n====================================\n"
)

cat(
  "Mean Runtime (minutes)\n"
)

cat(
  "====================================\n"
)

print(
  mean(
    Time1,
    na.rm = TRUE
  )
)


cat(
  "\n\n====================================\n"
)

cat(
  "SD Runtime (minutes)\n"
)

cat(
  "====================================\n"
)

print(
  sd(
    Time1,
    na.rm = TRUE
  )
)


# ============================================================
# 15. Inspect replication 1
# ============================================================

cat(
  "\n\n====================================\n"
)

cat(
  "Diagnostics: replication 1\n"
)

cat(
  "====================================\n"
)


z025 <-
  Results1[[1]]$true_vs_pred[["0.25"]]

z050 <-
  Results1[[1]]$true_vs_pred[["0.5"]]

z075 <-
  Results1[[1]]$true_vs_pred[["0.75"]]

z100 <-
  Results1[[1]]$true_vs_pred[["1"]]


cat(
  "\nLandmark 0.25:\n"
)

print(
  head(
    z025
  )
)


cat(
  "\nLandmark 0.50:\n"
)

print(
  head(
    z050
  )
)


cat(
  "\nLandmark 0.75:\n"
)

print(
  head(
    z075
  )
)


cat(
  "\nLandmark 1.00:\n"
)

print(
  head(
    z100
  )
)


# ============================================================
# 16. Check that criteria are different across landmarks
# ============================================================

cat(
  "\n\n====================================\n"
)

cat(
  "Criteria by landmark - replication 1\n"
)

cat(
  "====================================\n"
)

print(
  Results1[[1]]$est2
)


# ============================================================
# 17. Final compact table
# ============================================================

FinalTable <- data.frame(
  
  Landmark = c(
    0.25,
    0.50,
    0.75,
    1.00
  ),
  
  AUC = round(
    MeanCriteria[1:4],
    4
  ),
  
  AUC_SD = round(
    SDCriteria[1:4],
    4
  ),
  
  Brier = round(
    MeanCriteria[5:8],
    4
  ),
  
  Brier_SD = round(
    SDCriteria[5:8],
    4
  )
)

cat(
  "\n\n====================================\n"
)

cat(
  "FINAL TABLE\n"
)

cat(
  "====================================\n"
)

print(
  FinalTable
)



mean(Time1);sd(Time1)


 