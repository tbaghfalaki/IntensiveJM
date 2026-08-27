# =============================================================================
# Libraries
# =============================================================================

rm(list = ls())

library(mvtnorm)
library(DPCri)


# =============================================================================
# Simulation parameters
# =============================================================================

NN <- 100

Results1 <- vector("list", NN)


# =============================================================================
# Simulation loop
# =============================================================================

for (kkk in 1:NN) {
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
  # Each marker has:
  #   - random intercept
  #   - random slope
  #
  # Intercept/slope correlation = 0.6 within each marker.
  # Different markers are independent.
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
  #
  # lambda_i(s)
  # ---------------------------------------------------------------------------

  haz_fun_i <- function(s, x1i, x2i, w1i, w2i, ui) {
    # Marker 1
    mu1 <- Beta1[1] +
      Beta1[2] * x1i +
      Beta1[3] * x2i +
      Beta1[4] * s +
      ui[1] +
      ui[2] * s

    # Marker 2
    mu2 <- Beta2[1] +
      Beta2[2] * x1i +
      Beta2[3] * x2i +
      Beta2[4] * s +
      ui[3] +
      ui[4] * s

    # Marker 3
    mu3 <- Beta3[1] +
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


  # ---------------------------------------------------------------------------
  # Survival function
  # ---------------------------------------------------------------------------

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


  # ---------------------------------------------------------------------------
  # Simulate event time
  # ---------------------------------------------------------------------------

  sim_event_time_i <- function(
      x1i,
      x2i,
      w1i,
      w2i,
      ui,
      Tmax = followup) {
    U <- runif(1)

    S_Tmax <- surv_i(
      T = Tmax,
      x1i = x1i,
      x2i = x2i,
      w1i = w1i,
      w2i = w2i,
      ui = ui
    )

    # Event occurs after Tmax
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

  Y1 <- matrix(NA_real_, nsample, length(t))
  Y2 <- matrix(NA_real_, nsample, length(t))
  Y3 <- matrix(NA_real_, nsample, length(t))


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
    td <- rexp(1, rate = 1)

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


  # =============================================================================
  # Long format
  # =============================================================================

  long.list <- vector("list", nsample)

  for (i in 1:nsample) {
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


  # =============================================================================
  # Survival data
  # =============================================================================

  surv.data <- data.frame(
    id       = id,
    x1       = x1,
    x2       = x2,
    w1       = w1,
    w2       = w2,
    survtime = time,
    Event    = death
  )


  # =============================================================================
  # Train / validation split
  # =============================================================================

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


  # =============================================================================
  # Quick checks
  # =============================================================================

  mean(table(long.data$id))
  median(table(long.data$id))
  min(table(long.data$id))
  max(table(long.data$id))

  mean(surv.data$survtime)
  mean(surv.data$Event)


  # =============================================================================
  # TRUE / ORACLE DYNAMIC PREDICTION
  # =============================================================================
  #
  # For landmark s and prediction horizon Delta:
  #
  #   DP_i(s, Delta)
  #
  # = P(s < T_i <= s + Delta | T_i > s, X_i, u_i)
  #
  # = 1 - exp{
  #       - integral_s^(s+Delta) lambda_i(v) dv
  #   }.
  #
  # Here the TRUE parameters and TRUE random effects are used.
  #
  # This is therefore the ORACLE / TRUE dynamic prediction.
  # =============================================================================


  S <- c(0.25, 0.5, 0.75, 1)

  t_pred <- 0.25

  # Make sure every prediction interval is inside follow-up
  stopifnot(
    all(S + t_pred <= followup)
  )


  # Rows = landmark times
  # Columns = criteria returned by DPCri::Criteria()
  estr <- matrix(
    0,
    nrow = length(S),
    ncol = 2
  )

  sdr <- matrix(
    0,
    nrow = length(S),
    ncol = 2
  )

  dp_results <- vector(
    "list",
    length(S)
  )


  for (r in seq_along(S)) {
    s <- S[r]


    # -------------------------------------------------------------------------
    # Subjects who are event-free at landmark s
    # -------------------------------------------------------------------------

    valid_at_risk <- subset(
      surv_valid,
      survtime >= s
    )

    valid_ids_s <- valid_at_risk$id


    # -------------------------------------------------------------------------
    # TRUE conditional survival probability
    # -------------------------------------------------------------------------

    surv_pred <- numeric(
      length(valid_ids_s)
    )


    for (i in seq_along(valid_ids_s)) {
      idv <- valid_ids_s[i]


      # -----------------------------------------------------------------------
      # True cumulative hazard from s to s + t_pred
      # -----------------------------------------------------------------------

      cumhaz <- integrate(
        function(v) {
          haz_fun_i(
            s = v,
            x1i = x1[idv],
            x2i = x2[idv],
            w1i = w1[idv],
            w2i = w2[idv],
            ui = u[idv, ]
          )
        },
        lower = s,
        upper = s + t_pred,
        rel.tol = 1e-6
      )$value


      # -----------------------------------------------------------------------
      # TRUE conditional survival probability
      #
      # P(T > s+t_pred | T > s, X_i, u_i)
      # -----------------------------------------------------------------------

      surv_pred[i] <- exp(
        -cumhaz
      )
    }


    # -------------------------------------------------------------------------
    # TRUE conditional EVENT probability
    #
    # P(s < T <= s+t_pred | T > s, X_i, u_i)
    # -------------------------------------------------------------------------

    event_pred <- 1 - surv_pred


    # -------------------------------------------------------------------------
    # Store TRUE DP
    # -------------------------------------------------------------------------

    dp_results[[r]] <- data.frame(
      id         = valid_ids_s,
      landmark   = s,
      horizon    = s + t_pred,
      surv_pred  = surv_pred,
      event_pred = event_pred
    )


    # -------------------------------------------------------------------------
    # Evaluate dynamic prediction criteria
    # -------------------------------------------------------------------------

    Crit <- Criteria(
      s = s,
      t = t_pred,
      Survt = valid_at_risk$survtime,
      CR = valid_at_risk$Event,
      P = event_pred,
      cause = 1
    )


    # -------------------------------------------------------------------------
    # Store criteria
    # -------------------------------------------------------------------------

    estr[r, ] <- Crit$Cri[, 1]

    sdr[r, ] <- Crit$Cri[, 2]
  }


  # =============================================================================
  # Store results for this simulation
  # =============================================================================

  Results1[[kkk]] <- list(

    # True / oracle dynamic prediction criteria
    TrueDP = estr,

    # Standard deviations / second output returned by Criteria
    TrueDP_SD = sdr,

    # Subject-level true predictions at each landmark
    DP = dp_results
  )


  # Progress
  if (kkk %% 10 == 0) {
    print(kkk)
  }
}


# =============================================================================
# Summarise simulation results
# =============================================================================

# Dimensions:
#
# NN simulations
# x
# length(S) landmarks
# x
# 2 criteria
#
# ------------------------------------------------------------------------------

TrueDP_array <- array(
  NA_real_,
  dim = c(
    NN,
    length(S),
    2
  )
)


for (kkk in 1:NN) {
  TrueDP_array[kkk, , ] <-
    Results1[[kkk]]$TrueDP
}


# =============================================================================
# Mean criteria across simulations
# =============================================================================

TrueDP_mean <- apply(
  TrueDP_array,
  c(2, 3),
  mean,
  na.rm = TRUE
)


TrueDP_sd <- apply(
  TrueDP_array,
  c(2, 3),
  sd,
  na.rm = TRUE
)


# =============================================================================
# Results table
# =============================================================================

TrueDP_results <- data.frame(
  landmark = S,
  criterion1_mean = TrueDP_mean[, 1],
  criterion1_sd = TrueDP_sd[, 1],
  criterion2_mean = TrueDP_mean[, 2],
  criterion2_sd = TrueDP_sd[, 2]
)


print(TrueDP_results)


# =============================================================================
# Example: inspect true DP for the first simulation
# =============================================================================

Results1[[1]]$DP


# =============================================================================
# Example: inspect true DP at landmark s = 0.5
# =============================================================================

Results1[[1]]$DP[[which(S == 0.5)]]



Results1[[kkk]]$TrueDP



Cr0 <- matrix(0, NN, 8)
for (kkk in 1:NN) {
  Cr0[kkk, ] <- as.numeric(Results1[[kkk]]$TrueDP)
}

cbind(apply(Cr0, 2, mean), apply(Cr0, 2, sd))
xtable::xtable(cbind(apply(Cr0, 2, mean), apply(Cr0, 2, sd)), digits = 3)
 