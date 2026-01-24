library(rxode2)
library(ggplot2)
library(tidyverse)


#TODO set a variable with common initial estimates for all the models tested
set.seed(32)
rxSetSeed(32)


ini_vals=c(ka=2,cl=3)

mod_correct <- function() {
  ini({
    tka <- ini_vals$ka
    tcl <- 50
    tvc <- 7
    # between subject variability
    eta.ka ~ 0.12
    eta.cl + eta.vc ~ c(
      0.05,
      0.01, 0.01
    )
    add.err <- sqrt(0.25)
    prop.err <- sqrt(0.01)
  })

  model({
    ka <- tka * exp(eta.ka)
    cl <- (1 - exp(-0.0199 * time)) * tcl * exp(eta.cl)
    vc <- tvc * exp(eta.vc)

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - cl / vc * central

    conc <- central / vc
    DV <- conc
    DV ~ add(add.err) + prop(prop.err)
  })
}

mod_resid_miss <- function() {
  ini({
    tka <- 2.5
    tcl <- 50
    tvc <- 7
    # between subject variability
    eta.ka ~ 0.12
    eta.cl + eta.vc ~ c(
      0.05,
      0.01, 0.01
    )
    add.err <- sqrt(0.25)
    prop.err <- sqrt(0.01)
  })

  model({
    ka <- tka * exp(eta.ka)
    cl <- (1 - exp(-0.0199 * time)) * tcl * exp(eta.cl)
    vc <- tvc * exp(eta.vc)

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - cl / vc * central

    conc <- central / vc
    DV <- conc
    DV ~ prop(prop.err)
  })
}

mod_struct_miss<- function() {
  ini({
    tka <- 2.5
    tcl <- 50
    tvc <- 7
    # between subject variability
    eta.ka ~ 0.12
    eta.cl + eta.vc ~ c(
      0.05,
      0.01, 0.01
    )
    add.err <- sqrt(0.25)
    prop.err <- sqrt(0.01)
  })

  model({
    ka <- tka * exp(eta.ka)
    cl <- tcl * exp(eta.cl)
    vc <- tvc * exp(eta.vc)

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - cl / vc * central

    conc <- central / vc
    DV <- conc
    DV ~ add(add.err) + prop(prop.err)
  })
}

create_sim_data <- function(error_size, iiv_size, sampling) {
  mod <- function() {
    ini({
      tka <- 2.5
      tcl <- 50
      tvc <- 7
      # between subject variability
      eta.ka ~ 0.12
      eta.cl + eta.vc ~ c(
        0.05 * iiv_size,
        0.01, 0.01 * iiv_size
      )
      add.err <- sqrt(0.25) * error_size
      prop.err <- sqrt(0.01) * error_size
    })

    model({
      ka <- tka * exp(eta.ka)
      cl <- (1 - exp(-0.0199 * time)) * tcl * exp(eta.cl)
      vc <- tvc * exp(eta.vc)

      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl / vc * central

      conc <- central / vc
      DV <- conc
      DV ~ add(add.err) + prop(prop.err)
    })
  }

  n_sub <- 100
  n_doses <- 3


  # dose records
  dose_et <- et(
    amountUnits = "mg",
    timeUnits = "days"
  ) %>%
    et(
      id = 1:n_sub,
      dose = 400,
      nbr.doses = n_doses,
      dosing.interval = 1
    )

  # eventable
  et <- dose_et %>% sampling()

  sim <- rxSolve(mod, et)

  demo <- sim %>%
    distinct(id) %>%
    rename(ID = id)

  obs_dt <- sim %>%
    as_tibble() %>%
    select(id, time, DV) %>%
    rename(ID = id, TIME = time) %>%
    mutate(AMT = 0, CMT = 2, EVID = 0)

  dose_dt <- dose_et$expand() %>%
    as_tibble() %>%
    rename(ID = id, AMT = amt, TIME = time) %>%
    select(-c("ii", "evid")) %>%
    mutate(EVID = 101, CMT = 1, DV = 0)

  sim_dt <- rbind(obs_dt, dose_dt) %>%
    arrange(ID, TIME, -EVID) %>%
    left_join(demo, by = "ID") %>%
    mutate(TSLD = TIME %% 1)

  return(sim_dt)
}

extract_fit_features <- function(fit) {
  # Safety checks
  if (is.null(fit$residuals)) stop("No residuals found in fit object")
  if (is.null(fit$eta)) stop("No ETA estimates found in fit object")

  res <- fit$residuals
  eta <- fit$eta

  # Ensure required columns exist
  required_cols <- c("time", "cwres", "dv", "ipred")
  missing_cols <- setdiff(required_cols, names(res))
  if (length(missing_cols) > 0) {
    stop(paste("Missing residual columns:", paste(missing_cols, collapse = ", ")))
  }

  # ----------------------------
  # 1. OFV / information
  # ----------------------------
  ofv  <- fit$objf
  aic  <- tryCatch(AIC(fit), error = function(e) NA)
  bic  <- tryCatch(BIC(fit), error = function(e) NA)

  # ----------------------------
  # 2. CWRES summary
  # ----------------------------
  cwres_mean <- mean(res$cwres, na.rm = TRUE)
  cwres_sd   <- sd(res$cwres, na.rm = TRUE)
  cwres_skew <- mean((res$cwres - cwres_mean)^3, na.rm = TRUE) / cwres_sd^3

  # ----------------------------
  # 3. CWRES vs time
  # ----------------------------
  lm_cwres_time <- lm(cwres ~ time, data = res)
  cwres_slope   <- coef(lm_cwres_time)[["time"]]
  cwres_r2      <- summary(lm_cwres_time)$r.squared

  # ----------------------------
  # 4. CWRES autocorrelation
  # ----------------------------
  cwres_acf1 <- tryCatch(
    acf(res$cwres, plot = FALSE)$acf[2],
    error = function(e) NA
  )

  # ----------------------------
  # 5. Prediction error features
  # ----------------------------
  pred_err <- res$dv - res$ipred
  pred_rmse <- sqrt(mean(pred_err^2, na.rm = TRUE))
  pred_bias <- mean(pred_err, na.rm = TRUE)

  # Early vs late split
  t_med <- median(res$time, na.rm = TRUE)
  early <- res$time <= t_med
  late  <- res$time >  t_med

  delta_bias_early_late <- mean(pred_err[late], na.rm = TRUE) -
    mean(pred_err[early], na.rm = TRUE)

  # ----------------------------
  # 6. ETA features
  # ----------------------------
  eta_cl_name <- grep("CL", names(eta), value = TRUE)[1]

  eta_cl_sd   <- if (!is.na(eta_cl_name)) sd(eta[[eta_cl_name]], na.rm = TRUE) else NA
  eta_cl_mean <- if (!is.na(eta_cl_name)) mean(eta[[eta_cl_name]], na.rm = TRUE) else NA

  # Correlation between first two ETAs (if available)
  eta_cor <- if (ncol(eta) >= 2) {
    cor(eta[[1]], eta[[2]], use = "complete.obs")
  } else {
    NA
  }

  # ----------------------------
  # 7. Shrinkage
  # ----------------------------
  shrink_cl <- tryCatch(
    fit$shrink[grep("CL", names(fit$shrink))[1]],
    error = function(e) NA
  )

  # ----------------------------
  # 8. Convergence / stability
  # ----------------------------
  converged <- as.numeric(isTRUE(fit$convergence))
  n_iter    <- tryCatch(fit$iterations, error = function(e) NA)

  # ----------------------------
  # Output feature vector
  # ----------------------------
  features <- c(
    ofv                     = ofv,
    aic                     = aic,
    bic                     = bic,
    cwres_mean              = cwres_mean,
    cwres_sd                = cwres_sd,
    cwres_skew              = cwres_skew,
    cwres_slope_time        = cwres_slope,
    cwres_r2_time           = cwres_r2,
    cwres_acf1              = cwres_acf1,
    pred_rmse               = pred_rmse,
    pred_bias               = pred_bias,
    delta_pred_bias_late    = delta_bias_early_late,
    eta_cl_mean             = eta_cl_mean,
    eta_cl_sd               = eta_cl_sd,
    eta_cor_1_2             = eta_cor,
    shrink_cl               = shrink_cl,
    converged               = converged,
    n_iterations            = n_iter
  )

  return(features)
}



create_ML_samples<-function(mod_string,
                            error_size,
                            iiv_size,
                            sampling){

  sim_data=create_sim_data( error_size,
                   iiv_size,
                   sampling)


  if(mod_string=="correct"){
    fit <- nlmixr(mod_correct, sim_data, est="focei")
  }else if(mod_string=="struct_miss"){
    fit <- nlmixr(mod_struct_miss, sim_data, est="focei")
  }else if(mod_string=="resid_miss"){
    fit <- nlmixr(mod_resid_miss, sim_data, est="focei")
  }else {
    stop()
    }

  features=extract_fit_features(fit)

  obj=list(mod_name=mod_string,
           features=features,
           error_size=error_size,
           iiv_size=iiv_size,
           sampling=sampling)
return(obj)
}