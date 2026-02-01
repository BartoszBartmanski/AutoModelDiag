library(rxode2)
library(nlmixr2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(purrr)
library(tibble)

set.seed(32)
rxSetSeed(32)


ini_vals=c(ka=2,cl=3)

mod_correct <- function() {

  ini({
    tka <- 3
    tcl <- 13
    tvc <- 130

    eta.ka ~ 0.12
    eta.cl + eta.vc ~ c(
      0.05,
      0.01, 0.01
    )

    add.err  <- sqrt(0.25)
    prop.err <-  sqrt(0.01)
  })

  model({
  ka = tka * exp(eta.ka)
  cl = (1 - exp(-0.0199 * time)) * tcl * exp(eta.cl)
  vc = tvc * exp(eta.vc)
  d/dt(depot) = -ka * depot
  d/dt(central) = ka * depot - cl / vc * central
  conc = central / vc
  DV = conc
  DV ~ add(add.err) + prop(prop.err) })
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


create_event_table_dataframe <- function(
    frac_dense=0.8,
    jitter_min = 10,
    seed = 355232
) {

  set.seed(seed)

  jitter_hr <- jitter_min / 60

  ## ---- Nominal sampling times (hours) ----

  #Phase 1 subjects

  num_subj_phase1=70
  #Single dose 100 mg
  ph1_pk_sampling <- c(0.5, 1, 1.5, 2, 4, 6, 12, 24, 36, 48, 60, 72, 96)

  #Phase two subjects
  # 5 cycles, each 21 days, with 100 mg q.d.
 dense_times <-  c(
    7*24 + c(0.5, 1, 2, 3, 4, 6, 8, 9, 24, 48, 72, 96,120),
    1*24 + c(-0.5, 1, 2, 4),
    8*24 + c(-0.5, 1, 2, 4),
    15*24 + c(-0.5, 0.5, 1, 2, 3, 4, 6, 8, 9, 24)
    )  %>% sort() %>% unique()

  # predose on day 1 of cycles 1–5 and predose on cycle 7 day 1, cycle 8 day 1, and cycle 10 day 1
   sparse_times <- map(1+seq(0,4)*24, ~.x -0.5)%>% unlist()  %>% sort() %>% unique()

  # Phase 1 eventable
  pk_sampling1=tibble(id = seq(num_subj_phase1)) %>%
    mutate(
      time = map(
        id,
        ~ ph1_pk_sampling +
          runif(length(ph1_pk_sampling), -jitter_hr, jitter_hr)
      )
    )%>%
    unnest(time) %>%
    mutate(time=round(time,3)) %>%
    arrange(id, time)

  # ## Phase 2 eventable
  num_subj_phase2=276

   pk_sampling_ph2 = tibble(id=seq((num_subj_phase1+1),num_subj_phase1+num_subj_phase2),
                 pk_sampling_type = sample(c("dense", "sparse"),
                                     num_subj_phase2,
                                     replace = TRUE,
                                     prob = c(frac_dense, 1-frac_dense)
                                     )) %>%
    mutate(
      time = map(
        pk_sampling_type,
        ~ {
          nominal <- if (.x == "dense") dense_times else sparse_times
          nominal + runif(length(nominal), -jitter_hr, jitter_hr)
        }
      )
    ) %>%
    unnest(time) %>%
    mutate(time = round(time, 3)) %>%
    arrange(id, time) %>%
    select(-pk_sampling_type)



df_ev=NULL

i=1

 for (i in seq_len(num_subj_phase1)){

   # dose records
   ev  = et(amountUnits = "mg",timeUnits   = "hours")%>%
         et(amt = 100,nbr.doses = 1,start.time=0 )

  # adding the pk_sampling
  ev=ev %>% et(pk_sampling1%>% filter(id==i)%>% pull(time))

  # adding it to the event df
  df_ev=rbind(df_ev,tibble(ID=i,ev$get.EventTable()))
}

i=num_subj_phase1+1
for (i in seq(num_subj_phase1+1,num_subj_phase1+num_subj_phase2)){

  # dose records
  ev   = et(amountUnits = "mg",timeUnits   = "hours")%>%
  et(dose = 100, nbr.doses =15, start.time=0, dosing.interval = 24)
  ev_exp=etExpand(ev)

  # adding the pk_sampling
  ev2=ev_exp %>% et(pk_sampling_ph2 %>% filter(id==i)%>% pull(time))

  # adding it to the event df
  df_ev=rbind(df_ev,tibble(ID=i,ev2$get.EventTable())%>% select(-ii))
}

  return(df_ev)
}


create_sim_data <- function(error_scale, iiv_scale, frac_dense = 0.5) {
  mod <- function() {
    ini({
      tka <- 3.13
      tcl <- 14.472
      tvc <- 120
      # between subject variability
      eta.ka ~ 0.12 * iiv_scale
      eta.cl + eta.vc ~ c(
        0.05 * iiv_scale,
        0.01 * iiv_scale, 0.01 * iiv_scale
      )
      add.err <- sqrt(0.25) * error_scale
      prop.err <- sqrt(0.01) * error_scale
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


  df_ev <- create_event_table_dataframe(frac_dense=frac_dense)

  sim <- rxSolve(mod, df_ev)

  demo <- sim %>%
    distinct(id) %>%
    rename(ID = id)

  obs_dt <- sim %>%
    as_tibble() %>%
    select(id, time, DV) %>%
    rename(ID = id, TIME = time) %>%
    mutate(AMT = 0, CMT = 2, EVID = 0)

  dose_dt <- df_ev %>%
    as_tibble() %>%
    filter(evid==1)%>%
    select(-c( "evid")) %>%
    rename(AMT = amt, TIME = time) %>%
    mutate(EVID = 101, CMT = 1, DV = 0)

  sim_dt <- rbind(obs_dt, dose_dt) %>%
    arrange(ID, TIME, -EVID) %>%
    left_join(demo, by = "ID") %>%
    mutate(TSLD = TIME %% 24)

  return(sim_dt)
}

extract_fit_features <- function(fit) {
  res <- data.frame(time=fit$TIME, cwres=fit$CWRES, dv=fit$DV, ipred=fit$IPRED)
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
  ofv <- fit$objDf$OBJF
  aic <- fit$objDf$AIC
  bic <- fit$objDf$BIC

  # ----------------------------
  # 2. CWRES summary
  # ----------------------------
  # TODO: check this formula
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
  # TODO: why the second element of the auto-correlations array
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

  delta_bias_early_late <- (
    mean(pred_err[late], na.rm = TRUE) -
    mean(pred_err[early], na.rm = TRUE)
  )

  # ----------------------------
  # 6. ETA features
  # ----------------------------
  eta_cl_name <- grep("cl", names(eta), value = TRUE)[1]

  eta_cl_sd   <- if (!is.na(eta_cl_name)) sd(eta[[eta_cl_name]], na.rm = TRUE) else NA
  eta_cl_mean <- if (!is.na(eta_cl_name)) mean(eta[[eta_cl_name]], na.rm = TRUE) else NA

  # Correlation between first two ETAs (if available)
  eta_cor <- if (ncol(eta) >= 2) {
    cor(eta[, 2], eta[, 3], use = "complete.obs")
  } else {
    NA
  }

  # ----------------------------
  # 7. Shrinkage
  # ----------------------------
  shrink_cl <- tryCatch(
    fit$shrink[grep("cl", names(fit$shrink))[1]],
    error = function(e) NA
  )

  # ----------------------------
  # 8. Convergence / stability
  # ----------------------------
  # TODO: maybe try to get the following two parameters from fit object
  # it is possible to set maxIter -> not get it from fit obj
  # converged <- as.numeric(isTRUE(fit$convergence))
  # n_iter <- tryCatch(fit$iterations, error = function(e) NA)


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
    shrink_cl               = shrink_cl
  )

  return(features)
}

#'
#'
#' @param mod_name_string String with the name of a rxode 2 model. Currently accepts: "correct", "struct_miss", "resid_miss".
#' @param error_scale Double that will indicate the scale in which the error/additive will be increased.
#' @param iiv_scale Double that will indicate the scale in which the variance of the IIV of all variables will be increased.
#' @param frac_dense String that defines the PK sampling: "dense", "sparse"
#'
#' @return A list with the mod name string, fit features, error_scale, iiv_scale and pk_sampling
#'
create_ML_samples <- function(
    mod_string, error_scale, iiv_scale, frac_dense
) {
  sim_data <- create_sim_data(
    error_scale,
    iiv_scale,
    frac_dense
  )

  if(mod_string=="correct"){

    fit <- nlmixr(
      mod_correct,
      sim_data,
      "saem",
      # control=list(print=0),
      table=list(cwres=TRUE, npde=TRUE)
    )

  } else if (mod_string=="struct_miss") {

    fit <- nlmixr(
      mod_struct_miss,
      sim_data,
      "saem",
      # control=list(print=0),
      table=list(cwres=TRUE, npde=TRUE)
    )

  }else if(mod_string=="resid_miss"){

    fit <- nlmixr(
      mod_resid_miss,
      sim_data,
      "saem",
      # control=list(print=0),
      table=list(cwres=TRUE, npde=TRUE)
    )

  } else {
    stop("Wrong option for mod_string argument")
  }

  fit_features <- extract_fit_features(fit)

  obj=list(
    mod_name=mod_string,
    fit_features=fit_features,
    error_scale=error_scale,
    iiv_scale=iiv_scale,
    pk_sampling=pk_sampling
  )

  return(obj)
}

error_scale=1
iiv_scale=1
frac_dense=1

create_ML_samples("correct", error_scale, iiv_scale, frac_dense)
