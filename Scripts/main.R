library(rxode2)
library(ggplot2)
library(tidyverse)

set.seed(32)
rxSetSeed(32)

mod_correct <- function() {
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
    DV ~ add(add.err) + prop(prop.err)
  })
}

mod_just_prop_err <- function() {
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

mod_no_time_dep <- function() {
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

creat_sim_data <- function(error_size, iiv_size, sampling) {
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

  sigma <- lotri(add.sd ~ 0.05)

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

  sim <- rxSolve(mod, et) # , sigma = sigma)

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
