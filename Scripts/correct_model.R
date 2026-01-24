library(rxode2)
library(ggplot2)
library(tidyverse)

set.seed(32)
rxSetSeed(32)

mod <- function() {
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
  })

  model({
    ka <- tka * exp(eta.ka)
    cl <- tcl * exp(eta.cl)
    vc <- tvc * exp(eta.vc)

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - cl / vc * central

    conc <- central / vc
    DV <- conc + add.sd
  })
}

sigma <- lotri(add.sd ~ 0.05)

n_sub <- 100


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
et <- dose_et %>%
  et(
    map(
      0:(n_doses - 1),
      ~ .x + c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    ) %>%
      unlist()
  )

sim <- rxSolve(mod, et, sigma = sigma)

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

write_csv(sim_dt, here::here("Data/02_sim_dt.csv"))
