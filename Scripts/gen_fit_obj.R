library(rxode2)
library(nlmixr2)
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
      0.5,
      0.01, 0.1
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
    tka <- 3
    tcl <- 13
    tvc <- 130

    # between subject variability
    eta.ka ~ 0.12
    eta.cl + eta.vc ~ c(
      0.5,
      0.01, 0.1
    )
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
    tka <- 3
    tcl <- 13
    tvc <- 130

    # between subject variability
    eta.ka ~ 0.12
    eta.cl + eta.vc ~ c(
      0.5,
      0.01, 0.1
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

  # Nominal sampling times (hours)

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
        0.5 * iiv_scale,
        0.01, 0.1 * iiv_scale
      )
      add.err <- sqrt(0.25 * error_scale)
      prop.err <- sqrt(0.01 * error_scale)
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

get_model <-function(model_string) {

  if(mod_string=="correct"){
      return(mod_correct)
  } else if (mod_string=="struct_miss") {
      return(mod_struct_miss)
  } else if(mod_string=="resid_miss"){
      return(mod_resid_miss)
  } else {
    stop("wrong option for mod_string argument")
  }

}

# Get snakemake variables
error_scale <- as.numeric(snakemake@wildcards[["error_scale"]])
iiv_scale <- as.numeric(snakemake@wildcards[["iiv_scale"]])
frac_dense <- as.numeric(snakemake@wildcards[["frac_dense"]])
mod_string <- snakemake@wildcards[["mod_string"]]
output_path <- snakemake@output[[1]]

# Create simulation data
sim_data <- create_sim_data(
  error_scale,
  iiv_scale,
  frac_dense
)

# Fit to the data a model specified with `mod_string`
fit <- nlmixr(
  get_model(mod_string),
  sim_data,
  est="saem",
  control=saemControl(print=50, nBurn=200, nEm=300),
  table=tableControl(cwres=TRUE),
)

# Save all parameters and fit object in RDS
list(
  mod_name=mod_string,
  error_scale=error_scale,
  iiv_scale=iiv_scale,
  frac_dense=frac_dense,
  fit=fit
) %>% saveRDS(output_path)
