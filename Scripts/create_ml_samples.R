
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
      table=list(cwres=true, npde=true)
    )

  } else if (mod_string=="struct_miss") {

    fit <- nlmixr(
      mod_struct_miss,
      sim_data,
      "saem",
      # control=list(print=0),
      table=list(cwres=true, npde=true)
    )

  } else if(mod_string=="resid_miss"){

    fit <- nlmixr(
      mod_resid_miss,
      sim_data,
      "saem",
      # control=list(print=0),
      table=list(cwres=true, npde=true)
    )

  } else {
    stop("wrong option for mod_string argument")
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