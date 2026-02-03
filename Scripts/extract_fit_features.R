library(tidyverse)
library(nlmixr2)
library(jsonlite)


extract_fit_features <- function(fit) {
  res <- data.frame(time=fit$TIME, cwres=fit$CWRES, dv=fit$DV, ipred=fit$IPRED)
  eta <- fit$eta

  # Ensure required columns exist
  required_cols <- c("time", "cwres", "dv", "ipred")
  missing_cols <- setdiff(required_cols, names(res))
  if (length(missing_cols) > 0) {
    stop(paste("Missing residual columns:", paste(missing_cols, collapse = ", ")))
  }

  # 1. OFV / information
  ofv <- fit$objDf$OBJF
  aic <- fit$objDf$AIC
  bic <- fit$objDf$BIC

  # 2. CWRES summary
  # TODO: check this formula
  cwres_mean <- mean(res$cwres, na.rm = TRUE)
  cwres_sd   <- sd(res$cwres, na.rm = TRUE)
  cwres_skew <- mean((res$cwres - cwres_mean)^3, na.rm = TRUE) / cwres_sd^3

  # 3. CWRES vs time
  lm_cwres_time <- lm(cwres ~ time, data = res)
  cwres_slope   <- coef(lm_cwres_time)[["time"]]
  cwres_r2      <- summary(lm_cwres_time)$r.squared

  # 4. CWRES autocorrelation
  # TODO: why the second element of the auto-correlations array
  cwres_acf1 <- tryCatch(
    acf(res$cwres, plot = FALSE)$acf[2],
    error = function(e) NA
  )

  # 5. Prediction error features
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

  # 6. ETA features
  eta_cl_name <- grep("cl", names(eta), value = TRUE)[1]

  eta_cl_sd   <- if (!is.na(eta_cl_name)) sd(eta[[eta_cl_name]], na.rm = TRUE) else NA
  eta_cl_mean <- if (!is.na(eta_cl_name)) mean(eta[[eta_cl_name]], na.rm = TRUE) else NA

  # Correlation between first two ETAs (if available)
  eta_cor <- if (ncol(eta) >= 2) {
    cor(eta[, 2], eta[, 3], use = "complete.obs")
  } else {
    NA
  }

  # 7. Shrinkage
  shrink_cl <- tryCatch(
    fit$shrink[grep("cl", names(fit$shrink))[1]],
    error = function(e) NA
  )

  # Output feature vector
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

# Get snakemake variables
error_scale <- as.numeric(snakemake@wildcards[["error_scale"]])
iiv_scale <- as.numeric(snakemake@wildcards[["iiv_scale"]])
frac_dense <- as.numeric(snakemake@wildcards[["frac_dense"]])
mod_string <- snakemake@wildcards[["mod_string"]]
input_file <- snakemake@input[[1]]
output_path <- snakemake@output[[1]]

params_and_fit <- readRDS(input_file)

# extract features
features <- extract_fit_features(params_and_fit$fit)

# Append fit parameters
features <- c(
  features,
  error_scale=error_scale,
  iiv_scale=iiv_scale,
  frac_dense=frac_dense,
  mod_string=mod_string
)

# Save in json
jsonlite::write_json(features, output_path)
