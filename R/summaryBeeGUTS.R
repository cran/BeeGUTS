#' Summary of \code{beeSurvFit} objects
#'
#' @description This is the generic \code{summary} S3 method for the \code{beeSurvFit} class.
#' It shows the quantiles of priors and posteriors on parameters.
#'
#' @param object An object of class \code{beeSurvFit}
#' @param ... Additional arguments to be parsed to the generic \code{summary} method (not used)
#'
#' @return A summary of the \code{beeSurvFit} object
#' @export
#'
#' @examples
#' \donttest{
#' data(fitBetacyfluthrin_Chronic)
#' summary(fitBetacyfluthrin_Chronic)
#' }
summary.beeSurvFit <- function(object, ...) {

  cat("Computing summary can take some time. Please be patient...")

  # Prepare prior
  lsData_fit <- object$dataFit
  lsData_fit$nDatasets <- ifelse(is.null(lsData_fit$nDatasets), 1, lsData_fit$nDatasets)

  ## Common parameters
  hb <- 10^qnorm(p = c(0.5, 0.025, 0.975),
                 mean = lsData_fit$hbMean_log10,
                 sd = lsData_fit$hbSD_log10)

  kd <- 10^qnorm(p = c(0.5, 0.025, 0.975),
                 mean = lsData_fit$kdMean_log10,
                 sd = lsData_fit$kdSD_log10)

  ## Model specific parameters
 if (object$modelType == "SD") {

   zw <- 10^qnorm(p = c(0.5, 0.025, 0.975),
                  mean = lsData_fit$zwMean_log10,
                  sd = lsData_fit$zwSD_log10)

   bw <- 10^qnorm(p = c(0.5, 0.025, 0.975),
                  mean = lsData_fit$bwMean_log10,
                  sd = lsData_fit$bwSD_log10)

   outPrior <- data.frame(parameters = c("hb", "kd", "zw", "bw"),
                     median = c(hb[1], kd[1], zw[1], bw[1]),
                     Q2.5 = c(hb[2], kd[2], zw[2], bw[2]),
                     Q97.5 = c(hb[3], kd[3], zw[3], bw[3]))

 } else if (object$modelType == "IT") {

   mw <- 10^qnorm(p = c(0.5, 0.025, 0.975),
                  mean = lsData_fit$mwMean_log10,
                  sd = lsData_fit$mwSD_log10)

   beta <- 10^qunif(p = c(0.5, 0.025, 0.975),
                  min = lsData_fit$betaMin_log10,
                  max = lsData_fit$betaMax_log10)

   outPrior <- data.frame(parameters = c("hb", "kd", "mw", "beta"),
                          median = c(hb[1], kd[1], mw[1], beta[1]),
                          Q2.5 = c(hb[2], kd[2], mw[2], beta[2]),
                          Q97.5 = c(hb[3], kd[3], mw[3], beta[3]))
 }

  # Prepare posteriors

  tmpRes <- rstan::monitor(object$stanFit, print = FALSE)

  ## Common parameters
  hb_med <- c()
  hb_inf95 <- c()
  hb_sup95 <- c()
  for(i in 1:lsData_fit$nDatasets){
    parname <- ifelse(lsData_fit$nDatasets == 1, "hb_log10", paste0("hb_log10[",i,"]"))
    hb_med[i] <- 10^tmpRes[[parname, "50%"]]
    hb_inf95[i] <- 10^tmpRes[[parname, "2.5%"]]
    hb_sup95[i] <- 10^tmpRes[[parname, "97.5%"]]
  }
  kd_med <- 10^tmpRes[["kd_log10", "50%"]]
  kd_inf95 <- 10^tmpRes[["kd_log10", "2.5%"]]
  kd_sup95 <- 10^tmpRes[["kd_log10", "97.5%"]]

  ## Model specific parameters
  if (object$modelType == "SD") {

    zw_med <- 10^tmpRes[["zw_log10", "50%"]]
    zw_inf95 <- 10^tmpRes[["zw_log10", "2.5%"]]
    zw_sup95 <- 10^tmpRes[["zw_log10", "97.5%"]]

    bw_med <- 10^tmpRes[["bw_log10", "50%"]]
    bw_inf95 <- 10^tmpRes[["bw_log10", "2.5%"]]
    bw_sup95 <- 10^tmpRes[["bw_log10", "97.5%"]]

    outPost <- data.frame(parameters = c("kd", "zw", "bw"),
                      median = c(kd_med, zw_med, bw_med),
                      Q2.5 = c(kd_inf95, zw_inf95, bw_inf95),
                      Q97.5 = c(kd_sup95, zw_sup95, bw_sup95))


  } else if (object$modelType == "IT") {

    mw_med <- 10^tmpRes[["mw_log10", "50%"]]
    mw_inf95 <- 10^tmpRes[["mw_log10", "2.5%"]]
    mw_sup95 <- 10^tmpRes[["mw_log10", "97.5%"]]

    beta_med <- 10^tmpRes[["beta_log10", "50%"]]
    beta_inf95 <- 10^tmpRes[["beta_log10", "2.5%"]]
    beta_sup95 <- 10^tmpRes[["beta_log10", "97.5%"]]

    outPost <- data.frame(parameters = c("kd", "mw", "beta"),
                          median = c(kd_med, mw_med, beta_med),
                          Q2.5 = c(kd_inf95, mw_inf95, beta_inf95),
                          Q97.5 = c(kd_sup95, mw_sup95, beta_sup95))

  }

  hbNames <- c()
  for(i in 1:lsData_fit$nDatasets){
    hbNames[i] <- paste0("hb[",i,"]")
  }
  outPost_hb <- data.frame(parameters = hbNames,
                           median = hb_med,
                           Q2.5 = hb_inf95,
                           Q97.5 = hb_sup95)

  # Format and output
  outPrior <- format(data.frame(outPrior), scientific = TRUE, digit = 6)
  outPost <- format(data.frame(outPost), scientific = TRUE, digit = 6)
  outPost_hb <- format(data.frame(outPost_hb), scientific = TRUE, digit = 6)
  maxRhat <- max(rstan::summary(object$stanFit)$summary[,"Rhat"], na.rm= TRUE)
  minBulk_ESS <- min(tmpRes$Bulk_ESS)
  minTail_ESS <- min(tmpRes$Tail_ESS)

  cat("Summary: \n\n")
  cat("Bayesian Inference performed with Stan.\n",
      "Model type:", object$modelType, "\n",
      "Bee species:", object$data$beeSpecies, "\n\n",
      "MCMC sampling setup (select with '$setupMCMC')\n",
      "Iterations:", object$setupMCMC$nIter, "\n",
      "Warmup iterations:", object$setupMCMC$nWarmup, "\n",
      "Thinning interval:", object$setupMCMC$thinInterval, "\n",
      "Number of chains:", object$setupMCMC$nChains)
  cat("\n\nPriors of the parameters (quantiles) (select with '$Qpriors'):\n\n")
  print(outPrior, row.names = FALSE)
  cat("\nPosteriors of the parameters (quantiles) (select with '$Qposteriors'):\n\n")
  print(outPost_hb, row.names = FALSE)
  print(outPost, row.names = FALSE)
  cat("\n\n Maximum Rhat computed (na.rm = TRUE):", maxRhat, "\n",
      "Minimum Bulk_ESS:", minBulk_ESS, "\n",
      "Minimum Tail_ESS:", minTail_ESS, "\n",
      "Bulk_ESS and Tail_ESS are crude measures of effecting sampling size for
      bulk and tail quantities respectively. An ESS > 100 per chain can be
      considered as a good indicator. Rhat is an indicator of chains convergence.
      A Rhat <= 1.05 is a good indicator of convergence. For detail results,
      one can call 'rstan::monitor(YOUR_beeSurvFit_OBJECT$stanFit)")
  cat("\n\n EFSA Criteria (PPC, NRMSE, and SPPE) can be accessed via 'x$EFSA'")

  critEFSA <- criteriaCheck(object)

  invisible(list(
      setupMCMC = object$setupMCMC,
      Qpriors = outPrior,
      Qposteriors_hb = outPost_hb,
      Qposteriors = outPost,
      EFSA = critEFSA))
}



#' Summary of \code{LCx} objects
#'
#' @description This is the generic \code{summary} S3 method for the \code{LCx} class.
#' It shows the median and 95% credible interval of the calculated LCx.
#'
#' @param object An object of class \code{LCx}
#' @param ... Additional arguments to be parsed to the generic \code{summary} method (not used)
#'
#' @return A summary of the \code{LCx} object
#' @export
#'
#' @examples
#' \donttest{
#' data(fitBetacyfluthrin_Chronic)
#' out <- LCx(fitBetacyfluthrin_Chronic)
#' summary(out)
#' }
summary.LCx <- function(object, ...) {
  cat("Summary: \n\n")
  cat("LC",object$X_prop*100, " calculation. \n",
      "Time for which the LCx is calculated:", object$timeLCx, "\n",
      "Bee species:", object$beeSpecies, "\n",
      "Test type:", object$testType, "\n",
      "LCx:", "\n")
  print(object$dfLCx)
}
