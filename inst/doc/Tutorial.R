## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/Tutorial-",
  out.width = "100%"
)

## ----example------------------------------------------------------------------
library(BeeGUTS)
file_location <- system.file("extdata", "betacyfluthrin_chronic_ug.txt", package = "BeeGUTS") # Load the path to one of the example file
lsData <- dataGUTS(file_location = file_location, test_type = 'Chronic_Oral') # Read the example file
plot(lsData) # Plot the data
fit <- fitBeeGUTS(lsData, modelType = "SD", nIter = 2000, nChains = 1) # Fit a SD model. This can take some time...
traceplot(fit) # Produce a diagnostic plot of the fit
plot(fit) # Plot the fit results
summary(fit) # Gives a summary of the results
validation <- validate(fit, lsData) # produce a validation of the fit (here it uses the same dataset as calibration as an example, so not relevant…)
plot(validation) # plot the validation results
dataPredict <- data.frame(time = c(1:5, 1:15), conc = c(rep(5, 5), rep(15, 15)),  replicate = c(rep("rep1", 5), rep("rep3", 15))) # Prepare data for forwards prediction
prediction <- predict(fit, dataPredict) # Perform forwards prediction. At the moment, no concentration recalculation is performed in the forwards prediction. The concentrations are taken as in a chronic test
plot(prediction) # Plot of the prediction results

