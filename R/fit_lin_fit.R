#' Processes a flow dataset of induction timeseries
#' 
#' Fits a linear function to some data and returns a dataframe containing the slope and intercept
#' @param dat The dataframe containing columns corresponding to recovery times ("timepoint") and log ratios of cell types ("growth.ratio")
#' @export
#' @return dataframe with slope and intercept



fit_lin_fit <- function(dat) { 
  the_fit <- lm(growth.ratio ~ timepoint, dat)
  setNames(data.frame(t(coef(the_fit))), c("intercept", "slope"))
}