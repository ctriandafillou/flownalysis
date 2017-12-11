#' pHluorin Calibration Curve Analysis
#'
#' This function takes raw data from a calibration curve experiment (yeast expressing pHluorin, flow cytometry data) and does all processing steps to construct a calibration curve that maps fluorescence ratio to pH. Takes *raw* (untransformed) FITC and BV510 channels as inputs!
#' @param subset The dataframe containing the raw calibration curve data. Identifier column should be called 'exp' and formatted as 'cc_pH.fcs'
#' @param background The dataframe containing the raw background data. Identifier column should be called 'exp' and formatted as 'experiment_strain_background' with background being either 'buffer' or 'media'
#' @param id Name/date of experiment, will end up as part of the cc function name so use plaintext characters and no dashes
#' @param buffer.values A list of buffer values for that day, default is FALSE (inferred from filenames)
#' @param inputs A list of column entries that were used during data collection, default is FALSE (inferred from filenames)
#' @param xmin minimum pH predicted (depends on range of buffer values used)
#' @param xmax maximum pH predicted (depends on range of buffer values used)
#' @param FITC.thresh Value that untransformed FITC.A (area) signal must exceed; default is 800
#' @param BV.thresh Value that untransformed BV510.A (area) signal must exceed; default is 800
#' @param start.list Edit the starting fitting parameters; default is list(a=.5, b=2, c=7, d=0.25) (keep this form but change numbers)
#' @param return possible values: "plot" for full analysis, "value" for the value only
#' @export
#' @return convert.to.pH function, plot of calibration curve, analysis of fit quality



pH.pka.analysis <- function(subset, background, id, inputs=FALSE, buffer.values=FALSE, xmin=4.9, xmax=8.6, FITC.thresh = 800, BV.thresh = 800, start.list = list(a=.5, b=2, c=7, d=0.25), return = "plot"){
  bkgs <- background %>%
    separate(exp, c("experiment", "strain", "background"), convert=T, extra="drop") %>%
    group_by(strain, background, experiment) %>%
    summarise_all(median)

  FITC.bkg <- as.numeric(filter(bkgs, strain == "by4743" & background == "buffer")$FITC.A)
  BV.bkg <- as.numeric(filter(bkgs, strain == "by4743" & background == "buffer")$BV510.A)

  FITC.m.bkg <- as.numeric(filter(bkgs, strain == "by4743" & background == "media")$FITC.A)
  BV.m.bkg <- as.numeric(filter(bkgs, strain == "by4743" & background == "media")$BV510.A)

  if (inputs == FALSE){
    cc <- subset %>%
      filter(FITC.A > FITC.thresh & BV510.A > BV.thresh) %>% # High-quality points only
      separate(exp, c("experiment", "pH"), extra="drop") %>%
      mutate(pH = as.numeric(gsub("p", ".", pH)),
             pH.ratio = (BV510.A - BV.bkg) / (FITC.A - FITC.bkg))
  } else{
    cc <- subset %>%
      filter(FITC.A > FITC.thresh & BV510.A > BV.thresh) %>% # High-quality points only
      separate(exp, c("experiment", "pH"), extra="drop") %>%
      mutate(pH = as.numeric(plyr::mapvalues(pH, from=inputs, to=buffer.values)), pH.ratio = (BV510.A - BV.bkg) / (FITC.A - FITC.bkg))
  }

  cc.quality <- nrow(cc)/nrow(subset)*100


  cc <- cc %>%
    group_by(pH) %>%
    summarise(med.ratio = median(pH.ratio), high = quantile(pH.ratio, 0.75), low = quantile(pH.ratio, 0.25))


  sigmoid <- function(params, x) {
    (params[1] / (1 + exp(-params[2] * (x-params[3])))) + params[4]
  }

  x <- cc$pH
  y <- cc$med.ratio

  pka <- data.frame(x,y) %>%
    mutate(log.ratio = -log10((y-max(y))/(min(y)-y))) %>%
    filter(!is.na(log.ratio) & !is.infinite(log.ratio))
  colnames(pka) <- c("pH", "med.ratio", "log.ratio")
  pka.fit <- lm(pka$log.ratio~pka$pH)
  pka.model <- coef(pka.fit)

  if (return == "value") {
    return(-pka.model[[1]]/pka.model[[2]])
  } else {

  plot(pka$pH, pka$log.ratio, xlab="pH", ylab="log((med.ratio - max)/(min-med.ratio))", xlim = c(5.4,8.0), ylim=c(-1.5,0.6))
  abline(pka.model)
  abline(h=0, lty="dotted")
  title(paste("pKa plot", id, sep = " "))

  # Completion statement
  cat(id, "pKa =", -pka.model[[1]]/pka.model[[2]], sep = " ")
  }
}
