#' pHluorin Calibration Curve Analysis
#' 
#' This function takes raw data from a calibration curve experiment (yeast expressing pHluorin, flow cytometry data) and does all processing steps to construct a calibration curve that maps fluorescence ratio to pH. Takes *raw* (untransformed) FITC and BV510 channels as inputs!
#' @param subset The dataframe containing the raw calibration curve data. Identifier column should be called 'exp' and formatted as 'cc_pH.fcs'
#' @param background The dataframe containing the raw background data. Identifier column should be called 'exp' and formatted as 'experiment_strain_background' with background being either 'buffer' or 'media'
#' @param id Name/date of experiment
#' @param buffer.values A list of buffer values for that day
#' @param inputs A list of column entries that were used during data collection
#' @param xmin minimum pH predicted (depends on range of buffer values used)
#' @param xmax maximum pH predicted (depends on range of buffer values used)
#' @export
#' @return convert.to.pH function, plot of calibration curve, analysis of fit quality



cc.analysis <- function(subset, background, id, inputs, buffer.values, xmin=4.9, xmax=8.6){
  bkgs <- background %>%
    separate(exp, c("experiment", "strain", "background"), convert=T, extra="drop") %>%
    group_by(strain, background, experiment) %>%
    mutate(rel.red = PEDazzle594.A / FSC.A) %>%
    summarise_all(median)
  
  FITC.bkg <- as.numeric(filter(bkgs, strain == "by4743" & background == "buffer")$FITC.A)
  BV.bkg <- as.numeric(filter(bkgs, strain == "by4743" & background == "buffer")$BV510.A)
  
  FITC.m.bkg <- as.numeric(filter(bkgs, strain == "by4743" & background == "media")$FITC.A)
  BV.m.bkg <- as.numeric(filter(bkgs, strain == "by4743" & background == "media")$BV510.A)

  cc <- subset %>%
    filter(FITC.A > 800 & BV510.A > 800) %>% # High-quality points only
    separate(exp, c("experiment", "pH"), extra="drop") %>%
    mutate(pH = as.numeric(plyr::mapvalues(pH, from=inputs, to=buffer.values)), pH.ratio = (BV510.A - BV.bkg) / (FITC.A - FITC.bkg))
  
  cc.quality <- nrow(cc)/nrow(subset)
  
  cc <- cc %>%
    group_by(pH) %>%
    summarise(med.ratio = median(pH.ratio), high = quantile(pH.ratio, 0.75), low = quantile(pH.ratio, 0.25))
  
  sigmoid <- function(params, x) {
    (params[1] / (1 + exp(-params[2] * (x-params[3])))) + params[4]
  }
  
  x <- cc$pH
  y <- cc$med.ratio
  
  fitmodel <- nls(y~a/(1 + exp(-b * (x-c))) + d, start=list(a=.5, b=2, c=7, d=0.25))
  params = coef(fitmodel)
  
  x2 <- seq(xmin, xmax, 0.01)
  y2 <- sigmoid(params, x2)
  
  q <- ggplot(data=NULL) + geom_point(aes(x=x, y=y)) + geom_line(aes(x=x2, y=y2), size=0.5) 
  q <- q+ geom_errorbar(data=cc, aes(x=pH, ymin=low, ymax=high), size=0.5, width=0.05)
  q <- q + theme_minimal() + labs(
    x = "pH",
    y = "Ratio 405/488",
    title = paste("pHluorin calibration curve", id, sep = " ")) + geom_text(aes(x=5.8, y=0.55, label="ratio==frac(a,1+e^(-b(x-c))) + d"), parse = T)
  
  convert.to.pH <- function(FITC.value=NA, BV.value=NA, ratio.value=NA, p=params, lims=y2, ratio=F){
    # ratio = [a / (1 + exp(-b(pH-c)))] + d
    # pH = [log(a / (ratio - d) - 1) / -b] + c
    if (ratio == T) corrected.ratio = ratio.value
    else corrected.ratio <- (BV.value - BV.m.bkg) / (FITC.value - FITC.m.bkg)
    min.ratio = min(y2)
    max.ratio = max(y2)
    ifelse(corrected.ratio < min.ratio, NA, ifelse(corrected.ratio < max.ratio, (log((p[1]/(corrected.ratio-p[4])) - 1)/-p[2]) + p[3], NA))
  }
  fxn.name <- paste("convert.to.pH", as.character(id), sep=".")
  assign(fxn.name, convert.to.pH, envir = .GlobalEnv)
  
  # Completetion statement:
  cat("Calibration curve with experiment ID \'", id, "\' complete:\n", sep="")
  print(q)
  cat(cc.quality, "% of data was retained\n", sep="")
  cat("To convert with this calibration curve, use the function \'", fxn.name, "\'\n", sep="")
}
