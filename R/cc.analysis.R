#' pHluorin Calibration Curve Analysis
#'
#' This function takes raw data from a calibration curve experiment (yeast expressing pHluorin, flow cytometry data) and does all processing steps to construct a calibration curve that maps fluorescence ratio to pH. Takes *raw* (untransformed) FITC and BV510 channels as inputs.
#' 
#' NEW as of 01-2019: now the resulting conversion function can do adaptive background subtraction
#' 
#' WARNING proper background subtraction only achieved when the data are grouped first (see examples)
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
#' @param return.plot Boolean, whether or not to create a plot of the resulting calibration curve. Default is TRUE
#' 
#' @examples
#' # Given a dataframe 'df' that contains multiple samples with different backgrounds
#' # background is specified in the column 'recovery'
#' processed_df <- df %>%
#' group_by(recovery) %>%
#' mutate(pH = convert.to.pH(channel1, channel2, media.type = recovery))
#' 
#' @export
#' @return convert.to.pH function, plot of calibration curve, analysis of fit quality



cc.analysis <- function(subset, background, id, inputs=FALSE, buffer.values=FALSE, xmin=4.9, xmax=8.6, FITC.thresh = 800, BV.thresh = 800, start.list = list(a=.5, b=2, c=7, d=0.25), return.plot = TRUE){
  bkgs <- background %>%
    separate(exp, c("experiment", "strain", "background"), convert=T, extra="drop") %>%
    group_by(strain, background, experiment) %>%
    summarise_all(median)

  # Generate background fluorescence for subtraction
  ## Buffer
  FITC.bkg <- as.numeric(filter(bkgs, strain == "by4743" & grepl("buffer", background))$FITC.A)
  BV.bkg <- as.numeric(filter(bkgs, strain == "by4743" & grepl("buffer", background))$BV510.A)

  ## Media; new as of 2019-01, adapts to different types of media
  
  # Default values for background subtraction
  # Note that if multiple medias are present only the first will be returned, giving spurrious results
  # A better behavior would be to explicitly select for the correct background, but not sure exactly how to do that yet
  FITC.m.bkg <- as.numeric(filter(bkgs, strain == "by4743" & grepl("media", background))$FITC.A)
  BV.m.bkg <- as.numeric(filter(bkgs, strain == "by4743" & grepl("media", background))$BV510.A)
  
  # Values for SC buffered to pH 7.5
  FITC.m.7p5.bkg <- as.numeric(filter(bkgs, strain == "by4743" & grepl("7p5media", background))$FITC.A)
  BV.m.7p5.bkg <- as.numeric(filter(bkgs, strain == "by4743" & grepl("7p5media", background))$BV510.A)
  
  # Values for regular SC when explicitly called out
  FITC.m.4p0.bkg <- as.numeric(filter(bkgs, strain == "by4743" & grepl("4p0media", background))$FITC.A)
  BV.m.4p0.bkg <- as.numeric(filter(bkgs, strain == "by4743" & grepl("4p0media", background))$BV510.A)
  
  
  ## Fit calibration curve
  
  if (inputs == FALSE){
    cc <- subset %>%
      filter(FITC.A > FITC.thresh & BV510.A > BV.thresh) %>% # High-quality points only
      separate(exp, c("experiment", "pH"), extra="drop") %>%
      mutate(pH = as.numeric(gsub("p", ".", pH)),
             pH.ratio = (BV510.A - BV.bkg) / (FITC.A - FITC.bkg))
  } else {
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

  fitmodel <- nls(y~a/(1 + exp(-b * (x-c))) + d, start=start.list)
  params = coef(fitmodel)

  x2 <- seq(xmin, xmax, 0.01)
  y2 <- sigmoid(params, x2)

  
  
  
  
  # Add to function that will be exported to the main environment
  ## Still blows my mind that this is allowed
  ### Thanks R
  convert.to.pH <- function(FITC.value=NA, BV.value=NA, ratio.value=NA, p=params, lims=y2, ratio=FALSE,
                            media.type="unspecified", verbose=FALSE){
    # ratio = [a / (1 + exp(-b(pH-c)))] + d
    # pH = [log(a / (ratio - d) - 1) / -b] + c
    
    if (ratio == T) corrected.ratio = ratio.value
    else { # only calculate the corrected ratio if it's going to be used
      
      # Consider changing these to logical greps rather than exact values to make things more robust?
      if (media.type == "7p5") { # Condition where recovery was in SC buffered to pH 7.5
        corrected.ratio <- (BV.value - BV.m.7p5.bkg) / (FITC.value - FITC.m.7p5.bkg)
      } else if (media.type == "4p0") { # Recovery in 4.0 SC, explicitly stated
        corrected.ratio <- (BV.value - BV.m.4p0.bkg) / (FITC.value - FITC.m.4p0.bkg)
      } else if (media.type == "unspecified") { # Default condition; recovery in SC not explicitly stated
        if (length(FITC.m.bkg != 1)) {
          cat("WARNING: Multiple backgrounds in input backgrounds dataframe detected. Defaulting to first entry.\nIf this behavior is not expected, check number of background readings submitted to 'cc.analysis'.")
        }
        corrected.ratio <- (BV.value - BV.m.bkg) / (FITC.value - FITC.m.bkg)
      } else {
        # This won't print for every row if applied by mutate, which is nice
        cat('Unrecognized media type entered, not converting')
        return(NA)
      }
      
    }
    
    # Set min/max
    min.ratio = min(y2)
    max.ratio = max(y2)
    
    # Extra info printed if verbose is TRUE:
    if (verbose) {
      #cat("The default background value for the FITC channel is", FITC.m.bkg, "\n", sep = " ")
      #cat("The default background value for the BV510 channel is ", BV.m.bkg, "\n", sep = " ")
      cat("The background value for the FITC channel for 7.5 media is ", FITC.m.7p5.bkg, "\n", sep = "")
      cat("The background value for the BV510 channel for 7.5 media is ", BV.m.7p5.bkg, "\n", sep = "")
      cat("The background value for the FITC channel for 4.0 media is ", FITC.m.4p0.bkg, "\n", sep = "")
      cat("The background value for the BV510 channel for 4.0 media is ", BV.m.4p0.bkg, "\n", sep = "")
    }
    
    # Return statement
    return(ifelse(corrected.ratio < min.ratio, NA, ifelse(corrected.ratio < max.ratio, (log((p[1]/(corrected.ratio-p[4])) - 1)/-p[2]) + p[3], NA)))

  }
  
  fxn.name <- paste("convert.to.pH", as.character(id), sep=".")
  assign(fxn.name, convert.to.pH, envir = .GlobalEnv)

  if(return.plot) {
    q <- ggplot(data=NULL) + geom_point(aes(x=x, y=y)) + geom_line(aes(x=x2, y=y2), size=0.5)
    q <- q + geom_errorbar(data=cc, aes(x=pH, ymin=low, ymax=high), size=0.5, width=0.05)
    q <- q + theme_minimal() + labs(
      x = "pH",
      y = "Ratio 405/488",
      title = paste("pHluorin calibration curve", id, sep = " "))
    # Uncomment to add equation to plot
    #+ geom_text(aes(x=5.8, y=0.55, label="ratio==frac(a,1+e^(-b(x-c))) + d"), parse = T)

    print(q)
  }
  
  # Completion statement:
  cat("Calibration curve with experiment ID \'", id, "\' complete:\n", sep="")
  cat(cc.quality, "% of data was retained\n", sep="")
  cat("To convert with this calibration curve, use the function \'", fxn.name, "\'\n", sep="")
  
}
