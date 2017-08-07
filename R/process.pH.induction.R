#' Processes a flow dataset of induction timeseries
#' 
#' This function takes a dataframe produced by merge_flowSet and produces a processed dataframe with induction data; currently requires pre-classification of populations into a column called 'population'
#' Remaining to do: 1. generate calibration curve and convert to pH within function; 2. find a way to change the induction channel name;
#' @param df The dataframe containing the raw induction data. Identifier column should be called 'exp' and formatted as background data: 'controls_strain_background.fcs'; or induction data: 'timepoint_treatment_fitness_pH.fcs'
#' @param cc whether or not to generate a calibration curve and calculate intracellular pH; default is FALSE : CURRENTLY UNAVAILABLE
#' @param summarize function returns either a median-summarised timeseries (TRUE) or unsummarized data (FALSE; default)
#' @param channel induction channel name; default is "PEDazzle.A" : CURRENTLY UNAVAILABLE
#' @param ref.strain the reference strain; default is "ycgt028"
#' @param ref.background the background of the reference strain; default is "media"
#' @export
#' @return processed dataframe



process.pH.induction <- function(df, cc=FALSE, summarise=FALSE, ref.strain='ycgt028', ref.background='media') {
  bkgs <- filter(df, grepl("controls", exp)) %>%
    separate(exp, c("exp", "strain", "background"), extra = "drop")
  
  ind.bkg <- filter(bkgs, strain == ref.strain & background == ref.background) %>%
    mutate(rel.red = PEDazzle.A / FSC.A) %>%
    summarise(rel.red = median(rel.red, na.rm = T)) %>%
    as.numeric()
  
  ts <- filter(df, !grepl("controls", exp) & !grepl("cc", exp)) %>%
    separate(exp, c("timepoint", "treatment", "shock.pH", "fitness"), extra = "drop", convert = TRUE) %>%
    mutate(rel.red = (PEDazzle.A / FSC.A) / ind.bkg)
  
  if (summarise == TRUE) {
    ts <- ts %>%
      group_by(timepoint, treatment, shock.pH, fitness, population) %>%
      summarise(med.red = median(rel.red, na.rm = T),
                red.high = quantile(rel.red, 0.75, na.rm = T),
                red.low = quantile(rel.red, 0.25, na.rm = T))
    return(ts)
  } else {
    return(ts)
  }
}
