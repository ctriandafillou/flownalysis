#' Processes a flow dataset of induction timeseries
#' 
#' This function takes a dataframe produced by merge_flowSet and produces a processed dataframe with induction data
#' @param df The dataframe containing the raw induction data. Identifier column should be called 'exp' and formatted as background data: 'controls_strain_background.fcs'; or induction data: 'timepoint_treatment_fitness_pH.fcs'
#' @param cc whether or not to generate a calibration curve and calculate intracellular pH; default is FALSE
#' @param summarize function returns either a median-summarised timeseries (TRUE) or distributions (FALSE; default)
#' @param channel induction channel name; default is "PEDazzle.A"
#' @export
#' @return processed dataframe



process.pH.induction <- function(df, cc=FALSE, summarise=FALSE, channel='PEDazzle.A') {
  bkgs <- filter(df, grepl("controls", exp)) %>%
    separate(exp, c("exp", "strain", "background"), extra = "drop") %>%
    group_by(strain, background) %>%
    summarise(med.)
  
}
