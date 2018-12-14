#' Get background red fluorescence for induction timecourses
#'
#' From a dataframe representing a stress and recovery experiment (timecourse measuring induction of a red fluorophore-tagged protein), generate the background level of fluorescence in untreated cells.
#' @param df Name of the dataframe which contains the induction data. Should be an output of read.chicago.flowSet or similar formatting.
#' @param experiment.name The precise name of the experiment to use to generate the background. Will be in the 'exp' column. Default value is FALSE, where the experiment name will be inferred.
#' @param media Is the background acidic or buffered neutral SC? If so, specify 'sc' or 'buffered.sc'; otherwise returns the median of all rows which have 'backgrounds' and 'ycgt028' in the 'exp' column.
#' @export
#' @return A single numerical value which is the median of the distribution of size-normalized red fluroescence per cell (PEDazzle area / forward scatter area).



generate.red.bkg <- function(df, experiment.name = FALSE, media = "none"){
  # Will always look in 'exp' column by default (can change later but could be messy)
  if (experiment.name == FALSE) {
    if (media == "sc") {
      bkg <- filter(df, grepl("controls[[:graph:]]+ycgt028[[:graph:]]+4p", exp))
    } else if (media == "buffered.sc") {
      bkg <- filter(df, grepl("controls[[:graph:]]+ycgt028[[:graph:]]+7p", exp))
    }
    else {
      bkg <- filter(df, grepl("controls[[:graph:]]+ycgt028[[:graph:]]+", exp))
    }
  } else {
    bkg <- filter(df, grepl(experiment.name, exp))
  }
  
  if ("PEDazzle594.A" %in% colnames(df)) {
    bkg <- mutate(bkg, rel.red = PEDazzle594.A / FSC.A) %>%
      summarise(med.red = median(rel.red, na.rm = T)) %>%
      as.numeric()
  } else if ("PEDazzle.A" %in% colnames(df)) {
    bkg <- mutate(bkg, rel.red = PEDazzle.A / FSC.A) %>%
      summarise(med.red = median(rel.red, na.rm = T)) %>%
      as.numeric()
  } else {stop("Do you have a column for red fluorescence? If so, maybe time to make this function more general.")}
  
  return(bkg)
}

