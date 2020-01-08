#' Read in flow cytometry data generated at the UChicago Flow Cytometry Core
#'
#' Read a folder of FCS files generated on either the LSR Fortessa or the LSR Fortessa HTS at the flow cytometry core at the University of Chicago. Deals with the naming conventions of the individual instruments and produces a tidy dataframe with all specimen represented.
#' @param filename The name of the folder which contains the data. The function will default look in Cat Triandafillou's personal folder, so not useful for other users (sorry) but...
#' @param fullpath Should your data be located somewhere other than the default folder containing data, set this value to 'TRUE' and then put the full filepath instead of the folder name for the filename argument.
#' @param instrument Options are the HTS model ('HTS') or the other model ('nHTS'). 'nHTS' is default, specifying anything else will fail.
#' @export
#' @return tidy dataframe with all data from specified folder; individual FCS files appear in the 'exp' column, and all other data in the file appears in the dataframe.



read.chicago.flowSet <- function(filename, instrument = "nHTS", fullpath = FALSE) {
  # Insert a check if file exists here
  
  # Default behavior is to look in 'flow' folder
  if(fullpath){
    filename = filename
  }
  else {
    filename = paste0("/Users/Triandafillou/Dropbox (Drummond Lab)/cat/data/flow/", filename)
  }
  
  cat(paste0("Looking for files at: \n", filename, "\n"))
  
  if(instrument == "nHTS") {
    return(merge_flowSet(flowCore::read.flowSet(path = filename, alter.names = TRUE)))
  }
  else if (instrument == "HTS") {
    return(merge_flowSet(flowCore::read.flowSet(path = filename), method = "new"))
  }
  else (stop("instrument must be 'nHTS' or 'HTS'"))
}

