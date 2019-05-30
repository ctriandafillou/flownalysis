#' Build a dataframe out of a flowframe
#'
#' This function takes a flowFrame object from the flowCore package as its only argument, and returns a tidy dataframe with the experiment data in columns with names according to the cytometer channels, a relative time column called 'Time', and an absolute time column in standard strptime format called 'abstime'.
#' @param ff The flowFrame
#' @export
#' @return a tidy dataframe of flow cytometry data; rows are observations, columns are channels; the identifier column is by default the name of the fcs file and is called 'exp'


build_ff_df_w_abstime <- function(ff) {
  to.return <- as.data.frame(exprs(ff))
  datetime = strptime(paste(description(ff)$`$DATE`, description(ff)$`$BTIM`), format = "%d-%b-%Y %H:%M:%S")
  to.return$abstime = as.character(datetime)
  return(to.return)
}