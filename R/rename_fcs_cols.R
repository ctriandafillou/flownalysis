#' Rename the columns of a dataframe intelligently
#'
#' @param dataframe the name of the dataframe to be renamed
#' @param print.names whether or not to print the resulting column names; default is F
#' @export
#' @return a dataframe in the same format as merge_flowSet, but with shorter column names

rename_fcs_cols <- function(dataframe, print.names = F){
  default.column.names <- colnames(dataframe)[1:(ncol(dataframe)-2)]
  m <- gregexpr("[[:alpha:]]+[[:digit:]]*.[[:upper:]]{1}", default.column.names)
  colnames(dataframe)[1:(ncol(dataframe)-2)] <- regmatches(default.column.names, m)
  if (print.names == TRUE){
    print(regmatches(default.column.names, m))
  }
  return(dataframe)
}
