#' Rename the columns of a dataframe intelligently
#'
#' @param dataframe the name of the dataframe to be renamed
#' @param print.names whether or not to print the resulting column names; default is F
#' @param method method of column renaming; options are "old" (first iteration, better for non-HTS instrument, use with character swapped values), "new" (second version, better for HTS instrument, use with UNSWAPPED values); default is "new"
#' @export
#' @return a dataframe in the same format as merge_flowSet, but with shorter column names

rename_fcs_cols <- function(dataframe, print.names=F, method="new"){
  default.column.names <- colnames(dataframe)[1:(ncol(dataframe)-2)]
  if (method == "old"){
    m <- gregexpr("[[:alpha:]]+[[:digit:]]*.[[:upper:]]{1}", default.column.names)
    colnames(dataframe)[1:(ncol(dataframe)-2)] <- regmatches(default.column.names, m)
  }
  
  if (method == "new"){
    m <- gregexpr('\\(*[[:alnum:]]{3,5}\\)*\\-*[[:upper:]]{1}', default.column.names) #un-altered, gets everything but dazzle
    colnames(dataframe)[1:(ncol(dataframe)-2)] <- regmatches(default.column.names, m)
    colnames(dataframe) <- plyr::mapvalues(colnames(dataframe), from = c("594)-A", "594)-H", "594)-W"), to = c("PEDazzle.A", "PEDazzle.H", "PEDazzle.W"))
    colnames(dataframe) <- gsub("\\)", "", gsub("-", ".", gsub("\\(", "", colnames(dataframe))))
  }
  
  if (print.names == TRUE){
    print(regmatches(default.column.names, m))
  }
  return(dataframe)
}
