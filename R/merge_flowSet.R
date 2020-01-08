#' Merge a flowSet into a dataframe
#'
#' This function takes a flowSet object created by the read.flowSet function from the flowCore package
#' @param fs The flowset to be merged
#' @param rename.cols whether or not to rename the columns of the dataframe to short channel names; default is TRUE
#' @param method the method for renaming columns; will only be used if rename.cols = TRUE; default is "old"
#' @export
#' @return a tidy dataframe of flow cytometry data; rows are observations, columns are channels; the identifier column is by default the name of the fcs file and is called 'exp'

merge_flowSet <- function(fs, rename.cols = T, method = "old"){
  df <- plyr::rename(plyr::ldply(flowCore::fsApply(fs, function(x) as.data.frame(flowCore::exprs(x)), simplify = FALSE)), c(".id" = "exp")) %>%
    select(-exp, everything())
  if (rename.cols == T) {
    df <- rename_fcs_cols(df, method = method)
    return(df)
  } else {
    return(df)
  }
}
