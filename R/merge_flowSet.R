#' Merge a flowSet into a dataframe
#'
#' This function takes a flowSet object created by the read.flowSet function from the flowCore package
#' @param fs The flowset to be merged
#' @export
#' @return a tidy dataframe of flow cytometry data; rows are observations, columns are channels; the identifier column is by default the name of the fcs file and is called 'exp'

merge_flowSet <- function(fs, rename.cols = T){
  plyr::rename(plyr::ldply(fsApply(fs, function(x) as.data.frame(exprs(x)), simplify = FALSE)), c(".id" = "exp")) %>%
    select(-exp, everything())
  #if (rename.cols == T) {
  #  rename_fcs_cols(fs)
  #}
}
