#' Standard plotting theme
#' 
#' Adding this function to a ggplot changes the theme
#' @param base_size text size; default is 12
#' @param point.alpha if T sets the alpha of the legend to be 1; default is F
#' @export


theme_ct_greypanels <- function(base_size=12, point.alpha=F)  
  theme_grey(base_size=base_size) %+replace% 
  theme(#legend.position="none",
    panel.grid=element_blank(),panel.background=element_rect(colour="grey20"),
    axis.ticks=element_line(colour="grey20"), panel.border=element_rect(fill=NA)) 