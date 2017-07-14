#' Use kmeans clustering to detect subpopulations in a flow cytometry experiment
#'
#' @param flow.df the name of the flow experiment in the form of a dataframe (assumed output from merge_flowSet)
#' @param cluster.on the column names to cluster on; currently only clustering on subsets of the data is allowed
#' @param tol.u user-defined tolerance for clusters (goes to flowPeaks function); default is 1 (can go from 0 to 1)
#' @param h0.u user-defined h0 parameter; default is 1
#' @param h.u user-defined h parameter; default is 1.5
#' @param qc.exp optional qc experiment to generate a plot for to visually inspect clusters [currently not working]
#' @param qc.ax1 x axis of qc plot (must be a column name in flow.df)
#' @param qc.ax2 y axis of qc plot (must be a column name in flow.df)
#' @export
#' @return the input dataframe flow.df with an additional column called 'cluster' that contains the cluster assignment for each observation

cluster.flow <- function(flow.df, cluster.on, tol.u = 1, h0.u = 1, h.u = 1.5, qc.exp = FALSE, qc.ax1 = NULL, qc.ax2 = NULL){
  df <- select(flow.df, cluster.on)
  fp <- flowPeaks(df, tol = tol.u)
  all.the.clusters <- as.data.frame(fp$peaks.cluster)
  colnames(all.the.clusters) <- "cluster"
  to.return <- bind_cols(flow.df, all.the.clusters) %>%
    mutate(cluster = as.factor(cluster), cluster = plyr::mapvalues(cluster, from = c(1, 2, 3, 4, 5, 6, 7, 8), to = c("A", "B", "C", "D", "E", "F", "G", "H")))
  print(table(all.the.clusters))
  if (qc.exp != FALSE){
    if (qc.exp %in% unique(to.return$exp)){
      to.plot <- filter(to.return, exp == qc.exp)
      p <- ggplot(to.plot, aes(x=qc.ax1, y=qc.ax2, color=cluster)) +
        theme_ct_greypanels() + geom_point(alpha=0.1, size=0.5) + geom_density_2d()
      p
    }
    else {
      print("qc experiment not found")
    }
  }
  return(to.return)
}