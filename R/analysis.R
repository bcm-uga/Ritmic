#' Compute the down and up regulated genes in the two gene lists D and down 
#'
#' @param D_list 
#' @param U_list 
#'
#' @return
#' @export
generate_data_bypatient = function(D_list, U_list){
  down = colSums(D_list)
  up = colSums(U_list)
  total = down + up
  patient_names = colnames(D_list)
  patients = rep(factor(patient_names , levels = patient_names [order(total)]), 3)
  variable = c(rep("down", length(patient_names)),
               rep("up", length(patient_names)),
               rep("total", length(patient_names)))
  value = c(down, up, total)
  pc = c(down/nrow(D_list)*100, up/nrow(D_list)*100, total/nrow(D_list)*100)
  return(data.frame(patients = patients,
                    variable = variable,
                    value = value,
                    pc = round(pc,2)))
}


#' Title
#'
#' @param data_patients 
#'
#' @import ggplot2
#' @return
#' @export
plot_figure = function(data_patients){
  mytheme = theme(panel.background = element_blank(),
                  panel.grid.major = element_line(colour="black", size = (0.1)),
                  panel.grid.minor = element_blank())
  
  p1 = ggplot(data_patients, aes(x = patients, y = pc)) +
    geom_line(aes(group = variable), colour = "grey80") +
    mytheme +
    ylab("% of gene deregulation") + xlab("patients") +
    geom_point(aes(colour = variable), size = 0.5) +
    ylim(0, 80) +
    scale_x_discrete(breaks = NULL) +
    scale_colour_manual(
      name = "Gene deregulation per patient",
      values = c("blue", "black", "red"),
      labels = c("DOWN", "UP & DOWN", "UP")
    ) #+
  #  theme(legend.position = "none", axis.text.x = element_blank())
  return(p1)
}

#' Plot the heatmap of 
#'
#' @param data 
#'
#' @return
#' @export
plot_heatmap_hclust = function (data) {
  sum(apply(is.na(data), 1, any))
  data = data[!apply(is.na(data), 1, any), ]
  
  # clustering base on correlation for tissues
  tmp_d = data
  tmp_d = t(tmp_d) - apply(tmp_d, 2, mean)
  tmp_d = t(tmp_d)
  tmp_d = cor(tmp_d, method="pe")
  dim(tmp_d)
  hc_col = hclust(dist(1 - tmp_d), method="complete")
  
  Colv = as.dendrogram(hc_col)
  dendrogram="col"      
  
  # clustering base on eucl. dist. for genes
  d = dist(data)
  hc_row = hclust(d, method="complete")
  Rowv = as.dendrogram(hc_row)
  dendrogram="both"      
  
  # col
  colors=c("blue", "gray", "red")
  cols = colorRampPalette(colors)(20)
  
  foo = gplots::heatmap.2(data, Rowv=Rowv, Colv=Colv, dendrogram="col", trace="none", col=cols,
                          labRow = FALSE,labCol = FALSE,
                          main=paste0("Penda (", nrow(data), " genes x ", ncol(data), " samples)"), mar=c(10,5), useRaster=TRUE)
}