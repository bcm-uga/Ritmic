
#' Pretreatment of the Dcancer matrix
#' 
#' \code{pre_treat_Dcancer} pretreat the Dcancer matrix to remove genes with
#' a negative expression and samples with a proportion of cancer cells under
#' \code{prop_cancer_threshold}.
#'
#' @param Dcancer The matrix of gene expression in each cancer sample.
#' 
#' @param A_cancer A vector with the cancer cell proportion in each sample.
#' @param prop_cancer_threshold The minimum proportion of cancer cells needed to keep the sample.
#' 
#'
#' @return The Dcancer matrix without negative expression and samples with a 
#' cancer proportion under \code{prop_cancer_threshold}.
#' @export
pre_treat_Dcancer = function(Dcancer, A_cancer, prop_cancer_threshold = 0.1){
  colnames(Dcancer) = seq(ncol(Dcancer))
  print("Proportion of negative expression:")
  print(prop.table(table(Dcancer < 0 )))
  print(paste0("Number of samples with less than ",  prop_cancer_threshold, " cancer cells:"))
  print(table(A_cancer<= 0.1))
  Dcancer[Dcancer < 0] = 0
  Dcancer = Dcancer[, A_cancer > 0.1]
  return(Dcancer)
}
