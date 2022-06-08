#' Pre-treatment of the Penda results
#' 
#' \code{pre_treat_Penda} takes as input the Penda result lists and 
#' combines them to make only one matrix of size genes*samples, with 1 for 
#' deregulated genes and 0 for conserved genes. Then, the matrix is sorted to 
#' conserved only the genes with a different deregulation status in more than 
#' thres_p samples. 
#' 
#' @param penda_res List of two matrices, $up_genes and $down_genes with 
#'   1 or 0 for each gene and each sample obtained by the Penda method. 
#'   See: \code{\link[penda]{penda_test}}
#' @param thres_p Minimal number of samples in each deregulation group to 
#'  conserved the gene. 
#'
#' @return A binary matrix of deregulation 0/1 with the genes that can be used 
#'  for the statistical analysis and their status for each sample. 
#' 
#' @export
pre_treat_Penda = function(penda_res, thres_p){
  if(thres_p >= ncol(penda_res$up_genes)){
    stop("You can't reach a minimal number of samples higher than the total number of samples")
  } else {
    penda_res <- abs(penda_res$up_genes - penda_res$down_genes)
    
    #Genes with the same deregulation status in all the samples are removed
    g_same = which(apply(penda_res, 1, function(g) {
      length(table(g)) == 1
    }))
    penda_res = penda_res[-g_same, ]
    
    #Genes with less than thres_p samples in each deregulation group are removed
    g_sup = which(apply(penda_res, 1, function(g) {
      table(g)[1] > thres_p & table(g)[2] > thres_p
    }))
    penda_res = penda_res[g_sup, ]
    
    return(penda_res)
  }
}

#' Comparison between gene deregulation and cell type proportion.
#' 
#' \code{calc_dist} groups samples according to their deregulation status 
#' for a given gene. It then computes for each cell type different metrics to 
#' evaluate the distance between the two groups of cell proportions. An important 
#' distance suggests a link between the gene expression and the cell type proportion. \cr 
#' Metrics: 
#' \tabular{lll}{
#' \strong{Metric}          \tab \strong{Tested parameter} \tab \strong{Computed parameter}\cr 
#' Kantorovich distance     \tab Transportation           \tab Distance\cr
#' Student's t-test         \tab Mean comparison          \tab -log10(p-value)\cr
#' Kolmogorov-Smirnov test  \tab Distribution comparison  \tab -log10(p-value)\cr
#' }
#'
#' @param binary_penda The binary matrix of deregulation 0/1 for each gene and 
#'   each sample, can be obtain with the \code{pre_treat} function. 
#'   Dimension n (genes) * p (samples).
#' @param A The matrix of cell type proportions for each sample. 
#'   Dimension p (samples) * k (cell types).
#' 
#' @importFrom stats ks.test sd t.test
#' @importFrom ptlmapper kantorovich
#' @import progress
#' 
#' @return A data frame with for each gene and each cell type the results of the 3 metrics.
#' @export
#'
calc_dist = function(binary_penda, A){
  options(warn = -1)
  pb <- progress::progress_bar$new(format = "Running RiTMIC::calc_dist [:bar] :current/:total (:percent)", total = nrow(binary_penda), clear = FALSE, width = 80)
  
  res_dereg = c()
  for (g in rownames(binary_penda)) {
    pb$tick()
    gene = binary_penda[g, ]
    p_dereg = which(gene != 0)
    p_zero = which(gene == 0)
    if (length(p_dereg) != 0 & length(p_zero) != 0) {
      for (t in 1:nrow(A)) {
        x = A[t, p_dereg]
        y = A[t, p_zero]
        kanto_dist = ptlmapper::kantorovich(x, y)
        student = -log10(t.test(x, y)$p.value)
        statstud = t.test(x, y)$statistic
        ks = ks.test(x, y)
        res_dereg = rbind(res_dereg, c(g, length(p_dereg), 
                                       length(p_zero), t, kanto_dist, student, ks$statistic, 
                                       -log10(ks$p.value), statstud))
      }
    }
  }
  colnames(res_dereg) = c("Gene", "nb_dereg", "nb_zero", "cell_type", 
                          "kanto", "pval student", "ks d", "ks pvalue", "stat student")
  df = data.frame(genes = res_dereg[, 1], 
                  type = factor(res_dereg[, 4]), 
                  kanto = as.numeric(res_dereg[, 5]), 
                  st_pval = as.numeric(res_dereg[, 6]), 
                  ks_d = as.numeric(res_dereg[, 7]), 
                  ks_pval = as.numeric(res_dereg[, 8]), 
                  st_st = as.numeric(res_dereg[, 9]), stringsAsFactors = F)
  options(warn = 0)
  return(df)
}

#' Correlation between gene expression and cell type proportion.
#' 
#' \code{calc_corr} compute the correlation between the gene expression in 
#' tumors and the micro-environment proportion.
#'
#' @param D_cancer The matrix with the gene expression in each tumor.
#' @param A The matrix of cell type proportions for each sample.
#' 
#'
#' @return A data frame with for each gene and each cell type the
#'   correlation between the gene expression and the cell type proportion.
#' @export 
calc_corr = function(D_cancer, A) {
  
  options(warn = -1)
  res_cor = c()
  pb <- progress::progress_bar$new(
    format = "Running RiTMIC::calc_corr [:bar] :current/:total (:percent) in :elapsed",
    total = nrow(D_cancer), clear = FALSE, width = 80)
  for(g in rownames(D_cancer)){
    pb$tick()
    for(t in 1:nrow(A)){
      c = cor(D_cancer[g, ], A[t, ])
      res_cor = rbind(res_cor, c(g, t, c))
    }
  }
  options(warn = 0)
  df = data.frame(genes = res_cor[, 1],
                  type =  factor(res_cor[, 2]),
                  corr = as.numeric(res_cor[, 3]), 
                  stringsAsFactors = F)
  return(df)
}
