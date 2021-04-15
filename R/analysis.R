#' Pre-treatment of the Penda results
#' 
#' \code{pre_treat} takes as input the Penda result lists and 
#' combines them to make only one matrix of size genes*samples, with 1 for 
#' deregulated genes and 0 for conserved genes. Then, the matrix is sorted to 
#' conserved only the genes with a different deregulation status in more than 
#' thres_p samples. 
#' 
#' @param penda_res List of two matrices, $up_genes and $down_genes with 
#'   1 or 0 for each gene and each sample obtained by the Penda method. 
#'   See: \code{\link[penda]{penda_test}}
#' @param thres_p Minimum number of samples in each deregulation group to 
#'  conserved the gene. 
#'
#' @return A binary matrix of deregulation 0/1 with the genes that can be used 
#'  for the statistical analysis and their status for each sample. 
#' 
#' @export
pre_treat = function(penda_res, thres_p){
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
#' @return A data frame with for each gene and each cell type the result of the 3 metrics.
#' @export
#'
calc_dist = function(binary_penda, A){
  options(warn = -1)
  # Progress bar
  pb <- progress::progress_bar$new(
    format = "Running RiTMIC::calc_dist [:bar] :current/:total (:percent)",
    total = nrow(binary_penda), clear = FALSE, width = 80)
  
  res_dereg = c()
  #For each gene
  for(g in rownames(binary_penda)){
    pb$tick()
    gene = binary_penda[g, ]
    p_dereg = which(gene != 0)
    p_zero = which(gene == 0)
    
    #If different deregulation status exists
    if(length(p_dereg) != 0 & length(p_zero) != 0){
      #For each cell type
      for(t in 1:nrow(A)){
        x = A[t, p_dereg]
        y = A[t, p_zero]
      
        kanto_dist = ptlmapper::kantorovich(x, y)
        student = -log10(t.test(x, y)$p.value)
        ks = ks.test(x, y)
        
        res_dereg = rbind(res_dereg, c(g, length(p_dereg), length(p_zero), t, kanto_dist, student, ks$statistic, -log10(ks$p.value)))
      }
    }
  }
  
  colnames(res_dereg) = c("Gene", "nb_dereg", "nb_zero", "cell_type", "kanto", "student", "ks d", "ks pvalue")
  
  df =  data.frame(genes = res_dereg[, 1],
                   type =  factor(res_dereg[, 4]),
                   kanto = as.numeric(res_dereg[, 5]),
                   student = as.numeric(res_dereg[, 6]),
                   ks = as.numeric(res_dereg[, 8]),  stringsAsFactors = F)
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
#' @seealso \code{\link{compute_1_res}} ; \code{\link{plot_res}} 
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


#' For one p-value, test relevance of the results.
#' 
#' \code{eval_results} compares the computed results and the genes deregulated 
#' in the simulation to obtain the confusion matrix (TP, FP, TN, FN) necessary to 
#' plot ROC curves.
#'  
#' @param values The vector of values associated at each gene, can be p-values 
#'  of a statistical test, correlations, distances, etc. 
#'  See also: \code{\link[Ritmic]{calc_dist}}
#' @param genes The vector with the names of the genes. 
#' @param genes_dereg The vector with the names of the genes deregulated in the
#'  simulations.
#' @param pval The threshold of interest, adapted for the metric used in values. 
#'
#' @return A vector with the number of True Positive, False Positive, True 
#'   Negative and False Negative results, the False Positive Rate and the True 
#'   Positive Rate.
#' @export  
eval_results = function(values, genes, genes_dereg, pval){
  names(values) = genes
  values = values[!is.na(values)]
  
  #Values of the genes deregulated in simulations
  genes_dereg = values[genes_dereg]
  genes_dereg = genes_dereg[!is.na(genes_dereg)]
  
  TP = sum(genes_dereg > pval)
  
  FN = length(genes_dereg) - TP
  
  FP =  sum(values > pval) - TP
  
  TN =  length(values) - TP - FN - FP
  
  TPR  = TP / (TP + FN)
  FPR = FP / (TN + FP)
  return(c(TP, FP, TN, FN, FPR, TPR))
}


#' Plot the ROC curve of the detection of deregulated genes linked to micro-environment.
#' 
#' \code{plot_res} makes a ROC plot with the detection of genes deregulated in 
#' simulations. It makes one curve for the direct correlation between genes 
#' expression and micro-environment proportion, and three curves for the comparison
#' between PenDA results and micro-environment, with three metrics: \cr 
#' \itemize{
#'   \item{The pvalue of the Kolmogorov-Smirnov test} \item{The p-value of the Student's t-test} \item{The distance of Kantorovich}
#' } 
#'
#' @param corr_res The data frame output from \code{\link[Ritmic]{calc_corr}},
#'   with the correlation for each gene of the cell type deregulated.
#' @param dist_res The data frame output from \code{\link[Ritmic]{calc_dist}}, 
#'   with the metrics for each gene of the cell type deregulated.
#' @param genes_dereg The vector with the names of genes deregulated in the simulation.
#' @param pvalues The vector of thresholds to trace the ROC curve.
#' @param graph_title An optional title for the ROC curve plot.
#' 
#' @seealso \code{\link[Ritmic]{compute_1_res}}
#'
#' @return A ggplot object with the ROC curves. See \code{\link[ggplot2]{ggplot}}
#' @export
plot_res = function(corr_res, dist_res, genes_dereg, pvalues = c(0, 0.01, 0.05, 0.1, 1), graph_title = ""){
  
  genes_dereg = unique(genes_dereg[genes_dereg %in% corr_res$genes])
  res_ks = sapply(-log10(pvalues), function(x){
    Ritmic::eval_results(dist_res$ks, dist_res$genes, genes_dereg, x)
  })
  
  res_st = sapply(-log10(pvalues), function(x){
    Ritmic::eval_results(dist_res$student, dist_res$genes, genes_dereg, x)
  })
 
  res_kanto = sapply(pvalues, function(x){
    Ritmic::eval_results(dist_res$kanto, dist_res$genes, genes_dereg, x)
  })
  
  res_cor = sapply(pvalues, function(x){
    Ritmic::eval_results(abs(corr_res$corr), corr_res$genes, genes_dereg, x)
  })
  
  df = data.frame(pval = as.factor(rep(pvalues)),
                  FPR = c(res_ks[5, ], res_st[5, ], res_kanto[5, ], res_cor[5, ]),
                  TPR = c(res_ks[6, ], res_st[6, ], res_kanto[6, ], res_cor[6, ]),
                  metrique = rep(c("ks", "student", "kanto", "correlation"), each = length(pvalues)))
  
  plot = ggplot2::ggplot(df, aes(x = FPR, y = TPR, color = metrique, group = metrique)) +
    geom_point() + 
    geom_line() + 
    theme_minimal()+
    labs(title = graph_title) 
  
  return(plot)
}
