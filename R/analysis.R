#' Pre-treatment for the link between gene deregulation and microenvironment composition
#' @description Matrix are sorted to conserved only the genes with a different deregulation statut in more than thres_p samples  
#' @param penda_res The result of the Penda method with the deregulation of each gene for each sample  \code{\link[penda]{penda_test}}
#' @param thres_p Minimum number of samples for the deregulated and not deregulated groups. 
#'
#' @return A binary matrix with the genes that can be used for the statistical analysis and their deregulation status for each sample. 
#' 
#' @export
pre_treat = function(penda_res, thres_p){
  if (thres_p >= ncol(penda_res$up_genes)){
    stop("You can't reach a minimal number of samples higher than the total number of samples")
  } else {
    penda_res <- abs(penda_res$up_genes - penda_res$down_genes)
    
    #We removed the genes with the same deregulation status in all the samples
    g_egaux = which(apply(penda_res, 1, function(g) {
      length(table(g)) == 1
    }))
    penda_res = penda_res[-g_egaux, ]
    
    #We removed the genes with less than thres_p samples in each deregulation group
    g_sup = which(apply(penda_res, 1, function(g) {
      table(g)[1] > thres_p & table(g)[2] > thres_p
    }))
    penda_res = penda_res[g_sup, ]
    
    return (penda_res)
  }
}

#' Compute metrics between the two groups: samples with deregulated genes versus not deregulated genes
#' @description  For each cell type and each gene, 3 statistical tests are applied to evaluate the distances bewteen the cell type proportion of the two groups of deregulation: 
#' \tabular{lll}{
#' Statistical test   \tab Tested parameter               \tab Result type\cr 
#' Kantorovich        \tab Kantorovich distance           \tab distance\cr
#' Student t-test     \tab Mean comparison                \tab -log10(p-value)\cr
#' Kolmogorov-Smirnov \tab Repartition comparison         \tab -log(p-value)\cr
#' }
#'
#' @param binary_penda The binary matrix of deregulation, can be obtain with the \code{pre_treat} function
#' @param A The matrix of the different cell type proportion in each sample \cr \cr  Dimensions: k (cell types) * p (samples) 
#' 
#' @importFrom stats ks.test sd t.test
#' @importFrom ptlmapper kantorovich
#' @import progress
#' 
#' @return A data.frame with for each gene (genes) and each celle type (type), the result of the 3 metrics
#' @export
#'
calc_dist = function(binary_penda, A){
  options(warn = -1)
  res_dereg = c()
  genes = rownames(binary_penda)
  # Progress bar
  pb <- progress_bar$new(format = "  Running RiTMIC::calc_dist [:bar] :current/:total (:percent)",total = nrow(binary_penda), clear = FALSE, width= 80)
  
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
        cvx = sd(x) / mean(x)
        cvy = sd(y) / mean(y)
        
        res_dereg = rbind(res_dereg, c(g, length(p_dereg), length(p_zero), t, kanto_dist, student, ks$statistic, -log10(ks$p.value), cvx, cvy))
      }
    }
  }
  
  colnames(res_dereg) = c("Gene", "nb_dereg", "nb_zero", "cell_type", "kanto", "student", "ks d", "ks pvalue", "cvx", "cvy")
  
  df =  data.frame(genes = res_dereg[, 1],
                   type =  factor(res_dereg[, 4]),
                   kanto = as.numeric(res_dereg[, 5]),
                   student = as.numeric(res_dereg[, 6]),
                   ks = as.numeric(res_dereg[, 8]),  stringsAsFactors = F)
  options(warn = 0)
  return(df)
}

#' Compute the correlation between the gene expression in tumors and the micro-environment proportion
#'
#' @param D_cancer The matrix with the gene expression in each tumor
#' @param A The matrix of the different cell type proportion in each sample
#' 
#' @seealso \code{compute_1_res} \code{plot_res}
#'
#' @return A data.frame with for each gene and each cell type the correlation between the gene expression and the cell type proportion 
#' @export 
calc_corr <- function(D_cancer, A) {
  
  options(warn = -1)
  res_cor = c()
  pb <- progress_bar$new(format = "  Running RiTMIC::pre_plot_res [:bar] :current/:total (:percent) in :elapsed",
                         total = nrow(D_cancer), clear = FALSE, width= 80)
  for(g in rownames(D_cancer)){
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


#' For one p-value, test if the we find genes deregulated in the simulations  
#' 
#' @description Compute the confusion matrix(TP,FP,TN,FN) necessary to plot ROC curves 
#'  
#' @param values the result of the test for each gene (ex: distance value, p-value, etc.)
#' @param genes names of the genes analysed
#' @param genes_dereg names of the genes deregulated in the simulations
#' @param pval the threshold of interest, can be a -log10(pval) or a distance
#'
#' @return A vector with the number of TP, FP, TN, FN, the FPR and the TPR
#' @export  
eval_results = function(values, genes, genes_dereg, pval){
  # Recuperation of the deregulated gene list
  names(values) = genes
  values = values[!is.na(values)]
  #True deregulation simulated in corr_prop functions
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


#' Check the PenDA enrichment: Plot ROC curves  
#' @description Compare non-enriched data by PenDA with a correlation test versus the PenDA efficiency with tests: \cr \itemize{
#' \item{Kolmogorov-Smirnov} \item{Student} \item{Kantorovich}
#' } 
#'
#' @param corr_res Output from \code{calc_corr} function for the cell type deregulated
#' @param dist_res Output from \code{calc_dist} function for the cell type deregulated
#' @param genes_dereg The names of genes deregulated in the simulation
#' @param pvalues The threshold for the ROC curve
#' @param graph_title Title of the ROC curve graph 
#' 
#' @seealso \code{as_def_res} \code{pre_plot_res} \code{compute_1_res}
#'
#'@return ROC curves
#'@export
plot_res = function(corr_res, dist_res, genes_dereg, pvalues = c(0, 0.1, 0.5, 0.8), graph_title){
  
  genes_dereg = unique(genes_dereg[genes_dereg%in%corr_res$genes])
  res_ks = sapply(pvalues, function(x){
    Ritmic::eval_results(dist_res$ks, dist_res$genes, genes_dereg, x)
  })
  
  res_st = sapply(pvalues, function(x){
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
  
  plot = ggplot(df, aes(x = FPR, y = TPR, color = metrique, group = metrique)) +
    geom_point() + 
    geom_line() + 
    theme_minimal()+
    labs(title = graph_title) 
  
  return(plot)
}
