#' Pre-treatment for statistical analysis
#' @description Matrix are sorted to have the common patients with only genes that are deregulated in more than a minimum number of patients  
#' @param penda_res Output from \code{\link[penda]{penda_test_1ctrl}}
#' @param patientNumber Minimum number of patients with the gene deregulation
#'
#' @return Matrix without non deregulated genes 
#' 
#' @export
pre_treat = function(penda_res,patientNumber){
  if (patientNumber >= ncol(penda_res$up_genes)){
    stop("The minimum of patient number needs to be less than the sample dataset")
  } else {
    penda_res <- abs(penda_res$up_genes - penda_res$down_genes)
    g_egaux = which(apply(penda_res, 1, function(g) {
      length(table(g)) == 1
    }))
    penda_res = penda_res[-g_egaux, ]
    
    g_sup = which(apply(penda_res, 1, function(g) {
      table(g)[1] > patientNumber & table(g)[2] > patientNumber
    }))
    penda_res = penda_res[g_sup, ]
    
    return (penda_res)
  }
}

#' Compute distances between two groups: deregulated genes versus regulated genes
#' @description  For each cell line and gene, four statistical tests are applied to evaluate the distance: 
#' \tabular{lll}{
#' Statistical test   \tab Tested parameter               \tab Result type\cr 
#' Kantorovich metric \tab Wasserstein distance           \tab distance\cr
#' Student test       \tab Mean comparison                \tab -log10(p-value)\cr
#' Kolmogorov-Smirnov \tab Repartition comparison         \tab distance and   -log(p-value)\cr
#' Correlation        \tab Variability between two groups \tab Correlation coefficient 
#' }
#'
#' @param pre_treat_output Output from \code{pre_treat} function
#' @param matrix_A Matrix A of proportion of cell lines distribution per tumor \cr \cr  Dimensions: n(cell lines) x alpha(cell lines proportion per patient) 
#' 
#' @importFrom stats ks.test sd t.test
#' @importFrom ptlmapper kantorovich
#' @import progress
#' 
#' @return Matrix of dimension geneNumber*cellLineNumber\cr \cr Each element is filled out by the four statistical tests: kanto, t-test, ks and cor 
#' @export
#'
calc_dist = function(pre_treat_output,matrix_A){
  options(warn = -1)
  res_dereg = c()
  genes = rownames(pre_treat_output)
  # Progress bar
  pb <- progress_bar$new(format = "  Running RiTMIC::calc_dist [:bar] :current/:total (:percent)",total = nrow(pre_treat_output), clear = FALSE, width= 80)
  for(g in rownames(pre_treat_output)){
    pb$tick()
    gene = pre_treat_output[g, ]
    
    p_dereg = which(gene != 0)
    p_zero = which(gene == 0)
    
    if(length(p_dereg) != 0 & length(p_zero) != 0){
      for(t in 1:nrow(matrix_A)){
      
        x = matrix_A[t, p_dereg]
        y = matrix_A[t, p_zero]
      
        kanto_dist = kantorovich(x, y)
        student = -log10(t.test(x, y)$p.value)
        ks = ks.test(x, y)
        cvx = sd(x) / mean(x)
        cvy = sd(y) / mean(y)
        
        res_dereg = rbind(res_dereg, c(g, length(p_dereg), length(p_zero), t, kanto_dist, student, ks$statistic, -log10(ks$p.value), cvx, cvy))
      }
    }
  }
  
  colnames(res_dereg) = c("Gene", "nb_dereg", "nb_zero", "type", "kanto", "student", "ks d", "ks pvalue", "cvx", "cvy")
  
  df =  data.frame(genes = res_dereg[, 1],
                   kanto = as.numeric(res_dereg[, 5]),
                   type =  factor(res_dereg[, 4]),
                   student = as.numeric(res_dereg[, 6]),
                   ks = as.numeric(res_dereg[, 8]),  stringsAsFactors = F)
  options(warn = 0)
  return(list(df = df, genes = genes))
}


#Si on représente graphiquement ces matrices, on voit nettement les gènes qui donnent une plus grande distance entre les deux groupes pour chaque type : ces gènes doivent donc être propre au type cellulaire.

#Si pour un gène et un type cellulaire il y a une grande distance entre le groupe qui ne l'exprime pas et le groupe qui le surexprime, ça signifie que ce gène est probablement un marqueur de ce type cellulaire.

#Kanto : 75% des distances < 0.06

#Student : 75% des -log10(pvalue) < 0.8

#Ks : 75% des distances < 0.30

#' Compute the confusion matrix   
#' 
#' @description Compute the confusion matrix(TP,FP,TN,FN) necessary to plot ROC curves 
#'  
#' @param values statistical test used
#' @param genes gene names on the matrix D
#' @param genes_fibro gene expression in fibro cell line
#' @param genes_immune gene expression in immune cell line
#' @param pval list of p-values observed
#'
#' @return In a list: True positive, False positive, True negative and False negative with associated rates
#' @export  
compute_1_res = function(values, genes, genes_fibro, genes_immune, pval){
  # Recuperation of the deregulated gene list
  names(values) = stringr::str_c(genes, rep(1:4))
  values = values[!is.na(values)]
  #True deregulation simulated in corr_prop functions
  genes_dereg = c(values[stringr::str_c(genes_fibro, rep(2))], values[stringr::str_c(genes_immune, rep(3))])
  genes_dereg = genes_dereg[!is.na(genes_dereg)]
  
  TP = sum(genes_dereg > -log10(pval))
  
  FN = length(genes_dereg) - TP
  
  FP =  sum(values > -log10(pval)) - TP
  
  TN =  length(values) - TP - FN - FP
  
  TPR  = TP / (TP + FN)
  FPR = FP / (TN + FP)
  return(c(TP, FP, TN, FN, FPR, TPR))
}

#' Load environment required to \code{plot_res} function
#'
#' @param matrix_D matrix_D created with \code{simu_D}
#' @param matrix_A The matrix_A from \code{simu_A}
#' 
#' @seealso \code{compute_1_res} \code{plot_res}
#'
#' @return correlation values between gene expressions
#' @export 
calc_corr <- function(matrix_D,matrix_A) {
  matrix_T <- matrix_D$T$T
  options(warn = -1)
  cor_T60_c = c()
  pb <- progress_bar$new(format = "  Running RiTMIC::pre_plot_res [:bar] :current/:total (:percent) in :elapsed",
                         total = nrow(matrix_T), clear = FALSE, width= 80)
  for(g in rownames(matrix_T)){
    for(t in 1:nrow(matrix_A)){
      c = cor(matrix_T[g, ], matrix_A[t, ])
      
      cor_T60_c = rbind(cor_T60_c, c(g, t, c))
    }
  }
  options(warn = 0)
  return(list(cor_T = cor_T60_c, matrix_D = matrix_D))
}

#' Check the PenDA enrichment: Plot ROC curves  
#' @description Compare non-enriched data by PenDA with a correlation test versus the PenDA efficiency with tests: \cr \itemize{
#' \item{Kolmogorov-Smirnov} \item{Student} \item{Kantorovich}
#' } 
#'
#' @param calc_corr_output Output from \code{pre_plot_res} function
#' @param calc_dist_output Output from \code{calc_dist} function
#' @param graph_title Title of the ROC curve graph 
#' 
#' @seealso \code{as_def_res} \code{pre_plot_res} \code{compute_1_res}
#'
#'@return ROC curves
#' @export
plot_res = function(calc_corr_output, calc_dist_output, graph_title){
  "FPR" <- "TPR" <- "metrique" <- c()
  matrix_D <- calc_corr_output$matrix_D$matrix_D
  T <- calc_corr_output$matrix_D$T
  genes_c = c(rownames(T$T)[T$g_immune], rownames(T$T)[T$g_fibro])
  genes = genes_c[(genes_c %in% rownames(matrix_D))]
  cor_T <- calc_corr_output$cor_T
  pvalues <- c(0, 0.00005, 0.0001, 0.0005, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25,  0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
  g_fibro = unique(genes[genes%in%rownames(T$T)[T$g_fibro]])
  g_immune = unique(genes[genes%in%rownames(T$T)[T$g_immune]])
  df <- calc_dist_output$df
  res_ks = sapply(pvalues, function(x){
    compute_1_res(df$ks, df$genes, g_fibro, g_immune, x)
  })
  
  res_st = sapply(pvalues, function(x){
    compute_1_res(df$student, df$genes, g_fibro, g_immune, x)
  })
 
  res_kanto = sapply(pvalues, function(x){
    compute_1_res(df$kanto, df$genes, g_fibro, g_immune, x)
  })
  
  res_cor = sapply(pvalues, function(x){
    compute_1_res(abs(as.numeric(cor_T[, 3])), cor_T[, 1], rownames(T$T)[T$g_fibro], rownames(T$T)[T$g_immune], x)
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
