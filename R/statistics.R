
#' Pretreatment for statistical analysis
#'
#' @param penda_res Output from \code{\link[penda]{penda_test_1ctrl}}
#'
#' @return
#' @export
#'
pre_treat = function(penda_res){
  g_egaux = which(apply(penda_res, 1, function(g) {
    length(table(g)) == 1
  }))
  penda_res = penda_res[-g_egaux, ]
  
  g_sup = which(apply(penda_res, 1, function(g) {
    table(g)[1] > 5 & table(g)[2] > 5
  }))
  penda_res = penda_res[g_sup, ]
}


#On étudie la distance entre deux groupes, qui contiennent les gènes dérégulés VS les gènes inchangés, dans 100 individus chacun minimum.

#On passe chaque gène de Penda un par un. Pour chaque type cellulaire, on calcule la distance des deux groupes. On obtient pour chacun des trois paramètres une matrice de distances de taille g*k avec pour chaque gène et pour chaque type cellulaire les valeurs : 
  
#  * Distance de kantorovich
#* -log10(p-valeur) du test de Student
#* Distance et -log10(p-valeur) du test de Kolmogorov-Smirnov
#* cv de chaque groupe

#' Compute distances between two groups
#'
#' @param penda_res output from \code{pre_treat} function
#'
#' @return
#' @export
#'
#' @examples
calc_dist = function(penda_res){
  penda_1v1c = abs(penda_res_1v1c$up_genes - penda_res_1v1c$down_genes)
  options(warn = -1)
  res_dereg = c()
  
  #Pour chaque gène
  for(g in rownames(penda_res)){
    gene = penda_res[g, ]
    
    p_dereg = which(gene != 0)
    p_zero = which(gene == 0)
    
    if(length(p_dereg) != 0 & length(p_zero) != 0){
      #Pour chaque type cellulaire
      for(t in 1:nrow(A_60p)){
        #On fait les deux groupes
        x = A_60p[t, p_dereg]
        y = A_60p[t, p_zero]
        
        kanto_dist = kantorovich(x, y)
        student = -log10(t.test(x, y)$p.value)
        ks = ks.test(x, y)
        cvx = sd(x) / mean(x)
        cvy = sd(y) / mean(y)
        
        res_dereg = rbind(res_dereg, c(g, length(p_dereg), length(p_zero), t, kanto_dist, student, ks$statistic, -log10(ks$p.value), cvx, cvy))
      }
    }
  }
  options(warn = 0)
  colnames(res_dereg) = c("Gene", "nb_dereg", "nb_zero", "type", "kanto", "student", "ks d", "ks pvalue", "cvx", "cvy")
  return(res_dereg)
}


#Si on représente graphiquement ces matrices, on voit nettement les gènes qui donnent une plus grande distance entre les deux groupes pour chaque type : ces gènes doivent donc être propre au type cellulaire.

#Si pour un gène et un type cellulaire il y a une grande distance entre le groupe qui ne l'exprime pas et le groupe qui le surexprime, ça signifie que ce gène est probablement un marqueur de ce type cellulaire.

#Kanto : 75% des distances < 0.06

#Student : 75% des -log10(pvalue) < 0.8

#Ks : 75% des distances < 0.30

#' Observe deregulation between groups
#'
#' @param res_dereg the output from \code{calc_dist}
#'
#' @return Results 
#' @export
#'
#' @examples
as_df_res = function(res_dereg){
  df =  data.frame(genes = res_dereg[, 1],
                   kanto = as.numeric(res_dereg[, 5]),
                   type =  factor(res_dereg[, 4]),
                   student = as.numeric(res_dereg[, 6]),
                   ks = as.numeric(res_dereg[, 8]),  stringsAsFactors = F)
  return(df)
}

#' P values and ROC curves
#'
#' @param values 
#' @param genes 
#' @param genes_fibro 
#' @param genes_immune 
#' @param pval 
#'
#' @return
#' @export
#'
#' @examples
compute_1_res = function(values, genes, genes_fibro, genes_immune, pval){
  names(values) = stringr::str_c(genes, rep(1:4))
  values = values[!is.na(values)]
  
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

#' Plot the ROC curves of the ks, student, kantorovitch and correlation
#'
#' @param genes_f 
#' @param df_res 
#' @param cor_T 
#' @param titre_graph 
#'
#' @return
#' @export
#'
#' @examples
plot_res = function(genes_f, T, df_res, cor_T, titre_graph){
  g_fibro = unique(genes_f[genes_f%in%rownames(T$T)[T$g_fibro]])
  g_immune = unique(genes_f[genes_f%in%rownames(T$T)[T$g_immune]])
  res_ks = sapply(pvalues, function(x){
    compute_1_res(df_res$ks, df_res$genes, g_fibro, g_immune, x)
  })
  res_st = sapply(pvalues, function(x){
    compute_1_res(df_res$student, df_res$genes, g_fibro, g_immune, x)
  })
  res_kanto = sapply(pvalues, function(x){
    compute_1_res(df_res$kanto, df_res$genes, g_fibro, g_immune, x)
  })
  res_cor = sapply(pvalues, function(x){
    compute_1_res(abs(as.numeric(cor_T[, 3])), cor_T[, 1], rownames(T$T)[T$g_fibro], rownames(T$T)[T$g_immune], x)})
  df = data.frame(pval = as.factor(rep(pvalues)),
                  FPR = c(res_ks[5, ], res_st[5, ], res_kanto[5, ], res_cor[5, ]),
                  TPR = c(res_ks[6, ], res_st[6, ], res_kanto[6, ], res_cor[6, ]),
                  metrique = rep(c("ks", "student", "kanto", "correlation"), each = length(pvalues)))
  
  plot = ggplot(df, aes(x = FPR, y = TPR, color = metrique, group = metrique)) +
    geom_point() + 
    geom_line() + 
    theme_minimal()+
    labs(title = titre_graph) 
  
  return(plot)
}

#' Load environment required to \code{plot_res} function
#'
#' @param matrix_T The matrix_T output from corr_prop functions
#' @param matrix_A The matrix_A from \code{simu_A}
#' 
#' @seealso \code{plot_res}
#'
#' @return correlation values between gene expressions
#' @export 
#'
#' @examples
pre_plot_res <- function(matrix_T,matrix_A,res_deg_output) {
  genes = c(rownames(matrix_T$T)[matrix_T$g_immune], rownames(matrix_T$T)[matrix_T$g_fibro])
  genes_f = genes_c[(genes %in% rownames(res_deg_output))]
  cor_T60_c = c()
  options(warn = -1)
  for(g in rownames(matrix_T$T)){
    for(t in 1:nrow(matrix_A)){
      #On fait les deux groupes
      c = cor(matrix_T$T[g, ], matrix_A[t, ])
      cor_T60_c = rbind(cor_T60_c, c(g, t, c))
    }
    return(list(gene = genes_f, corr = cor_T60_c))
  }
  options(warn = 0)
  return(cor_T60_c)
}
