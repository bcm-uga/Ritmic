#' Simulate cell lines proportion distribution in tumors
#' 
#' @description \code{simu_A} generates the \strong{matrix A} : \cr Proportion of cell lines distribution per tumor \cr \cr  Dimensions: n x alpha 
#' @details Cell lines distributions are obtained with \code{\link[gtools]{rdirichlet}} method
#' @param n patients number
#' @param alpha cell type proportions of size k in a concatenated vector
#' 
#' @seealso \code{\link[gtools]{rdirichlet}}
#' 
#' @return Proportion of different cell lines in a matrix 
#' @export 
#'
#' @examples simu_A(20,c(1.5, 4.5, 3))
simu_A = function(n, alpha = c(1.5, 4.5, 1, 3)){
  tmp_mix_Lvar = t(gtools::rdirichlet(n = n, alpha = alpha))
}

#' Simulate a tumor micro-environment 
#' @description Generate a \strong{matrix T} with random \strong{gene expression} distributions in \strong{n} tumor cell lines \cr Matrix dimension : n x gene_name
#' @details gene_name is imported from the tumor_RNAseq data \cr Random gene expression distributions is performed by \code{\link{rnorm}} \cr Probability of gene expression in a cell line is computed by \code{\link[mixtools]{normalmixEM}}
#' @param tumor_RNAseq RNAseq dataset to simulate a random T matrix
#' @param n Supposed cell types number 
#' 
#' @importFrom stats rnorm
#' @importFrom mixtools normalmixEM
#'
#' @return Matrix A: Gene expressions distributions per cell line  
#' @export
#'
#' @examples
simu_T_cancer = function(tumor_RNAseq, n){
  # n_simu return the difference between the **n integer and the col number** of the tumor_RNAseq
  n_simu = n - ncol(tumor_RNAseq)
  # T_res is assigned to an empty matrix with in columns : the number of cell types n parameter and in rows : gene expressions from tumor_RNAseq
  T_res = matrix(NA, ncol = n, nrow = nrow(tumor_RNAseq))
  # Fill the boxes with those from @tumor_RNAseq
  T_res[, 1:ncol(tumor_RNAseq)] = tumor_RNAseq
  # Store the rownames from tumor_RNAseq in the new matrix 
  rownames(T_res) = rownames(tumor_RNAseq)
  t = 0
  
  for(g in rownames(tumor_RNAseq)){
    # Recuperation of gene expression
    dist_g = tumor_RNAseq[g, ]
    # Genes that don't work (à corriger un jour)
    # Compute formulas if at least, the half of gene expression is not null
    if(sum(dist_g == 0) < 0.5*length(dist_g) & !g%in%c("CELF3", "ZCCHC24", "DNAH17", "MMP10", "UPK2", "PCAT1", "DPH6-AS1", "PRAF2", "CNTF", "LINC00940", "SNHG7", "LY6G6C", "C2CD4A", "HIST3H2BB", "HIST1H2AI", "TMEM201", "OPTC", "PALM3", "ZPBP2", "MROH2A", "TMEM255B", "AIFM3", "BEGAIN", "DDN", "SERTAD2", "ZFAS1", "KCNJ4", "CXorf58", "PF4", "TMED6", "SMG1", "C18orf25", "FLI1", "HABP2", "CDKAL1", "YARS2", "KHK", "SPOCD1", "IPO8", "GCH1", "AOC3", "IDO1", "SH3BGRL", "VGF", "ZBP1", "FABP3", "CDC20", "SCN4A", "C2CD5", "PDGFRL", "MGAT4A", "UPB1", "ASCL1", "AKAP4", "RHAG", "IGSF1", "SLCO1A2", "LHFPL4", "PGLYRP3",  "PANX2", "KIF6", "CA1", "FNDC1", "TBX4", "KCNC2", "ASPG", "SOAT2", "CYP7B1", "CLCNKB", "ZYG11A", "SNORD121A", "PSG4", "CDH19", "BMP6", "CD300LB", "SOGA3", "SPEG", "TGM1", "VRK1",  "KRT23")){
      #Computing the bimodal distribution
      # Recuperation of parameters necessary to maximum likelihood methods 
      dists_sep <- normalmixEM(dist_g, maxrestarts=10000, k= 2)
      # Probability of gene expression in cell lines 
      probs <- dists_sep$lambda
      # Final mean parameters
      m <- dists_sep$mu
      # Final standard deviations
      s <- dists_sep$sigma
      N <- 1e5
      # Random sample : on all probs obtained, each gene are sampled many times to be randomly affected 
      grp <- sample(length(probs), N, replace=TRUE, prob=probs)
      # Random generation for the normal distribution with mean=m and sd = s for the gene in for loop  
      x <- rnorm(N, m[grp], s[grp])
      # All negative values from the normal distribution of gene expression equals to 0 
      x[x<0] = 0
    } else {
      x = dist_g
      t = t+1
    }
    # The values from the matrix T_res are replaced by x with n_simu reordered  
    T_res[g, (ncol(tumor_RNAseq)+1):n] = sample(x, n_simu, replace=TRUE)
  }
  return(T_res)
}

#' Simple deregulations in tumors micro-environment 
#' @description 
#' Generation of deregulations with a simple model in the matrix T generated from \code{simu_T_cancer} with cell lines proportions from the return of \code{simu_A} method \cr\cr
#' Step 1: G genes are drawn\cr \cr 
#' Step 2: Over-expression induced to these genes : \cr \cr
#' newGeneExpression = geneExpression + x * cellTypeProportion
#' @param T_cancer RNAseq data
#' @param G Gene number(integer type) to be sampled 
#' @param A Matrix A : distribution of cell lines per tumor \code{simu_A} 
#' @param x Coefficient, by default it's settled to 100
#' 
#' @return List with : Deregulated matrix $T and deregulated genes in $g_immune and $g_fibro
#' @export
#'
#' @examples
corr_prop_s = function(T_cancer, G, A, x = 100){
  # Sampling of the gene names from T_cancer 
  genes = sample(nrow(T_cancer), G, replace = F)
  # First half part of selected genes represents the immune cell lines
  g_immune = genes[1:(G/2)]
  # The other half represents the tumor cell lines
  g_fibro =  genes[(G/2 +1):G]
  T_f = T_cancer
  for(n in 1:ncol(T_cancer)){
    # Capturing immune and fibro cells proportions in the n tumor from matrix @parameter A 
    immune = A[3, n]
    fibro = A[2, n]
    # Deregulation of random genes for fibro cell types and immune cell types : Adding the coefficient x multiply by the cell type proportion to induce over-expression 
    T_f[g_immune, n] =  T_f[g_immune, n]+ x*immune
    T_f[g_fibro, n] =  T_f[g_fibro, n]+ x*fibro
  }
  return(list(T = T_f, g_immune = g_immune, g_fibro = g_fibro ))
}

#' Complex deregulations in tumors micro-environment 
#' @description 
#' Deregulations with a complex model in the matrix T generated from \code{simu_T_cancer} with cell lines proportions from the return of \code{simu_A} method \cr \cr
#' Step 1: G genes are drawn\cr \cr 
#' Step 2: Checking : proportions of cell type immune or fibro is respectively greater than \strong{tres_i} and \strong{tres_f} ? \cr \cr
#' Step 3: For the cell lines (fibro, immune) having a proportion greater than threshold: \cr \cr
#' newGeneExpression = geneExpression * y
#' 
#' @param T_cancer The matrix T, for more documentation : \code{simu_T_cancer}
#' @param G Number of random Genes to be drawn
#' @param A The matrix A, for more documentation \code{simu_A}
#' @param y Overexpression coefficient to multiply the gene expression, {default} = 2
#' @param thres_i Threshold to induce overexpression in immune cells, {default} = 0.1
#' @param thres_f Threshold to induce overexpression in fibro cells, {default} = 0.45
#'
#' @return List with : Deregulated matrix $T and deregulated genes in $g_immune and $g_fibro
#' @export
#'
#' @examples
corr_prop_c = function(T_cancer, G, A, y = 2, thres_i = 0.1, thres_f = 0.45){
  genes = sample(nrow(T_cancer), G, replace = F)
  g_immune = genes[1:(G/2)]
  g_fibro =  genes[(G/2 +1):G]
  T_f = T_cancer
  for(n in 1:ncol(T_cancer)){
    immune = A[3, n]
    fibro = A[2, n]
    if(immune > thres_i){
      T_f[g_immune, n] =  T_f[g_immune, n] * y
    }
    if(fibro > thres_f){
      T_f[g_fibro, n] =  T_f[g_fibro, n] * y
    }
  }
  return(list(T = T_f, g_immune = g_immune, g_fibro = g_fibro ))
}

#On pioche g gènes au hasard pour fibro et immune, puis pour chaque patient on fait expression*(1 + z*proportion(immune ou fibro))
#' New deregulations in tumors micro-environment 
#' @description 
#' Generation of deregulations from a new model in the matrix T generated from \code{simu_T_cancer} with cell lines proportions from the matrix A \code{simu_A} \cr \cr
#' Step 1: G genes are drawn\cr \cr 
#' Step 2: For the cell lines, over-expressions to every G genes is induced by: \cr \cr
#' newGeneExpression = geneExpression + [(1 + z)* cellTypeProportion)]
#'
#' @param T_cancer The matrix T, for more documentation \code{simu_T_cancer}
#' @param G Number of random Genes to be drawn
#' @param A The matrix A, for more documentation \code{simu_A}
#' @param z Overexpression coefficient to multiply the gene expression, by default z = 2
#'
#' @return A list with : Deregulated matrix $T and deregulated genes in $g_immune and $g_fibro
#' @export
#'
#' @examples
corr_prop_n = function(T_cancer, G, A, z = 2){
  genes = sample(nrow(T_cancer), G, replace = F)
  g_immune = genes[1:(G/2)]
  g_fibro =  genes[(G/2 +1):G]
  T_f = T_cancer
  for(n in 1:ncol(T_cancer)){
    immune = A[3, n]
    fibro = A[2, n]
    T_f[g_immune, n] =  T_f[g_immune, n]*(1 + z*immune)
    T_f[g_fibro, n] =  T_f[g_fibro, n]*(1 + z*fibro)
  }
  return(list(T = T_f, g_immune = g_immune, g_fibro = g_fibro ))
}

#' Add noise to the analysis or simulations of the tumor micro-environments gene expressions
#' @description Applications of adding noise : \cr Simulations : Better reflection of the biological observations \cr
#' Biological observation : improve the generalization error to avoid over-learning  
#' @param data matrix_D with or without deregulations obtained from the combinations of matrix_T * matrix_A
#' @param mean 
#' @param sd 
#' @param val_min 
#' @param val_max 
#' 
#' @seealso RiTMIC::simu_T_cancer 
#' @seealso RiTMIC::corr_prop functions...
#'
#' @return
#' @export
#'
#' @examples
add_noise = function(data, mean = 0, sd = 0.1, val_min = 0, val_max = 1){
  noise = matrix(rnorm(prod(dim(data)), mean = mean, sd = sd), nrow = nrow(data))
  datam = data + noise
  datam[datam < val_min] = data[datam < val_min]
  datam[datam > val_max] = data[datam > val_max]
  return(datam)
}