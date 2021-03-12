#' Simulation of the cell types proportion distribution (A matrix)
#' 
#' @description \code{simu_A} generates the \strong{matrix A} : \cr Proportion of cell types in samples \cr \cr  Dimensions: k*p
#' @details Cell lines distributions are obtained with \code{\link[gtools]{rdirichlet}} method
#' @param n Number of samples 
#' @param alpha Different proportion of cell types in a vector of size k (number of cell types)
#' 
#' @seealso \code{\link[gtools]{rdirichlet}}
#' 
#' @return This function return a matrix of size k*p with the cell type proportions.
#' @export 
#'
#' @examples simu_A(n = 20, alpha = c(1.5, 4.5, 3))
simu_A = function(n, alpha = c(1.5, 4.5, 1, 3)){
    tmp_mix_Lvar = t(gtools::rdirichlet(n = n, alpha = alpha))
    return(tmp_mix_Lvar)
}

#' Simulation of new tumors RNAseq profiles based on already known tumors expression
#' @description Generate a \strong{matrix T} with a \strong{gene expression} distributions for each \strong{n} tumor \cr Matrix dimension : n*p
#' @details Genes names are imported from the tumor_RNAseq data \cr Bimodal distribution of each gene expression is infered with \code{\link[mixtools]{normalmixEM}} \cr When the bimodal distribution doesn't work (too many 0, extreme value, etc.), the gene expression is sampled from existing values.
#' @param tumor_RNAseq Gene expression of tumors samples
#' @param n Final number of tumors
#' 
#' @importFrom stats rnorm
#' @importFrom mixtools normalmixEM
#' @import progress
#' @importFrom utils capture.output
#' 
#' @seealso \code{\link[mixtools]{normalmixEM}}
#'
#' @return This function return a matrix of size n*p with a different gene expression profile for each tumor
#' @export
simu_T_cancer = function(tumor_RNAseq, n){
  # Parameter checking 
  if(n <= ncol(tumor_RNAseq)) {
    stop("n : you already have the good number of tumors, please set n > ncol(tumor_RNAseq)")
  } else {
    
    # Number of tumors to simulate
    n_simu = n - ncol(tumor_RNAseq)
    # T_res is assigned to an empty matrix of size n*p
    T_res = matrix(NA, ncol = n, nrow = nrow(tumor_RNAseq))
    # Fill the already known values with those from @tumor_RNAseq
    T_res[, 1:ncol(tumor_RNAseq)] = tumor_RNAseq
    rownames(T_res) = rownames(tumor_RNAseq)
    t = 0
    
    # Progress bar
    pb <- progress_bar$new(format = "  Running RiTMIC::simu_T_cancer [:bar] :current/:total (:percent) in :elapsed",total = nrow(tumor_RNAseq), clear = FALSE, width= 80)
    
    #For each gene,
    for(g in rownames(tumor_RNAseq)){
      pb$tick()
      # Extraction of the gene expression
      dist_g = tumor_RNAseq[g, ]
      #if the gene is not expressed in any sample : we don't try to compute the normalmixtest
      if(sum(dist_g) == 0){
        x = dist_g
      } else {
        dists_sep = NA
        #Infering the bimodal distribution
        res <- try(capture.output(dists_sep <- normalmixEM(dist_g, maxrestarts=1000, k= 2)), silent = T)
        #If it doesn't work (too many 0s, extreme value, etc.), we sample in the existing values
        if(inherits(res, "try-error") || any(is.na(dists_sep$mu))){
          x = dist_g
          t = t+1    
          # Else we sample in the computed distribution
        } else { 
          # Probability of gene expression in cell types 
          probs <- dists_sep$lambda
          m <- dists_sep$mu
          s <- dists_sep$sigma
          N <- 1e5
          # Random sampling on all probs obtained : each gene are sampled many times to be randomly affected 
          grp <- sample(length(probs), N, replace=TRUE, prob=probs)
          # Random generation for the normal distribution with mean=m and sd = s for the gene in for loop  
          x <- rnorm(N, m[grp], s[grp])
          # Cutoff : All negative values from the normal distribution of gene expression equals to 0 
          x[x < 0] = 0
        }
      }
      T_res[g, (ncol(tumor_RNAseq)+1):n] = sample(x, n_simu, replace=TRUE)
    }
  cat("Number of genes failing the test normalmixEM:",t)
  return(T_res)
  }
}

#' Simple simulation of the correlation between tumor micro-environment proportion and tumor gene expression
#' @description 
#' Generation of deregulations with a simple model in the matrix T generated from \code{simu_T_cancer} with cell lines proportions from the return of \code{simu_A} method \cr\cr
#' Step 1: G genes are drawn\cr \cr 
#' Step 2: Over-expression induced to these genes : \cr \cr
#' newGeneExpression = geneExpression + x * microEnvironmentProportion
#' @param T_cancer The matrix of gene expression in each tumor, can be generated with: \code{simu_T_cancer}
#' @param G Number of genes correlated with the micro-environment composition
#' @param A_ME The vector with the proportion of the micro-environment in each tumor. Can be generated with a row of: \code{simu_A} 
#' @param x Deregulation coefficient
#' 
#' @return List with : \cr Deregulated matrix $T and deregulated genes in $g_dereg
#' @export
corr_prop_s = function(T_cancer, G, A_ME, x = 100){
  if (ncol(T_cancer) != length(A_ME)){
    stop("Patient number is not equal between the vector A_ME and the matrix T_cancer")    
  } else {
    # Sampling of the gene names from T_cancer 
    genes = sample(nrow(T_cancer), G, replace = F)
    T_cancer[genes, ] =  T_cancer[genes, ] + x*A_ME
    return(list(T = T_cancer, g_dereg = genes))
  }
}

#' Simulation of the correlation between tumor micro-environment proportion and tumor gene expression based on a Threshold
#' @description 
#' #' Generation of deregulations with a model based on threshold in the matrix T generated from \code{simu_T_cancer} with cell lines proportions from the return of \code{simu_A} method \cr\cr
#' Step 1: G genes are drawn\cr \cr 
#' Step 2: Checking if the proportion of the ME is greater than \strong{tres_ME} \cr \cr
#' Step 3: If yes, \cr \cr
#' newGeneExpression = geneExpression * y
#' 
#' @param T_cancer The matrix of gene expression in each tumor, can be generated with: \code{simu_T_cancer}
#' @param G Number of genes correlated with the micro-environment composition
#' @param A_ME The vector with the proportion of the micro-environment in each tumor. Can be generated with a row of: \code{simu_A} 
#' @param y Deregulation coefficient
#' @param thres_ME Threshold to induce the overexpression
#'
#' @return List with : Deregulated matrix T and deregulated genes in g_dereg
#' @export
corr_prop_t = function(T_cancer, G, A_ME, y = 2, thres_ME = 0.1){
  if (ncol(T_cancer) != length(A_ME)){
    stop("Patient number is not equal between the vector A_ME and the matrix T_cancer")    
  } else {
    genes = sample(nrow(T_cancer), G, replace = F)
    T_f = T_cancer
    for(n in 1:ncol(T_cancer)){
      if(A_ME[n] > thres_ME){
        T_f[genes, n] =  T_f[genes, n] * y
      }
    }
    return(list(T = T_f, g_dereg = genes))
  }
}

#' Simulation of the correlation between tumor micro-environment proportion and tumor gene expression based on a Factor
#' @description 
#' Generation of deregulations with a model based on factor in the matrix T generated from \code{simu_T_cancer} with cell lines proportions from the matrix A \code{simu_A} \cr \cr
#' Step 1: G genes are drawn\cr \cr 
#' Step 2: The over-expressions of G genes is induced by: \cr \cr
#' newGeneExpression = geneExpression x [1 + (z x cellTypeProportion)]
#' 
#' @param T_cancer The matrix of gene expression in each tumor, can be generated with: \code{simu_T_cancer}
#' @param G Number of genes correlated with the micro-environment composition
#' @param A_ME The vector with the proportion of the micro-environment in each tumor. Can be generated with a row of: \code{simu_A} 
#' @param z Deregulation coefficient
#'
#' @return A list with : \cr Deregulated matrix $T and deregulated genes in $g_dereg
#' @export
corr_prop_f = function(T_cancer, G, A_ME, z = 2){
  if (ncol(T_cancer) != length(A_ME)){
    stop("Patient number is not equal between the vector A_ME and the matrix T_cancer")    
  } else {
    genes = sample(nrow(T_cancer), G, replace = F)
    T_cancer[genes, ] =  T_cancer[genes, ]*(1 + z*A_ME)
    return(list(T = T_cancer, g_dereg = genes))
  }
}

#' Add a Gaussian noise to a matrix
#' 
#' @description Applications of \code{add_noise} \cr \cr Simulations : Better reflection of the biological observations \cr
#' @param matrix_D The matrix where to add a noise  
#' @param mean Mean of noise 
#' @param sd Standard deviation of the noise
#' @param val_min Minimum value of the matrix
#' @param val_max Maximum value of the matrix
#' 
#' @return The matrix_D with a gaussian noise of mean = mean, sd = sd, without exceeding the lower threshold val_min and the upper threshold val_max
#' @export
add_noise = function(matrix_D, mean = 0, sd = 0.1, val_min = 0, val_max = 1){
  noise = matrix(rnorm(prod(dim(matrix_D)), mean = mean, sd = sd), nrow = nrow(matrix_D))
  datam = matrix_D + noise
  datam[datam < val_min] = matrix_D[datam < val_min]
  datam[datam > val_max] = matrix_D[datam > val_max]
  return(datam)
}

#' Build the D matrix: gene expression per sample
#' @description Simulation of a complex RNAseq matrix, obtained with the matrix multiplication of the A and T matrices, with a different tumor profile for each sample
#'
#' @param A The matrix of the different cell type proportion in each sample, can be generated with: \code{simu_A}. Note: the last row will be used for the tumor.
#' @param T The matrix with the gene expression profile of each cell type, except tumors 
#' @param T_cancer The matrix of gene expression in each tumor, can be generated with: \code{simu_T_cancer}
#' @param noise TRUE for adding a Gaussian noise on the product matrix, FALSE for not (see \code{add_noise})
#' @param mean Mean of the noise (see \code{add_noise})
#' @param sd Standard deviation of the noise (see \code{add_noise})
#' @param val_min Minimum value of the matrix after the noise (see \code{add_noise})
#' @param val_max Maximum value of the matrix after the noise (see \code{add_noise})
#'
#' @return The D matrix of size n*p with the complex gene expression for each sample. If noise = T, return a list of the D matrix before the noise (D) and after the noise (D_noise)
#' @export
simu_D <- function(A, T, T_cancer, noise = F, mean = 0, sd = 0.1, val_min = 0, val_max = 1){
  matrix_D = c()
  for(p in 1:ncol(A)){
    tumor = T_cancer$T[, p]
    D_p = cbind(T, tumor) %*% A[, p]
    matrix_D = cbind(matrix_D, D_p)
  }
  colnames(matrix_D) = colnames(A)
  if (isTRUE(noise)){
    D_noise = add_noise(matrix_D, mean = mean, sd = sd, val_min = val_min, val_max = val_max)
    matrix_D = list(D = matrix_D, D_noise = D_noise)
  }
  return(matrix_D)
}

