#' Simulation of the A matrix
#' 
#' \code{simu_A} generates the \strong{A matrix}: \cr 
#' Different cell types proportions in each sample. \cr 
#' Dimensions: k*p
#' 
#' @details The cell type distribution in a sample is simulated by a Dirichlet 
#' distribution. See: \code{\link[gtools]{rdirichlet}}
#' 
#' @param n The number of samples.
#' @param alpha The vector of size kwith the different cell types proportion.
#' 
#' @seealso \code{\link[gtools]{rdirichlet}}
#' 
#' @return This function return a matrix of size k*p with the cell type proportions.
#' @export 
#'
simu_A = function(n, alpha = c(1.5, 4.5, 1, 3)){
    A_matrix = t(gtools::rdirichlet(n = n, alpha = alpha))
    return(A_matrix)
}

#' Simulation of new tumors RNAseq profiles based on already known tumors expression
#' 
#' \code{simu_T_cancer} generates a \strong{gene expression} matrix, with a different
#' expression for each of the \strong{n} tumors. \cr 
#' Matrix dimension : n*p
#' 
#' @details Genes names are imported from the tumor_RNAseq data \cr
#'  Bimodal distribution of each gene expression is infered with the
#'  \code{\link[mixtools]{normalmixEM}} function \cr 
#'  When the bimodal distribution doesn't work (too many 0, extreme value, etc.), 
#'  the gene expression is sampled from existing values.
#'  
#' @param tumor_RNAseq The matrix of gene expression for the already known tumors samples.
#' @param n Final number of tumors.
#' 
#' @seealso \code{\link[mixtools]{normalmixEM}}
#'
#' @return This function return a matrix of size n*p with a different gene 
#' expression profile for each tumor.
#' @export
simu_T_cancer = function(tumor_RNAseq, n){
  # Parameter checking 
  if(n <= ncol(tumor_RNAseq)) {
    stop("n : you already have the good number of tumors, 
         please set n > ncol(tumor_RNAseq)")
  } else {
    
    # Number of tumors to simulate
    n_simu = n - ncol(tumor_RNAseq)
    # T_res is assigned to an empty matrix of size n*p
    T_res = matrix(NA, ncol = n, nrow = nrow(tumor_RNAseq))
    # Fill the already known values with those from tumor_RNAseq
    T_res[, 1:ncol(tumor_RNAseq)] = tumor_RNAseq
    rownames(T_res) = rownames(tumor_RNAseq)
    t = 0
    
    # Progress bar
    pb <- progress::progress_bar$new(
      format = "Running RiTMIC::simu_T_cancer [:bar] :current/:total (:percent) in :elapsed",
      total = nrow(tumor_RNAseq), clear = FALSE, width = 80)
    
    #For each gene,
    for(g in rownames(tumor_RNAseq)){
      pb$tick()
      # Extraction of the gene expression
      dist_g = tumor_RNAseq[g, ]
      #if the gene is not expressed in any sample: we don't try to compute the normalmixtest
      if(sum(dist_g) == 0){
        x = dist_g
      } else {
        dists_sep = NA
        #Infering the bimodal distribution
        res = try(capture.output(
          dists_sep <- mixtools::normalmixEM(dist_g, maxrestarts = 1000, k = 2)), 
          silent = T)
        #If it doesn't work (too many 0s, extreme value, etc.), we sample in the existing values
        if(inherits(res, "try-error") || any(is.na(dists_sep$mu))){
          x = dist_g
          t = t + 1    
          # Else we sample in the computed distribution
        } else { 
          # Probability of gene expression in cell types 
          probs <- dists_sep$lambda
          m <- dists_sep$mu
          s <- dists_sep$sigma
          N <- 1e5
          # Random sampling on all probs obtained: each gene are sampled many times to be randomly affected 
          grp <- sample(length(probs), N, replace = TRUE, prob = probs)
          # Random generation for the normal distribution with mean = m and sd = s for the gene in for loop  
          x <- rnorm(N, m[grp], s[grp])
          # Cutoff : All negative values from the normal distribution of gene expression equals to 0 
          x[x < 0] = 0
        }
      }
      T_res[g, (ncol(tumor_RNAseq)+1):n] = sample(x, n_simu, replace=TRUE)
    }
  cat("Number of genes sampled directly in current expression values:",t)
  return(T_res)
  }
}


#' Simulation of the correlation between gene expression and tumor 
#' micro-environment proportion based on a threshold.
#' 
#' \code{corr_prop_t} simulates a link between the gene expression of \code{G} 
#' genes and the micro-environment (ME) proportion, only if this proportion is greater 
#' than \code{thres_ME}. 
#' 
#' \code{G} genes are drawn, then for each tumor in \code{T_cancer}: \cr
#' If the proportion of the ME in this tumor is greater than \code{thres_ME}: \cr
#' New gene expression = gene expression * \code{dereg_coeff}.
#' 
#' @param T_cancer The matrix of gene expression in each tumor. Can be 
#'   generated with: \code{simu_T_cancer}.
#' @param A_ME The vector with the proportion of the micro-environment in each tumor.
#'   Can be generated with a row of: \code{simu_A}. 
#' @param G The number of genes correlated with the micro-environment composition.
#' @param dereg_coeff The deregulation coefficient.
#' @param thres_ME The threshold to induce the overexpression.
#'
#' @return A list of $T the deregulated T matrix and $g_dereg the deregulated genes.
#' @export
corr_prop_t = function(T_cancer, A_ME, G = 100, dereg_coeff = 2, thres_ME = 0.1){
  if (ncol(T_cancer) != length(A_ME)){
    stop("Patient number is not equal between the vector A_ME and the matrix T_cancer")    
  } else {
    genes = sample(nrow(T_cancer), G, replace = F)
    T_f = T_cancer
    for(n in 1:ncol(T_cancer)){
      if(A_ME[n] > thres_ME){
        T_f[genes, n] = T_f[genes, n] * dereg_coeff
      }
    }
    return(list(T = T_f, g_dereg = rownames(T_cancer)[genes]))
  }
}


#' Simulation of the correlation between gene expression and tumor 
#' micro-environment proportion based on a factor.
#' 
#' \code{corr_prop_f} simulates a link between the gene expression of \code{G} 
#' genes and the micro-environment (ME) proportion, by multiplying the gene expression
#' by a factor of the ME proportion.
#' 
#' \code{G} genes are drawn, then for each tumor in \code{T_cancer}: \cr
#' New gene expression = gene expression * (1 + \code{dereg_coeff} * ME proportion)
#' 
#' @param T_cancer The matrix of gene expression in each tumor. Can be 
#'   generated with: \code{simu_T_cancer}.
#' @param A_ME The vector with the proportion of the micro-environment in each tumor.
#'   Can be generated with a row of: \code{simu_A}. 
#' @param G The number of genes correlated with the micro-environment composition.
#' @param dereg_coeff The deregulation coefficient.
#'
#' @return A list of $T the deregulated T matrix and $g_dereg the deregulated genes.
#' @export
corr_prop_f = function(T_cancer, A_ME, G = 100, dereg_coeff = 2){
  if (ncol(T_cancer) != length(A_ME)){
    stop("Patient number is not equal between the vector A_ME and the matrix T_cancer")    
  } else {
    genes = sample(nrow(T_cancer), G, replace = F)
    T_cancer[genes, ] =  T_cancer[genes, ] * (1 + dereg_coeff * A_ME)
    return(list(T = T_cancer, g_dereg = rownames(T_cancer)[genes]))
  }
}

#' Add a Gaussian noise to a matrix.
#' 
#' \code{add_noise} Apply a gaussian noise with a lower cap and an upper cap for 
#' the values.
#' 
#' @param D The matrix where to add a noise.
#' @param mean The mean of the noise.
#' @param sd The standard deviation of the noise.
#' @param val_min The minimum value for the final matrix.
#' @param val_max The maximum value of the final matrix.
#' 
#' @return A matrix corresponding to the D matrix with a gaussian noise.
#' @export
add_noise = function(D, mean = 0, sd = 0.1, val_min = 0, val_max = 1){
  noise = matrix(rnorm(prod(dim(D)), mean = mean, sd = sd), nrow = nrow(D))
  datam = D + noise
  datam[datam < val_min] = D[datam < val_min]
  datam[datam > val_max] = D[datam > val_max]
  return(datam)
}

#' Compute the D matrix
#' 
#' Simulation of a complex RNAseq matrix, obtained with the matrix multiplication 
#' of the A and T matrices, with a different tumor profile for each sample.
#'
#' @param A The matrix of the different cell type proportion in each sample.
#'   Can be generated with: \code{simu_A}. Note: the last row will be used for 
#'   the tumor cell type.
#' @param T The matrix with the gene expression profile of each cell type, except tumors.
#' @param T_cancer The matrix of gene expression in each tumor. Can be 
#'   generated with: \code{simu_T_cancer}.
#' @param noise A boolean: TRUE for adding a Gaussian noise on the product matrix, 
#'   FALSE for not. See \code{add_noise}).
#' @param mean If noise = TRUE, mean of the noise (see \code{add_noise}).
#' @param sd If noise = TRUE, standard deviation of the noise (see \code{add_noise}).
#' @param val_min If noise = TRUE, minimum value of the matrix after 
#'   the noise (see \code{add_noise}).
#' @param val_max If noise = TRUE, maximum value of the matrix after the noise
#'  (see \code{add_noise}).
#'
#' @return The D matrix of size n*p with the complex gene expression for each sample. 
#'  If noise = T, return a list of the D matrix before the noise (D) and after the noise (D_noise).
#' @export
simu_D <- function(A, T, T_cancer, noise = F, mean = 0, sd = 0.1, val_min = 0, val_max = 1){
  matrix_D = c()
  for(p in 1:ncol(A)){
    tumor = T_cancer$T[, p]
    D_p = cbind(T, tumor) %*% A[, p]
    matrix_D = cbind(matrix_D, D_p)
  }
  colnames(matrix_D) = colnames(A)
  if(isTRUE(noise)){
    D_noise = add_noise(matrix_D, mean = mean, sd = sd, val_min = val_min, val_max = val_max)
    matrix_D = list(D = matrix_D, D_noise = D_noise)
  }
  return(matrix_D)
}

