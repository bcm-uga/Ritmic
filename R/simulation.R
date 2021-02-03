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

#' Simulate an enhanced tumor micro-environment 
#' @description Generate a \strong{matrix T} with random \strong{gene expression} distributions in \strong{n} samples \cr Matrix dimension : n x gene_name
#' @details gene_name is imported from the tumor_RNAseq data \cr Random gene expression distributions is performed by \code{\link{rnorm}} \cr Probability of gene expression in a cell line is computed by \code{\link[mixtools]{normalmixEM}}
#' @param tumor_RNAseq RNAseq dataset to simulate a random T matrix
#' @param n Supposed sample number 
#' 
#' @importFrom stats rnorm
#' @importFrom mixtools normalmixEM
#' @import progress
#' @importFrom utils capture.output
#' 
#' @seealso \code{\link[mixtools]{normalmixEM}}
#'
#' @return Matrix A: \cr Gene expressions distributions per cell line  
#' @export
simu_T_cancer = function(tumor_RNAseq, n){
  # Parameter checking 
  if(n <= ncol(tumor_RNAseq)) {
    stop("n : the new simulation number needs to be greater than the previous dataset")
  } else {
    
    # Number of simulated tumors
    n_simu = n - ncol(tumor_RNAseq)
    # T_res is assigned to an empty matrix with in columns : the number of cell types @n parameter and in rows : gene expressions from tumor_RNAseq
    T_res = matrix(NA, ncol = n, nrow = nrow(tumor_RNAseq))
    # Fill the already known values with those from @tumor_RNAseq
    T_res[, 1:ncol(tumor_RNAseq)] = tumor_RNAseq
    rownames(T_res) = rownames(tumor_RNAseq)
    t = 0
    # Progress bar
    pb <- progress_bar$new(format = "  Running RiTMIC::simu_T_cancer [:bar] :current/:total (:percent) in :elapsed",total = nrow(tumor_RNAseq), clear = FALSE, width= 80)
    # For loop to create new simulations : gene per gene 
    for(g in rownames(tumor_RNAseq)){
      pb$tick()
      # Extraction of the gene expression
      dist_g = tumor_RNAseq[g, ]
      #if the gene is not expressed in any sample : we don't try to compute the normalmixtest
      if(sum(dist_g) == 0){
        x = dist_g
      } else {
        dists_sep = NA
        #Computing the bimodal distribution
        res <- try(capture.output(dists_sep <- normalmixEM(dist_g, maxrestarts=1000, k= 2)), silent = T)
        # If it doesn't work (too much 0s, extreme value, etc.), we sample in the gene distribution
        if((inherits(res, "try-error")) || (is.na(dists_sep))){
          x = dist_g
          t = t+1    
          # Else we sample in the computed distribution
        } else { 
          # Probability of gene expression in cell types 
          probs <- dists_sep$lambda
          # Final mean parameters
          m <- dists_sep$mu
          # Final standard deviations
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
      # Import the new generated values in the matrix T_res, following the previous dataset   
      T_res[g, (ncol(tumor_RNAseq)+1):n] = sample(x, n_simu, replace=TRUE)
    }
    cat("Number of genes failing the test normalmixEM:",t)
    return(T_res)
  }
}

#' Simple deregulations in tumors micro-environment 
#' @description 
#' Generation of deregulations with a simple model in the matrix T generated from \code{simu_T_cancer} with cell lines proportions from the return of \code{simu_A} method \cr\cr
#' Step 1: G genes are drawn\cr \cr 
#' Step 2: Over-expression induced to these genes : \cr \cr
#' newGeneExpression = geneExpression + x * cellTypeProportion
#' @param T_cancer The matrix T, for more documentation : \code{simu_T_cancer}
#' @param G Gene number to be sampled, this function requires an \strong{even} Gene number 
#' @param A Matrix A : distribution of cell lines per tumor \code{simu_A} 
#' @param x Coefficient, by default it's settled to 100
#' 
#' @return List with : \cr Deregulated matrix $T and deregulated genes in $g_immune and $g_fibro
#' @export
corr_prop_s = function(T_cancer, G, A, x = 100){
  if (G%%2 != 0) {
    stop("G the gene number needs to be an even number")
  } else if (ncol(T_cancer) != ncol(A)){
    stop("Patient number is not equal between the matrices A and T_cancer")    
  } else {
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
#' @param G Gene number to be sampled, this function requires an \strong{even} Gene number 
#' @param A The matrix A, for more documentation \code{simu_A}
#' @param y Overexpression coefficient to multiply the gene expression, {default} = 2
#' @param thres_i Threshold to induce overexpression in immune cells, {default} = 0.1
#' @param thres_f Threshold to induce overexpression in fibro cells, {default} = 0.45
#'
#' @return List with : Deregulated matrix T and deregulated genes in g_immune and g_fibro
#' @export
corr_prop_c = function(T_cancer, G, A, y = 2, thres_i = 0.1, thres_f = 0.45){
  if (G%%2 != 0) {
    stop("G the gene number needs to be an even number")
  } else if (ncol(T_cancer) != ncol(A)){
    stop("Patient number is not equal between the matrices A and T_cancer")
  } else {
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
}

#' New deregulations in tumors micro-environment 
#' @description 
#' Generation of deregulations from a new model in the matrix T generated from \code{simu_T_cancer} with cell lines proportions from the matrix A \code{simu_A} \cr \cr
#' Step 1: G genes are drawn\cr \cr 
#' Step 2: For the cell lines, over-expressions to every G genes is induced by: \cr \cr
#' newGeneExpression = geneExpression + [(1 + z)* cellTypeProportion)]
#'
#' @param T_cancer The matrix T, for more documentation \code{simu_T_cancer}
#' @param G Gene number to be sampled, this function requires an \strong{even} Gene number 
#' @param A The matrix A, for more documentation \code{simu_A}
#' @param z Overexpression coefficient to multiply the gene expression, by default z = 2
#'
#' @return A list with : \cr Deregulated matrix $T and deregulated genes in $g_immune and $g_fibro
#' @export
corr_prop_n = function(T_cancer, G, A, z = 2){
  if (G%%2 != 0) {
    stop("G the gene number needs to be an even number")
  } else if (ncol(T_cancer) != ncol(A)){
    stop("Patient number is not equal between the matrices A and T_cancer")    
  } else {
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
}

#' Add noise to the analysis or simulations of the tumor micro-environments gene expressions
#' @description Applications of \code{add_noise} \cr \cr Simulations : Better reflection of the biological observations \cr
#' Biological observation : improve the generalization error to avoid over-learning  
#' @param matrix_D matrix_D with or without deregulations obtained from the combinations of matrix_T * matrix_A
#' @param mean Mean of the data, {default} = 0 
#' @param sd Standard deviation of the data, {default} = 0.1 
#' @param val_min Minimum value, {default} = 0
#' @param val_max Maximum value, {default} = 1
#' 
#' @seealso \code{\link[RiTMIC]{simu_T_cancer}}, \code{\link[RiTMIC]{corr_prop_s}}, \code{\link[RiTMIC]{corr_prop_n}}, \code{\link[RiTMIC]{corr_prop_c}}
#' @return New matrix D with dimensions : \cr gene_name per tumor_number
#' @export
add_noise = function(matrix_D, mean = 0, sd = 0.1, val_min = 0, val_max = 1){
  noise = matrix(rnorm(prod(dim(matrix_D)), mean = mean, sd = sd), nrow = nrow(matrix_D))
  datam = matrix_D + noise
  datam[datam < val_min] = matrix_D[datam < val_min]
  datam[datam > val_max] = matrix_D[datam > val_max]
  return(datam)
}

#' Simulation of an RNAseq matrix, obtained with the matrices combinations of the matrix A and T for each sample
#'
#' @param matrix_A matrix obtained with \code{simu_A} 
#' @param matrix_T matrix obtained with \code{simu_T_cancer}
#' @param corr_matrix_T matrix obtained with corr_prop_functions
#' @param noise Boolean statement, {default} = T \cr\cr Do you want to use \code{add_noise} to your resulting matrix to increase the variance ? 
#' @param mean Mean of the noise, {default} = 0
#' @param sd Standard deviation of the noise, {default} = 0.1
#' @param val_min Minimum value of the noise, {default} = 0
#' @param val_max Maximum value of the noise, {default} = 1
#'
#' @seealso \code{\link[RiTMIC]{add_noise}}
#'
#' @return matrix D of gene expression per sample
#' @export
simu_D <- function(matrix_A, matrix_T, corr_matrix_T, noise = T, mean = 0, sd = 0.1, val_min = 0, val_max = 1){
  
  matrix_D = c()
  for(p in 1:dim(matrix_A)[2]){
    tumor = corr_matrix_T$T[, p]
    D_p = cbind(matrix_T, tumor) %*% matrix_A[, p]
    matrix_D = cbind(matrix_D, D_p)
  }
  colnames(matrix_D) <- paste(1:dim(matrix_D)[2])
  if (isTRUE(noise)){
    matrix_D <- add_noise(matrix_D, mean = mean, sd = sd, val_min = val_min, val_max = val_max)
  }
  return(matrix_D)
}

