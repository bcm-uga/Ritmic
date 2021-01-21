#' Simulate cell line distributions in tumors
#' 
#' @description 
#' Creates a new matrix with @n : the patient number and @alpha : the vector of cell lines proportion of size k
#'
#' @param n patients number
#' @param alpha cell type proportions of size k in a concatenated vector
#'
#' @return Proportion of different cell lines in a matrix 
#' @export 
#'
#' @examples simu_A(20,c(1.5, 4.5, 3))
simu_A = function(n, alpha = c(1.5, 4.5, 1, 3)){
  tmp_mix_Lvar = t(gtools::rdirichlet(n = n, alpha = alpha))
}

#' Simulate a tumoral environment 
#' @Description
#' Generate a matrix with random gene expression distribution in n tumoral cell lines with genes from RNA_seq dataset 
#'
#' @param tumor_RNAseq RNAseq dataset to simulate a random T matrix
#' @param n Supposed cell types number 
#'
#' @return A matrix : Gene expression distribution per cell type  
#' @export
#'
#' @examples
simu_T_cancer = function(tumor_RNAseq, n){
  # n_simu return the difference between the @n integer and the col number of the @tumor_RNAseq
  n_simu = n - ncol(tumor_RNAseq)
  # T_res is assigned to an empty matrix with in columns : the number of cell types @n parameter and in rows : gene expressions from tumor_RNAseq
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
      # Probability of gene expression in cell types 
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
      # Checking the wrong genes 
      x = dist_g
      t = t+1
    }
    # The values from the matrix T_res are replaced by x with n_simu reordered  
    T_res[g, (ncol(tumor_RNAseq)+1):n] = sample(x, n_simu, replace=TRUE)
  }
  return(T_res)
}


