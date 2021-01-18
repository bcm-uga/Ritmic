#' Title simu_A
#'
#' @param n number of patients
#' @param alpha proportions of different cell types
#'
#' @return the matrix of proportions
#' @export
#'
#' @examples simu_A(n = 30, alpha = c(1.5, 4.5, 1, 3))
simu_A = function(n, alpha = c(1.5, 4.5, 1, 3)){
  tmp_mix_Lvar = t(gtools::rdirichlet(n = n, alpha = alpha))
}
