#draws on http://www.jameskeirstead.ca/blog/gaussian-process-regression-with-r/

cov_fun_rethink <- function(d, etasq, rhosq){
  etasq * exp( -rhosq * d^2)
}

#For the GPL2 formulation from rethinking
#assumed mean of 0
predict_gpl2 <- function(xold, xnew, yold,
                       etasq = 1, rhosq = 1, 
                       cfun = cov_fun_rethink, 
                       out="sim"){

    dmat_old <- as.matrix(dist(xold))
  
  if(class(xold) == "matrix"){
    dmat_new_all <- as.matrix(dist(c(xold, xnew)))
  }else{
    dmat_new_all <- as.matrix(dist(c(xold, xnew)))
  }
  
  dmat_old_new <- dmat_new_all[-c(1:length(xold)), -c((length(xold)+1):ncol(dmat_new_all))]
  dmat_new_old <- dmat_new_all[-c((length(xold)+1):nrow(dmat_new_all)), -c(1:length(xold))]
  dmat_new_new <- dmat_new_all[-c(1:length(xold)), -c(1:length(xold))]
  
  k_old <- cfun(dmat_old, etasq = etasq, rhosq = rhosq)
  k_old_new <- cfun(dmat_old_new, etasq = etasq, rhosq = rhosq)
  
  #for cov
  k_new_old <- cfun(dmat_new_old, etasq = etasq, rhosq = rhosq)
  k_new_new <- cfun(dmat_new_new, etasq = etasq, rhosq = rhosq)
  
  diag(k_old) <- etasq
  diag(k_new_new) <- etasq
  
  mu <- k_old_new %*% solve(k_old) %*%yold
  cov_pred <- k_new_new - k_old_new %*%  solve(k_old) %*% k_new_old
  
  if(out=="mu") return(mu)
  
  if(out=="sim") return(as.vector(MASS::mvrnorm(1, mu, cov_pred)))
  
  return(data.frame(fit = mu, fit_se = sqrt(diag(cov_pred))))
  
}

#wrapper
predict_gpl2_fromsamp <- function(xold, xnew, yold_mat, etasq, rhosq, 
                                  n = length(etasq)){
  t(sapply(1:n,
         function(i) predict_gpl2(xold = xold, 
                                  xnew = xnew, 
                                  yold = yold_mat[i,],
                                  etasq = etasq[i],
                                  rhosq = rhosq[i])))
}