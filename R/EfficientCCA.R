library(foreach)
library(caret)
library(doParallel)
registerDoParallel(makeCluster(8))
library(SMUT)


soft_thresh <- function(A, lambda){
  result = sign(A) * pmax(abs(A) - lambda, 0)
  return(result)
}

soft_thresh2 <- function(A, lambda){
  result = A * pmax(1 - lambda/(sqrt(sum(A^2))), 0)
  return(result)
}



##' Fit lasso problem via ADMM
##' 
##' @param X:  n by p data matrix
##' @param Y:  n by q data matrix
##' @param lambda: penalty level
##' @param rho: parameter in augmented lagrangian
##' @param r: target rank of cca
##' @param niter: max number of iterations
##' 
##' @return C: U Lambda V^T
lasso_cca_admm <- function(X, Y, lambda = 1, rho = 1, r = 2, niter = 500, B_init = NULL, groups = NA){
  
  p = dim(X)[2]
  q = dim(Y)[2]
  n = dim(X)[1]
  
  Sigma_XY = eigenMapMatMult(t(X), Y)/ n
  Sigma_XX = eigenMapMatMult(t(X), X) / n
  Sigma_YY = eigenMapMatMult(t(Y), Y) / n
  
  if(is.null(B_init)){
    B = matrix(0, p, q)
  } else{
    B = B_init
  }
  U = matrix(0, p, q)
  Z = B
  
  eigen1 = eigen(Sigma_XX, symmetric = T)
  eigen2 = eigen(Sigma_YY, symmetric = T)
  
  b =outer(eigen1$values, eigen2$values) + rho
  
  B_1 = eigenMapMatMult(eigenMapMatMult(t(eigen1$vectors), Sigma_XY),  (eigen2$vectors) )
  
  for(i in 1:niter){
    
    # Update B
    
    
    oldB <- B
    B <- B_1 +   eigenMapMatMult(eigenMapMatMult(t(eigen1$vectors), ( - rho * U + rho * Z)),  (eigen2$vectors) ) 
    B <- B / b
    B <- eigenMapMatMult(eigenMapMatMult((eigen1$vectors), B),  t(eigen2$vectors) )
    
    
    # Update Z
    Z = B + U
    
    if (is.na(groups)){
      
      Z = soft_thresh(Z, lambda/rho)
    }
    else{
      for (g in length(groups)){
        Z[[g]] =  soft_thresh2(Z[[g]], lambda/rho)
      }
      
    }
    
    # Update U
    
    U = U + rho * (B - Z)
    
    if(max(abs(oldB - B)) < 1e-4){
      break
    }
  }
  
  if (i >= niter) {
    warning("ADMM does not converge")
  }
  
  

  return(list(result = Z, eigen1 = eigen1, eigen2 = eigen2))
  
  
}


##' Choose lasso hyperparameter via cross validation
##' 
##' @param X:  n by p data matrix
##' @param Y:  n by q data matrix
##' @param rho: parameter in augmented lagrangian
##' @param lambda_max: maximum lambda for cv
##' @param num_lambda: number of lambdas for cv
##' @param r: target rank of cca
##' @param niter: max number of iterations
##' @param nfold: number of folds in CV
##' 
##' @return C: U Lambda V^T
##' 
cv_admm_lasso <- function(X, Y, rho = 1, lambda_max = .2, num_lambda = 10, r = 2, niter = 500,
                          nfold = 8, B_init = NULL, groups = NA){
  p = dim(X)[2]
  q = dim(Y)[2]
  n = dim(X)[1]
  
  ## Create folds
  cv = createFolds(1:n, k = nfold, list = T )
  
  ## choose penalty lambda
  lambda_values <- rho * seq(from =0, to = 1, by = 1/num_lambda) * lambda_max
  cv_results <- data.frame(lambda = numeric(), mse = numeric(), std = numeric())
  
  ## Cross validation
  for(lambda in lambda_values) {
    mse_values = NULL
    
    ## Parallelly run cross validation
    results <- foreach(fold = cv, .export = c("lasso_cca_admm"), .packages = "SMUT") %dopar% {
      
      soft_thresh <- function(A, lambda){
        result = sign(A) * pmax(abs(A) - lambda, 0)
        return(result)
      }
      
      test_indices <- fold
      train_indices <- setdiff(1:n, test_indices)
      
      ## Fit lasso model
      B <- lasso_cca_admm(X[train_indices, ], Y[train_indices, ],  lambda, rho, r, niter, B_init, groups)$result
      
      (c(sum((X[test_indices, ] %*% B %*% t(Y[test_indices,]) / length(test_indices) -diag(length(test_indices))  )^2)))
      
      
    }
    
    # Store average MSE for this lambda
    cv_results <- rbind(cv_results, data.frame(lambda = lambda, mse = mean(unlist(results)),
                                               std = sd(unlist(results)) ))
    # Break the cv process if lambda is so large that B = 0
    if(all(unlist(results)  == n )){
      break
    }
    
  }
  return(cv_results)
  
}

##' Fit lasso problem via ADMM
##' 
##' @param X:  n by p data matrix
##' @param Y:  n by q data matrix
##' @param lambda: penalty level. (selected by cv if by default)
##' @param rho: parameter in augmented lagrangian
##' @param r: target rank of cca
##' @param niter: max number of iterations
##' @param lambda_max: maximum lambda for cv
##' @param num_lambda: number of lambdas for cv
##' @param nfold: number of folds in CV
##' 
##' @return C: U Lambda V^T
lasso_cca <- function(X, Y, lambda = NULL, rho = 1, r= 2, niter = 500,
                      lambda_max = .2, num_lambda = 10, nfold = 8, one_sd = F, B_init = NULL, groups = NA){
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(Y)[2]
  
  # select lambda via cv if not specified
  if(is.null(lambda)){
    cv_result = cv_admm_lasso(X, Y, rho,  lambda_max, num_lambda, r , niter = niter, nfold = nfold, 
                              B_init = B_init, groups = groups)
    min_ind = which(cv_result$mse == min(cv_result$mse) )
    if(one_sd){
      lambda = max(cv_result$lambda[which( cv_result$mse <= cv_result$mse[min_ind] + cv_result$std[min_ind] )  ] )
    } else{
      lambda = cv_result$lambda[min_ind]
    }
  }
  
  print(lambda)
  # Fit lasso
  B = lasso_cca_admm(X, Y, lambda, rho, r, niter, B_init, groups)

  # Convert result to CCA estimates
  
  
  Sigma_X_root = eigenMapMatMult(eigenMapMatMult(B$eigen1$vectors, diag(sqrt(B$eigen1$values + 1e-10))),  t(B$eigen1$vectors)) 
  Sigma_Y_root = eigenMapMatMult(eigenMapMatMult(B$eigen2$vectors, diag(sqrt(B$eigen2$values + 1e-10))),  t(B$eigen2$vectors)) 
  
  B = B$result
  C = eigenMapMatMult(eigenMapMatMult(Sigma_X_root, B ),  Sigma_Y_root) 

  result = svd(C)
  
  if(r > 1){
    U_0 = result$u[, 1:r]
    V_0 = result$v[, 1:r]
    Lambda_0 = result$d[1:r]
    U =  eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(B, Sigma_Y_root ),  V_0)  , diag(1/Lambda_0))
    V =  t(eigenMapMatMult( eigenMapMatMult(eigenMapMatMult( diag(1/Lambda_0), t(U_0)  ),  Sigma_X_root)  ,B))
  }
  else{
    U_0 = matrix( result$u[, 1], ncol = 1)
    V_0 = matrix( result$v[, 1], ncol = 1)
    Lambda_0 = result$d[1]
    
    U = B %*%  Sigma_Y_root %*% V_0  /Lambda_0
    V = t( (1/Lambda_0) %*% t(U_0) %*% Sigma_X_root %*% B)
    
  }
  
  
  
  
  return(list(Pi= B, U = U, V= V))
  
}

