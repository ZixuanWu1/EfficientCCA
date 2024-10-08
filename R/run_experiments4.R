library(MASS)
library(stats)
library(CVXR)
library(geigen)
library(pracma)
library(tidyverse)
library(CCA)
library(VGAM)
library(matlib)
library(PMA)
library(mvtnorm)
library(glmnet)
library(caret)

#Example of TGD on Sparse CCA
#n = 500, p1 = p2 = 100, s_u = s_v = 5
#k = 20, eta = 0.0025, lambda =0.01, T = 12000

source("EfficientCCA/R/GCA/utils.R")
source("EfficientCCA/R/GCA/gca_to_cca.R")
source("EfficientCCA/R/GCA/init_process.R")
source("EfficientCCA/R/GCA/sgca_init.R")
source("EfficientCCA/R/GCA/sgca_tgd.R")
source("EfficientCCA/R/GCA/subdistance.R")
source("EfficientCCA/R/GCA/adaptive_lasso.R")

source('EfficientCCA/R/alternative_methods/SAR.R')
source('EfficientCCA/R/alternative_methods/Parkhomenko.R')
source('EfficientCCA/R/alternative_methods/Witten_CrossValidation.R')
source('EfficientCCA/R/alternative_methods/Waaijenborg.R')

source("EfficientCCA/R/generate_examples.R")
source("EfficientCCA/R/evaluate_results.R")
source("EfficientCCA/R/EfficientCCA.R")
source("EfficientCCA/R/All_methods.R")
source("EfficientCCA/R/metrics.R")


args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
name_exp <- args[2]
set.seed(seed)
n <- as.numeric(args[3])
p <- as.numeric(args[4])
nnzeros <- as.numeric(args[5])
r <- as.numeric(args[6])

q = p

r_pca = 5


library(gradfps)

t1 = NULL
t2 = NULL

output = NULL
groups <- vector("list", length = (p / 5) * (q / 5))
m <- 1
indices_i <- rep(1:(p / 5), each = (q / 5))
indices_j <- rep(1:(q / 5), times = (p / 5))

# Loop through indices to fill the groups list
for (index in seq_along(indices_i)) {
  i <- indices_i[index]
  j <- indices_j[index]
  
  # Generate combinations of indices
  new_group <- as.matrix(expand.grid(5 * (i - 1) + 1:5, 5 * (j - 1) + 1:5))
  
  # Store result in the list
  groups[[m]] <- new_group
  
  m <- m + 1
}


for(strength_theta in c("strong", "medium", "weak")){
  if (strength_theta == "strong"){
    theta = diag(rep(.9, r))
  }
  if (strength_theta == "medium"){
    theta = diag(rep(.7, r))
  }
  if (strength_theta == "weak"){
    theta = diag(rep(.5, r))
  }
  
  
  
  data = generate_example_none_trivial_pca(n, p, q, r_pca = r_pca, nnzeros = nnzeros, 
                                           theta = theta, overlapping_amount = 0,
                                           lambda_pca = 1, r = r)
  
  
  
  tryCatch({
    temp1 = Sys.time()
    
    result_admm1 = lasso_cca(data$X, data$Y, lambda = sqrt(log(p + q)/n), r = r)
    
    temp2 = Sys.time()
    
    lasso_admm_dist1<- evaluate_results(Uhat= result_admm1$U, 
                                        Vhat = result_admm1$V, 
                                        example = data, 
                                        name_method="lasso_theory", 
                                        overlapping_amount=0,
                                        lambdax= NA,
                                        lambday = NA, 
                                        normalize_diagonal=T,
                                        criterion="prediction",
                                        r_pca = r_pca, nnz= nnzeros,
                                        signal_strength= strength_theta)
    
    t1 = c(t1, temp1 )
    
    t2 = c(t2, temp2)
    
    output = rbind(output,  lasso_admm_dist1)
    
    
  }, error = function(e) {
    # Print the error message
    cat("Error occurred in alternative methods", ":", conditionMessage(e), "\n")
    # Skip to the next iteration
  })
  
  
  tryCatch({
    temp1 = Sys.time()
    
    result_admm2 = lasso_cca(data$X, data$Y, lambda_max = .2, r = r, num_lambda = 20)
    
    temp2 = Sys.time()
    
    lasso_admm_dist2<- evaluate_results(Uhat= result_admm2$U, 
                                        Vhat = result_admm2$V, 
                                        example = data, 
                                        name_method="lasso_cv", 
                                        overlapping_amount=0,
                                        lambdax= NA,
                                        lambday = NA, 
                                        normalize_diagonal=T,
                                        criterion="prediction",
                                        r_pca = r_pca, nnz= nnzeros,
                                        signal_strength= strength_theta)
    
    t1 = c(t1, temp1 )
    
    t2 = c(t2, temp2)
    
    output = rbind(output,  lasso_admm_dist2)
    
    
  }, error = function(e) {
    # Print the error message
    cat("Error occurred in alternative methods", ":", conditionMessage(e), "\n")
    # Skip to the next iteration
  })
  
  
  tryCatch({
    temp1 = Sys.time()
    
    result_admm3 = lasso_cca(data$X, data$Y, lambda = 0.5 * sqrt(log(p + q)/n), groups = groups, r = r)
    
    temp2 = Sys.time()
    
    lasso_admm_dist3<- evaluate_results(Uhat= result_admm3$U, 
                                        Vhat = result_admm3$V, 
                                        example = data, 
                                        name_method="lasso_group_theory", 
                                        overlapping_amount=0,
                                        lambdax= NA,
                                        lambday = NA, 
                                        normalize_diagonal=T,
                                        criterion="prediction",
                                        r_pca = r_pca, nnz= nnzeros,
                                        signal_strength= strength_theta)
    
    t1 = c(t1, temp1 )
    
    t2 = c(t2, temp2)
    
    output = rbind(output,  lasso_admm_dist3)
    
    
  }, error = function(e) {
    # Print the error message
    cat("Error occurred in alternative methods", ":", conditionMessage(e), "\n")
    # Skip to the next iteration
  })
  
  
  
  
  tryCatch({
    temp1 = Sys.time()
    
    result_admm4 = lasso_cca(data$X, data$Y, groups = groups, r = r, lambda_max = .15, num_lambda= 30)
    temp2 = Sys.time()
    
    
    lasso_admm_dist4<- evaluate_results(Uhat= result_admm4$U, 
                                        Vhat = result_admm4$V, 
                                        example = data, 
                                        name_method="lasso_group_cv", 
                                        overlapping_amount=0,
                                        lambdax= NA,
                                        lambday = NA, 
                                        normalize_diagonal=T,
                                        criterion="prediction",
                                        r_pca = r_pca, nnz= nnzeros,
                                        signal_strength= strength_theta)
    output = rbind(output,  lasso_admm_dist4)
    
    
    
    t1 = c(t1, temp1 )
    
    t2 = c(t2, temp2)
    
  }, error = function(e) {
    # Print the error message
    cat("Error occurred in alternative methods", ":", conditionMessage(e), "\n")
    # Skip to the next iteration
  })

  

  

}



output$time = t2- t1

write_excel_csv(output, paste0("EfficientCCA/results2/", name_exp, "_", ".csv"))