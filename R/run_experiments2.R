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
  


  
  
  methods <- c("FIT_SAR_BIC", "FIT_SAR_CV", 
               "Witten_Perm", "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
               "SCCA_Parkhomenko", "Canonical Ridge-Author")
  
  for(i in 1:length(methods)){
    
    
    tryCatch({
      temp1 = Sys.time()
      
      result <- additional_checks(data$X, data$Y, method.type = methods[i], rank = r)
      temp2 = Sys.time()
      
      result <- evaluate_results(Uhat= result$u, 
                                 Vhat = result$v, 
                                 example = data, 
                                 name_method= methods[i], 
                                 overlapping_amount=0,
                                 lambdax= NA,
                                 lambday = NA, 
                                 normalize_diagonal=T,
                                 criterion="prediction",
                                 r_pca = r_pca, nnz= nnzeros,
                                 signal_strength= strength_theta )
      output = rbind(output,  result)
      
      t1 = c(t1, temp1)
      t2 = c(t2, temp2)
      
    }, error = function(e) {
      # Print the error message
      cat("Error occurred in alternative methods", ":", conditionMessage(e), "\n")
      # Skip to the next iteration
    })
    
    
    
    
  }
  
  
  tryCatch({
    result <- GCA_checks(data, r = r)
    t1 = c(t1, result$t1)
    t2 = c(t2, result$t2)
    
    result1 <- evaluate_results(Uhat= result$fantope$u, 
                                Vhat = result$fantope$v, 
                                example = data, 
                                name_method= "Fantope", 
                                overlapping_amount=0,
                                lambdax= NA,
                                lambday = NA, 
                                normalize_diagonal=T,
                                criterion="prediction",
                                r_pca = r_pca, nnz= nnzeros,
                                signal_strength= strength_theta )
    
    result2 <- evaluate_results(Uhat= result$gca$u, 
                                Vhat = result$gca$v, 
                                example = data, 
                                name_method= "GCA", 
                                overlapping_amount=0,
                                lambdax= NA,
                                lambday = NA, 
                                normalize_diagonal=T,
                                criterion="prediction",
                                r_pca = r_pca, nnz= nnzeros,
                                signal_strength= strength_theta )
    
    
    output = rbind(output,  result1)
    output = rbind(output, result2)
    
    
    
  }, error = function(e) {
    # Print the error message
    cat("Error occurred in alternative methods", ":", conditionMessage(e), "\n")
    # Skip to the next iteration
  })
  
  
  tryCatch({
    temp1 = Sys.time()
    temp2 = Sys.time()
    
    zero_u = matrix(0, p, r)
    zero_v = matrix(0, q, r)
    result <- evaluate_results(Uhat= zero_u,
                               Vhat = zero_v, 
                               example = data, 
                               name_method= "zero_benchmark", 
                               overlapping_amount=0,
                               lambdax= NA,
                               lambday = NA, 
                               normalize_diagonal=T,
                               criterion="prediction",
                               r_pca = r_pca, nnz= nnzeros,
                               signal_strength= strength_theta )
    output = rbind(output,  result)
    
    t1 = c(t1, temp1)
    t2 = c(t2, temp2)
    
  }, error = function(e) {
    # Print the error message
    cat("Error occurred in alternative methods", ":", conditionMessage(e), "\n")
    # Skip to the next iteration
  })

    tryCatch({
    temp1 = Sys.time()
    
    result_admm1 = scca_chao(data$X, data$Y, r = r)
    
    temp2 = Sys.time()
    
    lasso_admm_dist1<- evaluate_results(Uhat= result_admm1$u, 
                                        Vhat = result_admm1$v, 
                                        example = data, 
                                        name_method="scca_chao", 
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
  
  
  
}



output$time = t2- t1

write_excel_csv(output, paste0("EfficientCCA/results2/", name_exp, "_", ".csv"))

