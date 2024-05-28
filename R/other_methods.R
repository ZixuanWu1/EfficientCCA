
additional_checks <- function(X_train, Y_train, S=NULL, 
                              rank=2, kfolds=5, method.type="FIT_SAR_BIC",
                              lambdax = 10^seq(from=-3,to=2,length=100),
                              lambday = c(0, 1e-7, 1e-6, 1e-5)){
  
  X_train = as.matrix(data.frame(X_train) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))))
  Y_train = as.matrix(data.frame(Y_train) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))))
  p1 <- dim(X_train)[2]
  p2 <- dim(Y_train)[2]
  p <- p1 + p2;
  n <- nrow(X_train)
  pp <- c(p1,p2);
  if(is.null(S)){
    S = cov(cbind(X_train, Y_train))
  }
  
  if (method.type=="FIT_SAR_BIC"){
    method<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                      lambdaAseq=seq(from=0.2,to=0.02,length=10),
                      lambdaBseq=seq(from=0.2,to=0.02,length=10),
                      max.iter=100,conv=10^-2,
                      selection.criterion=1,n.cv=5)
    a_estimate = rbind(method$uhat, method$vhat)
    
  }
  if(method.type=="FIT_SAR_CV"){
    method<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                      lambdaAseq=seq(from=0.2,to=0.02,length=10),
                      lambdaBseq=seq(from=0.2,to=0.02,length=10),
                      max.iter=100,conv=10^-2, selection.criterion=2, n.cv=5)
    a_estimate = rbind(method$uhat, method$vhat)
    
  }
  if (method.type=="Witten_Perm"){
    Witten_Perm <- CCA.permute(x=X_train,z=Y_train,
                               typex="standard",typez="standard", 
                               penaltyxs =matrix(seq(from=0,to=1,length=20),nrow=1),
                               penaltyzs = matrix(seq(from=0,to=1,length=20),nrow=1),
                               nperms=50)
    method<-CCA(x=X_train, z=Y_train, typex="standard",typez="standard",K=rank,
                penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=FALSE)
    a_estimate = rbind(method$u, method$v)
  }
  if(method.type=="Witten.CV"){
    Witten_CV<-Witten.CV(X=X_train,Y=Y_train, n.cv=5,
                         lambdax=matrix(seq(from=0,to=1,length=20),nrow=1),
                         lambday=matrix(seq(from=0,to=1,length=20),nrow=1))
    method <-CCA(x=X_train,z=Y_train,typex="standard",typez="standard",
                 K=rank,penaltyx=Witten_CV$lambdax.opt,
                 penaltyz=Witten_CV$lambday.opt,trace=FALSE)
    a_estimate = rbind(method$u, method$v)
    
  }
  if(method.type=="Waaijenborg-Author"){
    method<-Waaijenborg(X=X_train,Y=Y_train,
                        lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                        lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                        rank=rank,selection.criterion=1)
    a_estimate = rbind(method$vhat, method$uhat)
    
  }
  if(method.type=="Waaijenborg-CV"){
    method<-Waaijenborg(X=X_train,
                        Y=Y_train,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                        lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                        rank=rank, selection.criterion=2)
    a_estimate = rbind(method$vhat, method$uhat)
    
  }
  if(method.type=="SCCA_Parkhomenko"){
    method<- SCCA_Parkhomenko(x.data=X_train, y.data=Y_train, Krank=rank,
                              lambda.v.seq = matrix(seq(from=0,to=2,length=20),nrow=1),
                              lambda.u.seq = matrix(seq(from=0,to=2,length=20),nrow=1))
    a_estimate = rbind(method$uhat, method$vhat)
    
  }
  if(method.type=="Canonical Ridge-Author"){
    RCC_cv<-estim.regul_crossvalidation(X_train,Y_train,n.cv=5, 
                                        lambda1grid=matrix(seq(from=0,to=1,length=20),nrow=1),
                                        lambda2grid=matrix(seq(from=0,to=1,length=20),nrow=1))
    method<-rcc(X_train,Y_train, RCC_cv$lambda1.optim, RCC_cv$lambda2.optim)
    a_estimate = rbind(method$xcoef[,1:rank], method$ycoef[,1:rank])
    
    
  }
  a_estimate <- gca_to_cca(a_estimate, S, pp)
  return(a_estimate)
}

