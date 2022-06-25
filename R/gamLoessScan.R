##Generalized Additive Model using LOESS



gamLoessScan=function(genotype,traits,U,cv_method="adaptive_cv",model_metric="RMSE",n_hyperparameter_search=10,...){

  
  set.seed(123)
  econtrol <- caret::trainControl(## 5-fold CV, repeat 5 times
    method = cv_method,
    number = 5,
    ## repeated 5 times
    repeats = 5,
    adaptive = list(min = 5, alpha = 0.05,method = "gls", complete = TRUE),
    search = "random")
  
  
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }
  
  
  #  genotype=as.data.frame(apply(genotype, 2, normalize))
  train_method="gamLoess"
  genotype_u=cbind(genotype,U)
  
  print(paste0("Target tait--",colnames(traits)," training"))
  Sys.time()
  print(paste0("Tuning hyperparameters --",colnames(traits),". _/Searching_",n_hyperparameter_search, "_parameter combinations from hyperparameter space"))
  Sys.time()
  #env=normalize(simdata[[j]][,2:11])
  # genotype_norm=as.data.frame(apply(sim_example,2,normalize))
  # simf=as.formula(paste(para[i],paste(names(genotype), collapse="+"),sep="~"))
  model_qt=caret::train(y=traits, x=as.data.frame(genotype_u),
                        method="gamLoess",
                        metric = model_metric,## "Accuracy", "RMSE"
                        #preProcess=c("scale"),
                        tuneLength = n_hyperparameter_search, ### search 10 combinations of parameters
                        # verbose=1 is reporting the progress,o is sclience
                        trControl = econtrol)
  print(paste0("Tuning hyperparameter finished"))
  Sys.time()
  #save(model_qt,file=paste0("Target_trait_",colnames(traits),"_",train_method,"_trained_model_COV.RData"))
  ###
  Model_imp=caret::varImp(model_qt)
  
  print("Exporting variant weights")
  Sys.time()
  write.csv(Model_imp$importance,file = paste0("Target_trait_",colnames(traits),"_",train_method,"_model_imp_COV.csv"))
  print("Exporting tuning results")
  write.csv(model_qt$results,file = paste0("Target_trait_",colnames(traits),"_",train_method,"_trained_model_tuning_results_COV.csv"))
  print(paste0("Target_trait_",colnames(traits),"_",train_method,"_training finished"))
  Sys.time()
  
  #################
  print("############################### Estimating p values ################")
  Sys.time()
  
  pvalue_norm=function(weight){
    
    normalize <- function(x) {
      return ((x - min(x)) / (max(x) - min(x)))
    }
    normdat <- apply(weight, 2, normalize)
    asindat=apply(normdat,2, function(x) {asin(sqrt(x))})
    normpval=pnorm(asindat, mean = mean(asindat,na.rm=TRUE), sd = sd(asindat,na.rm = TRUE), lower.tail = FALSE,log.p = FALSE)
    return(normpval)
  }
  
  Wpvalue_norm=pvalue_norm(Model_imp$importance)
  colnames(Wpvalue_norm)="p.value"
  
  pvalue_chi<-function(weight,K=1)
  {
    
    normalize <- function(x) {
      return ((x - min(x)) / (max(x) - min(x)))
    }
    normdat <- apply(weight, 2, normalize)
    asindat=apply(normdat,2, function(x) {asin(sqrt(x))})
    # resmaha <- covRob(asindat, distance = TRUE, na.action= na.omit, estim="donostah")$dist
    lambda <- median(asindat)/qchisq(0.5,df=K)
    reschi2test <- pchisq(asindat/lambda,K,lower.tail=FALSE)
    colnames(reschi2test)="p.value"
    
    padj <- p.adjust(reschi2test,method="bonferroni")
    return(data.frame(p.values=reschi2test, padj=padj,lambda=lambda))
  }
  
  Wpvalue_chi=pvalue_chi(Model_imp$importance)
  return(list(W=Model_imp$importance,p_norm=Wpvalue_norm,pvalue_chi=Wpvalue_chi))
  #registerDoSEQ()
  #stopCluster(cl)
  
}




