##  loess smoother, The $R^2$ statistic is calculated for this model against the intercept only null mode
##Local Polynomial Regression Fitting,Fit a polynomial surface determined by one or more numerical predictors, using local fitting.


VScan=function (x, y, U=NULL,methods = "loess",span = 0.65, family="gaussian",...) {
  
  if (methods=="lm") 
    nonpara="FALSE"
  else nonpara="TRUE"
  
  rocPerCol <- function(dat, cls){
    roc_auc <- ModelMetrics::auc(cls, dat)
    max(roc_auc, 1 - roc_auc)
  }
  
  if (is.null(U))
    x=x
  else
    x=as.data.frame(cbind(x,U))
  
  asNumeric <- function(data){
    fc <- sapply(data, is.factor)
    utils::modifyList(data, lapply(data[, fc], as.numeric))
  }
  
  ### now start the feature selection
  
  notNumber <- sapply(x, function(x) !is.numeric(x))
  x = asNumeric(x)
  
  if(is.factor(y)){
    classLevels <- levels(y)
    k <- length(classLevels)
    
    if(k > 2){
      
      Combs <- utils::combn(classLevels, 2)
      CombsN <- combn(1:k, 2)
      
      lStat <- lapply(1:ncol(Combs), FUN = function(cc){
        yLevs <- as.character(y) %in% Combs[,cc]
        tmpX <- x[yLevs,]
        tmpY <- as.numeric(y[yLevs] == Combs[,cc][2])
        apply(tmpX, 2, rocPerCol, cls = tmpY)
      })
      Stat = do.call("cbind", lStat)
      
      loutStat <- lapply(1:k, function(j){
        apply(Stat[,CombsN[,j]], 1, max)
      })
      
      outStat = do.call("cbind", loutStat)
      
    } else {
      tmp <- apply(x, 2, rocPerCol, cls = y)
      outStat <- cbind(tmp, tmp)
    }
    
    outStat <- as.data.frame(outStat, stringsAsFactors = FALSE)
    colnames(outStat) <- classLevels
    rownames(outStat) <- dimnames(x)[[2]]
    outStat <- data.frame(outStat)
  } else {
    
    paraFoo <- function(data, y) abs(stats::coef(summary(lm(y ~ data, na.action = na.omit)))[2, "t value"])
    nonparaFoo <- function(x, y, ...)
    {
      meanMod <- sum((y - mean(y, rm.na = TRUE))^2)
      nzv <- caret::nearZeroVar(x, saveMetrics = TRUE)
      
      if(nzv$zeroVar) return(NA)
      if(nzv$percentUnique < 20)
      {
        regMod <- lm(y~x, na.action = na.omit, ...)
      } else {
        regMod <- try(loess(y~x, na.action = na.omit, ...), silent = TRUE)
        
        if(inherits(regMod, "try-error") | any(is.nan(regMod$residuals))) try(regMod <- lm(y~x, ...))
        if(inherits(regMod, "try-error")) return(NA)
      }
      
      pR2 <- 1 - (sum(resid(regMod)^2)/meanMod)
      if(pR2 < 0) pR2 <- 0
      pR2
    }
    
    testFunc <- if(nonpara) nonparaFoo else paraFoo
    
    outStat <- apply(x, 2, testFunc, y = y)
    outStat <- data.frame(Weights = outStat)
  }
  
  
  
  pvalue_norm=function(weight){
    
    normalize <- function(x) {
      return ((x - min(x)) / (max(x) - min(x)))
    }
    normdat <- apply(weight, 2, normalize)
    asindat=apply(normdat,2, function(x) {asin(sqrt(x))})
    normpval=pnorm(asindat, mean = mean(asindat,na.rm=TRUE), sd = sd(asindat,na.rm = TRUE), lower.tail = FALSE,log.p = FALSE)
    return(normpval)
  }
  
  Wpvalue_norm=pvalue_norm(outStat)
  colnames(Wpvalue_norm)="p.value"
  
  if (is.null(U))
    Wpvalue_norm=Wpvalue_norm
  else {Xrep_norm=!rownames(Wpvalue_norm) %in% colnames(U) 
  Wpvalue_norm=data.frame(Wpvalue_norm[Xrep_norm,,drop=FALSE])}
  
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
    
    return(list(p=data.frame(p.values=reschi2test, padj=padj),lambda=lambda))
  }
  Wpvalue_chi=pvalue_chi(outStat)$p
  lambda=pvalue_chi(outStat)$lambda
  if (is.null(U))
    Wpvalue_chi=Wpvalue_chi
  else {
    Xrep=!rownames(Wpvalue_chi) %in% colnames(U) 
    Wpvalue_chi=Wpvalue_chi[Xrep,,drop=FALSE]
    
    #  Wpvalue_chi$p.values=Wpvalue_chi$p.values[Xrep,drop=FALSE]
    #  Wpvalue_chi$padj=Wpvalue_chi$padj[Xrep,drop=FALSE]
    
  }
  
  return(list(W=outStat,p_norm=Wpvalue_norm,pvalue_chi=Wpvalue_chi,lambda=lambda))
  
}
