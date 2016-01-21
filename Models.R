# Generic model object that returns the gradient vector and Hessian matrix
# Adapted from here:
# http://oddhypothesis.blogspot.com.au/2014/08/optimizing-with-r-expressions.html

# Generic model object 
ModelObject = setRefClass('ModelObject', 
                          fields = list(
                            name = 'character',
                            expr = 'expression'
                          ),
                          methods = list(
                            value = function(p, data){
                              eval(.self$expr, c(as.list(p), as.list(data)))
                            },
                            jacobian = function(p, data){
                              J = t(sapply(all.vars(.self$expr), function(v, p, data){
                                eval(D(.self$expr, v), c(as.list(p), as.list(data)))
                              }, p=p, data=data))
                              
                              return(J[names(p),,drop=F])
                            },
                            gradient = function(p, data){
                              r = data$y - value(p, data)
                              return(-jacobian(p, data) %*% r)
                            },
                            hessian = function(p, data){
                              J = jacobian(p, data)
                              return(J %*% t(J))
                            },
                            dev =function(p, data){
                              eval(.self$expr, c(as.list(p), as.list(data)))-data$y
                            }
                          )
)

# Function for getting best fit parameters and other relevant info
fit = function(mo,p,data){
  # Solve for best fit paramters
  Opt<-nlminb(p, function(p, data){
    r = data$y - mo$value(p,data)
    sum(dnorm(data$y,mean=mo$value(p,data),sd=SD,log=T))
    return(r %*% r)
  }, gradient = mo$gradient, hessian = mo$hessian, data=data)
  # Residuals
  Opt$residuals<-mo$dev(Opt$par,data)
  # Residual degrees of freedom
  Opt$df.residual<-dim(data)[1]-length(p)
  # Sigma
  Opt$sigma<-sqrt(sum(Opt$residuals^2)/Opt$df.residual)
  # Hessian matrix
  Opt$hessian<-mo$hessian(Opt$par,data)
  # Gradient
  Opt$gradient<-mo$gradient(Opt$par,data)
  # Variance covariance matrix
  Opt$vcov<-solve(0.5*(Opt$hessian + t(Opt$hessian)))*Opt$sigma^2
  # Paramter standard errors
  Opt$se<-sqrt(diag(Opt$vcov))
  # Organise parameters and standard errors
  Opt$parlist<-list()
    for(i in 1:length(Opt$par)){
      Opt$parlist[[i]]<-c(Opt$par[i],Opt$se[i])
      names(Opt$parlist)[i]<-paste(names(Opt$par[i]))
    }
  Opt$parlist<-data.frame(Opt$parlist)
  # 
  Opt$confInts<-predictNLMINB(Opt,data)
  Opt$data<-data
  return(Opt)
}


VB3 = ModelObject(
  name = 'VB3', 
  expr = expression( b2+(b1-b2)*(1-exp(-b3*x)) )
)

# Function for propagating uncertainty and getting 95% confidence intervals
#confInts<-Vectorize(function(x){
#  List<-Opt$parlist
#  List$x=c(x,0)
#  conf<-t(propagate(expr=mo$expr,data=List,type="stat",do.sim=FALSE)$prop[c(1,3,5,6)])
#  predLower<-conf[1]+sqrt(conf[2]^2+Opt$sigma^2)*qt(0.025,Opt$df.residual)
#  predUpper<-conf[1]+sqrt(conf[2]^2+Opt$sigma^2)*qt(0.975,Opt$df.residual)
#  conf<-c(conf,predLower,predUpper)
#  return(conf)
#})

# A hacked version of the predictNLS function to work with NLMINB output
predictNLMINB<-function(Opt,data){
  newdata<-data.frame(x=seq(min(data$x),max(data$x),length.out=50))
  ALPHA=0.05
  EXPR<-mo$expr
  VARS <- all.vars(EXPR)
  COEF <- Opt$par
  predVAR <- setdiff(VARS, names(COEF))
  if (!identical(names(newdata)[1:length(predVAR)], predVAR)) 
    stop("newdata should have name(s) ", predVAR, "!\n")
  VCOV <- Opt$vcov
  NR <- NROW(newdata)
  outMAT <- matrix(nrow = NR, ncol = 6)
  propLIST <- vector("list", length = NR)
  for (i in 1:NR) {
    tempDATA <- newdata[i, ]
    names(tempDATA) <- colnames(newdata)
    SEL <- which(names(tempDATA) == predVAR)
    predVEC <- c(COEF, as.numeric(tempDATA[SEL]))
    names(predVEC) <- c(names(COEF), predVAR)
    DF <- rbind(predVEC, 0)
    row.names(DF) <- NULL
    if (NCOL(tempDATA) == length(predVAR)){
      forCOV <- rep(0, length(predVAR))}else
        {forCOV <- tempDATA[, (length(predVAR) + 1):(2*length(predVAR))]}
    COV <- VCOV
    for (k in 1:length(forCOV)) {
      COV <- mixCov(COV, forCOV[k])
    }
    SEL <- tail(1:nrow(COV), length(predVAR))
    dimnames(COV)[[1]][SEL] <- dimnames(COV)[[2]][SEL] <- predVAR
    PROP <- propagate(expr = EXPR, data = DF, use.cov = COV, 
                      alpha = ALPHA, do.sim = FALSE)
    propLIST[[i]] <- PROP
    outPROP <- PROP$prop
    MEAN<-outPROP[2]
    SD<-outPROP[4]
    TQUAN <- qt(1 - ALPHA/2, Opt$df.residual)
        outPROP[5] <- MEAN - TQUAN * SD
        outPROP[6] <- MEAN + TQUAN * SD
        outPROP[7] <- MEAN - TQUAN * sqrt(SD^2 + Opt$sigma^2)
        outPROP[8] <- MEAN + TQUAN * sqrt(SD^2 + Opt$sigma^2)
    outMAT[i, ] <- outPROP[-c(1,3)]
  }
  outMAT <- as.data.frame(cbind(newdata,outMAT))
  colnames(outMAT) <- c("x","Mean","SD","cLower","cUpper","pLower","pUpper") 
  return(outMAT)
}