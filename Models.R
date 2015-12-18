# Generic model object that returns the gradient vector and Hessian matrix
# Stolen frmo here:
# http://oddhypothesis.blogspot.com.au/2014/08/optimizing-with-r-expressions.html

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
                            }
                          )
)

mo = ModelObject(
  name = 'gompertz', 
  expr = expression( y0*exp(A*exp(-exp((u*exp(1)/A)*(l-x)+1))) )
)
