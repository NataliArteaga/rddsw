#' @title Kernel function2
#' @description estimation of the kernel function using sample weights
#' @param kernel_ kernel type: triangular, epanechnikov, tricube, uniform, gaussian y cosine
#' @param weigths pik, NULL if there aren't sample weights
#' @param X continuous variable
#' @param c cutpoint
#' @param h bandwidth
#' @return matriz con pesos kernel incluyendo o no los pesos de muestreo
#' @examples
#'

kernel_sRD <-function(kernel_,weights,h,X,c){
  if(is.null(weights)){
    weights <- rep("null",length(X))
  }
  if(kernel_=="triangular"){
    W_ <- ifelse(X>=c-h & X<=c+h,ifelse(weights=="null",(1-abs((X-c)/h)),(1-abs((X-c)/h))/(h*weights)),0)
  }else{if(kernel_=="epanechnikov"){
    W_     <- ifelse(X>=c-h & X<=c+h,ifelse(weights=="null",((3/4)*(1-((X-c)/h)^2))/h,((3/4)*(1-((X-c)/h)^2))/(h*weights)),0)
  }else{if(kernel_=="tricube"){
    W_     <- ifelse(X>=c-h & X<=c+h,ifelse(weights=="null",((1-abs((X-c)/h)^3)^3)/h,((1-abs((X-c)/h)^3)^3)/(h*weights)),0)
  }else{if(kernel_=="uniform"){
    W_     <- ifelse(X>=c-h & X<=c+h,ifelse(weights=="null",1/2,1/(2*weights)),0)
  }else{if(kernel_=="gaussian"){
    W_     <- ifelse(weights=="null",(1/sqrt(2 * pi) * exp(-1/2 * ((X-c)/h)^2)),(1/sqrt(2 * pi) * exp(-1/2 * (X-c)^2))/(h*weights))
  }else{if(kernel_=="cosine"){
    W_     <- ifelse(X>=c-h & X<=c+h,ifelse(weights=="null",(pi/4 * cos(pi/2 * ((X-c)/h)))/h,((pi/4 * cos(pi/2 * ((X-c)/h)))/(h*weights))),0)
  }}}}}}
  W <- matrix(0,length(X),length(X))
  diag(W) <- W_
  return(W)
}
