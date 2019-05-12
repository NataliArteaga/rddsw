#' @title RDD local regression2
#' @description no parametric regression discontinuity using local regression estimator with and without sample weights
#' @param D variable with: 1 for treatment group and 0 for control group
#' @param Z continuous variable
#' @param Y outcome variable
#' @param c cutpoint
#' @param type sharp or fuzzy
#' @param TreatLeft T if the treatment is applied for below of the cutpoint
#' @param W sample weights 1/pik, NULL if there aren't sample weights
#' @param svydesign_ design sample using survey library
#' @param vcov_ variances and covariances type: HCO, HC1, HC2 y HC3
#' @param kernel_ kernel type: triangular, epanechnikov, tricube, uniform, gaussian y cosine
#' @param h bandwidth
#' @return Impact estimation, e.e. and t test
#' @examples
#'
RDker <- function(D = muestra1$D, Z = muestra1$x, Y = muestra1$y,c =0,type="sharp",TreatLeft=T,weights=NULL,svydesign_=muestra,vcov_,kernel_="uniforme",h=2){
  library(survey)
  if(!is.null(weights)){
    muestra        <- data.frame(D,Z,Y,weights)
  }else{
    muestra        <- data.frame(D,Z,Y)
  }
  muestra$intercepto  <- 1
  muestra$Z_C         <- muestra$Z-c
  muestra$ZD          <- muestra$Z*muestra$D
  if(kernel_!="gaussian"){
    muestra <- muestra[which(muestra$Z>=c-h & muestra$Z<=c+h),]
    if(!is.null(weights)){
      svydesign_ <- subset(svydesign_,Z>=c-h & Z<=c+h)
    }
  }
  Y        <- muestra$Y

  if(!is.null(weights)){
    W <- kernel_sRD(kernel_=kernel_,weights=1/muestra$weights,h=h,X=muestra$Z,c=c)
  }else{
    W <- kernel_sRD(kernel_=kernel_,weights=NULL,h=h,X=muestra$Z,c=c)
  }
  if(type!="sharp"){
    muestra$T_          <- ifelse(muestra$Z<c,1,0)
    if(TreatLeft==F){
      muestra$T_          <- 1-muestra$T_
    }
    muestra$ZT_          <- muestra$Z*muestra$T_
  }
  if(type=="sharp"){
    Xm1                 <- as.matrix(muestra[,c("intercepto","Z_C","D","ZD")])
    impacto.1_2         <- solve(t(Xm1)%*%W%*%Xm1)%*%(t(Xm1)%*%W%*%Y)
    beta                <- impacto.1_2[3,]
  }else{if(type=="fuzzy"){
    Zm1            <- as.matrix(muestra[,c("intercepto","Z_C","T_","ZT_")])
    Xm1            <- as.matrix(muestra[,c("intercepto","Z_C","D","ZT_")])
    beta1stg       <- solve(t(Zm1)%*%W%*%Zm1)%*%(t(Zm1)%*%W%*%muestra$D)
    muestra$Xest   <- Zm1%*%beta1stg
    Xest           <- as.matrix(muestra[,c("intercepto","Z_C","Xest","ZT_")])
    impacto.1_2    <- solve(t(Xest)%*%W%*%Xest)%*%(t(Xest)%*%W%*%Y)
    beta           <- impacto.1_2[3,]
  }
  }

  # Varianza sin pesos
  Yest                <- Xm1%*%impacto.1_2
  sigma2              <- sum((Y - Yest)^2)/(dim(muestra)[1]-length(impacto.1_2))
  M_                  <- (Y - Yest)^2
  M                   <- matrix(0,length(M_),length(M_))
  diag(M)             <- M_

  if(type=="sharp"){
    if(is.null(weights)){
      if(vcov_=="HC0"){
        var              <- solve(t(Xm1)%*%W%*%Xm1)%*%t(Xm1)%*%W%*%M%*%t(W)%*%Xm1%*%solve(t(Xm1)%*%W%*%Xm1)
        #t(Xm1)%*%W%*%M%*%t(W)%*%Xm1
        #j     <- cov(c((Y - Yest))*t(W)%*%Xm1)*(dim(Xm1)[1]-1)
        #var   <- solve(t(Xm1)%*%W%*%Xm1)%*%j%*%solve(t(Xm1)%*%W%*%Xm1)
      }else{if(vcov_=="HC1"){
        mat <- matrix(0,dim(Xm1)[1],dim(Xm1)[1])
        diag(mat) <- (dim(Xm1)[1]/(dim(Xm1)[1]-length(impacto.1_2)))
        var       <- solve(t(Xm1)%*%W%*%Xm1)%*%t(Xm1)%*%W%*%mat%*%M%*%t(W)%*%Xm1%*%solve(t(Xm1)%*%W%*%Xm1)
      }else{if(vcov_=="HC2"){
        mat <- matrix(0,dim(Xm1)[1],dim(Xm1)[1])
        diag(mat) <- 1/((1-diag(Xm1%*%solve(t(Xm1)%*%Xm1)%*%t(Xm1))))
        var       <- solve(t(Xm1)%*%W%*%Xm1)%*%t(Xm1)%*%W%*%mat%*%M%*%t(W)%*%Xm1%*%solve(t(Xm1)%*%W%*%Xm1)
      }else{if(vcov_=="HC3"){
        mat <- matrix(0,dim(Xm1)[1],dim(Xm1)[1])
        diag(mat) <- 1/((1-diag(Xm1%*%solve(t(Xm1)%*%Xm1)%*%t(Xm1)))^2)
        var       <- solve(t(Xm1)%*%W%*%Xm1)%*%t(Xm1)%*%W%*%mat%*%M%*%t(W)%*%Xm1%*%solve(t(Xm1)%*%W%*%Xm1)
      }}}}
    }else{
      var1            <- solve(t(Xm1)%*%W%*%Xm1) # revisar si dejo o pongo sigma2
      muestra$ek2     <- muestra$Y-Xm1%*%impacto.1_2
      diag(M)         <- muestra$ek2
      mat             <- matrix(0,dim(Xm1)[1],dim(Xm1)[1])
      if(vcov_=="HC0"){
        diag(mat) <- 1
      }else{if(vcov_=="HC1"){
        diag(mat) <- (dim(Xm1)[1]/(dim(Xm1)[1]-length(impacto.1_2)))
      }else{if(vcov_=="HC2"){
        diag(mat) <- 1/((1-diag(Xm1%*%solve(t(Xm1)%*%Xm1)%*%t(Xm1))))
      }else{if(vcov_=="HC3"){
        diag(mat) <- 1/((1-diag(Xm1%*%solve(t(Xm1)%*%Xm1)%*%t(Xm1)))^2)
      }}}}
      for(i in 1:dim(Xm1)[2]){
        muestra[,c(paste0("Uk",i))] <- (mat%*%M%*%Xm1[,i])
      }
      svydesign_$variables$Uk1<-muestra$Uk1
      svydesign_$variables$Uk2<-muestra$Uk2
      svydesign_$variables$Uk3<-muestra$Uk3
      svydesign_$variables$Uk4<-muestra$Uk4

      varmuestra           <- vcov(svytotal(~Uk1+Uk2+Uk3+Uk4,svydesign_))
      var                  <- var1%*%(varmuestra)%*%t(var1)
      
    }
  }else{if(type=="fuzzy"){
    if(is.null(weights)){
      if(vcov_=="HC0"){
        var       <- solve(t(Xest)%*%W%*%Xest)%*%t(Xest)%*%W%*%M%*%t(W)%*%Xest%*%solve(t(Xest)%*%W%*%Xest)
      }else{if(vcov_=="HC1"){
        mat <- matrix(0,dim(Xest)[1],dim(Xest)[1])
        diag(mat) <- (dim(Xest)[1]/(dim(Xest)[1]-length(impacto.1_2)))
        var       <- solve(t(Xest)%*%W%*%Xest)%*%t(Xest)%*%W%*%mat%*%M%*%t(W)%*%Xest%*%solve(t(Xest)%*%W%*%Xest)
        #var       <- (dim(Xest)[1]/(dim(Xest)[1]-2))*var
      }else{if(vcov_=="HC2"){
        mat <- matrix(0,dim(Xest)[1],dim(Xest)[1])
        diag(mat) <- 1/((1-diag(Xest%*%solve(t(Xest)%*%Xest)%*%t(Xest))))
        var       <- solve(t(Xest)%*%Xest)%*%t(Xest)%*%mat%*%M%*%Xest%*%solve(t(Xest)%*%Xest)
      }else{if(vcov_=="HC3"){
        mat <- matrix(0,dim(Xest)[1],dim(Xest)[1])
        diag(mat) <- 1/((1-diag(Xest%*%solve(t(Xest)%*%Xm1)%*%t(Xest)))^2)
        var       <- solve(t(Xest)%*%Xest)%*%t(Xest)%*%mat%*%M%*%Xest%*%solve(t(Xest)%*%Xest)
      }}}}
    }else{
      var1            <- solve(t(Xest)%*%W%*%Xest) # revisar si dejo o pongo sigma2
      muestra$ek2     <- muestra$Y-Xm1%*%impacto.1_2
      diag(M)         <- muestra$ek2
      mat <- matrix(0,dim(Xm1)[1],dim(Xm1)[1])
      if(vcov_=="HC0"){
        diag(mat) <- 1
      }else{if(vcov_=="HC1"){
        diag(mat) <- (dim(Xm1)[1]/(dim(Xm1)[1]-length(impacto.1_2)))
      }else{if(vcov_=="HC2"){
        diag(mat) <- 1/((1-diag(Xm1%*%solve(t(Xm1)%*%Xm1)%*%t(Xm1))))
      }else{if(vcov_=="HC3"){
        diag(mat) <- 1/((1-diag(Xm1%*%solve(t(Xm1)%*%Xm1)%*%t(Xm1)))^2)
      }}}}
      for(i in 1:dim(Xm1)[2]){
        muestra[,c(paste0("Uk",i))] <- mat%*%M%*%Xm1[,i]
      }
      svydesign_$variables$Uk1<-muestra$Uk1
      svydesign_$variables$Uk2<-muestra$Uk2
      svydesign_$variables$Uk3<-muestra$Uk3
      svydesign_$variables$Uk4<-muestra$Uk4
      varmuestra           <- vcov(svytotal(~Uk1+Uk2+Uk3+Uk4,svydesign_))
      var                  <- (var1)%*%varmuestra%*%(var1)
    }
  }}
  seBeta              <- sqrt(var[3,3])
  significativo       <- ifelse(abs(beta/(seBeta))>=qt(0.975,dim(muestra)[1]-dim(Xm1)[2]),"Significativo","No Significativo")
  impact              <- data.frame(round(beta,2),round(seBeta,5),pvalue=(1-pt(abs(beta/(seBeta)),df = dim(muestra)[1]-length(impacto.1_2)))*2)
  colnames(impact)    <- c("Impacto","S.E.","P-Value")
  impact$significativo<- ifelse(impact$`P-Value`<0.01,"***",ifelse(impact$`P-Value`<0.05,"**",ifelse(impact$`P-Value`<0.1,"*","")))
  colnames(impact)[4] <- ""
  print(impact)
}



