Hard_f    <- function(B,K){
  library(irlba)
  B_svd <- irlba(B,K)
  pos     <- 1:K
  B_th    <- B_svd$u[,pos]%*%diag(B_svd$d[pos],nrow=length(pos))%*%t(B_svd$v[,pos])
  return(B_th)  
}

Hard_f0    <- function(B,K){
  library(irlba)
  B_svd <- irlba(B,K)
  pos     <- 1:K
  B_th    <- B_svd$u[,pos]%*%diag(B_svd$d[pos],nrow=length(pos))%*%t(B_svd$v[,pos])
  return(list(Bu=B_svd$u[,pos]%*%diag(B_svd$d[pos],nrow=length(pos)),Bv=B_svd$v[,pos],B=B_th) ) 
}




Pard_f    <- function(B,K){
  library(irlba)
  B_svd <- irlba(B,K)
  pos     <- 1:K
  B_th    <- B_svd$u[,pos]%*%diag(sqrt(B_svd$d[pos]),nrow=length(pos))
  return(B_th)  
}

# BIC_Xf     <- function(X_array,FX,H_array,rankX,rankH){
#     n       <- dim(X_array)[1]
#     Tt      <- dim(X_array)[3]
#     FX_array<- array(FX,dim=c(n,n,Tt))
#     LX      <- sum((-X_array*(FX_array+H_array)+log(1+exp(FX_array+H_array))))
#     # BIC_v1  <- 2*LX+log(n^2*Tt)*(rankX*n)+log(n^2)*Tt*rankH*n
#     # BIC_v2  <- 2*LX+log(n^2*Tt)*(rankX*n)+2*Tt*rankH*n
#     # return(list(LX=LX,BIC_v1=BIC_v1,BIC_v2=BIC_v2))
#     BIC_vx    <- LX+log(n^2*Tt)*rankX*n
#     BIC_vh    <- vector()
#     for(t in 1:Tt){
#       BIC_vh[t] <- 2*sum((-X_array[,,t]*(FX+H_array[,,t])+log(1+exp(FX+H_array[,,t]))))+log(n^2/2)*rankH*n
#     }
#     return(list(LX=LX,BIC_vx=BIC_vx,BIC_vh=BIC_vh))
#      
# }



BIC_Xf0     <- function(X_array,FX,H_array,rankX,rankH){
  n       <- dim(X_array)[1]
  Tt      <- dim(X_array)[3]
  FX_array<- array(FX,dim=c(n,n,Tt))
  LX      <- sum((-X_array*(FX_array+H_array)+log(1+exp(FX_array+H_array))))
  BIC_v1  <- LX+2*log(n^2*Tt)*(rankX*n)+0.5*log(n^2)*Tt*rankH*n
  BIC_v2  <- LX+log(n^2*Tt)*(rankX*n)+2*Tt*rankH*n
  return(list(LX=LX,BIC_v1=BIC_v1,BIC_v2=BIC_v2))
  # BIC_vx    <- LX+log(n^2*Tt)*rankX*n
  # BIC_vh    <- vector()
  # for(t in 1:Tt){
  #   BIC_vh[t] <- 2*sum((-X_array[,,t]*(FX+H_array[,,t])+log(1+exp(FX+H_array[,,t]))))+log(n^2/2)*rankH*n
  # }
  # return(list(LX=LX,BIC_vx=BIC_vx,BIC_vh=BIC_vh))
  
}



BIC_Yf0     <- function(Y_array,FY,H_array,A,rankY,rankA){
  n       <- dim(Y_array)[1]
  Tt      <- dim(Y_array)[3]
  FY_array<- array(FY,dim=c(n,n,Tt))
  LY      <- 0
  LY0     <- 0
  for(t in 1:Tt){
    LY    <- LY+sum((-Y_array[,,t]*(FY+A%*%(H_array[,,t])%*%t(A))+log(1+exp(FY+A%*%(H_array[,,t])%*%t(A)))))
    LY0   <- LY0+sum((Y_array[,,t]-1/(1+exp(-FY-A%*%(H_array[,,t])%*%t(A))))^2)
  }
  

  
  BIC_v1  <- LY+2*log(n^2*Tt)*(rankY*n)+0.5*log(n^2*Tt)*rankA*n
  BIC_v2  <- LY+log(n^2*Tt)*(rankY*n)+2*rankA*n
  return(list(LY=LY,LY0=LY0,BIC_v1=BIC_v1,BIC_v2=BIC_v2))
  
}