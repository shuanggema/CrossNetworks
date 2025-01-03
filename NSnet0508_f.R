NSnet0508_f  <- function(Y_array,H_int_array,FY_int,a_int,rankY,rankA,eta0=1){
  nx      <- dim(Y_array)[1]
  ny      <- dim(Y_array)[1]
  Tt      <- dim(Y_array)[3]
  
  etaF   <- eta0/Tt
  etaA   <- eta0/Tt  # 300,150 5; 500,150 2; 50,500, 
  
  
  
  
  FY      <- FY_int
  H_array <- H_int_array
  Q_array <- array(0,dim=c(nx,nx,Tt))
  P_array <- array(0,dim=c(ny,ny,Tt))
  #A       <- diag(sqrt(max(a_int,0.01)),ny,nx)
  A       <- diag(sqrt(a_int),ny,nx)
  Au      <- A
  Av      <- diag(1,ny)
  for(t in 1:Tt){
    Theta_ty             = FY+Au%*%t(Av)%*%H_array[,,t]%*%Av%*%t(Au)
    P_array[,,t]         = 1/(1+exp(-Theta_ty)) 
    diag(P_array[,,t])   = 0
  }
  
  d0     <- nx*nx*Tt
  for(l in 1:100){
    etaF   <- etaF*0.8
    etaA   <- etaA*0.8
    FY0                    = FY
    A0                     = A
    d00                    = d0 
    P0_array               = P_array
    DFY                    = FY+etaF*(apply(Y_array-P_array,c(1,2),sum))
    DA0                     = matrix(0,nr=ny,nc=nx)
    for(t in 1:Tt){
      DA0                   = DA0+etaA/nx*(Y_array[,,t]-P_array[,,t])%*%Au%*%(t(Av)%*%H_array[,,t])
    }
    DA                      = A+DA0
    
    FY                           = Hard_f(DFY,rankY)
    FY                           = (FY+t(FY))/2
    
    A                            = Hard_f(DA,rankA)
    
    for(t in 1:Tt){
      Theta_ty             = FY+A%*%H_array[,,t]%*%t(A)
      P_array[,,t]         = 1/(1+exp(-Theta_ty))
      diag(P_array[,,t])   = 0
    }
    if(sum((FY-FY0)^2)/(sum(FY^2)+10^(-5))+sum((A-A0)^2)/(sum(A^2)+10^(-5))<5*10^(-5))break
    if(sum(A^2)>nx*ny*100)break
    d0       <- sum((Y_array-P_array)^2)
    if(d0>d00){
      FY                    = FY0
      A                     = A0
      d0                    = d00
      P_array               = P0_array
      print(l)
    }
  }
  print(d0)
  return(list(A=A,FY=FY))
}