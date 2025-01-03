XYnet0830_f  <- function(X_array,Y_array,FX_int,H_int_array,FY_int,a_int,rankX,rankH,rankY,rankA,eta0=1){
  nx      <- dim(X_array)[1]
  ny      <- dim(Y_array)[1]
  Tt      <- dim(X_array)[3]
  
  etaF   <- eta0/Tt
  etaH   <- eta0
  etaA   <- eta0/Tt# 300,150 5; 500,150 2; 50,500, 
  
  
  
  FX      <- FX_int
  FY      <- FY_int
  H_array <- H_int_array
  Q_array <- array(0,dim=c(nx,nx,Tt))
  P_array <- array(0,dim=c(ny,ny,Tt))
  #A       <- diag(sqrt(max(a_int,0.01)),ny,nx)
  A       <- diag(sqrt(a_int),ny,nx)
  Au      <- A
  Av      <- diag(1,ny)
  U_array     = array(0,dim=c(nx,rankH,Tt))
  V_array     = array(0,dim=c(nx,rankH,Tt))
  for(t in 1:Tt){
    fit_Ht     = Hard_f0(H_array[,,t],rankH)
    U_array[,,t]=fit_Ht$Bu
    V_array[,,t]=fit_Ht$Bv
  }
  
  for(t in 1:Tt){
    AU                   = A%*%U_array[,,t]
    AV                   = A%*%V_array[,,t]
    Theta_ty             = FY+AU%*%t(AV)
    P_array[,,t]         = 1/(1+exp(-Theta_ty)) 
    diag(P_array[,,t])   = 0
  }
  
  
  
  d0     <- nx*nx*Tt
  for(l in 1:100){
    etaF   <- etaF*0.8
    etaH   <- etaH*0.8
    etaA   <- etaA*0.8
    
    FY0                    = FY
    A0                     = A
    Au0                    = Au
    Av0                    = Av 
    d00                    = d0 
    P0_array               = P_array            
    
    
    DFY                    = FY+etaF*(apply(Y_array-P_array,c(1,2),sum))
    
    DA0                     = matrix(0,nr=ny,nc=nx)
    for(t in 1:Tt){
      S1                    = (Y_array[,,t]-P_array[,,t])%*%Au
      S2                    = t(Av)%*%U_array[,,t]
      DA0                   = DA0+etaA/nx*S1%*%S2%*%t(V_array[,,t])
    }
    DA                      = A+DA0
    
    
    FY                           = Hard_f(DFY,rankY)
    FY                           = (FY+t(FY))/2
    
    A_fit                        = Hard_f0(DA,rankA)
    A                            = A_fit$B
    Au                           = A_fit$Bu
    Av                           = A_fit$Bv
    
    for(t in 1:Tt){
      AV0                  =t(Av)%*%U_array[,,t]
      AU                   = Au%*%AV0
      AV1                  =t(Av)%*%V_array[,,t]
      AV                   = Au%*%AV1
      
      Theta_ty             = FY+AU%*%t(AV)
      P_array[,,t]         = 1/(1+exp(-Theta_ty)) 
      diag(P_array[,,t])   = 0
    }
    
    if(sum((FY-FY0)^2)/(sum(FY^2)+10^(-5))+sum((A-A0)^2)/(sum(A^2)+10^(-5))<5*10^(-5))break
    if(sum(A^2)>nx*ny*100)break
    
    
    d0       <- sum((Y_array-P_array)^2)
    if(d0>d00){
      FY                    = FY0
      A                     = A0
      Au                    = Au0
      Av                    = Av0 
      d0                    = d00
      P_array               = P0_array
      print(l)
    }
  }
  print(d0)
  return(list(FX=FX,H_array=H_array,A=A,FY=FY))
}