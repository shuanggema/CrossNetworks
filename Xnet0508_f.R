Xnet0508_f<- function(X_array,FX_int,H_int_array,rankX,rankH,eta0=10){
   # consider the diag 
  n      <- dim(X_array)[1]
  Tt     <- dim(X_array)[3]
  
  etaF   <- eta0/Tt #### etaF has been selected.
  etaH   <- 2*eta0# T=300 etaH=10
  
  
  FX      <- FX_int
  H_array <- H_int_array
  Q_array <- array(0,dim=c(n,n,Tt))
  
  for(t in 1:Tt){
    Theta_t              = FX+H_array[,,t]
    Q_array[,,t]         = 1/(1+exp(-Theta_t))  
    diag(Q_array[,,t])   = 0
  }
  
  
  for(l in 1:100){
    etaF   <- etaF*0.8
    etaH   <- etaH*0.8
    FX0                    = FX
    H0_array               = H_array
    DF                     = FX+etaF*(apply(X_array-Q_array,c(1,2),sum))
    DH_array               = H0_array+etaH*(X_array-Q_array)
    
    
    # soft-f H,G,F
    
    FX                           = Hard_f(DF,rankX)
    FX                           = (FX+t(FX))/2
    #FX                           = FXT
    for(t in 1:Tt){
      H_array[,,t]                 = Hard_f(DH_array[,,t],rankH)
      H_array[,,t]                 = (H_array[,,t]+t(H_array[,,t]))/2
    }
    
    for(t in 1:Tt){
      Theta_t              = FX+H_array[,,t]
      Q_array[,,t]         = 1/(1+exp(-Theta_t))  
      diag(Q_array[,,t])   = 0
    }
    
    if(sum((FX-FX0)^2)/(sum(FX^2)+10^(-5))+sum((H_array-H0_array)^2)/(sum(H_array^2)+10^(-5))<5*10^(-5) )break
  }
  print(sum((X_array-Q_array)^2))
  
  return(list(FX=FX,H_array=H_array))
}