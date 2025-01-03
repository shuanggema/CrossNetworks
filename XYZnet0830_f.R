XYZnet0830_f <- function(X_array,Y_array,Z_mat,W_mat,FX_int,H_int_array,FY_int,A_int,rankX,rankH,rankY,rankA,eta0,tau=1){
  nx      <- dim(X_array)[1]
  ny      <- dim(Y_array)[1]
  Tt      <- dim(X_array)[3]
  
  mu      <- mean(Z_mat)
  alpha   <- 0
  sig     <- 1
  
  
  etaF   <- eta0/Tt
  etaH   <- eta0
  etaA   <- eta0/Tt # 300,150 5; 500,150 2; 50,500, 
  
  FY       = FY_int
  FX       = FX_int
  H_array  = H_int_array
  Q_array  = array(0,dim=c(nx,nx,Tt))
  P_array  = array(0,dim=c(ny,ny,Tt))
  Z_vec    = as.vector(Z_mat)
  
  A_fit   <- Hard_f0(A_int,rankA)
  A       <- A_fit$B
  Au      <- A_fit$Bu
  Av      <- A_fit$Bv
  
  U_array     = array(0,dim=c(nx,rankH,Tt))
  V_array     = array(0,dim=c(nx,rankH,Tt))
  for(t in 1:Tt){
   Ht_fit     = Hard_f0(H_array[,,t],rankH)
   U_array[,,t]=Ht_fit$Bu
   V_array[,,t]=Ht_fit$Bv
  }
  
  
  D_mat                  = matrix(0,nr=ny,nc=Tt)
  for(t in 1:Tt){
    Theta_tx             = FX+H_array[,,t]
    Q_array[,,t]         = 1/(1+exp(-Theta_tx)) 
    diag(Q_array[,,t])   = 0
    
    AU                   = A%*%U_array[,,t]
    AV                   = A%*%V_array[,,t]
    Theta_ty             = FY+AU%*%t(AV)
    P_array[,,t]         = 1/(1+exp(-Theta_ty)) 
    diag(P_array[,,t])   = 0
    D_mat[,t]            = Z_mat[,t]-mu-alpha*(Theta_ty)%*%W_mat[,t]
  }
  
  
  d0     <- sum((Y_array-P_array)^2)
  for(l in 1:100){
    etaF                   = etaF*0.95
    etaH                   = etaH*0.95
    etaA                   = etaA*0.95
    
    FX0                    = FX
    H0_array               = H_array
    FY0                    = FY
    A0                     = A
    Au0                    = Au
    Av0                    = Av 
    alpha0                 = alpha
    mu0                    = mu
    d00                    = d0 
    Q0_array               = Q_array
    P0_array               = P_array
    D0_mat                 = D_mat
    
    
    
    DFX                    = FX+etaF/nx*(apply(X_array-Q_array,c(1,2),sum))
    DFY                    = FY+etaF/nx*(apply(Y_array-P_array,c(1,2),sum))
    
    DW_array               = array(0,dim=c(nx,nx,Tt)) 
    for(t in 1:Tt){
      DW_array[,,t]        = W_mat[,t]%*%t(D_mat[,t])
      dDW_v                = diag(DW_array[,,t])
      DFY                  = DFY+tau*etaF*alpha/sig/nx*(DW_array[,,t] +t(DW_array[,,t] )-diag(dDW_v))
    }
    
    DA0                    = matrix(0,nr=ny,nc=nx)
    for(t in 1:Tt){
      S1                    = (Y_array[,,t]-P_array[,,t])%*%Au
      S2                    = t(Av)%*%U_array[,,t]
      S3                    = S2%*%t(V_array[,,t]) 
      S4                    = Au%*%S3 ##S4S3=AHt;S1S3=(Y-P)AH
      DA0                   = (DA0+etaA/nx*S1%*%S3+
                                 tau*alpha*etaA/nx/sig*(D_mat[,t]%*%t(W_mat[,t])%*%S4+W_mat[,t]%*%t(D_mat[,t])%*%S4))
    }
    DA                     = A+DA0
    
    
    
    DH0_array              = H0_array+etaH/nx*(X_array-Q_array)
    DH_array               = DH0_array
    for(t in 1:Tt){
      Q10                  = t(Au)%*%(Y_array[,,t]-P_array[,,t])%*%Au
      Q1                   = Av%*%Q10%*%t(Av) #t(A)(Y-P)A
      Q1d                  = Q1-1/2*diag(diag(Q1))
      
      Q20                  = t(Au)%*%W_mat[,t]%*%t(D_mat[,t])%*%Au
      Q2                   = Av%*%Q20%*%t(Av)
      Q2d                  = Q2+t(Q2)-diag(Q2)
      DH_array[,,t]        = DH0_array[,,t]+etaH/nx*Q1d+
        tau*alpha*etaH/nx/sig*Q2d
    }
    
    ##### update H, A, and FY
    for(t in 1:Tt){
      fit_Ht                       = Hard_f0(DH_array[,,t],rankH)
      H_array[,,t]                 = (fit_Ht$B+t(fit_Ht$B))/2
      U_array[,,t]                 = fit_Ht$Bu
      V_array[,,t]                 = fit_Ht$Bv
    }
    
    FY                     = Hard_f(DFY,rankY)
    FY                     = (FY+t(FY))/2
    FX                     = Hard_f(DFX,rankX)
    FX                     = (FX+t(FX))/2
    
    A_fit                        = Hard_f0(DA,rankA)
    A                            = A_fit$B
    Au                           = A_fit$Bu
    Av                           = A_fit$Bv
    ##########
    
    
    G_vec= vector()
    for(t in 1:Tt){
      G_vec= c(G_vec,as.vector((FY+A%*%H_array[,,t]%*%t(A))%*%W_mat[,t]))
    }
    fit_lm   <- lm(Z_vec~G_vec)
    mu       <- fit_lm$coef[1]
    alpha    <- fit_lm$coef[2]
    sig      <- var(fit_lm$residual)
    
    
    
    
    if(sum((FY-FY0)^2)/(sum(FY^2)+10^(-5))+sum((A-A0)^2)/(sum(A^2)+10^(-5))<5*10^(-5) & l>2 )break
    
    
    D_mat                  = matrix(0,nr=ny,nc=Tt)
    for(t in 1:Tt){
      Theta_tx             = FX+H_array[,,t]
      Q_array[,,t]         = 1/(1+exp(-Theta_tx)) 
      diag(Q_array[,,t])   = 0
      
      AU                   = A%*%U_array[,,t]
      AV                   = A%*%V_array[,,t]
      Theta_ty             = FY+AU%*%t(AV)
      P_array[,,t]         = 1/(1+exp(-Theta_ty)) 
      diag(P_array[,,t])   = 0
      D_mat[,t]            = Z_mat[,t]-mu-alpha*(Theta_ty)%*%W_mat[,t]
    }
    d0       <- sum((Y_array-P_array)^2)
    # if(d0>d00){
    #   print(l)
    #   H_array               = H0_array
    #   FY                    = FY0
    #   FX                    = FX0
    #   A                     = A0
    #   Au                    = Au0
    #   Av                    = Av0 
    #   alpha                 = alpha0
    #   mu                    = mu0
    #   d0                    = d00
    #   D_mat                 = D0_mat
    #   P_array               = P0_array
    #   Q_array               = Q0_array
    # }
    
    
    #if(sum(A^2)>nx*ny*100)break
    # print(sum((Y_array-P_array)^2))
    # print(sum((A-AT)^2)/sum(AT^2))
    # print(sum((A-A0)^2)/(sum(A^2)+10^(-5)))
  }
  #print(sum((A-AT)^2)/sum(AT^2))
  print(d0)
  return(list(FX=FX,FY=FY,H_array=H_array,A=A,coef=c(mu,alpha),D_mat=D_mat))
}