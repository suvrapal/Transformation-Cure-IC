library(numDeriv)
library(knitr)
#x1 and x2 both linked to phi, PH structure
# generating interval censored data with two covariates from Box-Cox model for alpha in (0,1]
#special case where alpha = 0 -> poisson cure rate model
#using Weibull dist to generate t
#let baseline hazard be weibull, h=k*t^(k-1), k>0, lambda (other param)=1

LR_int=function(y1,len1,l1){
  if(y1>0 & y1<=l1){
    a = c(.Machine$double.eps,l1)
  }else{
    k = as.integer((y1-l1)/len1)+1
    a = c(l1+((k-1)*len1),l1+(k*len1))
  }
  return(a)
}

data_0_BC=function(n,b0,b1,b2,g1,g2,gam1,cenrate){ # gam1,gam2=1 baseline Weibull parameters
  
  L=rep(NA,n)
  R=rep(NA,n)
  d=rep(NA,n)
  
  x1 = rbinom(n=n, size=1, prob=0.5) # binary covariate
  x2=runif(n,min = 0.1,max = 20) # continuous covariate
  for(i in 1:n){
    if(x2[i]<=0){
      x2[i] = .Machine$double.eps
    }
  }
  phi = exp(b0 + (b1*x1) + (b2*x2))
  
  U = runif(n, min=0, max=1)
  C = rexp(n, rate=cenrate)
  
  p0 = exp(-phi)
  
  count.obs = 0
  count.cure = 0
  
  for(i in 1:n){
    if(U[i]<=p0[i]){
      L[i] = C[i]
      R[i] = Inf
      d[i] = 0
      count.cure=count.cure+1
    }else{
      U1 = runif(1,min=0,max=1)
      y = (-exp(-(g1*x1[i])-(g2*x2[i]))*log(1+(exp(-b0-(b1*x1[i])-(b2*x2[i]))*log(p0[i]+((1-p0[i])*U1)))))^gam1
      t = min(y,C[i])
      if(t==C[i]){    
        L[i] = C[i]
        R[i] = Inf
        d[i] = 0
      }else{
        len = runif(1,0.1,0.5)
        l = runif(1,0,1)
        ans = LR_int(t,len,l)
        L[i] = ans[1]
        R[i] = ans[2]
        d[i] = 1
        count.obs = count.obs + 1
      }# end of inner else
    }# end of outer else
  }# end of for
  
  # print("count.cure")
  # print(count.cure)
  # print("count.obs")
  # print(count.obs)
  # 
  return(data.frame(L,R,d,x1,x2))
  
}
data_gen_BC=function(n,alpha,b0,b1,b2,g1,g2,gam1,cenrate){ # alpha is the BC index parameter, lam is baselineexp parameter
  
  L=rep(NA,n)
  R=rep(NA,n)
  d=rep(NA,n)
  
  x1=rbinom(n=n,size=1,prob=0.5) # binary covariate
  x2=runif(n,min = 0.1,max = 20) # continuous covariate
  
  phi=exp(b0+(b1*x1)+(b2*x2))/(1+(alpha*exp(b0+(b1*x1)+(b2*x2))))
  
  U=runif(n,min=0,max=1)
  C=rexp(n,rate=cenrate)
  
  p0 = (1-(alpha*phi))^(1/alpha)
  
  count.obs=0
  count.cure=0
  
  for(i in 1:n){
    if(U[i]<=p0[i]){
      L[i]=C[i]
      R[i]=Inf
      d[i]=0
      count.cure=count.cure+1
    }else{
      U1 = runif(1,min=0,max=1)
      y=(-exp(-(g1*x1[i])-(g2*x2[i]))*log(((alpha*phi[i])+((p0[i]+((1-p0[i])*U1))^alpha)-1)/(alpha*phi[i])))^gam1
      t=min(y,C[i])
      if(t==C[i]){    
        L[i]=C[i]
        R[i]=Inf
        d[i]=0
      }else{
        len=runif(1,0.1,0.5)
        l=runif(1,0,1)
        ans=LR_int(t,len,l)
        L[i]=ans[1]
        R[i]=ans[2]
        d[i]=1
        count.obs = count.obs + 1
      }# end of inner else
    }# end of outer else
  }# end of for
  # 
  # print("count.cure")
  # print(count.cure)
  # print("count.obs")
  # print(count.obs)
  # 
  return(data.frame(L,R,d,x1,x2))
  
}
std0 = function(param){
  
  b0 = param[1]
  b1 = param[2]
  b2 = param[3]
  g1 = param[4]
  g2 = param[5]
  gam1 = param[6]
  
  phi.t = exp(b0+b1*x1t+b2*x2t)
  phi.c = exp(b0+b1*x1c+b2*x2c)
  
  S0.L.cens = exp(-L_cens^(1/gam1))
  S0.L.obs = exp(-L_obs^(1/gam1))
  S0.R.obs = exp(-R_obs^(1/gam1))
  
  F.c = 1 - S0.L.cens^(exp(g1*x1c+g2*x2c))
  F.tl = 1 - S0.L.obs^(exp(g1*x1t+g2*x2t))
  F.tr = 1 - S0.R.obs^(exp(g1*x1t+g2*x2t))
  Sp.c = exp(-phi.c*F.c)
  Sp.tl = exp(-phi.t*F.tl)
  Sp.tr = exp(-phi.t*F.tr)
  dum1 = rep(NA,length(Sp.tl))
  for(i in 1:length(dum1)){
    dum1[i] = max((Sp.tl[i]-Sp.tr[i]),.Machine$double.eps)
  }
  
  out1 = sum(log(Sp.c)) + sum(log(dum1))  #return log-likelihood function
  return(-out1)
}
stdgen = function(param){
  
  b0 = param[1]
  b1 = param[2]
  b2 = param[3]
  g1 = param[4]
  g2 = param[5]
  gam1 = param[6]
  
  phi.t = exp(b0+b1*x1t+b2*x2t)/(1+alpha*exp(b0+b1*x1t+b2*x2t))
  phi.c = exp(b0+b1*x1c+b2*x2c)/(1+alpha*exp(b0+b1*x1c+b2*x2c))
  
  S0.L.cens = exp(-L_cens^(1/gam1))
  S0.L.obs = exp(-L_obs^(1/gam1))
  S0.R.obs = exp(-R_obs^(1/gam1))
  
  F.c = 1 - S0.L.cens^(exp(g1*x1c+g2*x2c))
  F.tl = 1 - S0.L.obs^(exp(g1*x1t+g2*x2t))
  F.tr = 1 - S0.R.obs^(exp(g1*x1t+g2*x2t))
  Sp.c = (1-alpha*phi.c*F.c)^(1/alpha)
  Sp.tl = (1-alpha*phi.t*F.tl)^(1/alpha)
  Sp.tr = (1-alpha*phi.t*F.tr)^(1/alpha)
  dum1 = rep(NA,length(Sp.tl))
  for(i in 1:length(dum1)){
    dum1[i] = max((Sp.tl[i]-Sp.tr[i]),.Machine$double.eps)
  }
  
  out1 = sum(log(Sp.c)) + sum(log(dum1))  #return log-likelihood function
  return(-out1)
}
#EM algorithm for alpha=0
BC_0_EM_Wei=function(mydata,tol,maxit,increment,b0,b1,b2,g1,g2,gam1){
  
  data_obs = mydata[mydata$d==1,]
  data_cens = mydata[mydata$d==0,]
  L_obs = data_obs$L
  R_obs = data_obs$R
  L_cens = data_cens$L
  x1t = data_obs$x1
  x1c = data_cens$x1
  x2t = data_obs$x2
  x2c = data_cens$x2
  
  x1 = c(x1t,x1c)
  x2 = c(x2t,x2c)
  
  pnew = rep(0,6)
  pold = rep(0,6)
  
  
  pold[1] = runif(1,b0-increment*abs(b0), b0+increment*abs(b0))
  pold[2] = runif(1,b1-increment*abs(b1), b1+increment*abs(b1))
  pold[3] = runif(1,b2-increment*abs(b2), b2+increment*abs(b2))
  pold[4] = runif(1,g1-increment*abs(g1), g1+increment*abs(g1))
  pold[5] = runif(1,g2-increment*abs(g2), g2+increment*abs(g2))
  pold[6] = runif(1,gam1-increment*abs(gam1), gam1+increment*abs(gam1))
  
  continue = TRUE
  count = 0
  iter = 1
  
  while(continue){
    
    #print(iter)
    #Q function:E step: 
    #Initializing parameter estimates for current iteration 
    #and conditional expectation of cured status (e.W).
    #e.W only needed in Q for terms where d=0 (censored obs)
    
    
    S0.L.obs = exp(-L_obs^(1/pold[6])) #estimate of baseline hazard using Turnbull method
    S0.R.obs = exp(-R_obs^(1/pold[6]))
    S0.L.cens = exp(-L_cens^(1/pold[6]))
    phi.cens.r = exp(pold[1]+pold[2]*x1c+pold[3]*x2c)
    F.L.cens.r = 1-S0.L.cens^(exp(pold[4]*x1c+pold[5]*x2c))
    Sp.L.cens.r = exp(-phi.cens.r*F.L.cens.r)
    p0.r = exp(-phi.cens.r)  
    Su.r = (Sp.L.cens.r-p0.r)/(1-p0.r)
    e.W = (Sp.L.cens.r-p0.r)/Sp.L.cens.r
    
    
    #E step
    Q = function(par=c(b0,b1,b2,g1,g2,gam1)){
      S0.L.obs.est = exp(-L_obs^(1/par[6]))
      S0.R.obs.est = exp(-R_obs^(1/par[6]))
      S0.L.cens.est = exp(-L_cens^(1/par[6]))
      phi.obs.est = exp(par[1]+par[2]*x1t+par[3]*x2t)
      phi.cens.est = exp(par[1]+par[2]*x1c+par[3]*x2c)
      F.L.obs.est = 1-S0.L.obs.est^(exp(par[4]*x1t+par[5]*x2t))
      F.R.obs.est = 1-S0.R.obs.est^(exp(par[4]*x1t+par[5]*x2t))
      F.L.cens.est = 1-S0.L.cens.est^(exp(par[4]*x1c+par[5]*x2c))
      Sp.L.obs.est = exp(-phi.obs.est*F.L.obs.est)
      Sp.R.obs.est = exp(-phi.obs.est*F.R.obs.est)
      Sp.L.cens.est = exp(-phi.cens.est*F.L.cens.est)
      p0.est = exp(-phi.cens.est)  #only needed in Q for terms where d=0 (censored obs)
      Su.est = (Sp.L.cens.est-p0.est)/(1-p0.est) #only needed in Q for terms where d=0 (censored obs)
      for(i in 1:length(p0.est)){
        p0.est[i] = max(p0.est[i], .Machine$double.eps)
        Su.est[i] = max(Su.est[i], .Machine$double.eps)
      }
      dum = Sp.L.obs.est - Sp.R.obs.est
      
      for(i in 1:length(x1t)){
        dum[i] = max(dum[i], .Machine$double.eps)
      }
      
      dum1 = log(dum)
      dum2 = (1-e.W)*log(p0.est) + e.W*log((1-p0.est)*Su.est)
      result = sum(dum1)+sum(dum2)
      return(-result)
    } #end of Q Function
    
    #M step
    # pnew = tryCatch({optim(par=c(pold),fn=Q, method="Nelder-Mead",control=list(maxit=1))$par
    # },error=function(e){
    #   pnew = c(0,0,0,0,0,0)
    #   return(pnew)
    # }
    # )
    pnew = optim(par=c(pold),fn=Q, method="Nelder-Mead",control=list(maxit=100))$par
    #if error in estimation, stop checking for convergence and return zeros
    if(0 %in% pnew){ 
      
      continue = FALSE
      
      pnew = c(0,0,0,0,0,0)
      
      break
      
    }
    else{
      iter = iter+1
      continue = max(abs((pnew-pold)/pold))>tol 
      
      pold = pnew
    }#end of else
  }#end of while
  out = rep(NA,13)
  
  #if error in estimation, don't calculate SE and return all zeros
  if(0 %in% pnew){
    out = rep(0,13)
  }else{
    #updated "std" from parametric setting to use Turnbull estimate  
    std = function(param){
      
      b0 = param[1]
      b1 = param[2]
      b2 = param[3]
      g1 = param[4]
      g2 = param[5]
      gam1 = param[6]
      
      phi.t = exp(b0+b1*x1t+b2*x2t)
      phi.c = exp(b0+b1*x1c+b2*x2c)
      
      S0.L.cens = exp(-L_cens^(1/gam1))
      S0.L.obs = exp(-L_obs^(1/gam1))
      S0.R.obs = exp(-R_obs^(1/gam1))
      
      F.c = 1 - S0.L.cens^(exp(g1*x1c+g2*x2c))
      F.tl = 1 - S0.L.obs^(exp(g1*x1t+g2*x2t))
      F.tr = 1 - S0.R.obs^(exp(g1*x1t+g2*x2t))
      Sp.c = exp(-phi.c*F.c)
      Sp.tl = exp(-phi.t*F.tl)
      Sp.tr = exp(-phi.t*F.tr)
      dum1 = rep(NA,length(Sp.tl))
      for(i in 1:length(dum1)){
        dum1[i] = max((Sp.tl[i]-Sp.tr[i]),.Machine$double.eps)
      }
      
      out1 = sum(log(Sp.c)) + sum(log(dum1))  #return log-likelihood function
      return(out1)
    }
    
    hessmat1 = hessian(std,c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6]),method="Richardson")
    FI = solve(-1*hessmat1)
    if(any(is.na(FI))){
      print(FI)
      out = c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6],alpha,0,0,0,0,0,0,NA)
    }
    else if(FI[1,1]<0 | FI[2,2]<0 | FI[3,3]<0| FI[4,4]<0 | FI[5,5]<0 | FI[6,6]<0){ #added condition to handle if FI non-invertible
      out = c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6],alpha,0,0,0,0,0,0,NA)
    }
    else{
      std1 = sqrt(FI[1,1])
      std2 = sqrt(FI[2,2])
      std3 = sqrt(FI[3,3])
      std4 = sqrt(FI[4,4])
      std5 = sqrt(FI[5,5])
      std6 = sqrt(FI[6,6])
      
      std.mle = std(c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6]))
      
      out = c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6],0,std1,std2,std3,std4,std5,std6,std.mle)
      #print(out)
    }
  }#end of else
  return(out)
} 
#EM algorithm for alpha in (0,1]
BC_gen_EM_Wei=function(mydata,tol,maxit,increment,alpha,b0,b1,b2,g1,g2,gam1){
  data_obs = mydata[mydata$d==1,]
  data_cens = mydata[mydata$d==0,]
  L_obs = data_obs$L
  R_obs = data_obs$R
  L_cens = data_cens$L
  x1t = data_obs$x1
  x1c = data_cens$x1
  x2t = data_obs$x2
  x2c = data_cens$x2
  
  x1 = c(x1t,x1c)
  x2 = c(x2t,x2c)
  
  pnew = rep(0,6)
  pold = rep(0,6)
  
  pold[1] = runif(1,b0-increment*abs(b0), b0+increment*abs(b0))
  pold[2] = runif(1,b1-increment*abs(b1), b1+increment*abs(b1))
  pold[3] = runif(1,b2-increment*abs(b2), b2+increment*abs(b2))
  pold[4] = runif(1,g1-increment*abs(g1), g1+increment*abs(g1))
  pold[5] = runif(1,g2-increment*abs(g2), g2+increment*abs(g2))
  pold[6] = runif(1,gam1-increment*abs(gam1), gam1+increment*abs(gam1))
  
  
  continue = TRUE
  count = 0
  iter = 1
  
  while(continue){
    
    #print(iter)
    #Q function:E step: 
    #Initializing parameter estimates for current iteration 
    #and conditional expectation of cured status (e.W).
    #e.W only needed in Q for terms where d=0 (censored obs)
    
    S0.L.obs = exp(-L_obs^(1/pold[6])) #estimate of baseline hazard using Turnbull method
    S0.R.obs = exp(-R_obs^(1/pold[6]))
    S0.L.cens = exp(-L_cens^(1/pold[6]))
    phi.cens.r = exp(pold[1]+pold[2]*x1c+pold[3]*x2c)/(1+alpha*exp(pold[1]+pold[2]*x1c+pold[3]*x2c))
    F.L.cens.r = 1-S0.L.cens^(exp(pold[4]*x1c+pold[5]*x2c))
    Sp.L.cens.r = (1-alpha*phi.cens.r*(F.L.cens.r))^(1/alpha)
    p0.r= (1-alpha*phi.cens.r)^(1/alpha)  
    Su.r = (Sp.L.cens.r-p0.r)/(1-p0.r)
    e.W = (Sp.L.cens.r-p0.r)/Sp.L.cens.r
    
    
    #E step
    Q = function(par=c(b0,b1,b2,g1,g2,gam1)){
      S0.L.obs.est = exp(-L_obs^(1/par[6]))
      S0.R.obs.est = exp(-R_obs^(1/par[6]))
      S0.L.cens.est = exp(-L_cens^(1/par[6]))
      phi.obs.est = (exp(par[1]+par[2]*x1t+par[3]*x2t))/(1+alpha*exp(par[1]+par[2]*x1t+par[3]*x2t))
      phi.cens.est = (exp(par[1]+par[2]*x1c+par[3]*x2c))/(1+alpha*exp(par[1]+par[2]*x1c+par[3]*x2c))
      F.L.obs.est = 1-S0.L.obs.est^(exp(par[4]*x1t+par[5]*x2t))
      F.R.obs.est = 1-S0.R.obs.est^(exp(par[4]*x1t+par[5]*x2t))
      F.L.cens.est = 1-S0.L.cens.est^(exp(par[4]*x1c+par[5]*x2c))
      Sp.L.obs.est = (1-alpha*phi.obs.est*F.L.obs.est)^(1/alpha)
      Sp.R.obs.est = (1-alpha*phi.obs.est*F.R.obs.est)^(1/alpha)
      Sp.L.cens.est = (1-alpha*phi.cens.est*F.L.cens.est)^(1/alpha)
      p0.est = (1-alpha*phi.cens.est)^(1/alpha) #only needed in Q for terms where d=0 (censored obs)
      Su.est = (Sp.L.cens.est-p0.est)/(1-p0.est) #only needed in Q for terms where d=0 (censored obs)
      for(i in 1:length(p0.est)){
        p0.est[i] = max(p0.est[i], .Machine$double.eps)
        Su.est[i] = max(Su.est[i], .Machine$double.eps)
      }
      dum = Sp.L.obs.est - Sp.R.obs.est
      
      # for(i in 1:length(x1t)){
      #   dum[i] = max((Sp.L.obs.est - Sp.R.obs.est), .Machine$double.eps)
      # }
      
      dum1 = log(dum)
      dum2 = (1-e.W)*log(p0.est) + e.W*log((1-p0.est)*Su.est)
      result = sum(dum1)+sum(dum2)
      return(-result)
    } #end of Q Function
    
    #M step
    # pnew = tryCatch({optim(par=c(pold),fn=Q, method="Nelder-Mead",control=list(maxit=1))$par
    # },error=function(e){
    #   pnew = c(0,0,0,0,0,0)
    #   return(pnew)
    # }
    # )
    pnew = optim(par=c(pold),fn=Q, method="Nelder-Mead",control=list(maxit=100))$par
    #if error in estimation, stop checking for convergence and return zeros
    if(0 %in% pnew){ 
      
      continue = FALSE
      
      pnew = c(0,0,0,0,0,0)
      
      break
      
    }
    else{
      iter = iter+1
      continue = max(abs((pnew-pold)/pold))>tol 
      
      pold = pnew
    }#end of else
  }#end of while
  out = rep(NA,14)
  
  #if error in estimation, don't calculate SE and return all zeros
  if(0 %in% pnew){
    out = rep(0,14)
  }else{
    std = function(param){
      
      b0 = param[1]
      b1 = param[2]
      b2 = param[3]
      g1 = param[4]
      g2 = param[5]
      gam1 = param[6]
      
      phi.t = exp(b0+b1*x1t+b2*x2t)/(1+alpha*exp(b0+b1*x1t+b2*x2t))
      phi.c = exp(b0+b1*x1c+b2*x2c)/(1+alpha*exp(b0+b1*x1c+b2*x2c))
      
      S0.L.cens = exp(-L_cens^(1/gam1))
      S0.L.obs = exp(-L_obs^(1/gam1))
      S0.R.obs = exp(-R_obs^(1/gam1))
      
      F.c = 1 - S0.L.cens^(exp(g1*x1c+g2*x2c))
      F.tl = 1 - S0.L.obs^(exp(g1*x1t+g2*x2t))
      F.tr = 1 - S0.R.obs^(exp(g1*x1t+g2*x2t))
      Sp.c = (1-alpha*phi.c*F.c)^(1/alpha)
      Sp.tl = (1-alpha*phi.t*F.tl)^(1/alpha)
      Sp.tr = (1-alpha*phi.t*F.tr)^(1/alpha)
      dum1 = rep(NA,length(Sp.tl))
      for(i in 1:length(dum1)){
        dum1[i] = max((Sp.tl[i]-Sp.tr[i]),.Machine$double.eps)
      }
      
      out1 = sum(log(Sp.c)) + sum(log(dum1))  #return log-likelihood function
      return(out1)
    }
    
    hessmat1 = hessian(std,c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6]),method="Richardson")
    FI = solve(-1*hessmat1)
    
    if(any(is.na(FI))){
      print(FI)
      out = c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6],alpha,0,0,0,0,0,0,NA)
    }
    else if(FI[1,1]<0 | FI[2,2]<0 | FI[3,3]<0| FI[4,4]<0 | FI[5,5]<0 | FI[6,6]<0){ #added condition to handle if FI non-invertible
      out = c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6],alpha,0,0,0,0,0,0,NA)
    }
    else{
      std1 = sqrt(FI[1,1])
      std2 = sqrt(FI[2,2])
      std3 = sqrt(FI[3,3])
      std4 = sqrt(FI[4,4])
      std5 = sqrt(FI[5,5])
      std6 = sqrt(FI[6,6])
      
      std.mle = std(c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6]))
      out = c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6],alpha,std1,std2,std3,std4,std5,std6,std.mle)
    }
  }#end of else
  return(out)
} 

# calling the function
# g1 should be negative and g2 should be positive
# b1 is negative and b2 is positive
#Par Setting 1 (alpha = 0), n=200
b0 = 0.9
b1 = -0.6
b2 = -0.1
g1 = 1
g2 = -0.2
gam1 = 0.316
alpha.true=0.75
alpha=0.75
increment = 0.2
cenrate = 0.2
tol = 10^-4
maxit = 1000

M=250
mc.out.em = matrix(nrow=M,ncol=14)
mc.out.dm = matrix(nrow=M,ncol=14)
count.conv = 0
while(count.conv<M){
  alpha=alpha.true
  print(count.conv)
  # b0.init = b0+((-1)^rbinom(1,1,0.5))*runif(1,0,0.2)*abs(b0)
  # b1.init = b1+((-1)^rbinom(1,1,0.5))*runif(1,0,0.2)*abs(b1)
  # b2.init = b2+((-1)^rbinom(1,1,0.5))*runif(1,0,0.2)*abs(b2)
  # g1.init = g1+((-1)^rbinom(1,1,0.5))*runif(1,0,0.2)*abs(g1)
  # g2.init = g2+((-1)^rbinom(1,1,0.5))*runif(1,0,0.2)*abs(g2)
  # gam1.init = gam1+((-1)^rbinom(1,1,0.5))*runif(1,0,0.2)*abs(gam1)
  b0.init = b0+((-1)^rbinom(1,1,0.5))*runif(1,0.5,0.75)*abs(b0)
  b1.init = b1+((-1)^rbinom(1,1,0.5))*runif(1,0.5,0.75)*abs(b1)
  b2.init = b2+((-1)^rbinom(1,1,0.5))*runif(1,0.5,0.75)*abs(b2)
  g1.init = g1+((-1)^rbinom(1,1,0.5))*runif(1,0.5,0.75)*abs(g1)
  g2.init = g2+((-1)^rbinom(1,1,0.5))*runif(1,0.5,0.75)*abs(g2)
  gam1.init = gam1+((-1)^rbinom(1,1,0.5))*runif(1,0.5,0.75)*abs(gam1)
  
  if(alpha==0){
    mydata = data_0_BC(n=200,b0,b1,b2,g1,g2,gam1,cenrate) 
  }
  else{
    mydata = data_gen_BC(n=200,alpha,b0,b1,b2,g1,g2,gam1,cenrate) 
  }
  grid.alpha = seq(0,1,by=0.05)
  grid.out.em = matrix(NA,nrow=length(grid.alpha),ncol=14)
  grid.out.dm = matrix(NA,nrow=length(grid.alpha),ncol=14)
  
  data_obs = mydata[mydata$d==1,]
  data_cens = mydata[mydata$d==0,]
  L_obs = data_obs$L
  R_obs = data_obs$R
  L_cens = data_cens$L
  x1t = data_obs$x1
  x1c = data_cens$x1
  x2t = data_obs$x2
  x2c = data_cens$x2
  x1 = c(x1t,x1c)
  x2 = c(x2t,x2c)
  
  for(i in 1:length(grid.alpha)){

    if(grid.alpha[i]==0){
      alpha=0
      grid.out.em[i,] = BC_0_EM_Wei(mydata=mydata,tol,maxit,increment=0,b0=b0.init,b1=b1.init,
                                 b2=b2.init,g1=g1.init,g2=g2.init,gam1=gam1.init)
      dm.out = optim(par=c(b0.init,b1.init,b2.init,g1.init,g2.init,gam1.init),std0,method="Nelder-Mead")
      hessmat1 = hessian(std0,dm.out$par,method="Richardson")
      FI = solve(1*hessmat1)
      
        std1 = sqrt(FI[1,1])
        std2 = sqrt(FI[2,2])
        std3 = sqrt(FI[3,3])
        std4 = sqrt(FI[4,4])
        std5 = sqrt(FI[5,5])
        std6 = sqrt(FI[6,6])
        
        grid.out.dm[i,] = c(dm.out$par,alpha,std1,std2,std3,std4,std5,std6,-dm.out$value)
    }
    else{
      alpha=grid.alpha[i]
      grid.out.em[i,] = BC_gen_EM_Wei(mydata=mydata,tol,maxit,increment=0,grid.alpha[i],b0,b1,b2,g1,g2,gam1)
      dm.out = optim(par=c(b0.init,b1.init,b2.init,g1.init,g2.init,gam1.init),stdgen,method="Nelder-Mead")
      hessmat1 = hessian(stdgen,dm.out$par,method="Richardson")
      FI = solve(1*hessmat1)
      
      std1 = sqrt(FI[1,1])
      std2 = sqrt(FI[2,2])
      std3 = sqrt(FI[3,3])
      std4 = sqrt(FI[4,4])
      std5 = sqrt(FI[5,5])
      std6 = sqrt(FI[6,6])
      
      grid.out.dm[i,] = c(dm.out$par,alpha,std1,std2,std3,std4,std5,std6,-dm.out$value)
    }
  }
  alpha=alpha.true
  index.em=which.max(grid.out.em[,14])
  index.dm = which.max(grid.out.dm[,14])
  
  out.em = grid.out.em[index.em,]
  out.dm = grid.out.dm[index.dm,]
  print(out.em)
  print(out.dm)
  if(0 %in% out.em[1:6]==FALSE & any(is.na(out.dm))==FALSE){
    count.conv = count.conv +1
    mc.out.em[count.conv,]=out.em
    mc.out.dm[count.conv,]=out.dm
  }
  
}
alpha=1

em.count.b0 = 0
em.count.b1 = 0
em.count.b2 = 0
em.count.g1 = 0
em.count.g2 = 0
em.count.gam1 = 0
em.count.alpha = 0


dm.count.b0 = 0
dm.count.b1 = 0
dm.count.b2 = 0
dm.count.g1 = 0
dm.count.g2 = 0
dm.count.gam1 = 0
dm.count.alpha = 0


for(i in 1:M){
  #95 % Coverage probability of estimates if algorithm converged
  L.em = mc.out.em[i,1:6] - (qnorm(0.025,lower.tail=FALSE)*mc.out.em[i,8:13])
  L.dm = mc.out.dm[i,1:6] - (qnorm(0.025,lower.tail=FALSE)*mc.out.dm[i,8:13])
  U.em = mc.out.em[i,1:6] + (qnorm(0.025,lower.tail=FALSE)*mc.out.em[i,8:13])
  U.dm = mc.out.dm[i,1:6] + (qnorm(0.025,lower.tail=FALSE)*mc.out.dm[i,8:13])
  
  
  
  if((b0 > L.em[1]) & (b0 < U.em[1])){
    em.count.b0 = em.count.b0+1}
  
  if((b1 > L.em[2]) & (b1 < U.em[2])){
    em.count.b1 = em.count.b1+1}
  
  if((b2 > L.em[3]) & (b2 < U.em[3])){
    em.count.b2 = em.count.b2+1}
  
  if((g1 > L.em[4]) & (g1 < U.em[4])){
    em.count.g1 = em.count.g1+1}
  
  if((g2 > L.em[5]) & (g2 < U.em[5])){
    em.count.g2 = em.count.g2+1}
  
  if((gam1 > L.em[6]) & (gam1 < U.em[6])){
    em.count.gam1 = em.count.gam1+1}
  

  if((b0 > L.dm[1]) & (b0 < U.dm[1])){
    dm.count.b0 = dm.count.b0+1}
  
  if((b1 > L.dm[2]) & (b1 < U.dm[2])){
    dm.count.b1 = dm.count.b1+1}
  
  if((b2 > L.dm[3]) & (b2 < U.dm[3])){
    dm.count.b2 = dm.count.b2+1}
  
  if((g1 > L.dm[4]) & (g1 < U.dm[4])){
    dm.count.g1 = dm.count.g1+1}
  
  if((g2 > L.dm[5]) & (g2 < U.dm[5])){
    dm.count.g2 = dm.count.g2+1}
  
  if((gam1 > L.dm[6]) & (gam1 < U.dm[6])){
    dm.count.gam1 = dm.count.gam1+1}
  
  

}


em.avgest.b0 = sum(mc.out.em[,1])/M
em.avgest.b1 = sum(mc.out.em[,2])/M
em.avgest.b2 = sum(mc.out.em[,3])/M
em.avgest.g1 = sum(mc.out.em[,4])/M
em.avgest.g2 = sum(mc.out.em[,5])/M
em.avgest.gam1 = sum(mc.out.em[,6])/M
em.avgest.alpha = sum(mc.out.em[,7])/M

em.avgbias.b0 = sum(mc.out.em[,1]-b0)/M
em.avgbias.b1 = sum(mc.out.em[,2]-b1)/M
em.avgbias.b2 = sum(mc.out.em[,3]-b2)/M
em.avgbias.g1 = sum(mc.out.em[,4]-g1)/M
em.avgbias.g2 = sum(mc.out.em[,5]-g2)/M
em.avgbias.gam1 = sum(mc.out.em[,6]-gam1)/M
em.avgbias.alpha = sum(mc.out.em[,7]-alpha)/M

em.avgse.b0 = sum(mc.out.em[,8])/M
em.avgse.b1 = sum(mc.out.em[,9])/M
em.avgse.b2 = sum(mc.out.em[,10])/M
em.avgse.g1 = sum(mc.out.em[,11])/M
em.avgse.g2 = sum(mc.out.em[,12])/M
em.avgse.gam1 = sum(mc.out.em[,13])/M

em.avgrmse.b0 = sqrt(var(mc.out.em[,1]) + (em.avgbias.b0^2))
em.avgrmse.b1 = sqrt(var(mc.out.em[,2]) + (em.avgbias.b1^2))
em.avgrmse.b2 = sqrt(var(mc.out.em[,3]) + (em.avgbias.b2^2))
em.avgrmse.g1 = sqrt(var(mc.out.em[,4]) + (em.avgbias.g1^2))
em.avgrmse.g2 = sqrt(var(mc.out.em[,5]) + (em.avgbias.g2^2))
em.avgrmse.gam1 = sqrt(var(mc.out.em[,6]) + (em.avgbias.gam1^2))
em.avgrmse.alpha = sqrt(var(mc.out.em[,7]) + (em.avgbias.alpha^2))



dm.avgest.b0 = sum(mc.out.dm[,1])/M
dm.avgest.b1 = sum(mc.out.dm[,2])/M
dm.avgest.b2 = sum(mc.out.dm[,3])/M
dm.avgest.g1 = sum(mc.out.dm[,4])/M
dm.avgest.g2 = sum(mc.out.dm[,5])/M
dm.avgest.gam1 = sum(mc.out.dm[,6])/M
dm.avgest.alpha = sum(mc.out.dm[,7])/M

dm.avgbias.b0 = sum(mc.out.dm[,1]-b0)/M
dm.avgbias.b1 = sum(mc.out.dm[,2]-b1)/M
dm.avgbias.b2 = sum(mc.out.dm[,3]-b2)/M
dm.avgbias.g1 = sum(mc.out.dm[,4]-g1)/M
dm.avgbias.g2 = sum(mc.out.dm[,5]-g2)/M
dm.avgbias.gam1 = sum(mc.out.dm[,6]-gam1)/M
dm.avgbias.alpha = sum(mc.out.dm[,7]-alpha)/M

dm.avgse.b0 = sum(mc.out.dm[,8])/M
dm.avgse.b1 = sum(mc.out.dm[,9])/M
dm.avgse.b2 = sum(mc.out.dm[,10])/M
dm.avgse.g1 = sum(mc.out.dm[,11])/M
dm.avgse.g2 = sum(mc.out.dm[,12])/M
dm.avgse.gam1 = sum(mc.out.dm[,13])/M

dm.avgrmse.b0 = sqrt(var(mc.out.dm[,1]) + (dm.avgbias.b0^2))
dm.avgrmse.b1 = sqrt(var(mc.out.dm[,2]) + (dm.avgbias.b1^2))
dm.avgrmse.b2 = sqrt(var(mc.out.dm[,3]) + (dm.avgbias.b2^2))
dm.avgrmse.g1 = sqrt(var(mc.out.dm[,4]) + (dm.avgbias.g1^2))
dm.avgrmse.g2 = sqrt(var(mc.out.dm[,5]) + (dm.avgbias.g2^2))
dm.avgrmse.gam1 = sqrt(var(mc.out.dm[,6]) + (dm.avgbias.gam1^2))
dm.avgrmse.alpha = sqrt(var(mc.out.dm[,7]) + (dm.avgbias.alpha^2))

out = data.frame(Parameter=c("b0","b1","b2","G1","G2","GAM1","Alpha"),TV=c(b0,b1,b2,g1,g2,gam1,alpha),
                 MLE.em=round(c(em.avgest.b0,em.avgest.b1,em.avgest.b2,em.avgest.g1,em.avgest.g2,em.avgest.gam1,em.avgest.alpha),3),
                 SE.em=round(c(em.avgse.b0,em.avgse.b1, em.avgse.b2,em.avgse.g1,em.avgse.g2,em.avgse.gam1,NA),3),
                 Bias.em=round(c(em.avgbias.b0,em.avgbias.b1,em.avgbias.b2,em.avgbias.g1,
                 em.avgbias.g2,em.avgbias.gam1,em.avgbias.alpha),3), RMSE=round(c(em.avgrmse.b0,em.avgrmse.b1,
                 em.avgrmse.b2,em.avgrmse.g1,em.avgrmse.g2,em.avgrmse.gam1,em.avgrmse.alpha),3),
                 CP.em=c(em.count.b0/M,em.count.b1/M,em.count.b2/M,em.count.g1/M,em.count.g2/M,em.count.gam1/M,NA),
                 MLE.dm=round(c(dm.avgest.b0,dm.avgest.b1,dm.avgest.b2,dm.avgest.g1,dm.avgest.g2,dm.avgest.gam1,dm.avgest.alpha),3),
                 SE.dm=round(c(dm.avgse.b0,dm.avgse.b1, dm.avgse.b2,dm.avgse.g1,dm.avgse.g2,dm.avgse.gam1,NA),3),
                 Bias.dm=round(c(dm.avgbias.b0,dm.avgbias.b1,dm.avgbias.b2,dm.avgbias.g1,
                 dm.avgbias.g2,dm.avgbias.gam1,dm.avgbias.alpha),3), RMSE.dm=round(c(dm.avgrmse.b0,dm.avgrmse.b1,
                 dm.avgrmse.b2,dm.avgrmse.g1,dm.avgrmse.g2,dm.avgrmse.gam1,dm.avgrmse.alpha),3),
                 CP.dm=c(dm.count.b0/M,dm.count.b1/M,dm.count.b2/M,dm.count.g1/M,dm.count.g2/M,dm.count.gam1/M,NA))

kable(out, caption="EM vs DM Par Setting 2, N=200, alpha=0.75, dev 50-75")

#comparing log like
count=0
for(i in 1:M){
  if(mc.out.dm[i,14] < mc.out.em[i,14]){
    count=count+1
  }
}
out=data.frame(text=c("EM","DM"),avg=c(mean(mc.out.em[,14]),mean(mc.out.dm[,14])),prop=c(count/M, 1-count/M))
kable(out, caption="EM vs DM Max Log Like Par Setting 2, N=200, alpha=0.75, dev 50-75")


