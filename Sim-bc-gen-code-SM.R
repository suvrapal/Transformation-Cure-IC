
library(numDeriv)
library(tictoc)
library(knitr)

LR_int=function(y1,len1,l1){
  if(y1>0 & y1<=l1){
    a = c(.Machine$double.eps,l1)
  }else{
    k = as.integer((y1-l1)/len1)+1
    a = c(l1+((k-1)*len1),l1+(k*len1))
  }
  return(a)
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
  
  # print("count.cure")
  # print(count.cure)
  # print("count.obs")
  # print(count.obs)
  
  return(data.frame(L,R,d,x1,x2))
  
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
  
  pnew = rep(0,7)
  pold = rep(0,7)
  
  pold[1]=b0-runif(1,0.5,0.75)*abs(b0)
  pold[2]=b1-runif(1,0.5,0.75)*abs(b1)
  pold[3]=b2-runif(1,0.5,0.75)*abs(b2)
  pold[4]=g1-runif(1,0.5,0.75)*abs(g1)
  pold[5]=g2-runif(1,0.5,0.75)*abs(g2)
  pold[6]=gam1-runif(1,0.5,0.75)*abs(gam1)
  if(alpha==0){
    pold[7]=abs(runif(1,0.5,0.75))  
  }else{pold[7]=min(alpha-runif(1,0.5,0.75)*abs(alpha),1)}
  
  print(pold)
  
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
    phi.cens.r = exp(pold[1]+pold[2]*x1c+pold[3]*x2c)/(1+pold[7]*exp(pold[1]+pold[2]*x1c+pold[3]*x2c))
    F.L.cens.r = 1-S0.L.cens^(exp(pold[4]*x1c+pold[5]*x2c))
    Sp.L.cens.r = (1-pold[7]*phi.cens.r*(F.L.cens.r))^(1/pold[7])
    p0.r= (1-pold[7]*phi.cens.r)^(1/pold[7])  
    Su.r = (Sp.L.cens.r-p0.r)/(1-p0.r)
    e.W = (Sp.L.cens.r-p0.r)/Sp.L.cens.r
    
    
    #E step
    Q = function(par=c(b0,b1,b2,g1,g2,gam1,alpha)){
      S0.L.obs.est = exp(-L_obs^(1/par[6]))
      S0.R.obs.est = exp(-R_obs^(1/par[6]))
      S0.L.cens.est = exp(-L_cens^(1/par[6]))
      phi.obs.est = (exp(par[1]+par[2]*x1t+par[3]*x2t))/(1+par[7]*exp(par[1]+par[2]*x1t+par[3]*x2t))
      phi.cens.est = (exp(par[1]+par[2]*x1c+par[3]*x2c))/(1+par[7]*exp(par[1]+par[2]*x1c+par[3]*x2c))
      F.L.obs.est = 1-S0.L.obs.est^(exp(par[4]*x1t+par[5]*x2t))
      F.R.obs.est = 1-S0.R.obs.est^(exp(par[4]*x1t+par[5]*x2t))
      F.L.cens.est = 1-S0.L.cens.est^(exp(par[4]*x1c+par[5]*x2c))
      Sp.L.obs.est = (1-par[7]*phi.obs.est*F.L.obs.est)^(1/par[7])
      Sp.R.obs.est = (1-par[7]*phi.obs.est*F.R.obs.est)^(1/par[7])
      Sp.L.cens.est = (1-par[7]*phi.cens.est*F.L.cens.est)^(1/par[7])
      p0.est = (1-par[7]*phi.cens.est)^(1/par[7]) #only needed in Q for terms where d=0 (censored obs)
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
    #   pnew = c(0,0,0,0,0,0,0)
    #   return(pnew)
    # }
    # )
    pnew = optim(par=c(pold),fn=Q, method="Nelder-Mead",control=list(maxit=100))$par
    #if error in estimation, stop checking for convergence and return zeros
    if(0 %in% pnew){ 
      
      continue = FALSE
      
      pnew = c(0,0,0,0,0,0,0)
      
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
    
    std = function(par){
      
      b0 = par[1]
      b1 = par[2]
      b2 = par[3]
      g1 = par[4]
      g2 = par[5]
      gam1 = par[6]
      alpha = par[7]
      
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
    
    hessmat1 = hessian(std,c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6],pnew[7]),method="Richardson")
    FI = solve(-1*hessmat1)
    
    if(any(is.na(FI))){
      print(FI)
      out = c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6],pnew[7],0,0,0,0,0,0,0,iter)
    }
    else if(FI[1,1]<0 | FI[2,2]<0 | FI[3,3]<0| FI[4,4]<0 | FI[5,5]<0 | FI[6,6]<0 | FI[7,7]<0){ #added condition to handle if FI non-invertible
      out = c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6],pnew[7],0,0,0,0,0,0,0,iter)
    }
    else{
      std1 = sqrt(FI[1,1])
      std2 = sqrt(FI[2,2])
      std3 = sqrt(FI[3,3])
      std4 = sqrt(FI[4,4])
      std5 = sqrt(FI[5,5])
      std6 = sqrt(FI[6,6])
      std7 = sqrt(FI[7,7])
      
      out = c(pnew[1],pnew[2],pnew[3],pnew[4],pnew[5],pnew[6],pnew[7],std1,std2,std3,std4,std5,std6,std7,iter)
    }
  }#end of else
  return(out)
} 


std = function(par,mydata){
  
  b0 = par[1]
  b1 = par[2]
  b2 = par[3]
  g1 = par[4]
  g2 = par[5]
  gam1 = par[6]
  alpha = par[7]
  
  obs_data=mydata[mydata$d==1,]
  cens_data=mydata[mydata$d==0,]
  L_obs=obs_data$L
  R_obs=obs_data$R
  L_cens=cens_data$L
  x1t=obs_data$x1
  x1c=cens_data$x1
  x2t=obs_data$x2
  x2c=cens_data$x2
  
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

MC.Sim.log.like = function(n,b0,b1,b2,g1,g2,gam1,alpha,M,increment,cenrate){
  #defining null vectors of interest
  est.b0 = rep(NA,M)
  est.b1 = rep(NA,M)
  est.b2 = rep(NA,M)
  est.g1 = rep(NA,M)
  est.g2 = rep(NA,M)
  est.gam1 = rep(NA,M)
  est.alpha = rep(NA,M)
  count.div=0
  continue=TRUE
  count.conv = 0
  
  for(i in 1:M){ 
    print(i)
    
      b0.init = b0-runif(1,0.5,0.75)*abs(b0)
      b1.init = b1-runif(1,0.5,0.75)*abs(b1)
      b2.init = b2-runif(1,0.5,0.75)*abs(b2)
      g1.init = g1-runif(1,0.5,0.75)*abs(g1)
      g2.init = g2-runif(1,0.5,0.75)*abs(g2)
      gam1.init = gam1-runif(1,0.5,0.75)*abs(gam1)
      if(alpha==0){
        alpha.init = abs(runif(1,0.5,0.75))  
      }else{alpha.init = min(alpha-runif(1,0.5,0.75)*abs(alpha),1)}
      
    while(continue){ #to ensure we reach M successful MC runs
      mydata = data_gen_BC(n=200,alpha,b0,b1,b2,g1,g2,gam1,cenrate) 
      # result =optim(par=c(b0.init,b1.init,b2.init,g1.init,
      #                     g2.init,gam1.init,alpha.init),fn=std,mydata=mydata, 
      #               method="Nelder-Mead",control = list(maxit = 10000))
      
      result = tryCatch({optim(par=c(b0.init,b1.init,b2.init,g1.init,
                            g2.init,gam1.init,alpha.init),fn=std,mydata=mydata, 
                      method="Nelder-Mead",control = list(maxit = 10000))
        },error=function(e){
          result = c(0,0,0,0,0,0,0)
          return(result)
        }
        )
      print(result)
      if(result[1]==0){
        count.div=count.div+1
      }
      else if(result$convergence !=0 | result$par[6]<0 | result$par[7]<0 ){ 
        #g1, g2 Weibull parameters must be positive
        count.div=count.div+1
      }
      else{
        est.b0[i]=result$par[1]
        est.b1[i]=result$par[2]
        est.b2[i]=result$par[3]
        est.g1[i]=result$par[4]
        est.g2[i]=result$par[5]
        est.gam1[i]=result$par[6]
        est.alpha[i]=result$par[7]
        continue=FALSE
      }
    }#end of while loop
  }#end of for loop
  
  #extracting convergent estimates
  est.conv = rbind(est.b0,est.b1,est.b2,est.g1,est.g2,est.gam1,est.alpha)
  count.conv=ncol(est.conv)
  print(ncol(est.conv))
  print(count.div)

  
  #quantiles for all estimates
  Q = matrix(quantile(est.conv[1,]),byrow=TRUE,nrow=1)
  Q = rbind(Q,quantile(est.conv[2,]),quantile(est.conv[3,]),
            quantile(est.conv[4,]),quantile(est.conv[5,]),quantile(est.conv[6,]),quantile(est.conv[7,]))
  
  iqr.cov=rep(NA,7)
  
  for(i in 1:7){#counting number of times estimate was NOT in +/- 1.5 IQR range
    count.iqr.cov=0
    IQR = Q[i,4]-Q[i,2]
    for(j in 1:count.conv){
      if(est.conv[i,j]<(Q[i,2]- 1.5*IQR) | est.conv[i,j]>(Q[i,4]+ 1.5*IQR)){
        count.iqr.cov=count.iqr.cov+1
      }
      iqr.cov[i]=count.iqr.cov
    }
  }
  #returning median estimates and true coverage proportions
  
  out = data.frame(Parameter=c("B0","B1","B2","G1","G2","GAM1","Alpha"),TV=c(b0,b1,b2,g1,g2,gam1,alpha),
                   Med=c(Q[1,3],Q[2,3],Q[3,3],Q[4,3],Q[5,3],Q[6,3],Q[7,3]),
                   IQR.cov=c((iqr.cov[1]/count.conv)*100,(iqr.cov[2]/count.conv)*100,
                              (iqr.cov[3]/count.conv)*100,(iqr.cov[4]/count.conv)*100,
                              (iqr.cov[5]/count.conv)*100,(iqr.cov[6]/count.conv)*100,(iqr.cov[7]/count.conv)*100))
  return(out)
}#end of the main function

MC.Sim.BC.EM.IC = function(n,b0,b1,b2,g1,g2,gam1,alpha,M,maxit,tol,increment){
  #defining null vectors of interest
  est.b0 = rep(NA,M)
  est.b1 = rep(NA,M)
  est.b2 = rep(NA,M)
  est.g1 = rep(NA,M)
  est.g2 = rep(NA,M)
  est.gam1 = rep(NA,M)
  est.alpha = rep(NA,M)
  count.div=0
  
  for(i in 1:M){ 
    print(i)
    mydata = data_gen_BC(n=200,alpha,b0,b1,b2,g1,g2,gam1,cenrate) 
    result = BC_gen_EM_Wei(mydata=mydata,tol,maxit,increment,alpha,b0,b1,b2,g1,g2,gam1)
    if(result[1]==0){
      #set all estimates to zero
      est.b0[i] = 0
      est.b1[i] = 0
      est.b2[i] = 0
      est.g1[i] = 0
      est.g2[i] = 0
      est.gam1[i] = 0
      est.alpha[i] = 0
      
      #count number of divergent MC runs
      count.div = count.div + 1
    }
    
    else{
      est.b0[i] = result[1]
      est.b1[i] = result[2]
      est.b2[i] = result[3]
      est.g1[i] = result[4]
      est.g2[i] = result[5]
      est.gam1[i] = result[6]
      est.alpha[i] = result[7]
      print(result)
    }
  }#end of for loop
  
  
  #extracting convergent estimates
  est.conv = matrix(est.b0[est.b0 != 0],byrow=TRUE,nrow=1)
  est.conv = rbind(est.conv,est.b1[est.b1 != 0],est.b2[est.b2 != 0],
                   est.g1[est.g1 != 0],est.g2[est.g2 != 0],
                   est.gam1[est.gam1 != 0],est.alpha[est.alpha != 0])
  
  #quantiles for all estimates
  Q = matrix(quantile(est.conv[1,]),byrow=TRUE,nrow=1)
  Q = rbind(Q,quantile(est.conv[2,]),quantile(est.conv[3,]),quantile(est.conv[4,]),
            quantile(est.conv[5,]),quantile(est.conv[6,]),quantile(est.conv[7,]))
  
  iqr.cov=rep(NA,7)
  
  for(i in 1:7){#counting number of times estimate was NOT in +/- 1.5 IQR range
    count.iqr.cov=0
    IQR = Q[i,4]-Q[i,2]
    for(j in 1:(M-count.div)){
      if(est.conv[i,j]<(Q[i,2]- 1.5*IQR) | est.conv[i,j]>(Q[i,4]+ 1.5*IQR)){
        count.iqr.cov=count.iqr.cov+1
      }
      iqr.cov[i]=count.iqr.cov
    }
  }
  #returning median estimates and true coverage proportions
  
  out = data.frame(Parameter=c("B0","B1","B2","G1","G2","GAM1","Alpha"),TV=c(b0,b1,b2,g1,g2,gam1,alpha),
                   Med=c(Q[1,3],Q[2,3],Q[3,3],Q[4,3],Q[5,3],Q[6,3],Q[7,3]),
                   IQR.cov=c(iqr.cov[1]/(M-count.div),(iqr.cov[2]/(M-count.div)),
                              iqr.cov[3]/(M-count.div),iqr.cov[4]/(M-count.div),
                              iqr.cov[5]/(M-count.div),iqr.cov[6]/(M-count.div),iqr.cov[7]/(M-count.div))*100)
  print(count.div)
  return(out)
}#end of the main function
cenrate=0.2


#simulation part 2 - robustness study
tic("N=200, M=500, First parameter setting. 
    Robustness Setting 1 EM - Initial values= true values")
est.1a = MC.Sim.BC.EM.IC(n=200,b0=0.2,b1=-1.4,b2=0.1,g1=-1.2,g2=0.05,gam1=0.215,alpha=0.5,
                        M=500,tol=10^-4,maxit=1000,increment=0)
toc()


tic("N=200, M=500, First parameter setting.
    Robustness Setting 1 Log Like - Initial value= true value")
est.4a=MC.Sim.log.like(n=200,b0=0.2,b1=-1.4,b2=0.1,g1=-1.2,g2=0.05,gam1=0.215,alpha=0.5,
                       M=500,increment=0,cenrate=0.2)
toc()



tic("N=200, M=500, Second parameter setting. 
    Robustness Setting 1 EM - Initial values= true values")
est.1b = MC.Sim.BC.EM.IC(n=200,b0=0.9,b1=-0.6,b2=-0.1,g1=1,g2=-0.2,gam1=0.316,alpha=0.75,
                        M=500,tol=10^-4,maxit=1000,increment=0)
toc()


tic("N=200, M=500, First parameter setting. 
    Robustness Setting 2 EM - Initial value= true value + 0.5*U(0,1)*true_value")


est.2a = MC.Sim.BC.EM.IC(n=200,b0=0.2,b1=-1.4,b2=0.1,g1=-1.2,g2=0.05,gam1=0.215,alpha=0.5,
                         M=500,tol=10^-4,maxit=1000,increment=0.5)
toc()


tic("N=200, M=500, Second parameter setting. 
    Robustness Setting 2 EM - Initial value= true value + 0.5*U(0,1)*true_value")
est.2b = MC.Sim.BC.EM.IC(n=200,b0=0.9,b1=-0.6,b2=-0.1,g1=1,g2=-0.2,gam1=0.316,alpha=0.75,
                         M=500,tol=10^-4,maxit=1000,increment=0.5)
toc()


tic("N=200, M=500, First parameter setting. 
    Robustness Setting 3 EM - Initial value= true value - 0.5*U(0,1)*true_value")
est.3a = MC.Sim.BC.EM.IC(n=200,b0=0.2,b1=-1.4,b2=0.1,g1=-1.2,g2=0.05,gam1=0.215,alpha=0.5,
                         M=500,tol=10^-4,maxit=1000,increment= -0.5)
toc()


tic("N=200, M=500, Second parameter setting. 
    Robustness Setting 3 EM - Initial value= true value - 0.5*U(0,1)*true_value")
est.3b = MC.Sim.BC.EM.IC(n=200,b0=0.9,b1=-0.6,b2=-0.1,g1=1,g2=-0.2,gam1=0.316,alpha=0.75,
                         M=500,tol=10^-4,maxit=1000,increment=-0.5)
toc()


tic("N=200, M=500, Second parameter setting.
    Robustness Setting 1 Log Like - Initial value= true value")
est.4b=MC.Sim.log.like(n=200,b0=0.9,b1=-0.6,b2=-0.1,g1=1,g2=-0.2,gam1=0.316,alpha=0.75,
                       M=500,increment=0,cenrate=0.2)
toc()


tic("N=200, M=500, First parameter setting.
    Robustness Setting 2 Log Like - Initial value= true value + 0.5*U(0,1)*true_value")
est.5a=MC.Sim.log.like(n=200,b0=0.2,b1=-1.4,b2=0.1,g1=-1.2,g2=0.05,gam1=0.215,alpha=0.5,
                       M=500,increment=0.5,cenrate=0.2)
toc()


tic("N=200, M=500, Second parameter setting.
    Robustness Setting 2 Log Like - Initial value= true value + 0.5*U(0,1)*true_value")
est.5b=MC.Sim.log.like(n=200,b0=0.9,b1=-0.6,b2=-0.1,g1=1,g2=-0.2,gam1=0.316,alpha=0.75,
                       M=500,increment=0.5,cenrate=0.2)
toc()


tic("N=200, M=500, First parameter setting.
    Robustness Setting 3 Log Like - Initial value= true value -  0.5*U(0,1)*true_value")
est.6a=MC.Sim.log.like(n=200,b0=0.2,b1=-1.4,b2=0.1,g1=-1.2,g2=0.05,gam1=0.215,alpha=0.5,
                       M=500,increment=-0.5,cenrate=0.2)
toc()


tic("N=100, M=1000, Second parameter setting.
    Robustness Setting 3 Log Like - Initial value= true value -  0.5*U(0,1)*true_value")
est.6b=MC.Sim.log.like(n=200,b0=0.9,b1=-0.6,b2=-0.1,g1=1,g2=-0.2,gam1=0.316,alpha=0.75,
                       M=500,increment=-0.5,cenrate=0.2)
toc()


kable(est.1a, caption="N=200, M=500, Parameter setting 1, Robustness setting 1 EM")
kable(est.4a, caption="N=200, M=500, Parameter setting 1, Robustness setting 1 Log Like")


kable(est.1b, caption="N=200, M=500, Parameter setting 2, Robustness setting 1 EM")
kable(est.2a, caption="N=200, M=500, Parameter setting 1, Robustness setting 2 EM")
kable(est.2b, caption="N=200, M=500, Parameter setting 2, Robustness setting 2 EM")
kable(est.3a, caption="N=200, M=500, Parameter setting 1, Robustness setting 3 EM")
kable(est.3b, caption="N=200, M=500, Parameter setting 2, Robustness setting 3 EM")
kable(est.4b, caption="N=200, M=500, Parameter setting 2, Robustness setting 1 Log Like")
kable(est.5a, caption="N=200, M=500, Parameter setting 1, Robustness setting 2 Log Like")
kable(est.5b, caption="N=200, M=500, Parameter setting 2, Robustness setting 2 Log Like")
kable(est.6a, caption="N=200, M=500, Parameter setting 1, Robustness setting 3 Log Like")
kable(est.6b, caption="N=200, M=500, Parameter setting 2, Robustness setting 3 Log Like")







