library(survival)
library(parallel)
library(PWEALL)

set.seed(1000)
nsim = 200
censor_fun <- function(i, time){
  E = runif(nsim, min = 12, max = 60)
  L = c()
  for(i in 1:nsim){
    if( (60 - E[i]) < 24){
      L[i] = rexp(1, rate = 0.04)
    }else{
      L[i] = 1000
    }
  }
  X = rexp(nsim, rate = 0.04)
  dat = data.frame(dens = c(L, E, X), lines = rep(c("L", "E", "X"), each = nsim))
  C = c()
  for(i in 1:nsim){C[i] = min(L[i], E[i])}
  
  delta1 = c()
  Y = c()
  for(i in 1:nsim){
    Y[i] = min(C[i], X[i])
    if(X[i] < C[i]){delta1[i] = 1
    }else{delta1[i] = 0}
  }
  
  surv_obj1 = Surv(time = Y, event = (1-delta1))
  fit1 = survfit(surv_obj1~1)
  KM_c = summary(fit1,times=time,extend=TRUE)
  
  ti_c = time
  out_OKM = KM_c$surv
  OKM = out_OKM
  min_surv = max(OKM[which(OKM <0.5)])
  max_surv = min(OKM[which(OKM >=0.5)])
  
  both = c(min_surv, max_surv)
  diff = both - 0.5
  take_this = 2
  if(take_this == 1){
    closest = which(OKM == both[take_this])
    OKM_med = ti_c[closest[1]]
  }
  if(take_this == 2){
    closest = which(OKM == both[take_this])
    OKM_med = ti_c[closest[length(closest)]]
  }
  
  #===============================================
  
  #===============================================
  
  emp_E = c()
  for(i in 1:length(time)){
    emp_E[i] = length(which(E <= time[i]))/nsim
  }
  
  NKM = c()
  for(t in 1:(length(time))){
    ind = which(E > time[t])
    if(length(ind) > 0){
      Et = E[ind]
      Lt = L[ind]
      Xt = X[ind]
      
      delta = c()
      Y = c()
      for(i in 1:length(ind)){
        interm =  min(Et[i], Xt[i])
        Y[i] = min(Lt[i], interm)
        if(Lt[i] < interm){
          delta[i] = 1
        }else{
          delta[i] = 0
        }
      }
      
      surv_obj = Surv(time = Y, event = delta)
      fit = survfit(surv_obj~1)
      KM = summary(fit,times=time,extend=TRUE)
      ti=time
      interm = which(ti <= time[t])

      NKM[t] = KM$surv[interm[length(interm)]]*(1-emp_E[t])
     
    }
    else{
      NKM[t] = 0
    }
    
  }
  
  out_NKM = NKM
  min_surv = max(NKM[which(NKM <0.5)])
  max_surv = min(NKM[which(NKM >=0.5)])
  
  both = c(min_surv, max_surv)
  diff = both - 0.5
  take_this =2
  if(take_this == 1){
    closest = which(NKM == both[take_this])
    NKM_med = time[closest[1]]
  }
  if(take_this == 2){
    closest = which(NKM == both[take_this])
    NKM_med = time[closest[length(closest)]]
  }

  #=====================================
  #=====================================
  #this is the TRUE surv_object for C
  surv_objT = Surv(time = C, event = rep(1, length(X)))
  fitT = survfit(surv_objT~1)
  KM_T = summary(fitT,times=time,extend=TRUE)
  
  out_T = KM_T$surv
  min_surv = max(out_T[which(out_T  <0.5)])
  max_surv = min(out_T[which(out_T  >=0.5)])
  
  both = c(min_surv, max_surv)
  diff = both - 0.5
  take_this = 2
  if(take_this == 1){
    closest = which(out_T == both[take_this])
    T_med = time[closest[1]]
  }
  if(take_this == 2){
    closest = which(out_T == both[take_this])
    T_med = time[closest[length(closest)]]
    
  }
  
  out = matrix(c(out_NKM, out_OKM, out_T), nrow = length(out_NKM))
  out_quant = c(OKM_med,T_med, NKM_med )
  return(list(out, out_quant))
}

time = seq(0,56, by=0.1)

#odd column is augment, even column is conventional
out = lapply(c(1:5000), censor_fun, time = time)

surv = list()
med = list()
for(i in 1:length(out)){
  surv[[i]] = out[[i]][[1]]
  med[[i]] = out[[i]][[2]]
}

#### median
medi = do.call(rbind, med)
apply(medi, 2, mean, na.rm = T)
#col 1 = augmented, col 2 = conventional
quantile(medi[,1], c(0.025, 0.975))
quantile(medi[,2], c(0.025, 0.975))

