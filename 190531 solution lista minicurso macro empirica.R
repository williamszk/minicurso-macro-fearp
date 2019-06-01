#190531 solution lista minicurso macro empirica

library(optimx)
library(xts)

#search for optmin 


SMM_apply = 
  function(par){
  
  rho = par[1]
  beta = par[2]
  alpha = par[3]
  Tp = 50
  S = 1000
    
  hold_mu1 = double() #create holder vectors
  hold_mu2 = double() #
  hold_mu3 = double() #
  
  for (ii in 1:S) {
    
    set.seed(ii)
    e = rnorm(Tp,0,0.1)
    
    z = double()
    z[1] = 0
    for (t in 2:Tp) {
      z[t] = rho*z[t-1] + e[t]
    }
    
    K = double()
    K[1]=1
    for (t in 2:Tp) {
      K[t] = alpha*beta*exp(z[t-1])*K[t-1]^alpha
    }
    
    Y = double()
    Y = exp(z)*K^alpha
    
    C = double()
    C = Y - K[-1]
    
    rt_cons = diff(C)/C[-Tp]
    mu1 = mean(rt_cons[-1]) #eliminate the first observation because is infinite
    
    rt_out = diff(Y)/Y[-Tp]
    mu2 = mean(rt_out[-1])
    
    ratio_con_outp = C/Y
    mu3 = mean(ratio_con_outp); mu3
    
    hold_mu1[ii] = mu1
    hold_mu2[ii] = mu2
    hold_mu3[ii] = mu3
  }
  
  matrix1 = matrix(c(hold_mu1,hold_mu2,hold_mu3),ncol=3)
  
  W1 = cov(matrix1,matrix1)
  
  m_mu1 = mean(hold_mu1)
  m_mu2 = mean(hold_mu1)
  m_mu3 = mean(hold_mu1)
  
  theta1 = c(0.055-m_mu1, 0.0489-m_mu2, 0.9165-m_mu3)
  
  loss = t(theta1) %*% W1 %*% theta1
  as.numeric(loss)
}


optimx(par, SMM_apply , gr=NULL, hess=NULL, lower=-Inf, upper=Inf,
       method=c("Nelder-Mead","BFGS"))


#use the function

Tp = 50 #time periods of the data
S = 10 #number of simulations

rho = .55427
beta = 1.149524
alpha = 0.948
par = c(rho,beta,alpha)
SMM_apply(par)





#import data consuption brasil ----------------------------------------------

data_cons <- 
  read.csv2("C:/Users/willi/Desktop/working/RAW_DATA/IPEA DATA/consumo 190601 bens de consumo/ipeadata[01-06-2019-12-22] consumo bens de consumo.csv", 
            stringsAsFactors=FALSE)
dim(data_cons)
data_cons = data_cons[,-4] #delete the last column
names(data_cons) #see strange names
#indice de consumo
#indice de consumo desazonalisado
#ambos consumo de bens de consumo
names(data_cons) = c('date','cons','cons_deseason')
apply(data_cons,2,class) #all characters

date = data_cons$date
class(date)
date2 = paste(substr(date,1,4),substr(date,6,7),'01',sep='-') 
date3 = as.Date(date2,format = '%Y-%m-%d')
date = date3

cons = data_cons$cons
cons = as.numeric(cons)
class(cons)

cons_de = data_cons$cons_deseason
cons_de = as.numeric(cons_de)
class(cons_de)

cons_xts = xts(matrix(c(cons,cons_de),ncol=2),order.by = date)
class(cons_xts)
plot(cons_xts)
dim(cons_xts)
#290 time points

mean(diff(cons)/cons[-290]) #0.00424719
#very different from 0.055

aux_val1 = as.numeric(cons_xts['20120101'][1,1]) 
which(cons==aux_val1) #position of 20120101

#build the simulated consumption series ------------------------------------------------

rho = 0.55    #values found by SMM
beta = 1.33
alpha = 0.8
Tp = 290 #same time periods of the real database
S = 100

hold_mu1 = double() #create holder vectors
hold_mu2 = double() #
hold_mu3 = double() #

C_m = rep(0,Tp) #create series of consumption

for (ii in 1:S) {
  
  #set.seed(ii)
  e = rnorm(Tp,0,0.1)
  
  z = double()
  z[1] = 0
  for (t in 2:Tp) {
    z[t] = rho*z[t-1] + e[t]
  }
  
  K = double()
  K[1]=1
  for (t in 2:Tp) {
    K[t] = alpha*beta*exp(z[t-1])*K[t-1]^alpha
  }
  
  Y = double()
  Y = exp(z)*K^alpha
  
  C = double()
  C = Y - K
  
  C_m = C_m + C
  C_m = C_m/S
}

aux_val1 = C_m[205] #position 20120101
C_m2 = C_m*(100/aux_val1)
C_m = C_m2

plot(C_m,type='l')

C_xts1 = xts(C_m, order.by = date)

macro_xts1 = cbind(cons_xts,C_xts1)
class(macro_xts1)
plot(macro_xts1)









