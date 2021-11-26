#EM algorithm for K Gaussians
#EMKgauss function receives 1D data and a value for K
#plotEMK plots the gaussians and clustered points and the loglikelihood.

EMKgauss <- function(data,K, useKM = FALSE){
  #initial parameter values
  #initializes mu using random data from the sample
  if (useKM){
    km = kmeans(data,centers=K)
    mu = km$centers
    print('mu initialized using Kmeans')
  }
  else{
    mu = sample(data,K) 
    print('mum randomly initialized')
  }
  #sigsq as sample variance of the entire sample
  sigsq = rep(var(data),K)
  alpha = rep((1/K),K)
  N = length(data)
  #function to build responsability matrix
  getGamma <- function(data,mu,sigsq,alpha){
    aux = alpha[1]*dnorm(data,mean = mu[1],sd = sqrt(sigsq[1]))
    for (k in 2:K){
      aux = rbind(aux,alpha[k]*dnorm(data,mean = mu[k],sd = sqrt(sigsq[k])))
    }
    den = colSums(aux)
    wMtx = aux[1,]/den
    for (j in 2:K){
      wMtx = rbind(wMtx,aux[j,]/den)
    }
    return(wMtx)
  }
  
  getNk <- function(gamma){
    return(rowSums(gamma))
  }
  
  updateMu <- function(data,gamma,Nk){
    newMu = c()
    for (k in 1:K){
      aux = sum(gamma[k,]*data)
      newMu = append(newMu,aux/Nk[k])
    }
    return(newMu)
  }
  updateSigsq <- function(data,gamma,newMu,Nk){
    newSigsq = c()
    for (k in 1:K){
      aux = gamma[k,]*(data-newMu[k])**2
      newSigsq = append(newSigsq,sum(aux)/Nk[k])
    }
    return(newSigsq)
  }
  
  logVer <- function(data,gamma,alpha){
    aux = 0
    aux2 = 0
    for (k in 1:K){
      aux = aux + log(dnorm(data,mu[k],sqrt(sigsq[k])))*gamma[k,]
      aux2 = aux2 + log(alpha[k])*gamma[k,]
    }
    return(sum(aux+aux2))
  }
  
  
  
  count = 0
  gamma = sapply(data,getGamma,mu=mu,sigsq=sigsq,alpha = alpha)
  lvlist = c(logVer(data,gamma,alpha))
  
  conv = FALSE
  while (!conv) {
    #Expectation step
    gamma = sapply(data,getGamma,mu=mu,sigsq=sigsq,alpha = alpha)
    
    
    #maximization step
    Nkk = getNk(gamma)
    mu = updateMu(data,gamma,Nkk)
    sigsq = updateSigsq(data,gamma,mu,Nkk)
    alpha = Nkk/length(data)
    
    
    #computes loglikelihood
    lvlist = append(lvlist,logVer(data,gamma,alpha))
    print(class(lvlist[length(lvlist)]))
    
    if(length(lvlist)>13){
      condition = abs((lvlist[length(lvlist)]-lvlist[length(lvlist)-10])/lvlist[length(lvlist)-10])
      if(condition<0.0001){
        conv = TRUE
      }
    }
    count = count +1
    print(paste('Iteracion ',count))
  }
  toRet = list('mu'=mu,'sigsq' = sigsq,'gamma'=gamma,'data'=data,'K'=K,'lvlist'=lvlist)
  return(toRet)
}

plotEMK <- function(EMKans){
  #function that plots EMKgauss output--------------
  #---------------------------------------------------------------
  mu = EMKans$mu
  sigsq = EMKans$sigsq
  gamma = EMKans$gamma
  data = EMKans$data
  K = EMKans$K
  lvlist = EMKans$lvlist
  N = length(data)
  
  #loglikelihood plot
  plot(lvlist,type='l',lwd =3,col=4,main = 'Logverosimilitud',
       xlab = 'Iteraciones')
  
  #colour list
  cols = rep(0,N)
  for (i in 1:N){
    cols[i] = which.max(gamma[,i])
  }
  #computes gaussian maximum
  idx = which.min(sigsq)
  maxgauss =  dnorm(mu[idx],mu[idx],sqrt(sigsq[idx]))
  plot(data,rep(0,N),col=cols,ylim = c(0,1.1*maxgauss),
       main = paste('EM algorithm result with K=',K,' gaussians'),
       xlab = 'Sample observations')
  points(data,rep(maxgauss*1.1,N))
  x = seq(min(data),max(data),0.1)
  for (k in 1:K){
    lines(x,dnorm(x,mu[k],sqrt(sigsq[k])),col=k)
  }
}

#Generating sample data ''''''''''''''

#------Data #1
n1 = rnorm(100,-6,3)
n2 = rnorm(100,0,1)
n3 = rnorm(200,8,3)
data1 = c(n1,n2,n3)
ans = EMKgauss(data1,3,useKM = TRUE)
plotEMK(ans)


#------Data #2
# 10 000 weight measurements, half men half women
# obtained forom kaggle https://www.kaggle.com/mustafaali96/weight-height?select=weight-height.csv
df = read.csv('weight-height.csv')
data3 = df$Weight
ans3 = EMKgauss(data3,2,useKM = TRUE)
plotEMK(ans3)


