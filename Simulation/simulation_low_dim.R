source(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),'/Distributed_Huber.R'))
path = paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/simulation results/')

# simulations: Fixed n, increasing m with homogeneous or linear heterogenous model.

m_list = c(50,100,200,300,400)
n = 400
d = 20
heterogeneity = 'HetQua' # 'Hom', 'HetLin', 'HetQua'
total_runs = 500
robust.consts = c(1.25,2.16,2.64,2.64) 
beta.true = rep(1.5, d-1)
beta0.true = 1.5
dist_list = c('N','T','P','B')
set.seed(2020)
n_estimator = 5
inference = 'normal' # 'normal' or 'boot'
# output index: m, dist, estimator, runs, (beta, lCI, uCI), d
output = array(rep(0, length(m_list)*length(dist_list)*n_estimator*total_runs*3*d), dim = c(length(m_list), length(dist_list), n_estimator, total_runs, 3, d))

# tautable index: dist, runs*length(m_list), (m, tau_est, tau)
tautable = array(rep(0, length(dist_list)*total_runs*length(m_list)*3), dim = c(length(dist_list), total_runs*length(m_list), 3))
dimnames(tautable) = list(dist_list, 1:dim(tautable)[2], c('m','tau.est','tau.tf'))

for (i in 1:length(m_list)){
  m = m_list[i]
  N = n*m
  index = function(c,n) ((c-1)*n+1):(c*n)
  for (j in 1:total_runs){
    # X = matrix(runif(N*(d-1), -1.5,1.5),N,d-1)
    X = matrix(rnorm(N*(d-1)), N, d-1)
    for (k in 1:length(dist_list)){
      if (dist_list[k] == 'N') err = rnorm(N, 0, 1)
      if (dist_list[k] == 'T') err = rt(N,2)
      if (dist_list[k] == 'P') err = Pareto::rPareto(N, 4, 2) - 4*2/(2-1)
      if (dist_list[k] == 'L') err = rlnorm(N, 0, 1.5) - exp(0 + 1.5^2 / 2)
      if (dist_list[k] == 'B') {
        a=2;b=1;s=1  
        err = actuar::rburr(N,a,b,s) - a*beta(a-1/b,1+1/b)
      }
      if (heterogeneity == 'Hom'){
        y = beta0.true + X%*%beta.true + err
      }
      if (heterogeneity == 'HetLin'){
        y = beta0.true + X%*%beta.true + (X[,1]/2+1)*err
      }
      if (heterogeneity == 'HetQua'){
        y = beta0.true + X%*%beta.true + (X%*%beta.true)^2/sqrt(3)/sum(beta.true^2)*err
      }
      
      # 5 estimators global AHR, DC AHR, distributed OLS, dist AHR, DC OLS
      gAHR = Huber_GDBB(y,X)
      output[i,k,1,j,1,] = gAHR$coefficients
      output[i,k,2,j,1,] = rowMeans(sapply(1:m, function(x) Huber_GDBB(y[index(x,n)], X[index(x,n),])$coefficients))
      output[i,k,3,j,1,] = Distributed_OLS(y,X,m)
      dAHR = Distributed_Huber(y, X, m, robust_const = robust.consts[k])
      output[i,k,4,j,1,] = dAHR$coefficients
      tautable[k, j + total_runs*(i-1), ] = c(m, dAHR$tau/robust.consts[k], gAHR$tau)
      output[i,k,5,j,1,] = rowMeans(sapply(1:m, function(x) lm(y[index(x,n)]~as.matrix(X[index(x,n),]))$coefficients))
      
      # CIs
      X.1 = cbind(1, X)
      sigma1 = (sqrt(rowMeans(sapply(1:m, function(x) diag(n*solve(t(X.1[index(x,n),])%*%X.1[index(x,n),])/n))))* sqrt(sum((y - X.1%*%output[i,k,3,j,1,])^2)/m))[-1]
      output[i,k,3,j,2:3,-1] = rbind(output[i,k,3,j,1,][-1] - qnorm(0.975)/sqrt(N)*sigma1, output[i,k,3,j,1,][-1] + qnorm(0.975)/sqrt(N)*sigma1)
      if (inference == 'normal'){
        output[i,k,4,j,2:3,-1] = t(Distributed_Inference_Normal(y, X, m, beta=dAHR$coefficients, robust_const = robust.consts[k], tau = dAHR$tau, general = TRUE))
      }
      if (inference == 'boot'){
        output[i,k,4,j,2:3,-1] = Distributed_Inference_Boot(y, X, m, n.boot = 1000,  dAHR$coefficients, tau = dAHR$tau, robust_const = robust.consts[k], quantiles = c(0.025, 0.975)) 
      }
    }
    if (j %% 10 == 0) cat(paste0(j, ','))
  }
  cat('m =',m,'done.\n')
}


# save all outputs, beware not to overwrite previous outputs

# save(output, m_list, n_estimator, heterogeneity, dist_list, beta.true, beta0.true,robust.consts, total_runs,n, d, tautable, file = paste0(path, 'outputLow', heterogeneity, '.Rdata'))





# constant c tuning
par(mfrow = c(2,2))
for (ld in 1:(dim(tautable)[1])){
  plot(tautable[ld,,2], tautable[ld,,3], xlab = 'tau estimate (before muliply by c)', ylab = 'tau from global AHR',  main = paste0('error = ',dist_list[ld],', c = ', signif(median(tautable[ld,,3]/tautable[ld,,2]),3)))
}

