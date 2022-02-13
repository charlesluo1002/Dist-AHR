source(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),'/Distributed_Huber.R'))
path = paste0(dirname(rstudioapi::getActiveDocumentContext()$path), '/simulation results/')

# simulations: Fixed n, increasing m with homogeneous or linear heterogenous model.

m_list = c(10,20,30,40,50)
d = 1000
s = 5
n = 250
beta.true = c(rep(1.5,s-1),rep(0,d-s))
beta0.true = 1.5
heterogeneity = 'HetQua' # 'Hom', 'HetLin','HetQua'
total_runs = 100
robust.consts = c(1.19, 1.65, 1.87, 1.87)
dist_list = c('N','T','P','B')
validation_proportion = 0.25
set.seed(2020)

# output index: m, dist, estimator, runs,  d
output = array(rep(0, length(m_list)*length(dist_list)*4*total_runs*d), dim = c(length(m_list), length(dist_list), 4, total_runs, d))

# tautable index: dist, runs*length(m_list), (m, tau_est, tau)
tautable = array(rep(0, length(dist_list)*total_runs*length(m_list)*3), dim = c(length(dist_list), total_runs*length(m_list), 3))
dimnames(tautable) = list(dist_list, 1:dim(tautable)[2], c('m','tau.est','tau.tf'))

for (i in 1:length(m_list)){
  m = m_list[i]
  N = n*m
  N_valid = round(N*validation_proportion)
  index = function(c,n) ((c-1)*n+1):(c*n)
  for (j in 1:total_runs){
    # X = matrix(runif(N*(d-1), -1.5,1.5),N,d-1)
    # X_valid = matrix(runif(N_valid*(d-1), -1.5,1.5),N_valid,d-1)
    X = matrix(rnorm(N*(d-1)), N, d-1)
    X_valid = matrix(rnorm(N_valid*(d-1)), N_valid, d-1)
    for (k in 1:length(dist_list)){
      if (dist_list[k] == 'N') err = rnorm(N, 0, 1);err_valid = rnorm(N_valid, 0, 1)
      if (dist_list[k] == 'T') err = rt(N,2);err_valid = rt(N_valid,2)
      if (dist_list[k] == 'P') err = Pareto::rPareto(N, 4, 2) - 4*2/(2-1);err_valid = Pareto::rPareto(N_valid, 4, 2) - 4*2/(2-1)
      if (dist_list[k] == 'L') err = rlnorm(N, 0, 1.5) - exp(0 + 1.5^2 / 2); err_valid = rlnorm(N_valid, 0, 1.5) - exp(0 + 1.5^2 / 2)
      if (dist_list[k] == 'B') {
        a=2;b=1;s=1  
        err = actuar::rburr(N,a,b,s) - a*beta(a-1/b,1+1/b)
        err_valid = actuar::rburr(N_valid,a,b,s) - a*beta(a-1/b,1+1/b)
      }
      if (heterogeneity == 'Hom'){
        y = beta0.true + X%*%beta.true + err
        y_valid = beta0.true + X_valid%*%beta.true + err_valid
      }
      if (heterogeneity == 'HetLin'){
        y = beta0.true + X%*%beta.true + (X[,1]/2+1)*err
        y_valid = beta0.true + X_valid%*%beta.true + (X_valid[,1]/2+1)*err_valid
      }
      
      if (heterogeneity == 'HetQua'){
        y = beta0.true + X%*%beta.true + (X%*%beta.true)^2/sqrt(3)/sum(beta.true^2)*err
        y_valid = beta0.true + X_valid%*%beta.true + (X_valid%*%beta.true)^2/sqrt(3)/sum(beta.true^2)*err_valid
      }
      
      # 4 estimators global AHR, DC AHR, global Lasso, dist AHR
      output[i,k,2,j,] = rowMeans(sapply(1:m, function(x) ILAMM_validation(y[index(x,n)], X[index(x,n),], y_valid, X_valid, nlambda = 10, phi0 = 0.0001)$beta))
      lasso_out = glmnet(X,y,nlambda = 10, family = 'gaussian')
      output[i,k,3,j,] = as.numeric(coef(lasso_out)[,which.min(sapply(1:dim(lasso_out$beta)[2], function(x) sum((y_valid - cbind(1,X_valid)%*%as.numeric(coef(lasso_out)[,x]))^2)))])
      dAHR = Distributed_validation_ILAMM(y, X, y_valid, X_valid, m, nlambda = 10, robust_const = robust.consts[k], phi0 = 0.0001)
      gAHR = ILAMM_validation(y, X, y_valid, X_valid, lambda = dAHR$lambda, nlambda = 10, phi0 = 0.0001)
      output[i,k,1,j,] = gAHR$beta
      output[i,k,4,j,] = dAHR$beta
      tautable[k, j + total_runs*(i-1), ] = c(m, dAHR$tau/robust.consts[k], gAHR$tau)
    }
    if (j %% 5 == 0) cat(paste0(j, ','))
  }
  cat('m =',m,'done.\n')
}

# save all outputs, beware not to overwrite previous outputs

# save(output, m_list, heterogeneity, dist_list, beta.true, beta0.true,robust.consts, total_runs,n, d, s, tautable, validation_proportion, file = paste0(path, 'outputHigh', heterogeneity, '.Rdata'))


