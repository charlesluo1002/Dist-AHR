library(FarmTest)
library(glmnet)
library(ILAMM)
# implement adaptive huber regression using GDBB
Huber_GDBB <- function(y, X = NULL, delta = 1e-5, tau_const = 'adaptive', robust_const = 1, two.step = F, distributed = F, beta0, gradient_shift){
  # if distributed, then no rescale, with a given beta0 and an extra gradient_shift
  n = length(y);d = ncol(X) + 1- 1*distributed
  X_old = X
  # if (tau_const == 'adaptive') tau_const = robust_const*sqrt(n/(log(n) + d))
  if (distributed == F){
    # standardize y, X, add 1
    mean_X = colMeans(X)
    sd_X = apply(X, 2, sd)
    X = (X - matrix(rep(mean_X, n), byrow = T, nrow = n)) %*% diag(1/sd_X, nrow = length(sd_X))
    X = cbind(rep(1,n), X) 
  }
  # gradient function with the option of gradient shift for surrogate loss
  grad = function(y, X, beta, tau){
    n = length(y)
    res = as.vector(y - X %*% beta)
    w = res * (abs(res)<= tau) + sign(res) * tau * (abs(res)>tau)
    G = -t(X) %*% w / n
    if (distributed) return(G - gradient_shift)
    else return(G)
  }
  # tau update function
  tau_update_dis = function(y, X, beta, tau_const){
    res = y - X%*%beta
    return(tau_const*median(abs(res - median(res))))
  }
  
  tau_update = function(y, X, beta){
    n = nrow(X)
    d = ncol(X)
    res_square = as.vector(y - X%*%beta)^2
    x_max<-sum(res_square)*n/(n-d)
    x_min<-min(res_square)
    func <- function(x,r,n,d){
      sum(pmin(as.vector(r),x))/(x*(n-d))-(d+log(n*d))/n
    }
    return(sqrt(uniroot(func,c(x_min,x_max),r=res_square,n=n,d=d)$root))
  }
  
  # Algorithm loop starts, details see algorithm 2 of Smoothed quantile regression paper
  if (distributed == F){
    beta.t_1 = beta0 = lm(y~X[,-1])$coefficients; 
  }else{
    beta.t_1 = beta0
  }
  tau.t = tau0 = tau_update(y, X, beta0)
  grad.t = grad.0 = grad(y, X, beta.t_1, tau0)
  beta.t = beta.t_1 - grad.0
  count = 0
  while (max(abs(grad(y, X, beta.t, tau.t))) > delta & count <= 100){
    # if (distributed == T) {
    #   tau.t = tau_update_dis(y, X, beta.t, tau_const)
    # }else{
    #   tau.t = tau_update(y, X, beta.t)
    # }
    tau.t = tau_update(y, X, beta.t)
    delta.t = beta.t - beta.t_1
    grad.t = grad(y, X, beta.t, tau.t)
    g.t = grad.t - grad(y, X, beta.t_1, tau.t)
    # eta.1t = as.vector((delta.t %*% delta.t)/ (delta.t %*% g.t))
    # eta.2t = as.vector((delta.t %*% g.t)/ (g.t %*% g.t))
    # eta.t = ifelse(eta.1t > 0, min(eta.1t, eta.2t, 100), 1)
    eta.t = min(10, as.vector((t(delta.t)%*% delta.t)/ (t(delta.t) %*% g.t)))
    beta.t_1 = beta.t
    beta.t = beta.t_1 - eta.t*grad.t
    count = count + 1
  }
  if (distributed == F){
    # scale back beta
    beta.t[1] = beta.t[1] - sum(beta.t[-1]*mean_X/sd_X)
    beta.t[2:length(beta.t)] = beta.t[2:length(beta.t)]/sd_X
  }
  # re-estimate intercept by solving beta0 and tau successively
  if (two.step){
    tau.temp_update = function(x, beta){
      expr = function(x, beta, tau) -sum(sapply(x, function(xx) min(((xx - beta)/tau)^2, 1))) + log(n)
      bounds = c(min(abs(x - beta)), max(abs(x - beta))/sqrt(log(n)/n))
      while (diff(bounds) > delta){
        if (expr(x, beta, mean(bounds)) < 0){
          bounds[1] = mean(bounds)
        }else{
          bounds[2] = mean(bounds)
        }
      }
      return(mean(bounds))
    }
    beta0.temp_update = function(x, tau){
      expr = function(x, beta, tau) -sum(sapply(x, function(xx) min(abs(xx-beta), tau)*sign(xx-beta)))
      bounds = c(min(x), max(x))
      while (diff(bounds) > delta){
        if (expr(x, mean(bounds), tau) < 0){
          bounds[1] = mean(bounds)
        }else{
          bounds[2] = mean(bounds)
        }
      }
      return(mean(bounds))
    }
    delta_hat = as.vector(y_old - X_old%*%beta.t[2:(d+1)])
    tau.temp.last = tau_update(y, X, beta.t, sqrt(n/(log(n) + ncol(X))))
    beta0.temp.last = mean(delta_hat)
    tau.temp = tau.temp_update(delta_hat, beta0.temp.last)
    beta0.temp = beta0.temp_update(delta_hat, tau.temp)
    while (max(abs(tau.temp - tau.temp.last), abs(beta0.temp - beta0.temp.last)) > delta){
      tau.temp.last = tau.temp
      beta0.temp.last = beta0.temp
      tau.temp = tau.temp_update(delta_hat, beta0.temp.last)
      beta0.temp = beta0.temp_update(delta_hat, tau.temp)
    }
    beta.t[1] = beta0.temp
  }
  return(list(coefficients = as.vector(beta.t), tau = tau.t))
}

# implement distributed huber
Distributed_Huber <- function(y, X = NULL, m = 10, robust_const = 2){
  N = length(y)
  n = N/m
  n.iter = ceiling(max(2,log(m)))
  d = ncol(X) + 1
  index = function(c,n) ((c-1)*n+1):(c*n)
  tau_update_dd = function(y, X, beta){
    n = nrow(X)
    d = ncol(X)
    res_square = as.vector(y - X%*%beta)^2
    x_max<-sum(res_square)*n/(n-d)
    x_min<-min(res_square)
    func <- function(x,r,n,d){
      sum(pmin(as.vector(r),x))/(x*(n-d))-(d+log(n*d))/n
    }
    return(sqrt(uniroot(func,c(x_min,x_max),r=res_square,n=n,d=d)$root))
  }
  # gradient function for huber loss, this is not the surrogate loss
  grad = function(y, X, beta, tau){
    n = length(y)
    res = as.vector(y - as.matrix(X) %*% beta)
    w = res * (abs(res)<= tau) + sign(res) * tau * (abs(res)>tau)
    G = -t(X) %*% w / n
    return(G)
  }
  # standardize X using first machine data and then add 1 column of 1
  X_old = X
  mean_X = colMeans(as.matrix(X[index(1,n),]))
  sd_X = apply(as.matrix(X[index(1,n),]), 2, sd)
  X = (X - matrix(rep(mean_X, N), byrow = T, nrow = N)) %*% diag(1/sd_X, nrow = length(sd_X))
  X = cbind(rep(1,N), X)
  
  # initialize beta^(0) using DC OLS
  beta.list = list(rowMeans(sapply(1:m, function(x) lm(y[index(x,n)]~as.matrix(X[index(x,n),-1]))$coefficients)))
  
  kappas = sapply(1:m, function(x) tau_update_dd(y[index(x,n)], X[index(x,n),], beta.list[[1]]))
  kappa = kappas[1]
  tau = robust_const*sqrt(m*(d+log(n))/(d + log(N)))*median(kappas)
  
  # Do n.iter times of the one-step estimation with surrogate loss minimization.
  for (tt in 1:n.iter){
    # constant gradient shift that will be fed into Huber_GDBB
    global_grad = grad(y, X, beta.list[[tt]], tau)
    grad_shift = grad(y[index(1,n)], X[index(1,n),], beta.list[[tt]], kappa) - global_grad
    # Solve surrogate loss minimization using Huber_GDBB, needs input of constant gradient shift and first machine's data.
    beta.temp = Huber_GDBB(y[index(1,n)], as.matrix(X[index(1,n),]), distributed = T, beta0 = beta.list[[tt]], gradient_shift = grad_shift)$coefficients
    beta.list[[tt+1]] = beta.temp
    kappas = sapply(1:m, function(x) tau_update_dd(y[index(x,n)], X[index(x,n),], beta.temp))
    kappa = kappas[1]
    tau = robust_const*sqrt(m*(d+log(n))/(d + log(N)))*median(kappas)
    if (max(abs(beta.temp - beta.list[[tt]])) < 1e-5) break
  }
  # scale back beta
  for (l in 1:length(beta.list)){
    beta.list[[l]][1] = beta.list[[l]][1] - sum(beta.list[[l]][-1]*mean_X/sd_X)
    beta.list[[l]][2:length(beta.list[[l]])] = beta.list[[l]][2:length(beta.list[[l]])]/sd_X
  }
  
  return(list(coefficients = as.vector(beta.list[[length(beta.list)]]), tau = tau))
}


# implement distributed OLS
Distributed_OLS <- function(y, X = NULL, m = 10){
  N = length(y)
  n = N/m
  n.iter = ceiling(max(2,log(m)))
  d = ncol(X) + 1
  index = function(c,n) ((c-1)*n+1):(c*n)
  X_old = X
  X = cbind(rep(1,N), X)
  X1 = X[index(1,n),]
  y1 = y[index(1,n)]
  beta = rowMeans(sapply(1:m, function(x) lm(y[index(x,n)]~as.matrix(X_old[index(x,n),]))$coefficients))
  
  # Do n.iter times of the one-step estimation with surrogate loss minimization.
  XX_inv1 = solve(t(X1)%*%X1)
  local_beta = lm(y1~as.matrix(X1[,-1]))$coefficients
  for (tt in 1:n.iter){
    gradient_shift = t(X1)%*%(y1 - X1%*%beta)/n - t(X)%*%(y - X%*%beta)/N
    beta = local_beta - n*XX_inv1%*%gradient_shift 
  }
  return(beta)
}

# implement pointwise confidence interval of distributed huber using multiplier bootstrap.
Distributed_Inference_Boot <- function(y, X, m = 10, n.boot = 1000, beta, tau, robust_const, quantiles = c(0.025, 0.975)){
  # parameter and function setup
  N = length(y)
  d = ncol(X)
  n = N/m
  index = function(c,n) ((c-1)*n+1):(c*n)
  X = cbind(1, X)
  
  grad = function(y, X, beta, tau){
    n = length(y)
    if (is.vector(X)) X = t(X)
    res = as.vector(y - as.matrix(X) %*% beta) 
    w = res * (abs(res)<= tau) + sign(res) * tau * (abs(res)>tau)
    G = -t(X) %*% w / n
    return(G)
  }
  # calculate m - 1 local gradients, using Global adaptive tau
  grad.all.local = sapply(1:m, function(x) grad(y[index(x,n)], X[index(x,n),], beta, tau))
  grad.bar = rowMeans(grad.all.local)
  grad.local = grad.all.local[,-1] - grad.bar
  # calculate n first machine gradients, using Global adaptive tau
  grad.first = sapply(1:n, function(x) grad(as.vector(y[x]), X[x,], beta, tau)) - grad.bar
  gradients = cbind(grad.first, sqrt(n)*grad.local)[-1,]
  
  # estimate THETA (inverse hessian) using first machine data
  THETA = solve(cov(X[index(1,n),-1]))
  
  # compute multiplyee (-THETA %*% gradients/sqrt(n+m-1))
  multiplyee = t(-THETA %*% gradients/sqrt(n+m-1))
  
  # For 1:n.boot, sapply and produce a d x n.boot table (W.table) for bootstrapped W, no intercept.
  W.table = sapply(1:n.boot, function(x) colSums(multiplyee*rnorm(n+m-1)))
  
  # Get the corresponding quantiles for each column of W.table, output a d X 2 table with lower and upper CI.
  quants = t(apply(W.table, 1, function(x) quantile(x, sort(quantiles, decreasing = T))))/sqrt(N)
  CIs = beta[-1] - quants
  colnames(CIs) = c('2.5%', '97.5%')
  return(CIs)
}


# implement pointwise confidence interval of distributed huber using normal based method.
Distributed_Inference_Normal <- function(y, X, m = 10, beta, tau_const = NULL, robust_const = 2, tau = 0, alpha = 0.05, general = TRUE){
  # parameter and function setup
  N = length(y)
  d = ncol(X)
  n = N/m
  X = cbind(1, X)
  index = function(c,n) ((c-1)*n+1):(c*n)
  tau_update = function(y, X, beta, tau_const){
    res = y - X%*%beta
    return(tau_const*median(abs(res - median(res))))
  }
  tau_update_dd = function(y, X, beta){
    n = nrow(X)
    d = ncol(X)
    res_square = as.vector(y - X%*%beta)^2
    x_max<-sum(res_square)*n/(n-d)
    x_min<-min(res_square)
    func <- function(x,r,n,d){
      sum(pmin(as.vector(r),x))/(x*(n-d))-(d+log(n*d))/n
    }
    return(sqrt(uniroot(func,c(x_min,x_max),r=res_square,n=n,d=d)$root))
  }
  if (tau == 0){
    # tau = tau_update_dd(y,X,beta)
    kappa_median = median(sapply(1:m, function(x) tau_update_dd(y[index(x,n)], X[index(x,n),], beta)))
    tau = robust_const*sqrt(m*(d+log(n))/(d + log(N)))*kappa_median
  }
  res = as.vector(y - as.matrix(X) %*% beta)
  
  if (general == FALSE){
    # special independent case
    sigma_tau_hat = sqrt(sum(pmin(abs(res),tau)^2)/(N-d))
    Sigma_hat_inv = sqrt(rowMeans(sapply(1:m, function(x) diag(n*solve(t(X[index(x,n),])%*%X[index(x,n),]))[-1])))
    
    # Sigma_hat_inv = sqrt(diag(n*solve(t(X[index(1,n),])%*%X[index(1,n),])))[-1]
    width = qnorm(1 - alpha/2)*Sigma_hat_inv*sigma_tau_hat/sqrt(N)
  }
  else{
    # general
    # Sigma_hat1_inv = n*solve(t(X[index(1,n),])%*%X[index(1,n),])
    # Lambda_hat_1 = t(X[index(1,n),])%*%(pmin(abs(res[index(1,n)]),tau)^2*X[index(1,n),])/n
    # Sigma_hat1_inv = solve(cov(X))
    # Lambda_hat_1 = t(X)%*%(pmin(abs(res),tau)^2*X)/N
    # sigma_hat = sqrt(diag(Sigma_hat1_inv%*%Lambda_hat_1%*%Sigma_hat1_inv))[-1]
    local_sigma_hat_diag = function(X, r){
      Sigma_hat_i_inv = n*solve(t(X)%*%X)
      Lambda_hat_i = t(X)%*%(pmin(abs(r),tau)^2*X)/n
      return(diag(Sigma_hat_i_inv%*%Lambda_hat_i%*%Sigma_hat_i_inv)[-1])
    }
    sigma_hat = sqrt(rowMeans(sapply(1:m, function(x) local_sigma_hat_diag(X[index(x,n),], res[index(x,n)]))))
    width = qnorm(1 - alpha/2)*sigma_hat/sqrt(N) 
  }
  
  CI = cbind(beta[-1] - width, beta[-1] + width)
  return(CI)
}



# ILAMM for high-dimensional Adaptive Huber.
Huber_ILAMM <- function(y, X, tau = NULL, beta = NULL, lambda = NULL, epsilon = 1e-4, gamma_u = 1.05, phi0 = 0.0001){
  N = length(y)
  X1 = cbind(1,X)
  # helper functions
  tau_update_dd = function(y, X, beta){
    n = nrow(X)
    d = ncol(X)
    # s = ifelse(sum(beta[-1]!=0) < n/5, sum(beta[-1]!=0), 1)
    s = 1
    res_square = as.vector(y - X%*%beta)^2
    x_max<-sum(res_square)*n/(n-s)
    x_min<-min(res_square)
    func <- function(x,r,s,n,d){
      sum(pmin(as.vector(r),x))/(x*(n-s))-log(n*d)/n
    }
    return(sqrt(uniroot(func,c(x_min,x_max),r=res_square,s=s,n=n,d=d)$root))
  }
  soft_thresh <- function(x, lambda){
    sign(x)*pmax(abs(x) - lambda, 0)
  }
  T_lambda_phi <- function(beta, lambda, phi,Grad_last){
    unthresh_beta = beta - Grad_last/phi
    unthresh_beta[-1] = soft_thresh(unthresh_beta[-1], lambda/phi)
    return(as.numeric(unthresh_beta))
  }
  loss = function(y, X, beta, tau){
    res = as.vector(y - X %*% beta)
    return(mean(0.5*res^2*(abs(res)<=tau) + (tau*abs(res) - tau^2/2)*(abs(res) > tau)))
  }
  grad = function(y, X, beta, tau){
    n = length(y)
    res = as.vector(y - X %*% beta)
    w = res * (abs(res)<= tau) + sign(res) * tau * (abs(res)>tau)
    G = -t(X) %*% w / n
    return(G)
  }
  g_k <- function(beta, beta_last, phi, Loss_last, Grad_last){
    return(Loss_last + t(Grad_last)%*%(beta - beta_last) + phi/2*sum((beta - beta_last)^2))
  }
  lamm <- function(y, X, beta, tau, lambda, phi0=0.0001, gamma_u=1.05){
    phi = phi0
    for (k in 1:200){
      if (k != 1 && norm(beta - beta_last, type = '2')/sqrt(length(beta)) < epsilon) break
      beta_last = beta
      Grad_last = as.numeric(grad(y, X, beta, tau))
      Loss_last = as.numeric(loss(y, X, beta, tau))
      # find minimum phi, hence beta, that locally majorizes loss
      while (TRUE){
        beta = T_lambda_phi(beta_last, lambda, phi, Grad_last)
        if (g_k(beta, beta_last, phi, Loss_last, Grad_last) < as.numeric(loss(y, X, beta, tau))) {
          phi = phi*gamma_u
        }else{
          break
        }
      }
      # print(norm(beta - beta_last, type = '2'))
      phi = max(phi0,phi/gamma_u)
      # cat(phi,',')
      # if (phi <= phi0) cat(phi,',')
    }
    return(list(beta = beta, tau = tau, phi = phi))
  }
  
  # estimate tau using lasso with the specified lambda
  if (is.null(lambda)) lambda = cv.glmnet(X,y, nfolds = min(3, round(N/3)))$lambda.min
  if (is.null(beta)) beta = as.numeric(coef(glmnet(X,y, family = 'gaussian', lambda = lambda)))
  if (is.null(tau)) tau = tau_update_dd(y, X1, beta)
  for (i in 1:100){
    beta_last = beta
    output = lamm(y=y, X=X1, beta=beta, tau = tau, lambda=lambda ,phi0=phi0, gamma_u=gamma_u)
    beta = output$beta
    tau = tau_update_dd(y, X1, beta)
    if (max(abs(beta - beta_last))/sqrt(length(beta)) < epsilon) break
  }
  return(list(beta = beta, tau = tau, lambda = lambda))
}



ILAMM_validation <- function(y, X, y_valid, X_valid, lambda = NULL, nlambda = 10, phi0 = 0.0001){
  d = dim(X)[2]
  N = length(y)
  tau_update_dd = function(y, X, beta){
    n = nrow(X)
    d = ncol(X)
    # s = ifelse(sum(beta[-1]!=0) < n/5, sum(beta[-1]!=0), 1)
    s = 1
    res_square = as.vector(y - X%*%beta)^2
    x_max<-sum(res_square)*n/(n-s)
    x_min<-min(res_square)
    func <- function(x,r,s,n,d){
      sum(pmin(as.vector(r),x))/(x*(n-s))-log(n*d)/n
    }
    return(sqrt(uniroot(func,c(x_min,x_max),r=res_square,s=s,n=n,d=d)$root))
  }
  
  beta = as.numeric(coef(glmnet(X,y, family = 'gaussian', lambda = cv.glmnet(X,y, nfolds = min(3, round(n/3)))$lambda.min)))
  if (is.null(lambda)){
    tau = tau_update_dd(y, cbind(1,X), beta)
    lambda = tau*log(d)/N
  }
  lambda_seq = sort(unique(c(exp(seq(log(lambda/3), log(lambda), length.out = ceiling(nlambda/2 + 0.5))), exp(seq(log(lambda), log(lambda*3), length.out = ceiling(nlambda/2))))) , decreasing = T)  
  
  betalist = matrix(0, nlambda, d+1)
  taulist = SPEs = rep(0, nlambda)
  output1 = Huber_ILAMM(y, X, lambda = lambda_seq[1], phi0 = phi0)
  betalist[1,] = output1$beta
  taulist[1] = output1$tau
  SPEs[1] = sum((y_valid - cbind(1,X_valid)%*%betalist[1,])^2)
  for (i in 2:nlambda){
    output = Huber_ILAMM(y, X, beta = betalist[i-1,], lambda = lambda_seq[i], phi0 = phi0)
    betalist[i,] = output$beta
    taulist[i] = output$tau
    SPEs[i] = sum((y_valid - cbind(1,X_valid)%*%betalist[i,])^2)
    if (i !=1 && SPEs[i] > SPEs[i-1]*(1+1e-4)) break
  }
  index.min = which.min(SPEs[SPEs!=0])
  # Choose optimal lambda using Validation set
  return(list(beta = as.numeric(betalist[index.min,]), tau = taulist[index.min], lambda = lambda_seq[index.min], lambda_seq = lambda_seq))
}





Distributed_ILAMM <- function(y, X, m, epsilon = 1e-4, beta = NULL, tau = NULL, kappa = NULL, lambda = NULL, gamma_u = 1.05, phi0 = 0.0001, robust_const){
  N = length(y)
  n = N/m
  X1 = cbind(1,X)
  max_iter =  ceiling(max(2,log(m)))
  # helper functions
  index = function(c, n) ((c-1)*n+1):(c*n)
  tau_update_dd = function(y, X, beta){
    n = nrow(X)
    d = ncol(X)
    # s = ifelse(sum(beta[-1]!=0) < n/5, sum(beta[-1]!=0), 1)
    s = 1
    res_square = as.vector(y - X%*%beta)^2
    x_max<-sum(res_square)*n/(n-s)
    x_min<-min(res_square)
    func <- function(x,r,s,n,d){
      sum(pmin(as.vector(r),x))/(x*(n-s))-log(n*d)/n
    }
    return(sqrt(uniroot(func,c(x_min,x_max),r=res_square,s=s,n=n,d=d)$root))
  }
  soft_thresh <- function(x, lambda){
    sign(x)*pmax(abs(x) - lambda, 0)
  }
  T_lambda_phi <- function(beta, lambda, phi,Grad_last){
    unthresh_beta = beta - Grad_last/phi
    unthresh_beta[-1] = soft_thresh(unthresh_beta[-1], lambda/phi)
    return(as.numeric(unthresh_beta))
  }
  loss = function(y, X, beta, tau){
    res = as.vector(y - X %*% beta)
    return(mean(0.5*res^2*(abs(res)<=tau) + (tau*abs(res) - tau^2/2)*(abs(res) > tau)))
  }
  grad = function(y, X, beta, tau){
    n = length(y)
    res = as.vector(y - X %*% beta)
    w = res * (abs(res)<= tau) + sign(res) * tau * (abs(res)>tau)
    G = -t(X) %*% w / n
    return(G)
  }
  g_k <- function(beta, beta_last, phi, Loss_last, Grad_last){
    return(Loss_last + t(Grad_last)%*%(beta - beta_last) + phi/2*sum((beta - beta_last)^2))
  }
  lamm <- function(y, X, beta, kappa, gradient_shift, lambda, phi0=0.0001, gamma_u=1.05){
    # phi = max(phi0, phi_last/gamma_u^2)
    phi = phi0
    for (k in 1:200){
      if (k != 1 && norm(beta - beta_last, type = '2')/sqrt(length(beta)) < epsilon) break
      beta_last = beta
      Grad_last = as.numeric(as.numeric(grad(y, X, beta, kappa)) - gradient_shift)
      Loss_last = as.numeric(as.numeric(loss(y, X, beta, kappa)) - t(gradient_shift)%*%beta)
      # find minimum phi, hence beta, that locally majorizes loss
      while (TRUE){
        beta = T_lambda_phi(beta_last, lambda, phi, Grad_last)
        if (g_k(beta, beta_last, phi, Loss_last, Grad_last) < as.numeric(loss(y, X, beta, kappa) - t(gradient_shift)%*%beta)) {
          phi = phi*gamma_u
        }else{
          break
        }
      }
      # if (phi <= phi0) cat(phi,',')
      # print(norm(beta - beta_last, type = '2'))
      # cat(phi,',')
      phi = max(phi0,phi/gamma_u)
    }
    return(list(beta = beta,phi = phi))
  }
  
  # estimate tau using lasso with the specified lambda
  if (is.null(lambda)) lambda = cv.glmnet(X[index(1,n),],y[index(1,n)], nfolds = round(min(3, n/3)))$lambda.min
  if (is.null(beta)) beta = as.numeric(coef(glmnet(X[index(1,n),],y[index(1,n)], family = 'gaussian', lambda = lambda)))
  if (is.null(kappa)){
    kappas = sapply(1:m, function(x) tau_update_dd(y[index(x,n)], X1[index(x,n),], beta))
    kappa = kappas[1]
  }
  if (is.null(tau)) tau = robust_const*median(kappas)*sqrt(m*log(n*d)/log(N*d))
  
  # Surrogate loss alternating update between (kappa, tau) and beta using LAMM
  for (loop in 1:max_iter){
    beta_last = beta
    gradient_shift = grad(y[index(1,n)], X1[index(1,n),], beta, kappa) - grad(y, X1, beta, tau)
    # beta update
    output = lamm(y=y[index(1,n)], X=X1[index(1,n),], beta=beta, kappa=kappa, gradient_shift = gradient_shift, lambda=lambda ,phi0=phi0, gamma_u=gamma_u)
    beta = output$beta
    phi_last = output$phi
    
    # tau update
    kappas = sapply(1:m, function(x) tau_update_dd(y[index(x,n)], X1[index(x,n),], beta))
    kappa = kappas[1]
    tau = robust_const*median(kappas)*sqrt(m*log(n*d)/log(N*d))
    # kappa = tau_update_dd(y[index(1,n)], X1[index(1,n),], beta)
    # tau = tau_update_dd(y, X1, beta)
    
    # stop if max change < epsilon
    # cat(loop,max(abs(beta - beta_last)), '\n')
    if (max(abs(beta - beta_last))/sqrt(d+1) < epsilon) break
  }
  return(list(beta = beta, tau = tau, kappa = kappa, phi = phi_last, lambda = lambda))
}




Distributed_validation_ILAMM <- function(y, X, y_valid, X_valid, m, epsilon = 1e-4, tau = NULL, nlambda = 10, gamma_u = 1.05, phi0 = 0.0001, robust_const = 1){
  N = length(y)
  d = ncol(X)
  n = N/m
  
  # helper functions
  index = function(c, n) ((c-1)*n+1):(c*n)
  tau_update_dd = function(y, X, beta){
    n = nrow(X)
    d = ncol(X)
    # s = ifelse(sum(beta[-1]!=0) < n/5, sum(beta[-1]!=0), 1)
    s = 1
    res_square = as.vector(y - X%*%beta)^2
    x_max<-sum(res_square)*n/(n-s)
    x_min<-min(res_square)
    func <- function(x,r,s,n,d){
      sum(pmin(as.vector(r),x))/(x*(n-s))-log(n*d)/n
    }
    return(sqrt(uniroot(func,c(x_min,x_max),r=res_square,s=s,n=n,d=d)$root))
  }
  
  # find range of lambda using lasso on first machine
  beta = as.numeric(coef(glmnet(X[index(1,n),],y[index(1,n)], family = 'gaussian', lambda = cv.glmnet(X[index(1,n),],y[index(1,n)], nfolds = min(3, round(n/3)))$lambda.min)))
  # kappa = tau_update_dd(y[index(1,n)], cbind(1,X)[index(1,n),], beta)
  # tau = tau_update_dd(y, cbind(1,X), beta)
  
  kappas = sapply(1:m, function(x) tau_update_dd(y[index(x,n)], cbind(1,X)[index(x,n),], beta))
  kappa = kappas[1]
  tau = robust_const*median(kappas)*sqrt(m*log(d*n)/log(d*N))
  
  # lambda_max = max(abs(t(y)%*%X))/N
  lambda_approx = tau*log(d)/N
  lambda_seq = sort(exp(seq(log(lambda_approx/5), log(lambda_approx*5), length.out = nlambda)),decreasing = T)
  
  # estimate distributed beta for a sequence of tau with warm start
  betalist = matrix(0, nlambda, d+1)
  taulist = kappalist = SPEs = rep(0, nlambda)
  output1 = Distributed_ILAMM(y, X, m, beta = beta,epsilon = epsilon, tau = tau, kappa = kappa, lambda = lambda_seq[1], gamma_u = gamma_u, phi0 = phi0, robust_const = robust_const)
  betalist[1,] = output1$beta
  taulist[1] = output1$tau
  kappalist[1] = output1$kappa
  SPEs[1] = sum((y_valid - cbind(1,X_valid)%*%betalist[1,])^2)
  for (i in 2:nlambda){
    output = Distributed_ILAMM(y, X, m, beta = betalist[i-1,], epsilon = epsilon, tau = taulist[i-1], kappa = kappalist[i-1], lambda = lambda_seq[i], gamma_u = gamma_u, phi0 = phi0, robust_const = robust_const) 
    betalist[i,] = output$beta
    taulist[i] = output$tau
    kappalist[i] = output$kappa
    SPEs[i] = sum((y_valid - cbind(1,X_valid)%*%betalist[i,])^2)
    if (i !=1 && SPEs[i] > SPEs[i-1]*(1+1e-4)) break
  }
  index.min = which.min(SPEs[SPEs!=0])
  # Choose optimal lambda using Validation set
  return(list(beta = as.numeric(betalist[index.min,]), tau = taulist[index.min], kappa = kappalist[index.min], lambda = lambda_seq[index.min], lambda_seq = lambda_seq))
}


Distributed_validation_ncvxHuberReg <- function(y, X, y_valid, X_valid, lambda, tau = -1, nlambda = 10){
  lambda_seq = sort(unique(c(exp(seq(log(lambda/5), log(lambda), length.out = ceiling(nlambda/2 + 0.5))), exp(seq(log(lambda), log(lambda*5), length.out = ceiling(nlambda/2))))) , decreasing = T)
  betalist = matrix(0, nlambda, d+1)
  taulist = SPEs = rep(0, nlambda)
  output1 = ncvxHuberReg(X,y, lambda = lambda_seq[1], tau = tau, penalty = 'Lasso',intercept = T)
  betalist[1,] = output1$beta
  taulist[1] = output1$tau
  SPEs[1] = sum((y_valid - cbind(1,X_valid)%*%betalist[1,])^2)
  for (i in 2:nlambda){
    output = ncvxHuberReg(X,y, lambda = lambda_seq[i], tau = tau, penalty = 'Lasso',intercept = T)
    betalist[i,] = output$beta
    taulist[i] = output$tau
    SPEs[i] = sum((y_valid - cbind(1,X_valid)%*%betalist[i,])^2)
    if (i !=1 && SPEs[i] > SPEs[i-1]*(1+1e-4)) break
  }
  index.min = which.min(SPEs[SPEs!=0])
  # Choose optimal lambda using Validation set
  return(list(beta = as.numeric(betalist[index.min,]), tau = taulist[index.min], lambda = lambda_seq[index.min], lambda_seq = lambda_seq))
}
