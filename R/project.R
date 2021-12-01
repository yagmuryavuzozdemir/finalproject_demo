our_lm = function(resp, pred, alpha = 0.05){
      
  resp = as.vector(resp)
  pred = as.matrix(pred)
  
  n = length(resp)
  ones = rep(1, n)
  pred_ones = cbind(ones, pred)
  

  p = ncol(pred_ones)
  
  func = function(beta){
    
    sum((resp - pred_ones%*%beta)^2)
  
  }
  
  beta_0 = rep(NA, p)
  
  beta_0[1] = mean(resp)
  covar = cov(resp, pred_ones)
  
  for(i in 2:p){
    beta_0[i] = covar[i]/var(pred_ones[,i])
  }
  
  minimum = optim(beta_0, func, method = "L-BFGS-B")
  beta_hat = minimum$par
  
  resp_hat = pred_ones%*%beta_hat
  SST = sum((resp-mean(resp))^2)
  SSE = sum((resp-resp_hat)^2)
  R_square = 1-SSE/SST
      
  sigma2_hat = (1/(n-p))*t(resp-pred_ones%*%beta_hat)%*%(resp-pred_ones%*%beta_hat)
  C_p = SSE + 2*p*sigma2_hat
  
  stan_error = sqrt(as.numeric(sigma2_hat)*diag(solve(t(pred_ones)%*%pred_ones)))
  crit_value = qnorm(1-alpha/2)
  conf1 = beta_hat - crit_value*stan_error
  conf2 = beta_hat + crit_value*stan_error
  conf_int = cbind(conf1,conf2)
  
  dfm = p-1
  dfe = n-p
  SSM = sum((resp_hat-mean(resp))^2)
  MSM = SSM/dfm
  MSE = SSE/dfe
  F_star = MSM/MSE
  p_value = 1 - pf(F_star, dfm, dfe)
  residuals = resp - pred_ones%*%beta_hat
  fitted = pred_ones%*%beta_hat
  
  beta_hat = as.table(beta_hat)
  
  par(mar=c(1,1,1,1))
  par(mfrow = c(2, 2))
  plot1 = plot(residuals~fitted)
  abline(h=0, col = "pink")
  qqnorm(residuals)
  qqline(residuals, col = 'purple')
  hist(residuals)
  
  
  return(list(coefficients = beta_hat, confidence_interval = conf_int, R_square = R_square, C_p = C_p, F_stat = F_star, p = p_value, plot1))
}

