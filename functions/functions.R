EncoreDNM_single_disorder = function(y, N, m){
  M = 1000 # Monte Carlo sample size
  G = length(m) # Number of genes
  
  convergence = 0
  while(convergence<10){
    xi = matrix(rnorm(G*M), nrow = G, ncol = M)
    beta = 0
    sigma = 1
    theta = c(beta, sigma)
    re = gradient_descent_dispersion(theta, y, xi, N, m)
    theta1 = re[[1]]
    
    Ind = 0
    while(Ind<10){
      Ind = Ind + 1
      xi = matrix(rnorm(G*M), nrow = G, ncol = M)
      fn1 = fn_dispersion(theta1, y, xi, N, m)
      gr1 = gr_dispersion(theta1, y, xi, N, m)
      he1 = he_dispersion(theta1, y, xi, N, m)
      x = try(solve(-he_dispersion(theta1, y, xi, N, m)), silent = T)
      if(class(x)!="try-error"){
        if( min(eigen(x)$values) > 0 ){
          Ind = 100
          convergence = 100
          sd1 = sqrt(diag(x))
        }
      }
    }
  }
  l_alternative = fn1 - G*log(M)
  l_null = sum(y[y>0]*log(2*N*m[y>0]) - y[y>0] + y[y>0] * log(sum(y)/(2*N*sum(m))))
  
  beta_est = theta1[1]
  sigma_est = theta1[2]
  beta_sd = sd1[1]
  sigma_sd = sd1[2]
  beta_p = pnorm(-abs(beta_est/beta_sd))*2
  sigma_p = pnorm(-abs(sigma_est/sigma_sd))*2
  deviance_null = 2 * (l_alternative - l_null)
  p_null = pchisq(2 * deviance_null, df = 1, lower.tail = F)
  
  return(c(beta_est, beta_sd, beta_p, sigma_est, sigma_sd, sigma_p, deviance_null, p_null))
}


EncoreDNM_cross_disorders = function(y1, N1, y2, N2, m, theta0){
  M = 1000 # Monte Carlo sample size
  G = length(m) # Number of genes
  beta1_0 = theta0[1]
  beta2_0 = theta0[2]
  sigma1_0 = max(theta0[3], 0.2)
  sigma2_0 = max(theta0[4], 0.2)
  rho_0 = theta0[5]
  
  convergence = 0
  while(convergence<10){
    beta1 = runif(1, 0, abs(beta1_0))*(beta1_0 >= 0) - runif(1, 0, abs(beta1_0))*(beta1_0 < 0)
    beta2 = runif(1, 0, abs(beta2_0))*(beta2_0 >= 0) - runif(1, 0, abs(beta2_0))*(beta2_0 < 0)
    sigma1 = runif(1, 0.2, sigma1_0)
    sigma2 = runif(1, 0.2, sigma2_0)
    rho = rho_0
    theta = c(beta1, beta2, sigma1, sigma2, rho)
    
    convergence = convergence + 1
    xi1 = matrix(rnorm(G*M), nrow = G, ncol = M)
    xi2 = matrix(rnorm(G*M), nrow = G, ncol = M)
    re = gradient_descent_cor(theta, y1, y2, xi1, xi2, N1, N2, m)
    theta_est = re[[1]]
    
    Ind = 0
    while(Ind<5){
      Ind = Ind + 1
      xi1 = matrix(rnorm(G*M), nrow = G, ncol = M)
      xi2 = matrix(rnorm(G*M), nrow = G, ncol = M)
      x = try(solve(-he_cor(theta_est, y1, y2, xi1, xi2, N1, N2, m)), silent = T)
      if(class(x)!="try-error"){
        theta_sd = sqrt(diag(x))
        theta_p = sapply(-abs(theta_est/theta_sd), pnorm) * 2
        if( min(eigen(x)$values) > 0 ){
          Ind = 100
          convergence = 100
        }
      }
    }
  }
  
  return(c(theta_est, theta_sd, theta_p))
}




fn_dispersion = function(theta, y, xi, N, m){
  beta1 = theta[1]
  sigma1 = theta[2]
  lambda_ij = exp(beta1 + sigma1 * xi) * (2*N*m)
  temp = -lambda_ij + log(lambda_ij)*y
  diff_value = apply(temp, 1, max)
  Delta = exp(temp - diff_value)
  l = sum(log(apply(Delta, 1, sum)) + diff_value)
  return(l)
}

gr_dispersion = function(theta, y, xi, N, m){
  beta1 = theta[1]
  sigma1 = theta[2]
  lambda_ij = exp(beta1 + sigma1 * xi) * (2*N*m)
  temp = -lambda_ij + log(lambda_ij)*y
  diff_value = apply(temp, 1, max)
  Delta = exp(temp - diff_value)
  
  grad = rep(0, 2)
  denom = apply(Delta, 1, sum)
  numer = Delta*(-lambda_ij + y)
  grad[1] = sum(apply(numer, 1, sum)/denom)
  numer = Delta*(-lambda_ij + y)*xi
  grad[2] = sum(apply(numer, 1, sum)/denom)
  return(grad)
}

he_dispersion = function(theta, y, xi, N, m){
  beta1 = theta[1]
  sigma1 = theta[2]
  lambda_ij = exp(beta1 + sigma1 * xi) * (2*N*m)
  temp = -lambda_ij + log(lambda_ij)*y
  diff_value = apply(temp, 1, max)
  Delta = exp(temp - diff_value)
  
  denom = apply(Delta, 1, sum)
  hess = matrix(0, ncol = 2, nrow = 2)
  denom2 = denom^2
  numer1 = Delta * ((-lambda_ij + y)^2 - lambda_ij)
  numer2 = Delta * (-lambda_ij + y)
  numer = denom * apply(numer1, 1, sum) - apply(numer2, 1, sum)^2
  hess[1, 1] = sum(numer/denom2)
  numer1 = Delta * ((-lambda_ij + y)^2 - lambda_ij) * xi
  numer2 = Delta * (-lambda_ij + y)
  numer3 = Delta * (-lambda_ij + y) * xi
  numer = denom * apply(numer1, 1, sum) - apply(numer2, 1, sum) * apply(numer3, 1, sum)
  hess[1, 2] = sum(numer/denom2)
  hess[2, 1] = hess[1, 2]
  numer1 = Delta * ((-lambda_ij + y)^2 - lambda_ij) * xi^2
  numer2 = Delta * (-lambda_ij + y) * xi
  numer = denom * apply(numer1, 1, sum) - apply(numer2, 1, sum)^2
  hess[2, 2] = sum(numer/denom2)
  return(hess)
}

gradient_descent_dispersion = function(theta, y, xi, N, m){
  iter_max = 100
  tol = 1e-5
  
  theta0 = theta
  fn0 = fn_dispersion(theta0, y, xi, N, m)
  gr0 = gr_dispersion(theta0, y, xi, N, m)
  iter = 0
  Ind = T
  while(Ind & iter < iter_max){
    iter = iter + 1
    indic = T
    alpha = 0.01
    while(indic){
      theta1 = theta0 + alpha * gr0
      while(theta1[2]<0 | abs(theta1[1])>3 | abs(theta1[2])>3){
        alpha = alpha/2
        theta1 = theta0 + alpha * gr0
      }
      fn1 = fn_dispersion(theta1, y, xi, N, m)
      if(fn1<fn0){
        alpha = alpha/2
      }else{
        indic = F
      }
    }
    gr1 = gr_dispersion(theta1, y, xi, N, m)
    diff = sum((theta1 - theta0)^2)
    
    theta0 = theta1
    fn0 = fn1
    gr0 = gr1
    
    if(diff<tol){
      Ind = F
    }
  }
  he0 = he_dispersion(theta0, y, xi, N, m)
  result = list(theta0, fn0, gr0, he0)
  return(result)
}





fn_cor = function(theta, y1, y2, xi1, xi2, N1, N2, m){
  beta1 = theta[1]
  beta2 = theta[2]
  sigma1 = theta[3]
  sigma2 = theta[4]
  rho = theta[5]

  lambda_i1j = exp(beta1 + sigma1 * xi1) * (2*N1*m)
  lambda_i2j = exp(beta2 + sigma2 * (rho*xi1 + sqrt(1 - rho^2)*xi2) ) * (2*N2*m)
  
  temp = -lambda_i1j - lambda_i2j + log(lambda_i1j)*y1 + log(lambda_i2j)*y2
  diff_value = apply(temp, 1, max)
  Delta = exp(temp - diff_value)
  l = sum(log(apply(Delta, 1, sum)) + diff_value)
  return(l)
}

gr_cor = function(theta, y1, y2, xi1, xi2, N1, N2, m){
  beta1 = theta[1]
  beta2 = theta[2]
  sigma1 = theta[3]
  sigma2 = theta[4]
  rho = theta[5]

  lambda_i1j = exp(beta1 + sigma1 * xi1) * (2*N1*m)
  lambda_i2j = exp(beta2 + sigma2 * (rho*xi1 + sqrt(1 - rho^2)*xi2) ) * (2*N2*m)
  temp = -lambda_i1j - lambda_i2j + log(lambda_i1j)*y1 + log(lambda_i2j)*y2
  diff_value = apply(temp, 1, max)
  Delta = exp(temp - diff_value)
  
  temp1 = -lambda_i1j + y1
  temp2 = -lambda_i2j + y2
  temp3 = rho*xi1 + sqrt(1 - rho^2)*xi2
  temp4 = xi1 - rho/sqrt(1 - rho^2)*xi2

  grad = rep(0, 5)
  denom = apply(Delta, 1, sum)
  numer = Delta*temp1
  grad[1] = sum(apply(numer, 1, sum)/denom)
  numer = Delta*temp2
  grad[2] = sum(apply(numer, 1, sum)/denom)
  numer = Delta*temp1*xi1
  grad[3] = sum(apply(numer, 1, sum)/denom)
  numer = Delta*temp2*temp3
  grad[4] = sum(apply(numer, 1, sum)/denom)
  numer = Delta*temp2*temp4*sigma2
  grad[5] = sum(apply(numer, 1, sum)/denom)

  return(grad)
}

he_cor = function(theta, y1, y2, xi1, xi2, N1, N2, m){
  beta1 = theta[1]
  beta2 = theta[2]
  sigma1 = theta[3]
  sigma2 = theta[4]
  rho = theta[5]

  lambda_i1j = exp(beta1 + sigma1 * xi1) * (2*N1*m)
  lambda_i2j = exp(beta2 + sigma2 * (rho*xi1 + sqrt(1 - rho^2)*xi2) ) * (2*N2*m)
  temp = -lambda_i1j - lambda_i2j + log(lambda_i1j)*y1 + log(lambda_i2j)*y2
  diff_value = apply(temp, 1, max)
  Delta = exp(temp - diff_value)
  
  temp1 = -lambda_i1j + y1
  temp2 = -lambda_i2j + y2
  temp3 = rho*xi1 + sqrt(1 - rho^2)*xi2
  temp4 = xi1 - rho/sqrt(1 - rho^2)*xi2
  temp5 = Delta * temp1
  temp6 = Delta * temp2
  denom = apply(Delta, 1, sum)

  hess = matrix(0, ncol = 5, nrow = 5)
  denom2 = denom^2
  numer1 = Delta * (temp1^2 - lambda_i1j)
  numer = denom * apply(numer1, 1, sum) - apply(temp5, 1, sum)^2
  hess[1, 1] = sum(numer/denom2)
  numer1 = Delta * temp1 * temp2
  numer = denom * apply(numer1, 1, sum) - apply(temp5, 1, sum) * apply(temp6, 1, sum)
  hess[1, 2] = sum(numer/denom2)
  hess[2, 1] = hess[1, 2]
  numer1 = Delta * xi1 *(temp1^2 - lambda_i1j)
  numer2 = temp5 * xi1
  numer = denom * apply(numer1, 1, sum) - apply(temp5, 1, sum) * apply(numer2, 1, sum)
  hess[1, 3] = sum(numer/denom2)
  hess[3, 1] = hess[1, 3]
  numer1 = temp5 * temp2 * temp3
  numer2 = temp6 * temp3
  numer = denom * apply(numer1, 1, sum) - apply(temp5, 1, sum) * apply(numer2, 1, sum)
  hess[1, 4] = sum(numer/denom2)
  hess[4, 1] = hess[1, 4]
  numer1 = temp5 * temp2 * temp4 * sigma2
  numer2 = temp6 * temp4 * sigma2
  numer = denom * apply(numer1, 1, sum) - apply(temp5, 1, sum) * apply(numer2, 1, sum)
  hess[1, 5] = sum(numer/denom2)
  hess[5, 1] = hess[1, 5]
  numer1 = Delta * (temp2^2 - lambda_i2j)
  numer = denom * apply(numer1, 1, sum) - apply(temp6, 1, sum)^2
  hess[2, 2] = sum(numer/denom2)
  numer1 = temp5 * temp2 * xi1
  numer2 = temp5 * xi1
  numer = denom * apply(numer1, 1, sum) - apply(temp6, 1, sum) * apply(numer2, 1, sum)
  hess[2, 3] = sum(numer/denom2)
  hess[3, 2] = hess[2, 3]
  numer1 = Delta * temp3 * (temp2^2 - lambda_i2j)
  numer2 = temp6 * temp3
  numer = denom * apply(numer1, 1, sum) - apply(temp6, 1, sum) * apply(numer2, 1, sum)
  hess[2, 4] = sum(numer/denom2)
  hess[4, 2] = hess[2, 4]
  numer1 = Delta * temp4 * sigma2 * (temp2^2 - lambda_i2j)
  numer2 = temp6 * temp4 * sigma2
  numer = denom * apply(numer1, 1, sum) - apply(temp6, 1, sum) * apply(numer2, 1, sum)
  hess[2, 5] = sum(numer/denom2)
  hess[5, 2] = hess[2, 5]
  numer1 = Delta * xi1^2 * (temp1^2 - lambda_i1j)
  numer2 = temp5 * xi1
  numer = denom * apply(numer1, 1, sum) - apply(numer2, 1, sum)^2
  hess[3, 3] = sum(numer/denom2)
  numer1 = temp5 * xi1 * temp2 * temp3
  numer2 = temp5 * xi1
  numer3 = temp6 * temp3
  numer = denom * apply(numer1, 1, sum) - apply(numer2, 1, sum) * apply(numer3, 1, sum)
  hess[3, 4] = sum(numer/denom2)
  hess[4, 3] = hess[3, 4]
  numer1 = temp5 * xi1 * temp2 * temp4 * sigma2
  numer2 = temp5 * xi1
  numer3 = temp6 * temp4 * sigma2
  numer = denom * apply(numer1, 1, sum) - apply(numer2, 1, sum) * apply(numer3, 1, sum)
  hess[3, 5] = sum(numer/denom2)
  hess[5, 3] = hess[3, 5]
  numer1 = Delta * temp3^2 * (temp2^2 - lambda_i2j)
  numer2 = temp6 * temp3
  numer = denom * apply(numer1, 1, sum) - apply(numer2, 1, sum)^2
  hess[4, 4] = sum(numer/denom2)
  numer1 = Delta * temp4 * (temp2^2*sigma2*temp3 + (-lambda_i2j*sigma2)*temp3 + temp2)
  numer2 = temp6 * temp3
  numer3 = temp6 * temp4 * sigma2
  numer = denom * apply(numer1, 1, sum) - apply(numer2, 1, sum) * apply(numer3, 1, sum)
  hess[4, 5] = sum(numer/denom2)
  hess[5, 4] = hess[4, 5]
  numer1 = Delta * sigma2 * ( temp2^2*temp4^2*sigma2 + (-lambda_i2j)*temp4^2*sigma2 + temp2*(-xi2/(1-rho^2)^(3/2)) )
  numer2 = temp6 * temp4 * sigma2
  numer = denom * apply(numer1, 1, sum) - apply(numer2, 1, sum)^2
  hess[5, 5] = sum(numer/denom2)

  return(hess)
}

gradient_descent_cor = function(theta, y1, y2, xi1, xi2, N1, N2, m){
  iter_max = 100
  tol = 1e-5
  
  theta0 = theta
  fn0 = fn_cor(theta0, y1, y2, xi1, xi2, N1, N2, m)
  gr0 = gr_cor(theta0, y1, y2, xi1, xi2, N1, N2, m)
  iter = 0
  Ind = T
  while(Ind & iter < iter_max){
    iter = iter + 1
    indic = T
    alpha = 0.01
    
    while(indic){
      theta1 = theta0 + alpha * gr0
      while(abs(theta1[1])>3 | abs(theta1[2])>3 | abs(theta1[3])>3 | abs(theta1[4])>3 | theta1[3]<0 | theta1[4]<0 | theta1[5]>1 | theta1[5]< -1){
        alpha = alpha/2
        theta1 = theta0 + alpha * gr0
      }
      fn1 = fn_cor(theta1, y1, y2, xi1, xi2, N1, N2, m)
      if(fn1<fn0){
        alpha = alpha/2
      }else{
        indic = F
      }
    }
    gr1 = gr_cor(theta1, y1, y2, xi1, xi2, N1, N2, m)
    diff = sum((theta1 - theta0)^2)
    
    theta0 = theta1
    fn0 = fn1
    gr0 = gr1
    if(diff<tol){
      Ind = F
    }
  }
  he0 = he_cor(theta0, y1, y2, xi1, xi2, N1, N2, m)
  result = list(theta0, fn0, gr0, he0)
  return(result)
}




