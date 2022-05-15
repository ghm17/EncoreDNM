suppressMessages(require(mvtnorm))
suppressMessages(require(snowfall))
suppressMessages(require(optparse))

option_list = list(
  make_option('--dat1', action = 'store', default = NA, type = 'character'),
  make_option('--dat2', action = 'store', default = NA, type = 'character'),
  make_option('--N1', action = 'store', default = NA, type = 'numeric'),
  make_option('--N2', action = 'store', default = NA, type = 'numeric'),
  make_option('--dat1_encorednm_single', action = 'store', default = NA, type = 'character'),
  make_option('--dat2_encorednm_single', action = 'store', default = NA, type = 'character'),
  make_option('--mut', action = 'store', default = NA, type = 'character'),
  make_option('--out', action = 'store', default = NA, type = 'character'),
  make_option('--n_cores', action = 'store', default = NA, type = 'numeric')
)
opt = parse_args(OptionParser(option_list = option_list))
dat1 = read.table(opt$dat1, header = T, stringsAsFactors = F)
dat2 = read.table(opt$dat2, header = T, stringsAsFactors = F)
dat1_single = read.table(opt$dat1_encorednm_single, header = T, stringsAsFactors = F)
dat2_single = read.table(opt$dat2_encorednm_single, header = T, stringsAsFactors = F)
mutability = read.table(opt$mut, header = T, stringsAsFactors = F)
categ = colnames(mutability)[-c(1:2)]
N1 = opt$N1
N2 = opt$N2
file_out = opt$out
n_cores = opt$n_cores
n_jackknife = 100



##### some useful functions
EncoreDNM_cross_disorders = function(y1, N1, y2, N2, mut, theta0, n_jackknife, n_cores){
  M = 1000 # Monte Carlo sample size
  G = length(mut) # Number of genes
  theta = theta0
  xi1 = matrix(rnorm(G*M), nrow = G, ncol = M)
  xi2 = matrix(rnorm(G*M), nrow = G, ncol = M)
  
  re = gradient_descent_cross(theta, y1, y2, xi1, xi2, N1, N2, mut)
  theta_est = re[[1]]
  
  se_valid = 0
  ### inversion of Fisher information matrix to calculate standard error
  x = try(solve(-he_cross(theta_est, y1, y2, xi1, xi2, N1, N2, mut)), silent = T)
  if(class(x)[1] != 'try-error'){
    if( min(eigen(x)$values) > 0 & diag(x)[5] < 1){
      se_valid = 1
      theta_se = sqrt(diag(x))
      theta_p = sapply(-abs(theta_est/theta_se), pnorm) * 2
    }
  }
  if(se_valid == 0){ ### groupwise jackknife to calculate standard error
    jack_groups = list()
    remain = 1:G
    for(i in 1:n_jackknife){
      size = floor( length(remain) / (n_jackknife + 1 - i) )
      jack_groups[[i]] = sort(sample(remain, size = size, replace = F))
      remain = setdiff(remain, jack_groups[[i]])
    }
    
    apply.fun = function(jack_ind){
      y1_jack = y1[-jack_groups[[jack_ind]]]
      y2_jack = y2[-jack_groups[[jack_ind]]]
      mut_jack = mut[-jack_groups[[jack_ind]]]
      xi1_jack = xi1[-jack_groups[[jack_ind]], ]
      xi2_jack = xi2[-jack_groups[[jack_ind]], ]
      re = gradient_descent_cross(theta, y1_jack, y2_jack, xi1_jack, xi2_jack, N1, N2, mut_jack)
      theta_est = re[[1]]
      return(theta_est)
    }
    sfInit(parallel = TRUE, cpus = n_cores)
    sfLibrary(snowfall)
    sfExport('jack_groups', 'theta', 'y1', 'y2', 'mut', 'xi1', 'xi2', 'N1', 'N2')
    sfExport('fn_cross', 'gr_cross', 'he_cross', 'gradient_descent_cross')
    Result = sfLapply(1:n_jackknife, apply.fun)
    sfStop()
    result = matrix(0, ncol = 5, nrow = n_jackknife)
    for(jack_ind in 1:n_jackknife){
      result[jack_ind, ] = Result[[jack_ind]]
    }
    theta_se = sqrt(apply(result, 2, var) * (n_jackknife - 1)^2 / n_jackknife)
    theta_p = sapply(-abs(theta_est/theta_se), pnorm) * 2
  }
  
  return(c(theta_est, theta_se, theta_p))
}

fn_cross = function(theta, y1, y2, xi1, xi2, N1, N2, mut){
  beta1 = theta[1]
  beta2 = theta[2]
  sigma1 = theta[3]
  sigma2 = theta[4]
  rho = theta[5]
  
  lambda_i1j = exp(beta1 + sigma1 * xi1) * (2*N1*mut)
  lambda_i2j = exp(beta2 + sigma2 * (rho*xi1 + sqrt(1 - rho^2)*xi2) ) * (2*N2*mut)
  
  temp = -lambda_i1j - lambda_i2j + log(lambda_i1j)*y1 + log(lambda_i2j)*y2
  diff_value = apply(temp, 1, max)
  Delta = exp(temp - diff_value)
  l = sum(log(apply(Delta, 1, sum)) + diff_value)
  return(l)
}

gr_cross = function(theta, y1, y2, xi1, xi2, N1, N2, mut){
  beta1 = theta[1]
  beta2 = theta[2]
  sigma1 = theta[3]
  sigma2 = theta[4]
  rho = theta[5]
  
  lambda_i1j = exp(beta1 + sigma1 * xi1) * (2*N1*mut)
  lambda_i2j = exp(beta2 + sigma2 * (rho*xi1 + sqrt(1 - rho^2)*xi2) ) * (2*N2*mut)
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

he_cross = function(theta, y1, y2, xi1, xi2, N1, N2, mut){
  beta1 = theta[1]
  beta2 = theta[2]
  sigma1 = theta[3]
  sigma2 = theta[4]
  rho = theta[5]
  
  lambda_i1j = exp(beta1 + sigma1 * xi1) * (2*N1*mut)
  lambda_i2j = exp(beta2 + sigma2 * (rho*xi1 + sqrt(1 - rho^2)*xi2) ) * (2*N2*mut)
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

gradient_descent_cross = function(theta, y1, y2, xi1, xi2, N1, N2, mut){
  iter_max = 100
  tol = 1e-5
  
  theta0 = theta
  fn0 = fn_cross(theta0, y1, y2, xi1, xi2, N1, N2, mut)
  gr0 = gr_cross(theta0, y1, y2, xi1, xi2, N1, N2, mut)
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
      fn1 = fn_cross(theta1, y1, y2, xi1, xi2, N1, N2, mut)
      if(fn1<fn0){
        alpha = alpha/2
      }else{
        indic = F
      }
    }
    gr1 = gr_cross(theta1, y1, y2, xi1, xi2, N1, N2, mut)
    diff = sum((theta1 - theta0)^2)
    
    theta0 = theta1
    fn0 = fn1
    gr0 = gr1
    if(diff<tol){
      Ind = F
    }
  }
  he0 = he_cross(theta0, y1, y2, xi1, xi2, N1, N2, mut)
  result = list(theta0, fn0, gr0, he0)
  return(result)
}



##### main
result = data.frame(matrix(0, ncol = 16, nrow = length(categ)))
colnames(result) = c('variant_class', 
                     'beta1_est', 'beta2_est', 'sigma1_est', 'sigma2_est', 'rho_est', 
                     'beta1_se', 'beta2_se', 'sigma1_se', 'sigma2_se', 'rho_se', 
                     'beta1_p', 'beta2_p', 'sigma1_p', 'sigma2_p', 'rho_p')
result$variant_class = categ
for(K in 1:length(categ)){
  mut = mutability[, K + 2]
  y1 = dat1[, K + 2]
  y2 = dat2[, K + 2]
  
  beta1_0 = dat1_single$beta_est[K]
  beta2_0 = dat2_single$beta_est[K]
  sigma1_0 = dat1_single$sigma_est[K]
  sigma2_0 = dat2_single$sigma_est[K]
  rho_0 = 0
  theta0 = c(beta1_0, beta2_0, sigma1_0, sigma2_0, rho_0)
  
  result[K, 2:16] = EncoreDNM_cross_disorders(y1, N1, y2, N2, mut, theta0, n_jackknife, n_cores)
}
write.table(result, file_out, col.names = T, row.names = F, quote = F, sep = '\t')

