suppressMessages(require(mvtnorm))
suppressMessages(require(snowfall))
suppressMessages(require(optparse))

option_list = list(
  make_option('--dat', action = 'store', default = NA, type = 'character'),
  make_option('--N', action = 'store', default = NA, type = 'numeric'),
  make_option('--mut', action = 'store', default = NA, type = 'character'),
  make_option('--out', action = 'store', default = NA, type = 'character'),
  make_option('--n_cores', action = 'store', default = NA, type = 'numeric')
)
opt = parse_args(OptionParser(option_list = option_list))
dat = read.table(opt$dat, header = T, stringsAsFactors = F)
mutability = read.table(opt$mut, header = T, stringsAsFactors = F)
categ = colnames(mutability)[-c(1:2)]
N = opt$N
file_out = opt$out
n_cores = opt$n_cores
n_jackknife = 100



##### some useful functions
EncoreDNM_single_disorder = function(y, N, mut, theta0, n_jackknife, n_cores){
  M = 1000 # Monte Carlo sample size
  G = length(mut) # Number of genes
  theta = theta0
  xi = matrix(rnorm(G*M), nrow = G, ncol = M)
  
  re = gradient_descent_single(theta, y, xi, N, mut)
  theta_est = re[[1]]
  
  se_valid = 0
  ### inversion of Fisher information matrix to calculate standard error
  x = try(solve(-he_single(theta_est, y, xi, N, mut)), silent = T)
  if(class(x)[1] != 'try-error'){
    if( min(eigen(x)$values) > 0){
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
      y_jack = y[-jack_groups[[jack_ind]]]
      mut_jack = mut[-jack_groups[[jack_ind]]]
      xi_jack = xi[-jack_groups[[jack_ind]], ]
      re = gradient_descent_single(theta, y_jack, xi_jack, N, mut_jack)
      theta_est = re[[1]]
      return(theta_est)
    }
    sfInit(parallel = TRUE, cpus = n_cores)
    sfLibrary(snowfall)
    sfExport('jack_groups', 'theta', 'y', 'mut', 'xi', 'N')
    sfExport('fn_single', 'gr_single', 'he_single', 'gradient_descent_single')
    Result = sfLapply(1:n_jackknife, apply.fun)
    sfStop()
    result = matrix(0, ncol = 2, nrow = n_jackknife)
    for(jack_ind in 1:n_jackknife){
      result[jack_ind, ] = Result[[jack_ind]]
    }
    theta_se = sqrt(apply(result, 2, var) * (n_jackknife - 1)^2 / n_jackknife)
    theta_p = sapply(-abs(theta_est/theta_se), pnorm) * 2
  }
  
  beta_est = theta_est[1]
  sigma_est = theta_est[2]
  beta_se = theta_se[1]
  sigma_se = theta_se[2]
  beta_p = theta_p[1]
  sigma_p = theta_p[2]
  l_alternative = fn1 - G*log(M)
  l_null = sum(y[y>0]*log(2*N*mut[y>0]) - y[y>0] + y[y>0] * log(sum(y)/(2*N*sum(mut))))
  deviance_null = 2 * (l_alternative - l_null)
  p_null = pchisq(2 * deviance_null, df = 1, lower.tail = F)
  
  return(c(beta_est, beta_se, beta_p, sigma_est, sigma_se, sigma_p, deviance_null, p_null))
}

fn_single = function(theta, y, xi, N, mut){
  beta1 = theta[1]
  sigma1 = theta[2]
  lambda_ij = exp(beta1 + sigma1 * xi) * (2*N*mut)
  temp = -lambda_ij + log(lambda_ij)*y
  diff_value = apply(temp, 1, max)
  Delta = exp(temp - diff_value)
  l = sum(log(apply(Delta, 1, sum)) + diff_value)
  return(l)
}

gr_single = function(theta, y, xi, N, mut){
  beta1 = theta[1]
  sigma1 = theta[2]
  lambda_ij = exp(beta1 + sigma1 * xi) * (2*N*mut)
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

he_single = function(theta, y, xi, N, mut){
  beta1 = theta[1]
  sigma1 = theta[2]
  lambda_ij = exp(beta1 + sigma1 * xi) * (2*N*mut)
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

gradient_descent_single = function(theta, y, xi, N, mut){
  iter_max = 100
  tol = 1e-5
  
  theta0 = theta
  fn0 = fn_single(theta0, y, xi, N, mut)
  gr0 = gr_single(theta0, y, xi, N, mut)
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
      fn1 = fn_single(theta1, y, xi, N, mut)
      if(fn1<fn0){
        alpha = alpha/2
      }else{
        indic = F
      }
    }
    gr1 = gr_single(theta1, y, xi, N, mut)
    diff = sum((theta1 - theta0)^2)
    
    theta0 = theta1
    fn0 = fn1
    gr0 = gr1
    
    if(diff<tol){
      Ind = F
    }
  }
  he0 = he_single(theta0, y, xi, N, mut)
  result = list(theta0, fn0, gr0, he0)
  return(result)
}



##### main
result = data.frame(matrix(0, ncol = 9, nrow = length(categ)))
colnames(result) = c('variant_class', 'beta_est', 'beta_se', 'beta_p', 'sigma_est', 'sigma_se', 'sigma_p', 'deviance_null', 'p_null')
result$variant_class = categ
for(K in 1:length(categ)){
  mut = mutability[, K + 2]
  y = dat[, K + 2]
  theta0 = c(0, 1)
  result[K, 2:9] = EncoreDNM_single_disorder(y, N, mut, theta0, n_jackknife, n_cores)
}
write.table(result, file_out, col.names = T, row.names = F, quote = F, sep = '\t')

