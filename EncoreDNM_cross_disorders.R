library(mvtnorm)
source('./functions/functions.R')
args = commandArgs(trailingOnly=TRUE)
dir_dat1 = as.character(args[1])
dir_dat2 = as.character(args[2])
N1 = as.numeric(args[3])
N2 = as.numeric(args[4])
dir_dat1_single = as.character(args[5])
dir_dat2_single = as.character(args[6])
dir_mut = as.character(args[7])
dir_out = as.character(args[8])


dat1 = read.table(dir_dat1, header = T, stringsAsFactors = F)
dat2 = read.table(dir_dat2, header = T, stringsAsFactors = F)
dat1_single = read.table(dir_dat1_single, header = T, stringsAsFactors = F)
dat2_single = read.table(dir_dat2_single, header = T, stringsAsFactors = F)
mut = read.table(dir_mut, header = T, stringsAsFactors = F)
categ = colnames(mut)[-c(1:2)]


result = data.frame(matrix(0, ncol = 16, nrow = length(categ)))
colnames(result) = c('variant_class', 
                     'beta1_est', 'beta2_est', 'sigma1_est', 'sigma2_est', 'rho_est', 
                     'beta1_sd', 'beta2_sd', 'sigma1_sd', 'sigma2_sd', 'rho_sd', 
                     'beta1_p', 'beta2_p', 'sigma1_p', 'sigma2_p', 'rho_p')
result$variant_class = categ
for(K in 1:length(categ)){
  m = mut[, K + 2]
  y1 = dat1[, K + 2]
  y2 = dat2[, K + 2]
  
  beta1_0 = dat1_single$beta_est[K]
  beta2_0 = dat2_single$beta_est[K]
  sigma1_0 = dat1_single$sigma_est[K]
  sigma2_0 = dat2_single$sigma_est[K]
  rho_0 = 0
  theta_0 = c(beta1_0, beta2_0, sigma1_0, sigma2_0, rho_0)
  
  result[K, 2:16] = EncoreDNM_cross_disorders(y1, N1, y2, N2, m, theta_0)
}
write.table(result, dir_out, col.names = T, row.names = F, quote = F, sep = '\t')

