library(mvtnorm)
source('./functions/functions.R')
args = commandArgs(trailingOnly=TRUE)
dir_dat = as.character(args[1])
N = as.numeric(args[2])
dir_mut = as.character(args[3])
dir_out = as.character(args[4])


dat = read.table(dir_dat, header = T, stringsAsFactors = F)
mut = read.table(dir_mut, header = T, stringsAsFactors = F)
categ = colnames(mut)[-c(1:2)]


result = data.frame(matrix(0, ncol = 9, nrow = length(categ)))
colnames(result) = c('variant_class', 'beta_est', 'beta_sd', 'beta_p', 'sigma_est', 'sigma_sd', 'sigma_p', 'deviance_null', 'p_null')
result$variant_class = categ
for(K in 1:length(categ)){
  m = mut[, K + 2]
  y = dat[, K + 2]
  result[K, 2:9] = EncoreDNM_single_disorder(y, N, m)
}
write.table(result, dir_out, col.names = T, row.names = F, quote = F, sep = '\t')

