library(data.table)

setDTthreads(snakemake@threads)

sum_stats_dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

for(x in snakemake@params[['columns_to_drop']]) {
  if(x %in% names(sum_stats_dat)) {
    sum_stats_dat[, drop_col := NULL, env = list(drop_col = x)]
  }
}

setnames(sum_stats_dat,
         c(snakemake@params[['pval_col']], snakemake@params[['beta_col']], snakemake@params[['se_col']]),
         c('pval', 'beta', 'se'))

fwrite(sum_stats_dat, file = snakemake@output[[1]], sep = '\t')
