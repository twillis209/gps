library(data.table)

setDTthreads(snakemake@threads)

sum_stats_dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

sum_stats_dat[, p := 2*pnorm(abs(b/se), lower.tail = F), env = list(p = snakemake@params[['p_col']], b = snakemake@params[['beta_col']], se = snakemake@params[['se_col']])]
sum_stats_dat[, p := format(p, digits = 20), env = list(p = snakemake@params[['p_col']])]

fwrite(sum_stats_dat, file = snakemake@output[[1]], sep = '\t')
