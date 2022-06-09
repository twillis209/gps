library(data.table)

gwas_file <- snakemake@input[['gwas_file']]
log_file <- snakemake@input[['log_file']]
output_file <- snakemake@output[[1]]

gwas_dat <- fread(file = gwas_file, sep = '\t', header = T)

log_dat <- fread(file = log_file, sep = '\t', header = T)

gwas_dat <- gwas_dat[log_dat$Row]

fwrite(gwas_dat, sep = '\t', file = output_file)
