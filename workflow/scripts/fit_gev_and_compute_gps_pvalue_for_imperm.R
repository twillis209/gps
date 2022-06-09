library(argparse)
library(data.table)
library(evd)
library(fitdistrplus)

gps_file <- snakemake@input[['gps_file']]
coll_file <- snakemake@input[['collated_file']]
all_pvalues_file <- snakemake@output[['all_pvalues_file']]
median_pvalue_file <- snakemake@output[['median_pvalue_file']]
trait_name <- snakemake@params[['trait_name']]
gps_dat <- fread(gps_file, sep = '\t', header = T)

gps <- gps_dat[, GPS]

coll_dat <- fread(coll_file, sep = '\t', header = T, select = c('trait_A', 'trait_B', 'n', 'loc', 'loc.sd', 'scale', 'scale.sd', 'shape', 'shape.sd'))

for(i in 1:nrow(coll_dat)) {
  coll_dat[i, pval := pgev(gps, loc = coll_dat[i, loc], scale = coll_dat[i, scale], shape = coll_dat[i, shape], lower.tail = F)]
}

coll_dat[, gps := gps]

coll_dat <- coll_dat[, .(trait_A, trait_B, gps, n, loc, loc.sd, scale, scale.sd, shape, shape.sd, pval)]

fwrite(coll_dat, sep = '\t', file = all_pvalues_file)

coll_dat <- coll_dat[order(pval, decreasing = F)][ceiling(nrow(coll_dat)/2)]
coll_dat[, trait_B := trait_name]

fwrite(coll_dat[order(pval, decreasing = F)][ceiling(nrow(coll_dat)/2)], sep = '\t', file = median_pvalue_file)
