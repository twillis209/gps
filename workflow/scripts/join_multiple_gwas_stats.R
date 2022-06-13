library(cfdr)
library(data.table)
library(parallel)

principal_trait_gwas_file <- snakemake@input[['principal_trait_gwas_file']]
auxiliary_trait_gwas_files <- snakemake@input[['auxiliary_trait_gwas_files']]
pruned_auxiliary_trait_gwas_files <- snakemake@input[['pruned_auxiliary_trait_gwas_files']]
no_of_threads <- snakemake@threads
chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
p_col <- snakemake@params[['p_col']]
mhc <- snakemake@params[['mhc']]
output_file <- snakemake@output[[1]]

setDTthreads(no_of_threads)

save.image('join_multiple_gwas_stats.RData')

prin_dat <- fread(principal_trait_gwas_file, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col))

prin_dat <- prin_dat[!is.na(p_col), env = list(p_col = p_col)]

prin_dat[, (chr_col) := as.character(chrom), env = list(chrom = chr_col)]

for(i in seq_along(auxiliary_trait_gwas_files)) {
  aux_file <- auxiliary_trait_gwas_files[[i]]
  pruned_aux_file <- pruned_auxiliary_trait_gwas_files[[i]]

  aux_dat <- fread(aux_file, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col))

  aux_p_col <- sprintf('%s.%d', p_col, i)

  setnames(aux_dat, p_col, aux_p_col)

  aux_dat <- aux_dat[!is.na(p_col), env = list(p_col = aux_p_col)]

  aux_dat[, (chr_col) := as.character(chrom), env = list(chrom = chr_col)]

  aux_dat[, (ref_col) := toupper(ref), env = list(ref = ref_col)]
  aux_dat[, (alt_col) := toupper(alt), env = list(alt = alt_col)]

  pruned_dat <- fread(pruned_aux_file, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col))

  pruned_dat[, (chr_col) := as.character(chrom), env = list(chrom = chr_col)]

  pruned_dat[, prune_in := T]

  aux_dat <- merge(aux_dat, pruned_dat, all.x = T, by = c(chr_col, bp_col, ref_col, alt_col))

  aux_dat[is.na(prune_in), prune_in := F]

  setnames(aux_dat, 'prune_in', sprintf('prune_in.%d', i))

  aux_dat <- na.omit(aux_dat, cols = c(chr_col, bp_col))

  prin_dat <- merge(prin_dat, aux_dat, all.x = T, by = c(chr_col, bp_col, ref_col, alt_col))

  prin_dat <- na.omit(prin_dat, cols = c(chr_col, bp_col))
}

if(!mhc) {
  prin_dat <- prin_dat[!(chr_col == 6 & bp_col %between% c(24e6, 45e6)), env = list(chr_col = chr_col, bp_col = bp_col)]
}

fwrite(prin_dat, file = output_file, sep = '\t')
