library(cfdr)
library(data.table)
library(parallel)

gwas_file <- snakemake@input[['gwas_file']]
pruned_gwas_file <- snakemake@input[['pruned_gwas_file']]
no_of_threads <- snakemake@threads
chrom_col <- snakemake@params[['chrom_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
prin_col <- snakemake@params[['prin_col']]
aux_col <- snakemake@params[['aux_col']]
p_threshold <- snakemake@params[['p_threshold']]
v_col <- snakemake@params[['v_col']]
output_file <- snakemake@output[['results_file']]
vl_file <- snakemake@output[['vl_file']]

setDTthreads(no_of_threads)

gwas_dat <- fread(gwas_file, sep = '\t', header = T)

gwas_dat <- gwas_dat[!is.na(prin_col) & !is.na(aux_col), env = list(prin_col = prin_col, aux_col = aux_col)]

# Transforming p-values of 0 to Z-scores yields Inf values
gwas_dat[prin_col < 1e-300, prin_col := 1e-300, env = list(prin_col = prin_col)]
gwas_dat[aux_col < 1e-300, aux_col := 1e-300, env = list(aux_col = aux_col)]

pruned_gwas_dat <- fread(gwas_file, sep = '\t', header = T, select = c(chrom_col, bp_col, ref_col, alt_col))
pruned_gwas_dat[, prune_in := T]

gwas_dat <- merge(gwas_dat, pruned_gwas_dat, all.x = T, by = c(chrom_col, bp_col, ref_col, alt_col))

gwas_dat[is.na(prune_in), prune_in := F]

est_q0_pars <- fit.2g(P = gwas_dat[prune_in == T & prin_col > 0.5, aux_col, env = list(prin_col = prin_col, aux_col = aux_col)])$pars

folds <- mclapply(unique(gwas_dat[[chrom_col]]), function(x) which(gwas_dat[[chrom_col]]==x), mc.cores = no_of_threads)

candidate_indices <- gwas_dat[prin_col <= p_threshold, which = T, env = list(prin_col = prin_col)]

# Organise the candidate indices by fold
ind <- mclapply(folds, function(x) intersect(candidate_indices, x), mc.cores = no_of_threads)

# Exclude ind, folds with no data points
non_empty_indices <- ind[sapply(ind, function(x) length(x) > 0)]

folds_with_indices <- folds[sapply(ind, function(x) length(x) > 0)]

# Compute L-regions
v <- mcmapply(function(x,y) vl(gwas_dat[[prin_col]], gwas_dat[[aux_col]], indices = x, mode = 2, fold = y), non_empty_indices, folds_with_indices, mc.cores = no_of_threads, SIMPLIFY = F)

#if(vl_file) {
#  save(v, vl_file)
#}

# il calls are fast enough not to justify their being parallelised
# Integrate over L-regions to obtain v-values
for(j in 1:length(non_empty_indices)) {
    gwas_dat[non_empty_indices[[j]], v_col := il(v[[j]], pi0_null = est_q0_pars[1], sigma_null = est_q0_pars[2], distx = "norm"), env = list(v_col = v_col)];
}

fwrite(gwas_dat, file = output_file, sep = '\t')
