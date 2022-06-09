library(cfdr)
library(data.table)
library(parallel)

save.image('iterative_cfdr.RData')

gwas_file <- snakemake@input[[1]]
chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
prin_col <- snakemake@params[['prin_col']]
aux_cols <- snakemake@params[['aux_cols']]
prune_cols <- snakemake@params[['prune_cols']]
v_cols <- snakemake@params[['v_cols']]
v_cols <- c(prin_col, v_cols)
no_of_threads <- snakemake@threads
p_threshold <- snakemake@params[['p_threshold']]
output_file <- snakemake@output[[1]]

setDTthreads(no_of_threads)

gwas_dat <- fread(gwas_file, sep = '\t', header = T)

gwas_dat <- na.omit(gwas_dat, cols = c(chr_col, bp_col, prin_col))

gwas_dat[, (chr_col) := as.integer(chrom_col), env = list(chrom_col = chr_col)]

gwas_dat <- na.omit(gwas_dat, cols = chr_col)

gwas_dat <- unique(gwas_dat, by = c(chr_col, bp_col))

# Transforming p-values of 0 to Z-scores yields Inf values
gwas_dat[prin_p < 1e-300, (prin_col) := 1e-300, env = list(prin_p = prin_col)]

# v_cols is p and initial v_cols
for(i in seq_along(aux_cols)) {
  gwas_dat[aux_p < 1e-300, (aux_cols[i]) := 1e-300, env = list(aux_p = aux_cols[i])]

  # TODO subset columns
  sub_dat <- gwas_dat[!is.na(aux_p), env = list(aux_p = aux_cols[i])]

  # Estimate joint null distribution
  est_q0_pars <- fit.2g(P = sub_dat[!is.na(aux_p) & prune_in == T & prin_p > 0.5, aux_p, env = list(prune_in = prune_cols[i], prin_p = v_cols[i], aux_p = aux_cols[i])])$pars

  # Organise points into chromosome folds
  folds <- lapply(sort(unique(as.integer(sub_dat[[chr_col]]))),
                    function(x) sub_dat[chrom_col == x, which = T, env = list(chrom_col = chr_col)])

  # Not interested in any variants with a ur-principal p-value above the threshold
  candidate_indices <- sub_dat[prin_p <= p_threshold, which = T, env = list(prin_p = prin_col)]

  # List candidate indices by fold
  candidate_indices_by_fold <- lapply(folds, function(x) intersect(candidate_indices, x))

  # Exclude folds with no candidate indices
  non_empty_indices_by_fold <- candidate_indices_by_fold[sapply(candidate_indices_by_fold, function(x) length(x) > 0)]

  # Indices of points in each fold, omitting empty folds (note, all points, not just candidate indices with p < threshold)
  all_indices_by_fold <- folds[sapply(candidate_indices_by_fold, function(x) length(x) > 0)]

  # Compute L-regions
  v <- mcmapply(function(x,y) vl(
                                sub_dat[[v_cols[i]]],
                                 sub_dat[[aux_cols[i]]],
                                indices = x, mode = 2, fold = y),
                non_empty_indices_by_fold, all_indices_by_fold, mc.cores = no_of_threads, SIMPLIFY = F)

  save(v, file = sprintf("v-%d.RData", i))
  # il calls are fast enough not to justify their being parallelised
  # Integrate over L-regions to obtain v-values
  for(j in 1:length(non_empty_indices_by_fold)) {
    sub_dat[non_empty_indices_by_fold[[j]], (v_cols[i+1]) := il(v[[j]], pi0_null = est_q0_pars[1], sigma_null = est_q0_pars[2], distx = "norm")];
  }

  sub_cols <- c(chr_col, bp_col, ref_col, alt_col, v_cols[i+1])

  gwas_dat <- merge(gwas_dat, sub_dat[, ..sub_cols], all.x = T, suffixes = c('', ''), by = c(chr_col, bp_col, ref_col, alt_col))

  gwas_dat[is.na(v_curr), (v_cols[i+1]) := v_prev, env = list(v_curr = v_cols[i+1], v_prev = v_cols[i])]

  fwrite(gwas_dat, file = output_file, sep = '\t')
}

fwrite(gwas_dat, file = output_file, sep = '\t')
