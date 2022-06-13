library(data.table)

save.image('add_ld_block_column_to_merged_gwas.RData')

gwas_file <- snakemake@input[['gwas_file']]
block_file <- snakemake@input[['block_file']]
chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
output_file <- snakemake@output[[1]]

setDTthreads(snakemake@threads)

gwas_dat <- fread(gwas_file, sep = '\t', header = T)

block_dat <- fread(block_file, sep = '\t', header = T)

for(k in 1L:nrow(gwas_dat)) {
  chrom <- gwas_dat[k, chrom_col, env = list(chrom_col = chr_col)]
  bp <- gwas_dat[k, base_col, env = list(base_col = bp_col)]
  block <- block_dat[chr == chrom & bp %between% .(start, stop), block_no]

  if(length(block) == 0) {
    # If SNP lies outside LD blocks
    if(bp < block_dat[chr == chrom, min(start)]) {
      block <- block_dat[chr == chrom][start == min(start), block_no]
    } else if(bp > block_dat[chr == chrom, max(stop)]) {
      block <- block_dat[chr == chrom][stop == max(stop), block_no]
    } else {
      block <- NA
    }
  } else if(length(block) == 2) {
    # If SNP lies on boundary between LD blocks
    block <- max(block)
  } else if(length(block) != 1) {
    stop(sprintf("Invalid no. of blocks for row %d", k))
  }

  set(gwas_dat, i = k, j = 'block', value = block)
}

fwrite(gwas_dat, sep = '\t', file = output_file)
