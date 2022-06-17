library(data.table)

gwas_file_a <- snakemake@input[['A']]
gwas_file_b <- snakemake@input[['B']]
chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
p_col <- snakemake@params[['p_col']]
beta_col <- snakemake@params[['beta_col']]
se_col <- snakemake@params[['se_col']]
output_file <- snakemake@output[['AB']]
recalculate_p <- snakemake@params[['recalculate_p']]
mhc <- snakemake@params[['mhc']]
join <- snakemake@params[['join']]

setDTthreads(snakemake@threads)

dat_a <- fread(gwas_file_a, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col, beta_col, se_col))
dat_b <- fread(gwas_file_b, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col, beta_col, se_col))

dat_a[ , c(ref_col, alt_col) := list(toupper(get(ref_col)), toupper(get(alt_col)))]
dat_a <- dat_a[get(ref_col) %in% c('A','T','C','G') & get(alt_col) %in% c('A','T','C','G')]
dat_a[, (chr_col) := as.character(get(chr_col))]
dat_a <- na.omit(dat_a)

dat_b[ , c(ref_col, alt_col) := list(toupper(get(ref_col)), toupper(get(alt_col)))]
dat_b <- dat_b[get(ref_col) %in% c('A','T','C','G') & get(alt_col) %in% c('A','T','C','G')]
dat_b[, (chr_col) := as.character(get(chr_col))]
dat_b <- na.omit(dat_b)

if(join == 'inner') {
  merged_dat <- merge(dat_a, dat_b, by = c(chr_col, bp_col), suffixes = c('.A', '.B'))
} else if(join == 'left') {
  merged_dat <- merge(dat_a, dat_b, all.x = T, by = c(chr_col, bp_col), suffixes = c('.A', '.B'))
} else if(join == 'right') {
  merged_dat <- merge(dat_a, dat_b, all.y = T, by = c(chr_col, bp_col), suffixes = c('.A', '.B'))
} else {
  stop(sprintf("Unrecognised join param: %s", join))
}

# Removes the MHC
if(!mhc) {
  merged_dat <- merged_dat[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]
}

ref_a <- paste0(ref_col, '.A')
ref_b <- paste0(ref_col, '.B')
alt_a <- paste0(alt_col, '.A')
alt_b <- paste0(alt_col, '.B')

merged_dat <- merged_dat[
(ref_a == ref_b & alt_a == alt_b) |
(ref_a == alt_b & alt_a == ref_b),
env = list(ref_a = ref_a, ref_b = ref_b, alt_a = alt_a, alt_b = alt_b)
]

merged_dat[, c(ref_b, alt_b) := NULL]
setnames(merged_dat, c(ref_a, alt_a), c(ref_col, alt_col))

merged_dat <- unique(merged_dat, by = c(chr_col, bp_col))

if(recalculate_p) {
  p_a <- paste0(p_col, '.A')
  p_b <- paste0(p_col, '.B')
  beta_a <- paste0(beta_col, '.A')
  beta_b <- paste0(beta_col, '.B')
  se_a <- paste0(se_col, '.A')
  se_b <- paste0(se_col, '.B')

  merged_dat[, p_a := 2*pnorm(abs(beta_a/se_a), lower.tail = F),
             env = list(p_a = p_a, beta_a = beta_a, se_a = se_a)]
  merged_dat[, p_a := format(p_a, digits = 20), env = list(p_a = p_a)]

  merged_dat[, p_b := 2*pnorm(abs(beta_b/se_b), lower.tail = F),
             env = list(p_b = p_b, beta_b = beta_b, se_b = se_b)]
  merged_dat[, p_b := format(p_b, digits = 20), env = list(p_b = p_b)]
}

fwrite(merged_dat, file = output_file, sep = '\t', row.names = F)
