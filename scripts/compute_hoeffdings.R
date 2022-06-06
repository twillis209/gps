library(data.table)
library(independence)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
no_of_threads <- snakemake@threads
trait_A_col <- snakemake@params[['trait_A_col']]
trait_B_col <- snakemake@params[['trait_B_col']]

setDTthreads(no_of_threads)

dat <- fread(input_file, sep = '\t', header = T)

dat <- unique(unique(dat, by = trait_A_col), by = trait_B_col)

hoeffding_res <- hoeffding.D.test(xs = dat[[trait_A_col]], ys = dat[[trait_B_col]])

res_dat <- data.table(t(unlist(hoeffding_res[c('n', 'Dn', 'scaled', 'p.value')])))

res_dat[, `:=` (trait_A = trait_A_col, trait_B = trait_B_col)]

res_dat <- res_dat[, c('trait_A', 'trait_B', 'n', 'Dn', 'scaled', 'p.value')]

fwrite(res_dat, file = output_file, sep = '\t', col.names = T, row.names = F)
