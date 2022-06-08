library(data.table)

setDTthreads(snakemake@threads)

trait_A <- snakemake@params[['trait_A']]
trait_B <- snakemake@params[['trait_B']]

gwas_dat <- fread(snakemake@input[['gwas_file']], sep = '\t', header = T)

# The SumHer documentation states that A1 is the 'test allele' and A2 is the 'other allele', so I take this to mean A1 is the minor allele and A2 the major allele

# chr:bp:major:minor in 'variant' column of UKBB sum stats files
dat[, c('chr', 'bp', 'A2', 'A1') := tstrsplit(variant, split = ':', keep = 1:4)]

dat <- dat[!(chr %in% c('X', 'Y', 'MT'))]

dat[, `:=` (len.A1 = nchar(A1), len.A2 = nchar(A2))]

dat <- dat[len.A1 == 1 & len.A2 == 1]

# 'Predictor' has format chr:bp:A2:A1 in my simgwas file, not sure it is correct, though
dat[, Predictor := paste(chr, bp, A2, A1, sep = ':')]

setnames(dat, c(z_colname, sample_size_colname), c('Z', 'n'))

dat <- dat[!duplicated(dat, by = 'Predictor')]

dat <- na.omit(dat)

dat <- dat[, .(Predictor, A1, A2, n, Z)]

fwrite(dat, file = snakemake@output[[1]], sep = '\t', na = 'NA', quote = F)
