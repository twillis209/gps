library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Evaluate GPS statistic on GWAS p-values.')
parser$add_argument('-i', '--gwas_file', type = 'character', help = 'Path to merged GWAS summary statistics file')
parser$add_argument('-chr', type = 'character', help = 'Label of chromosome column in GWAS file', default = 'CHR38')
parser$add_argument('-bp', type = 'character', help = 'Label of BP column in GWAS file', default = 'BP38')
parser$add_argument('-ref', type = 'character', help = 'Label of reference allele column in GWAS file', default = 'REF')
parser$add_argument('-alt', type = 'character', help = 'Label of alternative allele column in GWAS file', default = 'ALT')
parser$add_argument('-prin', type = 'character', help = 'Label of principal p-value column', default = 'P.A')
parser$add_argument('-aux', type = 'character', help = 'Label of auxiliary p-value column', default = 'P.B')
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

#args_vec <- c('-i', 'gwas/pid_aster.tsv.gz', '-o', 'gwas/pid_aster/pruned/pid_aster.tsv.gz', '-nt', 8)
#args<-parser$parse_args(args_vec)

args<-parser$parse_args()

setDTthreads(threads = args$no_of_threads)
RcppParallel::setThreadOptions(numThreads = args$no_of_threads)

gwas_dat <- fread(args$gwas_file, sep = '\t', header = T, select = c(args$chr, args$bp, args$ref, args$alt, args$prin, args$aux))

gps_stat <- pidProjCode::gps_test_stat(gwas_dat[[args$prin]], gwas[[args$aux]])

results_dat <- data.table(gps = gps_stat, no_snps = nrow(gwas_dat))

fwrite(results_dat, file = args$output_file, row.names = F, sep = '\t', col.names = F, quote = F)
