library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Recalculate p-values to obtain higher precision')
parser$add_argument('-i', '--input_file', type = 'character', help = 'Path to GWAS summary statistics file for trait A')
parser$add_argument('-p_a', type = 'character', help = 'Label of p-value column for trait A', default = 'P.A')
parser$add_argument('-p_b', type = 'character', help = 'Label of p-value column for trait B', default = 'P.B')
parser$add_argument('-beta_a', type = 'character', help = 'Label of beta column for trait A', default = 'BETA.A')
parser$add_argument('-beta_b', type = 'character', help = 'Label of beta column for trait B', default = 'BETA.B')
parser$add_argument('-se_a', type = 'character', help = 'Label of SE column for trait A', default = 'SE.A')
parser$add_argument('-se_b', type = 'character', help = 'Label of SE column for trait B', default = 'SE.B')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$input_file, sep = '\t', header = T)

sum_stats_dat[, args$p_a := 2*pnorm(abs(get(args$beta_a)/get(args$se_a)), lower.tail = F)]
sum_stats_dat[, args$p_b := 2*pnorm(abs(get(args$beta_b)/get(args$se_b)), lower.tail = F)]

fwrite(sum_stats_dat, file = args$output_path, sep = '\t')
