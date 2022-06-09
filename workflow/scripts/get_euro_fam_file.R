library(argparse)

parser <- ArgumentParser(description = 'Produce PLINK-compatible fam file from 1000G population metadata')
parser$add_argument('-pe', '--ped_file', type = 'character', help = 'Path to ped file', required = T)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to fam file', required = T)

args <- parser$parse_args()

daf <- read.table(args$ped_file, header = T)

euro_daf <- subset(daf, Superpopulation == 'EUR' & FatherID == 0 & MotherID == 0)

# Get unrelated European samples
euro_fam <- data.frame(euro_daf[c('SampleID', 'SampleID', 'FatherID', 'MotherID', 'Sex')], Phenotype = -9)

write.table(euro_fam, file = args$output_file, sep = ' ', col.names = F, row.names = F, quote = F)
