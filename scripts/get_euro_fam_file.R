daf <- openxlsx::read.xlsx('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx')

LDAK_ROOT <- Sys.getenv('ldakRoot')

old_euro_fam <- read.table(file.path(LDAK_ROOT, '1000G_vcf/hg19/euro.fam'))

euro_pops <- c('GBR', 'FIN', 'IBS', 'CEU', 'TSI')
#euro_pops <- c('GBR', 'CEU', 'IBS')

euro_daf <- subset(daf, Population %in% euro_pops)

# TODO overwrite euro.fam
