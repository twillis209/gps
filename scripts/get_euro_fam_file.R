daf <- read.table('1000g/20130606_g1k_3202_samples_ped_population.txt', header = T)

euro_daf <- subset(daf, Superpopulation == 'EUR' & FatherID == 0 & MotherID == 0)

# Get unrelated European samples
euro_fam <- data.frame(euro_daf[c('SampleID', 'SampleID', 'FatherID', 'MotherID', 'Sex')], Phenotype = -9)

write.table(euro_fam, file = '1000g/euro.fam', sep = ' ', col.names = F, row.names = F, quote = F)
