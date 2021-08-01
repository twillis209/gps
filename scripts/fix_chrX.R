library(data.table)

bim_dat <- fread('1000g/chrX.bim', sep = '\t', header = F, col.names = c('CHR38', 'ID', 'Cm', 'BP38', 'ALT', 'REF'))

bim_dat[, ref_short := ifelse(nchar(REF) > 10, substr(REF, 1, 10), REF)]
bim_dat[, alt_short := ifelse(nchar(ALT) > 10, substr(ALT, 1, 10), ALT)]

bim_dat[, ID := paste(CHR38, BP38, alt_short, ref_short, sep = ':')]

bim_dat[, c('ref_short', 'alt_short') := NULL]

fwrite(bim_dat, file = '1000g/chrX_fixed.bim', sep = '\t', col.names = F, row.names = F, quote = F)
