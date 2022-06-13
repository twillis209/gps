library(pidProjCode)
library(data.table)
library(karyoploteR)
library(GenomicRanges)
library(magrittr)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
chrom_col <- snakemake@params[['chrom_col']]
bp_col <- snakemake@params[['bp_col']]
prin_col <- snakemake@params[['prin_col']]
prin_label <- snakemake@params[['prin_label']]

setDTthreads(snakemake@threads)

plot_height <- 5

plot_width <- 6

main_title_cex <- 1.1

chrom_names_cex <- 0.6

axis_label_cex <- 0.7

axis_label_offset <- 1

axis_tick_cex <- 0.5

points.cex <- 0.5

axis_label_margin <- 0.06

plot_params <- getDefaultPlotParams(plot.type = 4)

plot_params$leftmargin <- 0.08

plot_params$topmargin <- 40

gwas_dat <- fread(input_file, sep = '\t', header = T, select = c(chrom_col, bp_col, prin_col))

gwas_dat[!is.na(prin_col) & !is.na(chrom_col) & !is.na(bp_col), .(chrom_col, bp_col, prin_col), env = list(chrom_col = chrom_col, bp_col = bp_col, prin_col = prin_col)] %>%
  setNames(. , c('CHR38', 'BP38', 'P')) %>% gwas_to_granges(., min_p = 1e-15) -> prin_granges

png(output_file, width = plot_width, height = plot_height, unit = 'in', res = 300)

multitrack_manhattan(list(prin_granges), axis_labels = list('-log10(p)'), main = prin_label, chrom_names_cex = chrom_names_cex,  main_title_cex = main_title_cex, axis_label_cex = axis_label_cex, axis_label_offset = axis_label_offset, axis_tick_cex = axis_tick_cex, chrom_tick_dist = 5e7, axis_label_margin = axis_label_margin, plot_params = plot_params, points.cex = points.cex, points.col = 'brewer.set3')

dev.off()
