library(Rcpp)
library(data.table)
setDTthreads(snakemake@threads)

# TODO not using at the moment
sourceCpp(code = '#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double nextAfter(double x, double y) {
   return nextafter(x, y);
}
')

set.seed(snakemake@params[['seed']])

dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

#dat[, BETA.A_pert := BETA.A+runif(.N, min = 0, max = 1e-8)]
#dat[, BETA.B_pert := BETA.B+runif(.N, min = 0, max = 1e-8)]

dat[, P.A_pert := 2*pnorm(abs(BETA.A/SE.A), lower.tail = F)]
dat[, P.B_pert := 2*pnorm(abs(BETA.B/SE.B), lower.tail = F)]

dat[, P.A_pert := format(P.A_pert, digits = 20)]
dat[, P.B_pert := format(P.B_pert, digits = 20)]

fwrite(dat, snakemake@output[[1]], sep = '\t')
