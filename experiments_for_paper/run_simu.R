# Trick to check if we're on sherlock
sherlock = nchar(Sys.getenv("GROUP_SCRATCH")) > 0

if (!sherlock) {
  print('running local')
  setup = 7
  n = 400
  p = 3
  sigma = 1
  k = 2
  NREP = 2
  print(paste('running local', setup, n ,p, k, sigma, NREP))
} else {
  set.seed(as.numeric(Sys.getenv("SLURM_JOB_ID")))
  setwd("~/amlinear/experiments_for_paper/")
  setup = sample(1:7, 1)
  n = sample(c(200, 400, 800, 1600), 1)
  p = sample(c(6, 12), 1)
  sigma = sample(c(1), 1)
  k = sample(c(3, 6), 1)
  NREP = 10
  print(paste('running on sherlock', setup, n ,p, k, sigma, NREP))
}

# results filename, which uses a random suffix bit
# check for existence of results file
uniqstr = paste0(c(Sys.getenv(c('SLURM_JOB_ID', 'SLURM_JOB_NAME'))), collapse="")
fnm = paste("results/", uniqstr, ".csv",  sep="")
if(file.exists(fnm)) { stop('done'); }

source("utils2.R")
source("simulators.R")

results.list = lapply(1:NREP, function(iter) {
    X = matrix(rnorm(n*p), n, p)
    params = simulators[[setup]](X, k, sigma)
    W = params$w.fun()
    Y = params$y.fun(W)
    cmp = comparison(X, Y, W, gamma.oracle = params$oracle.fun(W))
    cmp$psi = params$psi
    cmp$cpsi = params$cpsi(W)
    cmp$estimand = params$estimand
    print(cmp)
    cmp
})

results = Reduce(rbind, results.list)
results = cbind(results, setup, n, p, sigma, k)
write.csv(results, file=fnm)

