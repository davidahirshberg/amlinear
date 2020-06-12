source("utils2.R")
source("simulators.R")

grid = expand.grid(setup = 1:length(simulators),
		   n = 400,
		   p = c(6,12),
		   sigma = 1,
		   k = 3,
		   rep=1:4)

results = Reduce(rbind, lapply(1:nrow(grid), function(ii) {
	setup = grid$setup[ii]
	n =     grid$n[ii]
	p =     grid$p[ii]
	sigma = grid$sigma[ii]
	k     = grid$k[ii]
	rep   = grid$rep[ii]

	X = matrix(rnorm(n*p), n, p)
	params = simulators[[setup]](X, k)
	W = params$w.fun()
	Y = params$mu + W * params$tau + sigma * params$sigma.mult * rnorm(n)
	btau = mean(params$tau)
	gamma.oracle = (W - params$w.mean) / params$w.var
	cmp = ape.comparison(X, Y, W, gamma.oracle = gamma.oracle)
	cmp$ape = params$ape
	cmp$cape = btau
	print(cmp)
	cmp
}))

