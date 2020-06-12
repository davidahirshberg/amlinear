simulators = list(
    function(X, k) {
        a = rowMeans(X[,1:k]) * sqrt(k)
        x = sign(a) * a^2
        tau = -0.2
        w.mean = pmax(0.05, pmin(1/(1 + exp(-x)), 0.95))
        mu = x - tau * (w.mean - 0.5)
        alpha = w.mean
        beta = (1 - w.mean)
        w.var = alpha * beta / (alpha + beta)^2 / (1 + alpha + beta)
        w.fun = function() rbeta(nrow(X), alpha, beta)
        list(mu=mu, tau=tau, w.mean=w.mean, w.var=w.var, w.fun=w.fun, sigma.mult = 1, ape = -0.2)
    },
    function(X, k) {
	k = k - 2
        pr = apply(2 * X[,1:k,drop=FALSE], 1, prod)
        mu = sign(pr) * sqrt(abs(pr)) / 2
        w.mean = 0.1 * sign(mu) + mu
        w.var = w.mean^2
        w.fun = function() w.mean * (1 + rnorm(nrow(X)))
        tau = pmax(X[,1] + X[,2], 0)/ 2
        list(mu=mu, tau=tau, w.mean=w.mean, w.var=w.var, w.fun=w.fun, sigma.mult = 1, ape = 0.5 / sqrt(pi))
    },
    function(X, k) {
        tau = rowMeans(cos(pi * X[,1:k] / 3))
        prob = 0.2 + tau^2
        mu = 4 * rowMeans(X) + 2 * prob
        w.fun = function() rpois(nrow(X), lambda = prob)
        list(mu=mu, tau=tau, w.mean=prob, w.var=prob, w.fun=w.fun, sigma.mult = 1, ape = 0.58)
    },
    function(X, k) {
	k = k + 1
        tau = sin(2 * pi * X[,1])
        rowm = rowMeans(X[,1:k]* sqrt(k))
        ln.mu = 1 / (1 + exp(-sign(rowm) * rowm^2))
        ln.sigma = 1/3
        w.mean = exp(ln.mu + ln.sigma^2/2)
        w.var = (exp(ln.sigma^2) - 1) * w.mean^2
        w.fun = function() exp(ln.mu + ln.sigma * rnorm(nrow(X)))
        mu = pmax(0, 2 * rowm)
        list(mu=mu, tau=tau, w.mean=w.mean, w.var=w.var, w.fun=w.fun, sigma.mult = 1, ape = 0)
    }
)

