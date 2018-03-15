library(EQL) # for Hermite polynomials
library(rpart)

generate.basis = function(X, order=3) {
    H = lapply(1:ncol(X), function(j) {
        sapply(1:order, function(k) hermite(X[,j], k, prob = TRUE) / sqrt(factorial(k)))
    })
    polys = lapply(1:order, function(r) {
        partitions = combn(r + ncol(X) -1, ncol(X) - 1,
                           function(vec) c(vec, r + ncol(X)) - c(0, vec) - 1)
        elems = sapply(1:ncol(partitions), function(iter) {
            part = partitions[,iter]
            idx = which(part > 0)
            elem = H[[idx[1]]][,part[idx[1]]]
            if (length(idx) > 1) {
                for (id in idx[-1]) {
                    elem = elem * H[[id]][,part[id]]
                }
            }
            elem
        })
        scale(elems) / sqrt(ncol(elems)) / r
    })
    Reduce(cbind, polys)
}

dyadic.basis = function(A, YA, min.n = 200) {
	 indicator = matrix(1, nrow(A), 1)
	 
	 # stop splitting
	 if (nrow(A) <= min.n / 2) {
	 	return(indicator)
	 } else if (nrow(A) <= min.n) {
	 # just one more greedy split
	 	improvement = sapply(1:ncol(A), function(j) {
	 		med.j = median(A[,j])
	 		if(max(A[,j]) == med.j) return (0)
	 		(mean(YA[A[,j] <= med.j]) - mean(YA[A[,j] > med.j]))^2
	 	})
	 	if (max(improvement) == 0) {
	 		return(indicator)
	 	}
	 	# do the greedy split
	 	j.max = which.max(improvement)
	 	med.j.max = median(A[,j.max])
	 	is.left = (A[,j.max] <= med.j.max)
	 	return(cbind(indicator, as.numeric(is.left), as.numeric(!is.left)))
	 }
	 sub.basis = lapply(1:ncol(A), function(j) {
	 	med.j = median(A[,j])
	 	if(max(A[,j]) == med.j) return (rep(1, nrow(A)))
	 	is.left = (A[,j] <= med.j)
	 	left = dyadic.basis(A[is.left,], YA[is.left], min.n)
	 	right = dyadic.basis(A[!is.left,], YA[is.left],  min.n)
	 	all = matrix(0, nrow(A), ncol(left) + ncol(right))
	 	all[is.left, 1:ncol(left)] = left
	 	all[!is.left, ncol(left) + 1:ncol(right)] = right
	 	all
	 })
	 cbind(indicator, Reduce(cbind, sub.basis))
}

make.tree.basis = function(X, W, num.tree = 5) {
	nobs = length(W)
	
	Reduce(cbind, lapply(1:num.tree, function(rep) {
		if (rep == 1) {
			features = 1:ncol(X)
		} else {
			features=sample(1:ncol(X), ceiling(3 * ncol(X)/4))
		}
		
		DF=data.frame(W=W, X[,features])
		tree = rpart(W~., data = DF, control=rpart.control(cp=0))
		nodes = matrix(0, nobs, nrow(tree$frame))
		for(iter in 1:nobs) {
			nodes[iter, tree$where[iter]] = 1
		}
		
		node.names = as.numeric(rownames(tree$frame))
		for(idx in nrow(tree$frame):2) {
			parent.name = floor(node.names[idx]/2)
			parent.idx = which(node.names == parent.name)
			nodes[,parent.idx] = nodes[,parent.idx] + nodes[,idx]
		}
		nodes
	}))
}

ape.comparison = function(X, Y, W, order=NULL, gamma.oracle = rep(1, length(W))) {
	
    if(is.null(order)) { 
		order = 1
		dd = min(1000, nrow(X) - ncol(X))
		while (dd > 0 & order < 10) {
			order = order + 1
			dd = dd - choose(ncol(X) - 1 + order, order)
		}
	}
	B = generate.basis(X, order)
	
	# pre-compute regression adjustment once, then re-use it for different methods
	lasso.out = rlasso(B, Y, W, standardize = FALSE)
	tau.hat = lasso.out$tau.hat
	w.hat = lasso.out$w.hat
	m.hat = lasso.out$y.hat + (W - w.hat) * lasso.out$tau.hat
	
	# make some more things to balance on
	w.rank = rank(w.hat) / (1 + length(w.hat))
	propensity.strata.5 = model.matrix(~cut(w.rank, breaks = (0:5) / 5))
	propensity.strata.10 = model.matrix(~cut(w.rank, breaks = (0:10) / 10))
	propensity.strata.20 = model.matrix(~cut(w.rank, breaks = (0:20) / 20))
	B.tree = make.tree.basis(X, Y, num.tree = 5)
	B.tree.nonpure = B.tree[, colSums(B.tree) >= 40]
	B.dyad = dyadic.basis(X, W, 1.1 * nrow(X) / 2^(floor(log(nrow(X)) / log(2 * ncol(X)))))
	B.plus = cbind(B, B.dyad/2, B.tree.nonpure/2, propensity.strata.5, propensity.strata.10, propensity.strata.20)
	
	# oracle estimate
	point.estimate = mean(tau.hat + gamma.oracle * (Y - m.hat))
	Vhat.oracle = mean((tau.hat - point.estimate)^2 + gamma.oracle^2 * (Y - m.hat)^2)
	standard.error.estimate = sqrt(Vhat.oracle / length(W))
	linear.point.estimate = mean(gamma.oracle * Y)
	ape.oracle = c(point.estimate=point.estimate, standard.error.estimate=standard.error.estimate, linear.point.estimate=linear.point.estimate)
	
	# feasible estimates
	ape.plugin = average_partial_effect(B, Y, W, balance.method = "plugin", fitted.model = lasso.out, verbose = FALSE, standardize = FALSE)
	ape.minimax = average_partial_effect(B, Y, W, balance.method = "minimax", fitted.model = lasso.out, verbose = TRUE, standardize = FALSE)
	ape.minimax.plus = average_partial_effect(B.plus, Y, W, balance.method = "minimax", fitted.model = lasso.out, verbose = TRUE, standardize = FALSE)
	
	c(oracle=ape.oracle, plugin=ape.plugin, minimax=ape.minimax, minimax.plus=ape.minimax.plus)
}
