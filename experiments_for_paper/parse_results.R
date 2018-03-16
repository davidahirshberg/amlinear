rm(list = ls())

library(xtable)
filenames = list.files("results", pattern="*", full.names=TRUE)

nms0 = c("raw", "oracle IPW", "LRB", "DML", "TMLE", "DML.ROB", "TMLE.ROB", "LRB.SIMPLE")
nms = c(paste("L1", nms0), paste("BART", nms0))
param.names = c("setup", "n", "p", "sigma", "spar")

nms = c("plugin", "minimax", "minimax.plus", "oracle")

raw = data.frame(t(sapply(filenames, function(fnm) {
	
	output = read.csv(fnm)[,-1]
	params = strsplit(fnm, "-")[[1]][2:6]
	
	ape = mean(output[,2])
	err = output[, c("plugin.point.estimate",
	              "minimax.point.estimate",
	              "minimax.plus.point.estimate",
	              "oracle.point.estimate")] - ape
	colnames(err) = nms
	
	se = output[, c("plugin.standard.error.estimate",
	             "minimax.standard.error.estimate",
	             "minimax.plus.standard.error.estimate",
	             "oracle.standard.error.estimate")]
	colnames(se) = nms
	
	err.lin = output[, c("plugin.linear.point.estimate",
	              "minimax.linear.point.estimate",
	              "minimax.plus.linear.point.estimate",
	              "oracle.linear.point.estimate")] - ape
	colnames(err.lin) = nms

	rmse = sqrt(colMeans(err^2))
	bias = colMeans(err)
	covg = colMeans(abs(err)/se < qnorm(0.975))
	rmse.lin = sqrt(colMeans(err.lin^2))
	bias.lin = colMeans(err.lin)
	
	c(params, rmse=sprintf("%.2f", round(rmse, 2)),
	  bias=sprintf("%.2f", round(bias, 2)),
	  covg=sprintf("%.2f", round(covg, 2)),
	  rlin=sprintf("%.2f", round(rmse.lin, 2)),
	  blin=sprintf("%.2f", round(bias.lin, 2)))
})))


rownames(raw) = 1:nrow(raw)


names(raw) = c(param.names,
               paste0(nms, ".rmse"),
               paste0(nms, ".bias"),
               paste0(nms, ".covg"),
               paste0(nms, ".rlin"),
               paste0(nms, ".blin"))

options(stringsAsFactors = FALSE)
raw = data.frame(apply(raw, 1:2, as.character))

raw = raw[order(as.numeric(raw$p)),]
raw = raw[order(as.numeric(raw$n)),]
raw = raw[order(raw$setup),]
rownames(raw) = 1:nrow(raw)

write.csv(raw, file="output.csv")

tab.all = raw[,c(2, 3, 5, 6, 10, 14, 7, 11, 15, 8, 12, 16, 9, 13, 17)]
rmse.idx = c(5, 8, 11) - 1
for(iter in 1:nrow(tab.all)) {
	best.idx = rmse.idx[which(as.numeric(tab.all[iter,rmse.idx]) == min(as.numeric(tab.all[iter,rmse.idx])))]
	tab.all[iter,best.idx] = paste("\\bf", tab.all[iter,best.idx])
}

xtab.all = xtable(tab.all, align=c("r", "|", "|", rep("c", 3), "|", "|", rep("c", 3), "|", rep("c", 3), "|", rep("c", 3), "|", rep("c", 3)))
print(xtab.all, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, hline.after = c(-1, -1, 0, 0, 4, 8, 8, 12, 16, 16, 20, 24, 24, 28, 32, 32), file = "simulation_results.tex")
