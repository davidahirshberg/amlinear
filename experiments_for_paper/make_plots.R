rm(list = ls())

library(RColorBrewer)

filenames = list.files("results", pattern="*", full.names=TRUE)

new.setup.order = c(1, 2, 3, 4)
param.names = c("setup.neworder", "setup", "n", "p", "sigma", "spar")
nms = c("plugin", "minimax", "minimax.plus", "oracle")

raw = data.frame(t(sapply(filenames, function(fnm) {
    
    output = read.csv(fnm)[,-1]
    params = strsplit(fnm, "-")[[1]][2:6]
    params = c(new.setup.order[as.numeric(params[1])], params)
    
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
    
    c(as.numeric(params), rmse=rmse, bias=bias, covg=covg, rlin=rmse.lin, blin=bias.lin)
})))


rownames(raw) = 1:nrow(raw)


names(raw) = c(param.names,
               paste0(nms, ".rmse"),
               paste0(nms, ".bias"),
               paste0(nms, ".covg"),
               paste0(nms, ".rlin"),
               paste0(nms, ".blin"))

cols = brewer.pal(4, "Set1")

pdf("augmented_cmp.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(raw$minimax.rmse, raw$minimax.rlin, log = "xy",
     col = cols[raw$setup.neworder], pch = 16, cex = 1.5,
     xlab = "augmented minimax linear rmse",
     ylab = "minimax linear rmse",
     xlim = range(raw$minimax.rmse, raw$minimax.rlin),
     ylim = range(raw$minimax.rmse, raw$minimax.rlin),
     xaxt="n", yaxt="n")
axis(1, at=c(0.1, 0.15, 0.2, 0.3, 0.4), labels=c(0.1, 0.15, 0.2, 0.3, 0.4))
axis(2, at=c(0.1, 0.15, 0.2, 0.3, 0.4), labels=c(0.1, 0.15, 0.2, 0.3, 0.4))
abline(0, 1, untf = TRUE, lwd = 2)
abline(0, sqrt(1.25), untf = TRUE, lwd = 2, lty = 2)
legend("bottomright", sapply(1:4, function(xx)paste("setup", xx)), pch = 16, cex = 1.5, col = cols)
par=pardef
dev.off()

pdf("augmented_plus_cmp.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(raw$minimax.plus.rmse, raw$minimax.plus.rlin, log = "xy",
     col = cols[raw$setup.neworder], pch = 16, cex = 1.5,
     xlab = "augmented minimax linear+ rmse",
     ylab = "minimax linear+ rmse",
     xlim = range(raw$minimax.plus.rmse, raw$minimax.plus.rlin),
     ylim = range(raw$minimax.plus.rmse, raw$minimax.plus.rlin),
     xaxt="n", yaxt="n")
axis(1, at=c(0.07, 0.1, 0.15, 0.2, 0.3, 0.4), labels=c(0.07, 0.1, 0.15, 0.2, 0.3, 0.4))
axis(2, at=c(0.07, 0.1, 0.15, 0.2, 0.3, 0.4), labels=c(0.07, 0.1, 0.15, 0.2, 0.3, 0.4))
abline(0, 1, untf = TRUE, lwd = 2)
abline(0, sqrt(1.25), untf = TRUE, lwd = 2, lty = 2)
legend("bottomright", sapply(1:4, function(xx)paste("setup", xx)), pch = 16, cex = 1.5, col = cols)
par=pardef
dev.off()
