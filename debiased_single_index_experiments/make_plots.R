rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
results = read.csv("results.csv")[,-1]

cols = RColorBrewer::brewer.pal(3, "Set2")

n.fact = factor(results$n)

BPDF = data.frame(confound=c(results$confound, results$confound),
                n=c(2 * as.numeric(n.fact) - 1, 2 * as.numeric(n.fact)),
                estimate=c(results$WZ, results$DEB))

ATmat = rbind(1:length(levels(n.fact)), 1:length(levels(n.fact)))
ATmat = ATmat + c(-0.15, 0.15)
ATvec = c(ATmat)

pdf("no_confound.pdf")
pardef = par(xpd = FALSE, mar = c(4.5, 5, 3, 3) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(estimate ~ n, data = subset(BPDF, confound == 0), col = cols[2:3],
        at = ATvec,  pars = list(boxwex = 0.2), xlim = range(ATvec),
        names = rep("", 8), xlab = "n", ylab = "estimate", xaxt = "n")
axis(1, at = 1:length(levels(n.fact)), labels = levels(n.fact))
abline(h=mean(results$oracle[results$confound == 0]), lwd = 2, lty = 3)

par(xpd = TRUE)
legend("bottomright", c("Wooldridge-Zhu", "Augmented Minimax Linear"), fill = cols[2:3], cex = 1.5)
par = pardef
dev.off()

pdf("confound.pdf")
pardef = par(xpd = FALSE, mar = c(4.5, 5, 3, 3) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(estimate ~ n, data = subset(BPDF, confound == 1), col = cols[2:3],
        at = ATvec,  pars = list(boxwex = 0.2), xlim = range(ATvec),
        names = rep("", 8), xlab = "n", ylab = "estimate", xaxt = "n")
axis(1, at = 1:length(levels(n.fact)), labels = levels(n.fact))
abline(h=mean(results$oracle[results$confound == 1]), lwd = 2, lty = 3)

par(xpd = TRUE)
legend("bottomright", c("Wooldridge-Zhu", "Augmented Minimax Linear"), fill = cols[2:3], cex = 1.5)
par = pardef
dev.off()