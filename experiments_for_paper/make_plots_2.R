library(RColorBrewer)
library(ggplot2)
library(plyr)

data = read.csv('merged.csv', header=T)
ggplot(data[data$k == unique(data$k)[2] & data$p == unique(data$p)[2],]) + 
    geom_boxplot(aes(x=estimator, y=estimate-cape)) + facet_grid(setup ~ n, scales='free_y', shrink=T) + geom_hline(yintercept=0) +
    theme(axis.text.x=element_text(angle=90,hjust=1))

plot.data = ddply(data, c('setup', 'estimator', 'n','p','k','sigma'), function(rows) {
    data.frame(bias = mean(rows$estimate - rows$cape, na.rm=T), 
	       rmse = sqrt(mean((rows$estimate - rows$cape)^2, na.rm=T)),
	       cov90 = mean(abs(rows$estimate - rows$cape)/rows$se.sample < qnorm(.95), na.rm=T))
})

ggplot(plot.data[plot.data$k==unique(plot.data$k)[1] & plot.data$p==unique(plot.data$p)[1], ]) +  
    geom_line(aes(x=n, y=rmse, color=estimator)) + facet_wrap(.~setup , scales='free_y', shrink=T)
