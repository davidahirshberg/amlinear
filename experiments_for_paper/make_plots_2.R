library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(plyr)

load('merged.Rdata') # loads results data.frame

ggplot(data[data$k == unique(data$k)[1] & data$p == unique(data$p)[1],]) + 
    geom_boxplot(aes(x=interaction(weights,outcome), y=estimate-cpsi)) + facet_grid(setup ~ n, scales='free_y', shrink=T) + geom_hline(yintercept=0) +
    theme(axis.text.x=element_text(angle=90,hjust=1))

plot.data = ddply(data, c('setup', 'weights', 'outcome', 'n','p','k','sigma'), function(rows) {
    data.frame(bias = mean(rows$estimate - rows$cpsi, na.rm=T), 
	       rmse = sqrt(mean((rows$estimate - rows$cpsi)^2, na.rm=T)),
	       cov90 = mean(abs(rows$estimate - rows$cpsi)/rows$se.sample < qnorm(.95), na.rm=T))
})

ggplot(plot.data[plot.data$k==unique(plot.data$k)[1] & plot.data$p==unique(plot.data$p)[1], ]) +  
    geom_line(aes(x=n, y=rmse, color=weights, linetype=outcome)) + facet_wrap(.~setup , scales='free_y', shrink=T)

ggplot(plot.data[plot.data$k==unique(plot.data$k)[2] & plot.data$p==unique(plot.data$p)[2], ]) +  
    geom_line(aes(x=n, y=rmse, color=weights, linetype=outcome)) + facet_wrap(.~setup , scales='free_y', shrink=T)
