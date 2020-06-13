library(data.table)
library(tidyverse)

results =
    list.files(path = "results") %>% 
    map_df(~fread(paste("results", ., sep='/'), stringsAsFactors = TRUE))
save(results, file="merged.Rdata")

