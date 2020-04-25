setwd("~/Desktop")
library(qtl)
map <- read.csv("Comparative genomics database of 3545 mapped SNPs2.csv",header = TRUE)

pos <- split(map[,4], map[,3])
mar <- split(map[,1], map[,3])
for(i in seq(along=pos))
  names(pos[[i]]) <- mar[[i]]
map <- pos
class(map) <- "map"
plot(map, main="")
legend("bottomright", bty = 'n', 
       c("Total map length: 1471.55 cM", "Total No. of markers: 3545", "Avg. inter-marker spacing: 0.42 cM"))


