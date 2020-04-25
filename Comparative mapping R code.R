setwd("~/OneDrive - University of Georgia/Bermudagrass")

library(qtl)
library(circlize)

Alignment <- read.csv("Comparative genomics database of 3545 mapped SNPs.csv", header = TRUE)
Bermuda.LG.max = matrix(NA, 18, 3)
for (i in 1:18){
  Bermuda.LG.max[i,1] <- i
  Bermuda.LG.max[i,2] <- max(Alignment[c(Alignment[,3]==i),4])
  Bermuda.LG.max[i,3] <- sum(Bermuda.LG.max[1:i,2])
}
colnames(Bermuda.LG.max) <- c("LG", "LG max", "Cumulative max")

# Oropetium
Othomaeum.Chr <- read.table("O. thomaeum Chr length.txt", header = TRUE)
Othomaeum.Chr.max = matrix(NA, 10, 3)
for (i in 1:10){
  Othomaeum.Chr.max[i,1] <- i
  Othomaeum.Chr.max[i,2] <- max(Othomaeum.Chr[Othomaeum.Chr[,1]==i,2])
  Othomaeum.Chr.max[i,3] <- sum(Othomaeum.Chr.max[1:i,2])
}
colnames(Othomaeum.Chr.max) <- c("GI", "GI max", "Cumulative max")


## Alignment between bermudagrass LGs and O. thomaeum genome
Bermuda.Thomaeum <- matrix(NA, 3568, 5)
colnames(Bermuda.Thomaeum) <- c("marker", "LG", "cum.pos", "Othomaeum.Chr", "Othomaeum.cum.pos")

for (i in 1:3568){
  Bermuda.Thomaeum[i,1] <- Alignment[i,1]
  Bermuda.Thomaeum[i,2] <- Alignment[i,3]
  Bermuda.Thomaeum[i,4] <- Alignment[i,5]
  Bermuda.Thomaeum[i,3] <- ifelse(Alignment[i,3]==1, Alignment[i,4],
                                  ifelse(Alignment[i,3]==2, Alignment[i,4]+Bermuda.LG.max[1,3],
                                         ifelse(Alignment[i,3]==3, Alignment[i,4]+Bermuda.LG.max[2,3],
                                                ifelse(Alignment[i,3]==4, Alignment[i,4]+Bermuda.LG.max[3,3],
                                                       ifelse(Alignment[i,3]==5, Alignment[i,4]+Bermuda.LG.max[4,3],
                                                              ifelse(Alignment[i,3]==6, Alignment[i,4]+Bermuda.LG.max[5,3],
                                                                     ifelse(Alignment[i,3]==7, Alignment[i,4]+Bermuda.LG.max[6,3],
                                                                            ifelse(Alignment[i,3]==8, Alignment[i,4]+Bermuda.LG.max[7,3],
                                                                                   ifelse(Alignment[i,3]==9, Alignment[i,4]+Bermuda.LG.max[8,3],
                                                                                          ifelse(Alignment[i,3]==10, Alignment[i,4]+Bermuda.LG.max[9,3],
                                                                                                 ifelse(Alignment[i,3]==11, Alignment[i,4]+Bermuda.LG.max[10,3],
                                                                                                        ifelse(Alignment[i,3]==12, Alignment[i,4]+Bermuda.LG.max[11,3],
                                                                                                               ifelse(Alignment[i,3]==13, Alignment[i,4]+Bermuda.LG.max[12,3],
                                                                                                                      ifelse(Alignment[i,3]==14, Alignment[i,4]+Bermuda.LG.max[13,3],
                                                                                                                             ifelse(Alignment[i,3]==15, Alignment[i,4]+Bermuda.LG.max[14,3],
                                                                                                                                    ifelse(Alignment[i,3]==16, Alignment[i,4]+Bermuda.LG.max[15,3],
                                                                                                                                           ifelse(Alignment[i,3]==17, Alignment[i,4]+Bermuda.LG.max[16,3], Alignment[i,4]+Bermuda.LG.max[17,3])))))))))))))))))
  Bermuda.Thomaeum[i,5] <- ifelse(Alignment[i,5]=="gi|1", Alignment[i,6],
                                  ifelse(Alignment[i,5]=="gi|2", Alignment[i,6]+Othomaeum.Chr.max[1,3],
                                         ifelse(Alignment[i,5]=="gi|3", Alignment[i,6]+Othomaeum.Chr.max[2,3],
                                                ifelse(Alignment[i,5]=="gi|4", Alignment[i,6]+Othomaeum.Chr.max[3,3],
                                                       ifelse(Alignment[i,5]=="gi|5", Alignment[i,6]+Othomaeum.Chr.max[4,3],
                                                              ifelse(Alignment[i,5]=="gi|6", Alignment[i,6]+Othomaeum.Chr.max[5,3],
                                                                     ifelse(Alignment[i,5]=="gi|7", Alignment[i,6]+Othomaeum.Chr.max[6,3],
                                                                            ifelse(Alignment[i,5]=="gi|8", Alignment[i,6]+Othomaeum.Chr.max[7,3],
                                                                                   ifelse(Alignment[i,5]=="gi|9", Alignment[i,6]+Othomaeum.Chr.max[8,3], Alignment[i,6]+Othomaeum.Chr.max[9,3])))))))))
}

## genetic map
map <- Alignment[,c(1,3,4)]
pos <- split(map[,3], map[,2])
mar <- split(map[,1], map[,2])
for(i in seq(along=pos))
  names(pos[[i]]) <- mar[[i]]
map <- pos
class(map) <- "map"

tiff("Figure 2. Genetic map and comparative mapping against oropetium.tiff", 
     units = "in", width = 8.5, height = 12.5, res = 300)

par(mfrow=c(2,1), mar=c(4.1, 4.1, 1.1, 2.1))

plot(map, main="", xlab=italic("Cynodon dactylon")~"linkage groups")
legend("bottomright", bty = 'n', 
       c("Total map length: 1471.55 cM", "Total No. of markers: 3544", "Avg. inter-marker spacing: 0.42 cM"))
fig_label("A", region = "figure", pos = "topleft", cex = 1.5)

plot(Bermuda.Thomaeum[,3], Bermuda.Thomaeum[,5], pch=20, col="blue", cex=0.5, 
     xlab = italic("Cynodon dactylon")~"linkage groups", ylab = expression(italic("Oropetium thomaeum")~"chromosomes"), xaxt='n', yaxt='n',
     #main = expression("Genomic synteny between"~italic("Cynodon dactylon")~"and"~italic("Zoysia japonica")),
     ylim = c(0,226591066),xaxs="i", yaxs="i")
fig_label("B", region = "figure", pos = "topleft", cex = 1.5)

abline(h=c(Othomaeum.Chr.max[,3][1:9]), lty=2, col="gray50")
abline(v=c(Bermuda.LG.max[,3][1:17]), lty=2, col="gray50")

#make lables for y axis
ylable.pos = matrix(NA, 10, 2)
for (i in 2:10){
  ylable.pos[i,1] <- i
  ylable.pos[i,2] <- Othomaeum.Chr.max[i-1,3]+(Othomaeum.Chr.max[i,2]/2)
}
ylable.pos[1,1] <- 1
ylable.pos[1,2] <- (34701904/2)


#make lables for x axis
xlable.pos = matrix(NA, 18, 2)
for (i in 2:18){
  xlable.pos[i,1] <- i
  xlable.pos[i,2] <- Bermuda.LG.max[i-1,3]+(Bermuda.LG.max[i,2]/2)
}
xlable.pos[1,1] <- 1
xlable.pos[1,2] <- (100.439/2)



axis(side=2, at=c(ylable.pos[,2]), cex.axis=1,
     labels=c("1","2","3","4","5","6","7","8","9","10"), 
     las=0)

axis(side=1, at=c(xlable.pos[,2]), cex.axis=1,
     labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"), 
     las=0)
axis(side=1, at=xlable.pos[,2][16], cex.axis=1, labels=c("16"), las=0)
dev.off()



# Load zoysia grass genome chromosome lengths
Zoysia.Chr <- read.table("Z. japonica Chr length.txt", header = TRUE)

Zoysia.Chr.max = matrix(NA, 20, 3)
for (i in 1:20){
  Zoysia.Chr.max[i,1] <- i
  Zoysia.Chr.max[i,2] <- max(Zoysia.Chr[Zoysia.Chr[,1]==i,2])
  Zoysia.Chr.max[i,3] <- sum(Zoysia.Chr.max[1:i,2])
}
colnames(Zoysia.Chr.max) <- c("Chr", "Chr max", "Cumulative max")


## Alignment between bermudagrass LGs and O. thomaeum genome
Bermuda.Zoysia <- matrix(NA, 3545, 5)
colnames(Bermuda.Zoysia) <- c("marker", "LG", "cum.pos", "Zoysia.Chr", "Zoysia.cum.pos")

for (i in 1:3545){
  Bermuda.Zoysia[i,1] <- as.character(Alignment[i,1])
  Bermuda.Zoysia[i,2] <- as.numeric(Alignment[i,3])
  Bermuda.Zoysia[i,4] <- as.character(Alignment[i,7])
  Bermuda.Zoysia[i,3] <- ifelse(Alignment[i,3]==1, Alignment[i,4],
                                ifelse(Alignment[i,3]==2, Alignment[i,4]+Bermuda.LG.max[1,3],
                                       ifelse(Alignment[i,3]==3, Alignment[i,4]+Bermuda.LG.max[2,3],
                                              ifelse(Alignment[i,3]==4, Alignment[i,4]+Bermuda.LG.max[3,3],
                                                     ifelse(Alignment[i,3]==5, Alignment[i,4]+Bermuda.LG.max[4,3],
                                                            ifelse(Alignment[i,3]==6, Alignment[i,4]+Bermuda.LG.max[5,3],
                                                                   ifelse(Alignment[i,3]==7, Alignment[i,4]+Bermuda.LG.max[6,3],
                                                                          ifelse(Alignment[i,3]==8, Alignment[i,4]+Bermuda.LG.max[7,3],
                                                                                 ifelse(Alignment[i,3]==9, Alignment[i,4]+Bermuda.LG.max[8,3],
                                                                                        ifelse(Alignment[i,3]==10, Alignment[i,4]+Bermuda.LG.max[9,3],
                                                                                               ifelse(Alignment[i,3]==11, Alignment[i,4]+Bermuda.LG.max[10,3],
                                                                                                      ifelse(Alignment[i,3]==12, Alignment[i,4]+Bermuda.LG.max[11,3],
                                                                                                             ifelse(Alignment[i,3]==13, Alignment[i,4]+Bermuda.LG.max[12,3],
                                                                                                                    ifelse(Alignment[i,3]==14, Alignment[i,4]+Bermuda.LG.max[13,3],
                                                                                                                           ifelse(Alignment[i,3]==15, Alignment[i,4]+Bermuda.LG.max[14,3],
                                                                                                                                  ifelse(Alignment[i,3]==16, Alignment[i,4]+Bermuda.LG.max[15,3],
                                                                                                                                         ifelse(Alignment[i,3]==17, Alignment[i,4]+Bermuda.LG.max[16,3], Alignment[i,4]+Bermuda.LG.max[17,3])))))))))))))))))
  Bermuda.Zoysia[i,5] <- ifelse(Alignment[i,7]=="Zjn_chr01", Alignment[i,8],
                                ifelse(Alignment[i,7]=="Zjn_chr02", Alignment[i,8]+Zoysia.Chr.max[1,3],
                                       ifelse(Alignment[i,7]=="Zjn_chr03", Alignment[i,8]+Zoysia.Chr.max[2,3],
                                              ifelse(Alignment[i,7]=="Zjn_chr04", Alignment[i,8]+Zoysia.Chr.max[3,3],
                                                     ifelse(Alignment[i,7]=="Zjn_chr05", Alignment[i,8]+Zoysia.Chr.max[4,3],
                                                            ifelse(Alignment[i,7]=="Zjn_chr06", Alignment[i,8]+Zoysia.Chr.max[5,3],
                                                                   ifelse(Alignment[i,7]=="Zjn_chr07", Alignment[i,8]+Zoysia.Chr.max[6,3],
                                                                          ifelse(Alignment[i,7]=="Zjn_chr08", Alignment[i,8]+Zoysia.Chr.max[7,3],
                                                                                 ifelse(Alignment[i,7]=="Zjn_chr09", Alignment[i,8]+Zoysia.Chr.max[8,3],
                                                                                        ifelse(Alignment[i,7]=="Zjn_chr10", Alignment[i,8]+Zoysia.Chr.max[9,3],
                                                                                               ifelse(Alignment[i,7]=="Zjn_chr11", Alignment[i,8]+Zoysia.Chr.max[10,3],
                                                                                                      ifelse(Alignment[i,7]=="Zjn_chr12", Alignment[i,8]+Zoysia.Chr.max[11,3],
                                                                                                             ifelse(Alignment[i,7]=="Zjn_chr13", Alignment[i,8]+Zoysia.Chr.max[12,3],
                                                                                                                    ifelse(Alignment[i,7]=="Zjn_chr14", Alignment[i,8]+Zoysia.Chr.max[13,3],
                                                                                                                           ifelse(Alignment[i,7]=="Zjn_chr15", Alignment[i,8]+Zoysia.Chr.max[14,3],
                                                                                                                                  ifelse(Alignment[i,7]=="Zjn_chr16", Alignment[i,8]+Zoysia.Chr.max[15,3],
                                                                                                                                         ifelse(Alignment[i,7]=="Zjn_chr17", Alignment[i,8]+Zoysia.Chr.max[16,3],
                                                                                                                                                ifelse(Alignment[i,7]=="Zjn_chr18", Alignment[i,8]+Zoysia.Chr.max[17,3],
                                                                                                                                                       ifelse(Alignment[i,7]=="Zjn_chr19", Alignment[i,8]+Zoysia.Chr.max[18,3], Alignment[i,8]+Zoysia.Chr.max[19,3])))))))))))))))))))
}

tiff("Figure S1. Genomic synteny between bermudagrass and zoysiagrass.tiff", 
     units = "in", width = 10, height = 8, res = 300)
plot(Bermuda.Zoysia[,3], Bermuda.Zoysia[,5], pch=20, col="blue", cex=0.5, 
     xlab = italic("Cynodon dactylon")~"linkage groups", ylab = expression(italic("Zoysia japonica")~"chromosomes"), xaxt='n', yaxt='n',
     main = expression("Genomic synteny between"~italic("Cynodon dactylon")~"and"~italic("Zoysia japonica")),
     xaxs="i")
abline(h=c(Zoysia.Chr.max[,3][1:19]), lty=2, col="gray50")
abline(v=c(Bermuda.LG.max[,3][1:17]), lty=2, col="gray50")

#make lables for y axis
ylable.pos = matrix(NA, 20, 2)
for (i in 2:20){
  ylable.pos[i,1] <- i
  ylable.pos[i,2] <- Zoysia.Chr.max[i-1,3]+(Zoysia.Chr.max[i,2]/2)
}
ylable.pos[1,1] <- 1
ylable.pos[1,2] <- (18983185/2)


#make lables for x axis
xlable.pos = matrix(NA, 18, 2)
for (i in 2:18){
  xlable.pos[i,1] <- i
  xlable.pos[i,2] <- Bermuda.LG.max[i-1,3]+(Bermuda.LG.max[i,2]/2)
}
xlable.pos[1,1] <- 1
xlable.pos[1,2] <- (100.439/2)



axis(side=2, at=c(ylable.pos[,2]), cex.axis=1,
     labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"), 
     las=0)

axis(side=1, at=c(xlable.pos[,2]), cex.axis=1,
     labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"), 
     las=0)
dev.off()



### Setaria
Sitalica.Chr <- read.table("S. italica chr length.txt", header = TRUE)

Sitalica.Chr.max = matrix(NA, 9, 3)
for (i in 1:9){
  Sitalica.Chr.max[i,1] <- i
  Sitalica.Chr.max[i,2] <- max(Sitalica.Chr[Sitalica.Chr[,1]==i,2])
  Sitalica.Chr.max[i,3] <- sum(Sitalica.Chr.max[1:i,2])
}
colnames(Bermuda.LG.max) <- c("Chr", "Chr max", "Cumulative max")


## Alignment between bermudagrass LGs and O. thomaeum genome
Bermuda.Sitalica <- matrix(NA, 3545, 5)
colnames(Bermuda.Sitalica) <- c("marker", "LG", "cum.pos", "Sitalica.Chr", "Sitalica.cum.pos")

for (i in 1:3545){
  Bermuda.Sitalica[i,1] <- Alignment[i,1]
  Bermuda.Sitalica[i,2] <- Alignment[i,3]
  Bermuda.Sitalica[i,4] <- Alignment[i,9]
  Bermuda.Sitalica[i,3] <- ifelse(Alignment[i,3]==1, Alignment[i,4],
                                  ifelse(Alignment[i,3]==2, Alignment[i,4]+Bermuda.LG.max[1,3],
                                         ifelse(Alignment[i,3]==3, Alignment[i,4]+Bermuda.LG.max[2,3],
                                                ifelse(Alignment[i,3]==4, Alignment[i,4]+Bermuda.LG.max[3,3],
                                                       ifelse(Alignment[i,3]==5, Alignment[i,4]+Bermuda.LG.max[4,3],
                                                              ifelse(Alignment[i,3]==6, Alignment[i,4]+Bermuda.LG.max[5,3],
                                                                     ifelse(Alignment[i,3]==7, Alignment[i,4]+Bermuda.LG.max[6,3],
                                                                            ifelse(Alignment[i,3]==8, Alignment[i,4]+Bermuda.LG.max[7,3],
                                                                                   ifelse(Alignment[i,3]==9, Alignment[i,4]+Bermuda.LG.max[8,3],
                                                                                          ifelse(Alignment[i,3]==10, Alignment[i,4]+Bermuda.LG.max[9,3],
                                                                                                 ifelse(Alignment[i,3]==11, Alignment[i,4]+Bermuda.LG.max[10,3],
                                                                                                        ifelse(Alignment[i,3]==12, Alignment[i,4]+Bermuda.LG.max[11,3],
                                                                                                               ifelse(Alignment[i,3]==13, Alignment[i,4]+Bermuda.LG.max[12,3],
                                                                                                                      ifelse(Alignment[i,3]==14, Alignment[i,4]+Bermuda.LG.max[13,3],
                                                                                                                             ifelse(Alignment[i,3]==15, Alignment[i,4]+Bermuda.LG.max[14,3],
                                                                                                                                    ifelse(Alignment[i,3]==16, Alignment[i,4]+Bermuda.LG.max[15,3],
                                                                                                                                           ifelse(Alignment[i,3]==17, Alignment[i,4]+Bermuda.LG.max[16,3], Alignment[i,4]+Bermuda.LG.max[17,3])))))))))))))))))
  Bermuda.Sitalica[i,5] <- ifelse(Alignment[i,9]=="scaffold_1", Alignment[i,10],
                                  ifelse(Alignment[i,9]=="scaffold_2", Alignment[i,10]+Sitalica.Chr.max[1,3],
                                         ifelse(Alignment[i,9]=="scaffold_3", Alignment[i,10]+Sitalica.Chr.max[2,3],
                                                ifelse(Alignment[i,9]=="scaffold_4", Alignment[i,10]+Sitalica.Chr.max[3,3],
                                                       ifelse(Alignment[i,9]=="scaffold_5", Alignment[i,10]+Sitalica.Chr.max[4,3],
                                                              ifelse(Alignment[i,9]=="scaffold_6", Alignment[i,10]+Sitalica.Chr.max[5,3],
                                                                     ifelse(Alignment[i,9]=="scaffold_7", Alignment[i,10]+Sitalica.Chr.max[6,3],
                                                                            ifelse(Alignment[i,9]=="scaffold_8", Alignment[i,10]+Sitalica.Chr.max[7,3], Alignment[i,10]+Sitalica.Chr.max[8,3]))))))))
}



tiff("Figure S2. Genomic synteny between bermudagrass and setaria.tiff", units = "in",
     width = 10, height = 8, res = 300)

plot(Bermuda.Sitalica[,3], Bermuda.Sitalica[,5], pch=20, col="blue", cex=0.5, 
     xlab = italic("Cynodon dactylon")~"linkage groups",  ylab = expression(italic("Setaria italica")~"chromosomes"), xaxt='n', yaxt='n',
     main = expression("Genomic synteny between"~italic("Cynodon dactylon")~"and"~italic("Setaria italica")),
     xaxs="i", yaxs="i")
abline(h=c(Sitalica.Chr.max[,3][1:8]), lty=2, col="gray50")
abline(v=c(Bermuda.LG.max[,3][1:17]), lty=2, col="gray50")

#make lables for y axis
ylable.pos = matrix(NA, 9, 2)
for (i in 2:9){
  ylable.pos[i,1] <- i
  ylable.pos[i,2] <- Sitalica.Chr.max[i-1,3]+(Sitalica.Chr.max[i,2]/2)
}
ylable.pos[1,1] <- 1
ylable.pos[1,2] <- (42145699/2)


#make lables for x axis
xlable.pos = matrix(NA, 18, 2)
for (i in 2:18){
  xlable.pos[i,1] <- i
  xlable.pos[i,2] <- Bermuda.LG.max[i-1,3]+(Bermuda.LG.max[i,2]/2)
}
xlable.pos[1,1] <- 1
xlable.pos[1,2] <- (100.439/2)



axis(side=2, at=c(ylable.pos[,2]), cex.axis=0.9,
     labels=c("1","2","3","4","5","6","7","8","9"), 
     las=0)

axis(side=1, at=c(xlable.pos[,2]), cex.axis=0.9,
     labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"), 
     las=0)
dev.off()


## Sorghum
Sbicolor.Chr <- read.table("S. bicolor Chr length.txt", header = TRUE)
Sbicolor.Chr.max = matrix(NA, 10, 3)
for (i in 1:10){
  Sbicolor.Chr.max[i,1] <- i
  Sbicolor.Chr.max[i,2] <- max(Sbicolor.Chr[Sbicolor.Chr[,1]==i,2])
  Sbicolor.Chr.max[i,3] <- sum(Sbicolor.Chr.max[1:i,2])
}
colnames(Bermuda.LG.max) <- c("Chr", "Chr max", "Cumulative max")


## Alignment between bermudagrass LGs and O. thomaeum genome
Bermuda.Sbicolor <- matrix(NA, 3545, 5)
colnames(Bermuda.Sbicolor) <- c("marker", "LG", "cum.pos", "Sbicolor.Chr", "Sbicolor.cum.pos")

for (i in 1:3545){
  Bermuda.Sbicolor[i,1] <- Alignment[i,1]
  Bermuda.Sbicolor[i,2] <- Alignment[i,3]
  Bermuda.Sbicolor[i,4] <- Alignment[i,11]
  Bermuda.Sbicolor[i,3] <- ifelse(Alignment[i,3]==1, Alignment[i,4],
                                  ifelse(Alignment[i,3]==2, Alignment[i,4]+Bermuda.LG.max[1,3],
                                         ifelse(Alignment[i,3]==3, Alignment[i,4]+Bermuda.LG.max[2,3],
                                                ifelse(Alignment[i,3]==4, Alignment[i,4]+Bermuda.LG.max[3,3],
                                                       ifelse(Alignment[i,3]==5, Alignment[i,4]+Bermuda.LG.max[4,3],
                                                              ifelse(Alignment[i,3]==6, Alignment[i,4]+Bermuda.LG.max[5,3],
                                                                     ifelse(Alignment[i,3]==7, Alignment[i,4]+Bermuda.LG.max[6,3],
                                                                            ifelse(Alignment[i,3]==8, Alignment[i,4]+Bermuda.LG.max[7,3],
                                                                                   ifelse(Alignment[i,3]==9, Alignment[i,4]+Bermuda.LG.max[8,3],
                                                                                          ifelse(Alignment[i,3]==10, Alignment[i,4]+Bermuda.LG.max[9,3],
                                                                                                 ifelse(Alignment[i,3]==11, Alignment[i,4]+Bermuda.LG.max[10,3],
                                                                                                        ifelse(Alignment[i,3]==12, Alignment[i,4]+Bermuda.LG.max[11,3],
                                                                                                               ifelse(Alignment[i,3]==13, Alignment[i,4]+Bermuda.LG.max[12,3],
                                                                                                                      ifelse(Alignment[i,3]==14, Alignment[i,4]+Bermuda.LG.max[13,3],
                                                                                                                             ifelse(Alignment[i,3]==15, Alignment[i,4]+Bermuda.LG.max[14,3],
                                                                                                                                    ifelse(Alignment[i,3]==16, Alignment[i,4]+Bermuda.LG.max[15,3],
                                                                                                                                           ifelse(Alignment[i,3]==17, Alignment[i,4]+Bermuda.LG.max[16,3], Alignment[i,4]+Bermuda.LG.max[17,3])))))))))))))))))
  Bermuda.Sbicolor[i,5] <- ifelse(Alignment[i,11]=="Chr01", Alignment[i,12],
                                  ifelse(Alignment[i,11]=="Chr02", Alignment[i,12]+Sbicolor.Chr.max[1,3],
                                         ifelse(Alignment[i,11]=="Chr03", Alignment[i,12]+Sbicolor.Chr.max[2,3],
                                                ifelse(Alignment[i,11]=="Chr04", Alignment[i,12]+Sbicolor.Chr.max[3,3],
                                                       ifelse(Alignment[i,11]=="Chr05", Alignment[i,12]+Sbicolor.Chr.max[4,3],
                                                              ifelse(Alignment[i,11]=="Chr06", Alignment[i,12]+Sbicolor.Chr.max[5,3],
                                                                     ifelse(Alignment[i,11]=="Chr07", Alignment[i,12]+Sbicolor.Chr.max[6,3],
                                                                            ifelse(Alignment[i,11]=="Chr08", Alignment[i,12]+Sbicolor.Chr.max[7,3],
                                                                                   ifelse(Alignment[i,11]=="Chr09", Alignment[i,12]+Sbicolor.Chr.max[8,3], Alignment[i,12]+Sbicolor.Chr.max[9,3])))))))))
}

tiff("Figure S3. Genomic synteny between bermudagrass and sorghum.tiff", 
     units = "in", width = 10, height = 8, res = 300)

plot(Bermuda.Sbicolor[,3], Bermuda.Sbicolor[,5], pch=20, col="blue", cex=0.5, 
     xlab = italic("Cynodon dactylon")~"linkage groups", ylab = expression(italic("Sorghum bicolor")~"chromosomes"), xaxt='n', yaxt='n',
     main = expression("Genomic synteny between"~italic("Cynodon dactylon")~"and"~italic("Sorghum bicolor")),
     xaxs="i", yaxs="i")
abline(h=c(Sbicolor.Chr.max[,3][1:9]), lty=2, col="gray50")
abline(v=c(Bermuda.LG.max[,3][1:17]), lty=2, col="gray50")

#make lables for y axis
ylable.pos = matrix(NA, 10, 2)
for (i in 2:10){
  ylable.pos[i,1] <- i
  ylable.pos[i,2] <- Sbicolor.Chr.max[i-1,3]+(Sbicolor.Chr.max[i,2]/2)
}
ylable.pos[1,1] <- 1
ylable.pos[1,2] <- (73727935/2)


#make lables for x axis
xlable.pos = matrix(NA, 18, 2)
for (i in 2:18){
  xlable.pos[i,1] <- i
  xlable.pos[i,2] <- Bermuda.LG.max[i-1,3]+(Bermuda.LG.max[i,2]/2)
}
xlable.pos[1,1] <- 1
xlable.pos[1,2] <- (100.439/2)



axis(side=2, at=c(ylable.pos[,2]), cex.axis=0.9,
     labels=c("1","2","3","4","5","6","7","8","9","10"), 
     las=0)

axis(side=1, at=c(xlable.pos[,2]), cex.axis=0.9,
     labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"), 
     las=0)

dev.off()


## Miscanthus
Msinensis.Chr <- read.table("M. sinensis chr length.txt", header = TRUE)

Msinensis.Chr.max = matrix(NA, 19, 3)
for (i in 1:19){
  Msinensis.Chr.max[i,1] <- i
  Msinensis.Chr.max[i,2] <- max(Msinensis.Chr[Msinensis.Chr[,1]==i,2])
  Msinensis.Chr.max[i,3] <- sum(Msinensis.Chr.max[1:i,2])
}
colnames(Bermuda.LG.max) <- c("Chr", "Chr max", "Cumulative max")


## Alignment between bermudagrass LGs and O. thomaeum genome
Bermuda.Msinensis <- matrix(NA, 3545, 5)
colnames(Bermuda.Msinensis) <- c("marker", "LG", "cum.pos", "Msinensis.Chr", "Msinensis.cum.pos")

for (i in 1:3545){
  Bermuda.Msinensis[i,1] <- Alignment[i,1]
  Bermuda.Msinensis[i,2] <- Alignment[i,3]
  Bermuda.Msinensis[i,4] <- Alignment[i,15]
  Bermuda.Msinensis[i,3] <- ifelse(Alignment[i,3]==1, Alignment[i,4],
                                   ifelse(Alignment[i,3]==2, Alignment[i,4]+Bermuda.LG.max[1,3],
                                          ifelse(Alignment[i,3]==3, Alignment[i,4]+Bermuda.LG.max[2,3],
                                                 ifelse(Alignment[i,3]==4, Alignment[i,4]+Bermuda.LG.max[3,3],
                                                        ifelse(Alignment[i,3]==5, Alignment[i,4]+Bermuda.LG.max[4,3],
                                                               ifelse(Alignment[i,3]==6, Alignment[i,4]+Bermuda.LG.max[5,3],
                                                                      ifelse(Alignment[i,3]==7, Alignment[i,4]+Bermuda.LG.max[6,3],
                                                                             ifelse(Alignment[i,3]==8, Alignment[i,4]+Bermuda.LG.max[7,3],
                                                                                    ifelse(Alignment[i,3]==9, Alignment[i,4]+Bermuda.LG.max[8,3],
                                                                                           ifelse(Alignment[i,3]==10, Alignment[i,4]+Bermuda.LG.max[9,3],
                                                                                                  ifelse(Alignment[i,3]==11, Alignment[i,4]+Bermuda.LG.max[10,3],
                                                                                                         ifelse(Alignment[i,3]==12, Alignment[i,4]+Bermuda.LG.max[11,3],
                                                                                                                ifelse(Alignment[i,3]==13, Alignment[i,4]+Bermuda.LG.max[12,3],
                                                                                                                       ifelse(Alignment[i,3]==14, Alignment[i,4]+Bermuda.LG.max[13,3],
                                                                                                                              ifelse(Alignment[i,3]==15, Alignment[i,4]+Bermuda.LG.max[14,3],
                                                                                                                                     ifelse(Alignment[i,3]==16, Alignment[i,4]+Bermuda.LG.max[15,3],
                                                                                                                                            ifelse(Alignment[i,3]==17, Alignment[i,4]+Bermuda.LG.max[16,3], Alignment[i,4]+Bermuda.LG.max[17,3])))))))))))))))))
  Bermuda.Msinensis[i,5] <- ifelse(Alignment[i,15]=="Chr01", Alignment[i,16],
                                   ifelse(Alignment[i,15]=="Chr02", Alignment[i,16]+Msinensis.Chr.max[1,3],
                                          ifelse(Alignment[i,15]=="Chr03", Alignment[i,16]+Msinensis.Chr.max[2,3],
                                                 ifelse(Alignment[i,15]=="Chr04", Alignment[i,16]+Msinensis.Chr.max[3,3],
                                                        ifelse(Alignment[i,15]=="Chr05", Alignment[i,16]+Msinensis.Chr.max[4,3],
                                                               ifelse(Alignment[i,15]=="Chr06", Alignment[i,16]+Msinensis.Chr.max[5,3],
                                                                      ifelse(Alignment[i,15]=="Chr07", Alignment[i,16]+Msinensis.Chr.max[6,3],
                                                                             ifelse(Alignment[i,15]=="Chr08", Alignment[i,16]+Msinensis.Chr.max[7,3],
                                                                                    ifelse(Alignment[i,15]=="Chr09", Alignment[i,16]+Msinensis.Chr.max[8,3],
                                                                                           ifelse(Alignment[i,15]=="Chr10", Alignment[i,16]+Msinensis.Chr.max[9,3],
                                                                                                  ifelse(Alignment[i,15]=="Chr11", Alignment[i,16]+Msinensis.Chr.max[10,3],
                                                                                                         ifelse(Alignment[i,15]=="Chr12", Alignment[i,16]+Msinensis.Chr.max[11,3],
                                                                                                                ifelse(Alignment[i,15]=="Chr13", Alignment[i,16]+Msinensis.Chr.max[12,3],
                                                                                                                       ifelse(Alignment[i,15]=="Chr14", Alignment[i,16]+Msinensis.Chr.max[13,3],
                                                                                                                              ifelse(Alignment[i,15]=="Chr15", Alignment[i,16]+Msinensis.Chr.max[14,3],
                                                                                                                                     ifelse(Alignment[i,15]=="Chr16", Alignment[i,16]+Msinensis.Chr.max[15,3],
                                                                                                                                            ifelse(Alignment[i,15]=="Chr17", Alignment[i,16]+Msinensis.Chr.max[16,3],
                                                                                                                                                   ifelse(Alignment[i,15]=="Chr18", Alignment[i,16]+Msinensis.Chr.max[17,3], Alignment[i,16]+Msinensis.Chr.max[18,3]))))))))))))))))))
}


tiff("Figure S4. Genomic synteny between bermudagrass and miscanthus.tiff", 
     units = "in", width = 10, height = 8, res = 300)
plot(Bermuda.Msinensis[,3], Bermuda.Msinensis[,5], pch=20, col="blue", cex=0.5, 
     xlab = italic("Cynodon dactylon")~"linkage groups", ylab = expression(italic("Miscanthus sinensis")~"chromosomes"), xaxt='n', yaxt='n',
     main = expression("Genomic synteny between"~italic("Cynodon dactylon")~"and"~italic("Miscanthus sinensis")),
     xaxs="i", yaxs="i")
abline(h=c(Msinensis.Chr.max[,3][1:18]), lty=2, col="gray50")
abline(v=c(Bermuda.LG.max[,3][1:17]), lty=2, col="gray50")

#make lables for y axis
ylable.pos = matrix(NA, 19, 2)
for (i in 2:19){
  ylable.pos[i,1] <- i
  ylable.pos[i,2] <- Msinensis.Chr.max[i-1,3]+(Msinensis.Chr.max[i,2]/2)
}
ylable.pos[1,1] <- 1
ylable.pos[1,2] <- (150798994/2)


#make lables for x axis
xlable.pos = matrix(NA, 18, 2)
for (i in 2:18){
  xlable.pos[i,1] <- i
  xlable.pos[i,2] <- Bermuda.LG.max[i-1,3]+(Bermuda.LG.max[i,2]/2)
}
xlable.pos[1,1] <- 1
xlable.pos[1,2] <- (100.439/2)



axis(side=2, at=c(ylable.pos[,2]), cex.axis=0.9,
     labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"), 
     las=0)

axis(side=1, at=c(xlable.pos[,2]), cex.axis=0.9,
     labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"), 
     las=0)
dev.off()


# Rice
Osativa.Chr <- read.table("O. sativa Chr length.txt", header = TRUE)
Osativa.Chr.max = matrix(NA, 12, 3)
for (i in 1:12){
  Osativa.Chr.max[i,1] <- i
  Osativa.Chr.max[i,2] <- max(Osativa.Chr[Osativa.Chr[,1]==i,2])
  Osativa.Chr.max[i,3] <- sum(Osativa.Chr.max[1:i,2])
}
colnames(Bermuda.LG.max) <- c("Chr", "Chr max", "Cumulative max")


## Alignment between bermudagrass LGs and O. thomaeum genome
Bermuda.Osativa <- matrix(NA, 3545, 5)
colnames(Bermuda.Osativa) <- c("marker", "LG", "cum.pos", "Osativa.Chr", "Osativa.cum.pos")

for (i in 1:3545){
  Bermuda.Osativa[i,1] <- Alignment[i,1]
  Bermuda.Osativa[i,2] <- Alignment[i,3]
  Bermuda.Osativa[i,4] <- Alignment[i,13]
  Bermuda.Osativa[i,3] <- ifelse(Alignment[i,3]==1, Alignment[i,4],
                                 ifelse(Alignment[i,3]==2, Alignment[i,4]+Bermuda.LG.max[1,3],
                                        ifelse(Alignment[i,3]==3, Alignment[i,4]+Bermuda.LG.max[2,3],
                                               ifelse(Alignment[i,3]==4, Alignment[i,4]+Bermuda.LG.max[3,3],
                                                      ifelse(Alignment[i,3]==5, Alignment[i,4]+Bermuda.LG.max[4,3],
                                                             ifelse(Alignment[i,3]==6, Alignment[i,4]+Bermuda.LG.max[5,3],
                                                                    ifelse(Alignment[i,3]==7, Alignment[i,4]+Bermuda.LG.max[6,3],
                                                                           ifelse(Alignment[i,3]==8, Alignment[i,4]+Bermuda.LG.max[7,3],
                                                                                  ifelse(Alignment[i,3]==9, Alignment[i,4]+Bermuda.LG.max[8,3],
                                                                                         ifelse(Alignment[i,3]==10, Alignment[i,4]+Bermuda.LG.max[9,3],
                                                                                                ifelse(Alignment[i,3]==11, Alignment[i,4]+Bermuda.LG.max[10,3],
                                                                                                       ifelse(Alignment[i,3]==12, Alignment[i,4]+Bermuda.LG.max[11,3],
                                                                                                              ifelse(Alignment[i,3]==13, Alignment[i,4]+Bermuda.LG.max[12,3],
                                                                                                                     ifelse(Alignment[i,3]==14, Alignment[i,4]+Bermuda.LG.max[13,3],
                                                                                                                            ifelse(Alignment[i,3]==15, Alignment[i,4]+Bermuda.LG.max[14,3],
                                                                                                                                   ifelse(Alignment[i,3]==16, Alignment[i,4]+Bermuda.LG.max[15,3],
                                                                                                                                          ifelse(Alignment[i,3]==17, Alignment[i,4]+Bermuda.LG.max[16,3], Alignment[i,4]+Bermuda.LG.max[17,3])))))))))))))))))
  Bermuda.Osativa[i,5] <- ifelse(Alignment[i,13]=="Chr1", Alignment[i,14],
                                 ifelse(Alignment[i,13]=="Chr2", Alignment[i,14]+Osativa.Chr.max[1,3],
                                        ifelse(Alignment[i,13]=="Chr3", Alignment[i,14]+Osativa.Chr.max[2,3],
                                               ifelse(Alignment[i,13]=="Chr4", Alignment[i,14]+Osativa.Chr.max[3,3],
                                                      ifelse(Alignment[i,13]=="Chr5", Alignment[i,14]+Osativa.Chr.max[4,3],
                                                             ifelse(Alignment[i,13]=="Chr6", Alignment[i,14]+Osativa.Chr.max[5,3],
                                                                    ifelse(Alignment[i,13]=="Chr7", Alignment[i,14]+Osativa.Chr.max[6,3],
                                                                           ifelse(Alignment[i,13]=="Chr8", Alignment[i,14]+Osativa.Chr.max[7,3],
                                                                                  ifelse(Alignment[i,13]=="Chr9", Alignment[i,14]+Osativa.Chr.max[8,3],
                                                                                         ifelse(Alignment[i,13]=="Chr10", Alignment[i,14]+Osativa.Chr.max[9,3],
                                                                                                ifelse(Alignment[i,13]=="Chr11", Alignment[i,14]+Osativa.Chr.max[10,3], Alignment[i,14]+Osativa.Chr.max[11,3])))))))))))
}



# circos plot for os6_os9_os6 (LG3_4), os2_os10_os2 (LG1_2), os5_os1_os5 (LG9_10)
os6_os9_os6 = Alignment[Alignment$Hongxu_Group == 1 | Alignment$Hongxu_Group == 2 |
                          Alignment$Hongxu_Group == 3 | Alignment$Hongxu_Group == 4,
                          #Alignment$Hongxu_Group == 9 | Alignment$Hongxu_Group == 10, 
                        c(3,4,13,14)]

C.dactylon_LG3_4 = os6_os9_os6[,1:2]
colnames(C.dactylon_LG3_4) <- c("chr", "pos")
O.sative_chr6_9_6 = os6_os9_os6[!is.na(os6_os9_os6$O..sativa),3:4]
O.sative_chr6_9_6 = O.sative_chr6_9_6[O.sative_chr6_9_6$O..sativa=="Chr6" | O.sative_chr6_9_6$O..sativa=="Chr9" |
                                        O.sative_chr6_9_6$O..sativa=="Chr2" | O.sative_chr6_9_6$O..sativa=="Chr10",]
O.sativa=data.frame("chr" = c("Chr2", "Chr2","Chr6", "Chr6","Chr9", "Chr9","Chr10","Chr10"),
                    "pos" = c(0,35937250,0,31248787,0,23012720,0, 23207287))

O.sativa[,2] <- O.sativa[,2]/500000
df <- rbind(C.dactylon_LG3_4, O.sativa)




tiff("Figure 4. Comparative mapping between bermudagrass and rice.tiff", 
     units = "in", width = 8.5, height = 11.5, res = 300)

par(mfrow=c(2,1), mar=c(4.1, 4.1, 1.1, 2.1))

plot(Bermuda.Osativa[,3], Bermuda.Osativa[,5], pch=20, col="blue", cex=0.5, 
     xlab = italic("Cynodon dactylon")~"LGs", ylab = expression(italic("Oryza sativa")~"chromosomes"), xaxt='n', yaxt='n',
     #main = expression("Genomic synteny between"~italic("Cynodon dactylon")~"and"~italic("Sativa japonica")),
     xaxs="i")
fig_label("A", region = "figure", pos = "topleft", cex = 1.5)

abline(h=c(Osativa.Chr.max[,3][1:11]), lty=2, col="gray50")
abline(v=c(Bermuda.LG.max[,3][1:17]), lty=2, col="gray50")

#make lables for y axis
ylable.pos = matrix(NA, 12, 2)
for (i in 2:12){
  ylable.pos[i,1] <- i
  ylable.pos[i,2] <- Osativa.Chr.max[i-1,3]+(Osativa.Chr.max[i,2]/2)
}
ylable.pos[1,1] <- 1
ylable.pos[1,2] <- (43270923/2)


#make lables for x axis
xlable.pos = matrix(NA, 18, 2)
for (i in 2:18){
  xlable.pos[i,1] <- i
  xlable.pos[i,2] <- Bermuda.LG.max[i-1,3]+(Bermuda.LG.max[i,2]/2)
}
xlable.pos[1,1] <- 1
xlable.pos[1,2] <- (100.439/2)



axis(side=2, at=c(ylable.pos[,2]), cex.axis=0.9,
     labels=c("1","2","3","4","5","6","7","8","9","10","11","12"), 
     las=0)

axis(side=1, at=c(xlable.pos[,2]), cex.axis=0.9,
     labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"), 
     las=0)



circos.par("track.height" = 0.1) #10% of circle radius
circos.initialize(factors = df$chr, x = df$pos)
circos.track(factors = df$chr, y = runif(nrow(df)), 
             bg.col=c(rep("grey",4), "cyan", "green", "orange", "yellow", "pink", "blue"),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(3, "mm"), 
                           CELL_META$sector.index, cex = 1)
               
             })

#Lines
O.sative_chr6_9_6_lines <- os6_os9_os6[!is.na(os6_os9_os6$O..sativa),]
O.sative_chr6_9_6_lines <- O.sative_chr6_9_6_lines[O.sative_chr6_9_6_lines$O..sativa=="Chr6"|
                                                     O.sative_chr6_9_6_lines$O..sativa=="Chr9" |
                                                     O.sative_chr6_9_6_lines$O..sativa=="Chr2" |
                                                     O.sative_chr6_9_6_lines$O..sativa=="Chr10",]

bed1 <- O.sative_chr6_9_6_lines[,1:2]
colnames(bed1) <- c("chr", "start")
bed1$end <- bed1$start
bed1$value1 <- runif(nrow(bed1))

bed2 <- O.sative_chr6_9_6_lines[,3:4]
colnames(bed2) <- c("chr", "start")
bed2$start <- bed2$start/500000
bed2$end <- bed2$start
bed2$value1 <- runif(nrow(bed2))

circos.genomicLink(bed1, bed2, 
                   col = ifelse(bed2$chr=="Chr6", "orange", 
                                ifelse(bed2$chr=="Chr9", "yellow", 
                                       ifelse(bed2$chr=="Chr2", "green","cyan"))),
                   lwd = 4,
                   border = NA)

fig_label("B", region = "figure", pos = "topleft", cex = 1.5)

dev.off()






# Oropetium chr 10:
tiff("Figure 5. Oropetium chr10 LGs.tiff", 
     units = "in", width = 8, height = 8, res = 300)
Oropetium.Chr.map <- read.csv("Chr10 markers map.csv", header = TRUE)
# circos plot
C.dactylon_LG = Oropetium.Chr.map[,1:2]
colnames(C.dactylon_LG) <- c("chr", "pos")
Oropetium_chr10 = Oropetium.Chr.map[Oropetium.Chr.map$Oropetium.Chr=="Chr10",]
Oropetium = data.frame("chr" = c("Chr10", "Chr10"),
                    "pos" = c(0,11017231))

Oropetium[,2] <- Oropetium[,2]/50000
df <- rbind(C.dactylon_LG, Oropetium)


circos.par("track.height" = 0.1) #10% of circle radius
circos.initialize(factors = df$chr, x = df$pos)
circos.track(factors = df$chr, y = runif(nrow(df)), 
             bg.col=c("cyan","blue","orange","pink", "grey"),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(3, "mm"), 
                           CELL_META$sector.index, cex = 1)
               
             })

#Lines
bed1 <- Oropetium_chr10[,1:2]
colnames(bed1) <- c("chr", "start")
bed1$end <- bed1$start
bed1$value1 <- runif(nrow(bed1))

bed2 <- Oropetium_chr10[,3:4]
colnames(bed2) <- c("chr", "start")
bed2$start <- bed2$start/50000
bed2$end <- bed2$start
bed2$value1 <- runif(nrow(bed2))

circos.genomicLink(bed1, bed2, 
                   col = ifelse(bed1$chr=="3_1", "cyan", 
                                ifelse(bed1$chr=="4_1", "blue",
                                       ifelse(bed1$chr=="5_1", "orange", "pink"))),
                   lwd = 2,
                   border = NA)

dev.off()


