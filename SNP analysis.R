setwd("~/OneDrive - University of Georgia/Bermudagrass")
Bermuda = read.table("Bermuda_hapmap.hmp.txt", header=T,sep="\t",comment.char="",stringsAsFactors=F )
Bermuda_Numeric=read.table("Bermuda_Numeric.txt", header = TRUE)
Bermuda_Numeric=t(Bermuda_Numeric)

#load(".RData")
# import genotype depth file
Depth=read.table("bermuda_depth.txt",header = T)
# Look at average deapth per genotype
table(unlist(Depth)>100) # 1430142 data points below 100x, 15671 data points above 100x, good, not many paralogs.
mean(unlist(Depth))
hist(unlist(Depth), breaks = 100, xlab="Read coverage per genotype", main="")

# Look at percentage of genotypes with read depth less than 4.
Depth.less.than.4.rate = apply(Depth[,-c(1,2)],
                               2, function(x) sum(x <= 4, na.rm=TRUE)/sum(!is.na(x)))
hist(Depth.less.than.4.rate)
table(Depth.less.than.4.rate<0.4) # Retain 172 inds, discard 20 inds.
table(Depth.less.than.4.rate<0.3) # Retain 191 inds, discard 1 inds.
# Convert genotypes that have read depth less than 5 to NAs.
Count.new=Depth[,-c(1:2)]
Bermuda_Numeric[Count.new<5]=NA
Bermuda_Num=cbind(Bermuda[,1:4],Bermuda_Numeric)
colnames(Bermuda_Num)


missing_by_SNP=apply(Bermuda_Num[,-c(1:4)], 1, function(x) sum(is.na(x))/131)
hist(missing_by_SNP, breaks=100, xlab = "Proportion of individuals with missing data", ylab = "No. of SNPs", main = "")
abline(v=0.2, lty=2, col="blue")
table(missing_by_SNP<=0.1) # 39013 SNPs with missing data less than 10%; 193153 with SNPS more than 10%

#Subset SNPs with missing data less than 10%
Bermuda_Num2=Bermuda_Num[missing_by_SNP<=0.1,]

missing_by_inds=apply(Bermuda_Num2[,-c(1:4)], 2, function(x) sum(is.na(x))/39013)
hist(missing_by_inds, breaks=100, xlab = "Proportion of SNPs with missing data", ylab = "No. of Individuals", main = "")



############################################################################################
# Look at allele freq across 66686 SNPs:
# Create my own function to calculate alternate (alt) allele freq
#Chi-square test for segregation distortion under 1:2:1 assumption
chisq_pval_intercross=c()
for (i in 1:nrow(Bermuda_Num2)){
  int_homo_ref_sum = sum(Bermuda_Num2[i,-c(1:4)]==0, na.rm=T)
  int_hetero_sum = sum(Bermuda_Num2[i,-c(1:4)]==0.5, na.rm=T)
  int_homo_alt_sum = sum(Bermuda_Num2[i,-c(1:4)]==1, na.rm=T)
  observed_geno_freqs = c(c(int_homo_ref_sum, int_hetero_sum, int_homo_alt_sum))
  chisq_pval_intercross[i] = chisq.test(observed_geno_freqs, p=c(0.25,0.5,0.25))$p.value
}
hist(chisq_pval_intercross, breaks=50, xlab = "P value", main = "")
abline(v=0.01, lty=2, col="blue")
table(chisq_pval_intercross > 0.01) # 7443 SNPs non-distorted at alpha=0.01, 31570 distorted SNPs
table(chisq_pval_intercross > 0.05) # 6674 SNPs non-distorted at alpha=0.05, 32339 distorted SNPs
table(chisq_pval_intercross > 0.001) # 7938 SNPs non-distorted at alpha=0.001
table(chisq_pval_intercross > 0.0001) # 8239 SNPs non-distorted at alpha=0.0001


# Subset mendelian segregating SNPs
Bermuda_Num3=Bermuda_Num2[chisq_pval_intercross > 0.0001,]
write.csv(Bermuda_Num3, "SNP dataset1 6674 SNPs x 131 inds (depth 6x, chi-square p=0.05).csv")

# Subset mendelian segregating SNPs
Bermuda_Num4=Bermuda_Num2[chisq_pval_intercross <= 0.0001,]
write.csv(Bermuda_Num4, "30774 distorted SNPs x 131 inds (depth 6x, chi-square p=0.0001).csv")


Bermuda_SNPs <- Bermuda_Num3[,-c(3,4)]
Bermuda_SNPs[Bermuda_SNPs==0.5]='hk'
Bermuda_SNPs[Bermuda_SNPs==0]='hh'
Bermuda_SNPs[Bermuda_SNPs==1]='kk'
Bermuda_SNPs[is.na(Bermuda_SNPs)]='--'

write.csv(Bermuda_SNPs, "Bermudagrass SNP dataset4 (8239 SNPs x 131 inds, depth 6x, chi-square p=0.0001).csv")


Bermuda_distorted_SNPs <- Bermuda_Num4[,-c(3,4)]
Bermuda_distorted_SNPs[Bermuda_distorted_SNPs==0.5]='hk'
Bermuda_distorted_SNPs[Bermuda_distorted_SNPs==0]='hh'
Bermuda_distorted_SNPs[Bermuda_distorted_SNPs==1]='kk'
Bermuda_distorted_SNPs[is.na(Bermuda_distorted_SNPs)]='--'

write.csv(Bermuda_distorted_SNPs, "30774 distorted SNPs x 131 inds (depth 6x, chi-square p=0.0001).csv")

# Subset genotypic dataset for LG 8, 17-18
dataset5.map <- read.csv("Comparative genomics database of 12947 filtered SNPs.csv", header = TRUE)
LG8.SNP.List <- dataset5.map[c(dataset5.map[,2]==9),1:2]
LG17_18.SNP.List <- dataset5.map[c(dataset5.map[,2]==8),1:2]

undistorted.SNP <- read.csv("Bermudagrass SNP dataset4 (8239 SNPs x 131 inds, depth 6x, chi-square p=0.0001).csv", header = TRUE)
distorted.SNP <- read.csv("30774 distorted SNPs x 131 inds (depth 6x, chi-square p=0.0001).csv", header = TRUE)
All.SNP <- rbind(undistorted.SNP, distorted.SNP)

LG8.SNP <- All.SNP[c(All.SNP[,1]%in%LG8.SNP.List[,1]),]
LG17_18.SNP <- All.SNP[c(All.SNP[,1]%in%LG17_18.SNP.List[,1]),]

write.csv(LG8.SNP, "LG8 SNP.csv", row.names = FALSE)
write.csv(LG17_18.SNP, "LG17_18 SNP.csv", row.names = FALSE)

###################### Extract tag seqeunces ######################
TOPM <- read.table("UNEAK.topm.txt", header = FALSE, sep = "\t")
TOPM <- TOPM[,1:8]
#
TOPM <- TOPM[c(TOPM[,6]%in%TagtoKeep[,1]),]

# Subset query seqeunces
tmp <- seq(from = 1, to = 464332, by=2)
TOPM_query <- TOPM[tmp,]
tmp <- seq(from = 2, to = 464332, by = 2)
TOPM_hit <- TOPM[tmp,]
#
TagName <- Bermuda[,1]

TOPM_query1 <- cbind(as.matrix(TagName), TOPM_query)
TOPM_hit1 <- cbind(as.matrix(TagName), TOPM_hit)
TOPM_query1 <- TOPM_query1[,1:2]
colnames(TOPM_query1) <- c("Marker name","Tag sequence 0")

TOPM_hit1 <- TOPM_hit1[,1:2]
colnames(TOPM_hit1) <- c("Marker name","Tag sequence 1")

Tags.in.columns <- cbind(TOPM_query1, TOPM_hit1[,2])
colnames(Tags.in.columns) <- c("Marker name","Tag sequence 0", "Tag sequence 1")
write.csv(Tags.in.columns, "Bermuda_232166_Tag_sequence.csv", row.names = FALSE)

##################################################################################
#OSU.Bermuda.SNP.Tag <- read.csv("Bermuda_merged_232166_Tag_sequence.csv", header = TRUE)
SNP.list <- read.csv("Bermudagrass SNP dataset4 (8239 SNPs x 131 inds, depth 6x, chi-square p=0.0001).csv",header = TRUE)
Bermuda_8239_SNPs_Tag_Sequence <- OSU.Bermuda.SNP.Tag[c(OSU.Bermuda.SNP.Tag[,3]%in%SNP.list[,1]),]
write.csv(Bermuda_8239_SNPs_Tag_Sequence, "Tag sequence of mapped 8239 SNPs.csv", row.names = FALSE)



### 05-09-2019
Final.Map <- read.csv("Comparative genomics database of 3545 mapped SNPs.csv", header = TRUE)
Final.Map <- Final.Map[,c(1,3,4)]
Final.Map.Geno <- Bermuda_Numeric[c(rownames(Bermuda_Numeric)%in%Final.Map[,1]), ]

# inter-marker spacing
inter.marker.dist=c()

for (i in 1:(nrow(Final.Map)-1)){
  inter.marker.dist[i] = Final.Map[i+1,3]-Final.Map[i,3]
}
  
inter.marker.dist[inter.marker.dist>=0]

LG.Length=c()
for (i in 1:18) {
  print(max(Final.Map[Final.Map[,2]==i,3]))
  LG.Length[i] = max(Final.Map[Final.Map[,2]==i,3])
}
sum(LG.Length)

# Look at the population principal coordinates plot:
d <- dist(Final.Map.Geno)
fit <- cmdscale(d, eig = TRUE, k=2)

# Principal coordinates plot for individuals
d1 <- dist(t(Final.Map.Geno))
fit1 <- cmdscale(d1)
load("PCoA.RData")
tiff("Figure 1. Population PCoA plot.tiff", units = "in", width = 7, height = 6, res = 300)
plot(fit1[,1], fit1[,2], pch=19, 
     xlab="Coordinate 1", ylab="Coordinate 2", col=ifelse(rownames(fit1)=="p118", "red", "blue"))
legend("bottomright", legend = c("Parent 'A12359'", "Selfed progeny"), col = c("red","blue"), pch = c(19,19), cex = 1)
dev.off()

