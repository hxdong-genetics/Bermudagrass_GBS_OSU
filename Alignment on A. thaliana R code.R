#### marker alignment
setwd("~/OneDrive - University of Georgia/Bermudagrass/Comparative_Genomics")
# read sam file generated by Bowtie, aligning HapMap fasta to Sorghum
mySam <- readLines("Arabidopsis_output.sam")
mySam <- mySam[-(1:9)] # get rid of header rows

# find which rows have successful alignments
aligned <- c(grep("Chr",mySam))
length(aligned)  # 3100
mySam[aligned[1:20]]
alignedC <- grep("Chr",mySam)
length(alignedC) # 3100

# convert into table
temp <- strsplit(mySam[alignedC], "\t")
postable <- data.frame(marker=sapply(temp,
                                     function(x) x[1]),
                       reference=sapply(temp, function(x) x[3]),
                       position=as.integer(sapply(temp, function(x) x[4])))
dim(postable)
postable[1:20,]
write.csv(postable, "Aligned 3100 tags on Arabidopsis genome.csv", row.names = FALSE)
