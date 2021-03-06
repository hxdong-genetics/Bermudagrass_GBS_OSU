#### marker alignment
setwd("~/OneDrive - University of Georgia/Bermudagrass/Comparative_Genomics")
# read sam file generated by Bowtie, aligning HapMap fasta to Sorghum
mySam <- readLines("Othomaeum_output.sam")
mySam <- mySam[-(1:224)] # get rid of header rows

# find which rows have successful alignments
aligned <- c(grep("gi",mySam),grep("lcl",mySam))
length(aligned)  # 53775
mySam[aligned[1:20]]
alignedC <- grep("gi",mySam)
length(alignedC) # 51421

# convert into table
temp <- strsplit(mySam[alignedC], "\t")
postable <- data.frame(marker=sapply(temp,
                                     function(x) x[1]),
                       reference=sapply(temp, function(x) x[3]),
                       position=as.integer(sapply(temp, function(x) x[4])))
dim(postable)
postable[1:20,]
write.csv(postable, "Aligned 51421 tags on O. thomaeum genome.csv", row.names = FALSE)
