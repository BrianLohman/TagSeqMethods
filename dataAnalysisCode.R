### TagSeq vs NEB Next Data Analysis

setwd("c:/Users/Lohman/Documents/GitHub/TagSeq-Improvements-and-benchmarking-paper")
library(ggplot2)
source('summarySE.R')
library(viridis)
source("RedingPlot.R")

counts <- read.csv("ERCCTruTagCounts.csv")
head(counts)

counts <- na.omit(counts)

logs <- log(counts[,c(4:ncol(counts))]+1,10) #log base 10 of seq counts +1
logs$ercc_ref <- log(counts$ercc_ref,10) #log base 10 of ercc concentrations
table(is.na(logs))

min(logs$ercc_ref)
logs$bin <- cut(logs$ercc_ref,seq(min(logs$ercc_ref)-0.1,max(logs$ercc_ref)+0.1,(max(logs$ercc_ref)-min(logs$ercc_ref))/4)) #divides data into 4 abundance classes and adds last column to tell which bin it's in
logs$bin[1] <- logs$bin[2] #small correction to remove NA
logs$bin <- as.numeric(logs$bin) #converts the bin range to bin number
str(logs)

#MS Figure 1
windows()
par(mai = c(0.4,0.3,0.3,0.3))
adj_r <- c()
res <- c() #makes empty vector
par(mfrow=c(2,4)) #plotting window
n <- "Gosling_1_TagSeq" #sets first sample followed by loop for plotting others
for (i in 1:8){
  n <- names(logs)[i]
  if (i<5) { col1 <- viridis(3)[2] } else { col1 <- viridis(3)[1]}
  col <- paste(substr(col1, 1, 7), "FE", sep = "")
  plot(counts[,n] ~ counts[,"ercc_ref"],log="xy",pch=16,col=col,cex=1.3,main=n,xlab="Spike Concentration",ylab="Count",mgp=c(2.3,1,0), cex.lab = 1.5, ylim = c(1,1000000))
  sub <- logs[logs[,n]>1,]
  l <- lm(sub[,n]~sub$ercc_ref) #linear model which regresses observed counts onto expected
  l <- lm(logs[,n]~logs$ercc_ref)	#same linear model as above
  abline(l,lwd = 1.5) #but plotted this time
  temp_r <- summary(l)
  adj_r <- c(adj_r, temp_r[9])
  res <- cbind(res,l$residuals)
}

adj_r #a vector which contains the adjusted R squared for each sample

adj_r <- unlist(as.vector(adj_r))

tag <- adj_r[1:4]
tru <- adj_r[5:8]

mean(tag)
mean(tru)

tResult1 <- t.test(tag,tru, paired = TRUE)
tResult1

#Plot to show difference in adjusted R squared between tru and tag
windows()
simple(list(tru,tag), ylim = c(0.77,0.9), ylab = "Rho", xlab = "Library Construction Method", lab = c("TotalRNAseq", "TagSeq"))

#----------------------------------------------------------
#calculate the spearman rho by library construction method
tag_result1 <- c()
for (i in 1:4){
  redData <- logs[,i]
  temp <- cor.test(redData, logs$ercc_ref, method = "spearman")
  tag_result1 <- c(tag_result1, temp$estimate)
}
tag_result1
mean(tag_result1)

tru_result1 <- c()
for (i in 5:8){
  redData <- logs[,i]
  temp <- cor.test(redData, logs$ercc_ref, method = "spearman")
  tru_result1 <- c(tru_result1, temp$estimate)
}
tru_result1
mean(tru_result1)

tResult2 <- t.test(tag_result1, tru_result1, paired = TRUE)
tResult2

#plotting results
#MS Figure 2
par(mfrow=c(1,1))
simple(list(tru_result1,tag_result1), ylim = c(0.8,0.95), ylab = "Rho", xlab = "Library Construction Method", lab = c("TotalRNAseq", "TagSeq"))

#-----------------------------------------------------------
#calculate spearman rho by rank abundance class
#separate out the tag and tru samples
tag_logs <- logs[,c(1:4,9,10)] 
tru_logs <- logs[,c(5:8,9,10)]

tag_result <- c()
for (i in 1:4){ #this loops through the abundance class
  redData <- tag_logs[tag_logs$bin == i ,]
  for (j in 1:4){ #this loops through samples
    temp <- cor.test(redData[,j], redData$ercc_ref, method = "spearman")
    tag_result <- c(tag_result, temp$estimate)  
  }
}
tag_result #first 4 elements are lowest abundance class

tru_result <- c()
for (i in 1:4){ #this loops through the abundance class
  redData <- tru_logs[tru_logs$bin == i ,]
  for (j in 1:4){ #this loops through samples
    temp <- cor.test(redData[,j], redData$ercc_ref, method = "spearman")
    tru_result <- c(tru_result, temp$estimate)  
  }
}
tru_result #first 4 elements are lowest abundance class

bin <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
all_results <- cbind(tag_result, tru_result, bin)
all_results

#plotting rho for each method by abundance class
#MS Figure 3
windows()
bar(list(all_results[1:4,2], all_results[1:4,1],
         all_results[5:8,2], all_results[5:8,1],
         all_results[9:12,2], all_results[9:12,1],
         all_results[13:16,2], all_results[13:16,1]
         ), jitter =F, SE= T,sample_size=F, bar_color=rep(viridis(3)[1:2],2), y_limits = c(-0.25,1), ylab = "Rho", xlab = "Expression Quartile")

FourthQuartile <- t.test(all_results[13:16, 1], all_results[13:16, 2], paired = TRUE)

#--------------------------------------------------------------
#Analysis of stickleback transcripts
#How similar are technial replicates of the same stickleback tissue?
tagCounts <- read.csv("stickleTagCounts.csv")
head(tagCounts)

rep1 <- cor.test(tagCounts$X1HA, tagCounts$X1HB, method = "spearman")
rep2 <- cor.test(tagCounts$X2HA, tagCounts$X2HB, method = "spearman")
rep3 <- cor.test(tagCounts$X3HA, tagCounts$X3HB, method = "spearman")
rep6 <- cor.test(tagCounts$X6HA, tagCounts$X6HB, method = "spearman")

rep7a <- cor.test(tagCounts$X8HA, tagCounts$X8HB, method = "spearman")
#rep7b <- cor.test(tagCounts$X8HA, tagCounts$X8HC, method = "spearman")
#rep7c <- cor.test(tagCounts$X8HB, tagCounts$X8HC, method = "spearman")

replicateMeanRho <- mean(rep1$estimate, rep2$estimate, rep3$estimate, rep6$estimate, rep7a$estimate)
replicateMeanRho


#----------------------------------------------------------------
#Spearman rank of stickleback transcripts, what is the correlation between tru and tag?
#This correlation is limited to the biolgoical samples that were prepped with both methods
tagCounts <- read.csv("stickleTagCounts.csv", header = TRUE, row.names = 1)
head(tagCounts)

truCounts <- read.csv("stickleTruCounts.csv", header = TRUE, row.names = 1)
head(truCounts)
#remove rows in truCounts where every entry is 0
truCounts <- truCounts[apply(truCounts[,-1], 1, function(x) !all(x==0)),]

#merge data frames to ensure that row names are the same between data sets
allCounts <- merge(tagCounts, truCounts, by = "row.names", all.y = TRUE)
dim(allCounts)

#Average technical replicats of TagSeq counts
tag1 <- (allCounts$X1HA + allCounts$X1HB)/2
tag3 <- (allCounts$X3HA + allCounts$X3HB)/2
tag6 <- (allCounts$X6HA + allCounts$X6HB)/2
tag7 <- allCounts$X7H

tagMeanCounts <- cbind(tag1, tag3, tag6, tag7)
head(tagMeanCounts)
dim(tagMeanCounts)

TagTruStickleRho <- cor.test(as.matrix(tagMeanCounts), as.matrix(allCounts[, 14:17]), method = "spearman")
TagTruStickleRho

#----------------------------------------------------------------
#Diversity of transcripts, tru vs tag
#What fraction of tru seq counts are recovered by tag seq?
tagCounts <- read.csv("stickleTagCounts.csv", header = TRUE, row.names = 1)
head(tagCounts)
truCounts <- read.csv("stickleTruCounts.csv", header = TRUE, row.names = 1)
head(truCounts)
#remove rows in truCounts where every entry is 0
truCounts <- truCounts[apply(truCounts[,-1], 1, function(x) !all(x==0)),]

dim(truCounts) #20678     4
dim(tagCounts) #19145    29

fraction <- nrow(tagCounts)/nrow(truCounts)
fraction #0.9258632

#----------------------------------------------------------------
#Barcode analyiss.
library(stringi)
barcodes <- read.csv("trimBarcodes.csv")
head(barcodes)

probability <- c(rep((1/64),nrow(barcodes)))
length(probability)

barcodes <- cbind(barcodes, probability)
head(barcodes)

gcContent <- stri_count_regex(barcodes$Barcode, c("G|C"))
gcContent <- gcContent/4

barcodes <- cbind(barcodes, gcContent)
head(barcodes)

write.csv(barcodes, "InlineBarcodeTesting.csv")

barcodes <- read.csv("InlineBarcodeTesting.csv", header = T)

result <- glm(Count~ gcContent, data = barcodes, family = poisson)
summary(result)
exp(-0.1222944) #0.8848878, a 1 unit increase in GC content reduces counts by a factor of that number
(1-0.8848878)/4 #the expted reduction in reads with the addition of a single G or C

#is the incorporation of barcodes random? No. 
chiresult <- chisq.test(barcodes[,1])
chiresult



