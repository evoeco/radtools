#!/usr/bin/Rscript

#Supplemental Script 1. R script for the simulation of the DBR success

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  write("Use as follows:", stderr())
  write("DBR_simulation.R <number of DBRs> <targe coverage> <number of permutations> <output png file>", stderr())
  quit(status=1)
}


#n = Number of small fragments (DBRs) to be ligated to genomic libraries
#n<-122888
n<-as.integer(args[1])

#C = Target coverage for the RAD libraries 
#C<-20
C<-as.integer(args[2])

#perm = Number of different permutations used for the calculation
#perm<-200
perm<-as.integer(args[3])

#Output filename
png<-args[4]

#Distribution of the different recovered DBRs
distribution<-data.frame()

#Calculation of the distribution
for (i in 1:perm){
  
  #Calculation of a single drawing
  drawing<-list()
  
  for (j in 1:C){
    drawing[j]<-sample(1:n,1)
  }
  distribution[i,1]<-length(unique(drawing))
  rm(drawing)
  cat(i, "\r", file = stderr())
}
cat("\n", file = stderr())

png(paste(png,".png",sep=""))

#Result tables 
sor <- rle(sort(distribution$V1))
result<-data.frame(number=sor$values, n=sor$lengths)

#Proportion of the different recovered DBRs
proportion<-data.frame(number=sor$values, prop=(sor$lengths/perm))

#Plot of the probability distribution to gain i different DBRs
plot(proportion, type="o", xlab="number of recovered DBRs", ylab="probability", font=1, 
     xlim=c(1,C) , ylim=c(0,max(proportion$prop)+0.05))
Cov = "Coverage = "
DBRs = ", n = "
text=paste(Cov, C, DBRs, n, sep="")
title(text)

#Probability to recover at least i different DBRs if n different small fragments are used
cumfreq = data.frame(number=sor$values, n=cumsum((sor$lengths/perm)))
cumfreq

#Plot of the proportion to recover at least i different DBRs if n different small fragments are used
plot(cumfreq, type="o", xlab="number of recovered DBRs", ylab="cumulative probability", font=1, xlim=c(1,C))
Cov = "Coverage = "
DBRs= ", n = "
text=paste(Cov, C, DBRs, n, sep="")
title(text)

