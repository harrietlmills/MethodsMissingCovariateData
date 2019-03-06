################### clean code to release for dividing CpG sites into the random bins
#### Author: Harriet L Mills
## Steps:
# choose size of bins
# state number of CpG sites
# divide into bins and save

rm(list = ls())

BIN_SIZE = 95 # between 10/1 or 3/1 ratio of samples to variables - can be a vector if want to create sets of different bin sizes

N_CpGs = 482739 # how many CpG sites looking at

for (binsize in BIN_SIZE){
  print(binsize)
  
  bins = matrix((1:N_CpGs)[order(runif(N_CpGs))], nrow=binsize) # we order the CpG sites randomly, and bin into binsize groups
  # will give warning about multiple of rows, so we correct for that with NAs in the last bin
  lastcol = N_CpGs%%binsize
  if (lastcol!=0){
    bins[(lastcol+1):binsize, dim(bins)[2]] = NA # we don't divide evenly by binsize, so we correct for this
  }
  
  saveRDS(bins, file=paste("RandomBins_size", binsize, ".rds", sep=""))
  
}
