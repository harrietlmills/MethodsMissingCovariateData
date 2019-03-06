################### download data from Tsaprouni paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4623553/, make some missing
#### Author: Harriet L Mills
## this data is used to illustrate the methods recommended in the paper
## Steps:
# (1) download data from Tsaprouni from the NCBI Gene Expression Omnibus (http://www.ncbi.nlm.nih.gov/geo/) accession number GSE50660
# (2) read in samples data from Tsaprouni
# make some samples have missing data
# save modified samples data
# (3) read in methylation data
# standardise for each CpG
################### 

#### (1) #### download data
# do this outside of R from the NCBI Gene Expression Omnibus (http://www.ncbi.nlm.nih.gov/geo/) accession number GSE50660
# name data files as:
# 1) methylation data: Tsaprouni_methylationdata.csv
# 2) samples data: Tsaprounisamples.csv

#### (2) #### deal with samples data
##### read in sample data
samples <- read.csv("Tsaprounisamples.csv")

##### create age, sex and smoking columns from the text-entry characteristics columns
samples$characteristics_ch1 <- as.character(samples$characteristics_ch1)
samples$characteristics_ch1.1 <- as.character(samples$characteristics_ch1.1)
samples$characteristics_ch1.2 <- as.character(samples$characteristics_ch1.2)

samples$smoking <- as.numeric(sapply(1:dim(samples)[1], function(e) substr(samples[e, "characteristics_ch1"], nchar(samples[e, "characteristics_ch1"]), nchar(samples[e, "characteristics_ch1"]))))
samples$age <- as.numeric(sapply(1:dim(samples)[1], function(e) substr(samples[e, "characteristics_ch1.1"], nchar(samples[e, "characteristics_ch1.1"])-1, nchar(samples[e, "characteristics_ch1.1"]))))
samples$sex <- sapply(1:dim(samples)[1], function(e) strsplit(samples[e, "characteristics_ch1.2"], " ")[[1]][2])

##### make some samples have missing data
# missingness mechanism: Missing with probability 75% for males aged 57 years and over. (MM1 in our paper)
# select the missing
age_limit = 57
gender = "Male"
which_poss = which(samples$age >= age_limit & samples$sex%in%gender)
prob_missing = 0.75
which_missing = which_poss[which(rbinom(length(which_poss), 1, prob=prob_missing)==1)]

# check proportion missing
100*length(which_missing) / dim(samples)[1]

## make it missing
samples$smoking[which_missing] = NA

##### save modified data
saveRDS(samples, file="samplesdata.rds")

#### (3) #### deal with methylation data
##### NB - this file is very large, best to do on an HPC, or read in and process in chunks of a manageable size
##### load the methylation data
methylation_all <- read.csv("Tsaprouni_methylationdata.csv", header=TRUE, row.names=1)

##### standardise for each CpG
for (cpg in 1:dim(methylation_all)[1]){
  tmp = as.numeric(methylation_all[cpg, ])
  m_cpg = mean(tmp)
  sd_cpg = sd(tmp)
  methylation_all[cpg, ] = (tmp - m_cpg) / sd_cpg
}

##### save it
saveRDS(methylation_all, "methylationdata_standardised.rds")

