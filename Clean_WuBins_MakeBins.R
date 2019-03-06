################### clean code to release for Wu Bins method - make the WuBins
#### Author: Harriet L Mills
## Steps:
# load the CC analysis
# pull out the top significant CpG sites
# select model with forward stepwise selection with BIC
# create bins always including these selected sites
################### 

rm(list = ls())

#### load packages required
library(survey)
library(mice)
library(nnet)

#### set details
binsize = 95 # between 10/1 or 3/1 ratio of samples to variables

#### load the samples data
samples = readRDS(file="samplesdata.rds")
table(samples$smoking, useNA="ifany")
# one row for every sample, with covariates as columns

#### load the methylation dataset
methylationdata <- readRDS(file="methylationdata_standardised.rds")
# should contain one row for each CpG site, one column for each sample
CpGs = rownames(methylationdata)

#### set bonferroni corrected p value that we use for the analysis
WuLimit = 0.05/dim(methylationdata)[1] # i.e. 0.05 / number of CpG

#### load results of the complete case analysis on the missing data
load("EWAS_results_CC.Rdata") 
# loads a matrix called res_methCC, one row for every CpG site and a column called "pval" for the pvalue for the analysis at each site

### sorting the p-values column from the CC into ascending order
res_methCC = res_methCC[order(res_methCC[, "pval"]), ]  

#### select the subset of CpGs to use for the imputation
## in the Wu paper they take the CpG with p value below 10^-7 but up to 100
res_methCC_inlimit = res_methCC[res_methCC[, "pval"]<=WuLimit, ]
if (dim(res_methCC_inlimit)[1]>100){
  CpG_selected = rownames(res_methCC_inlimit)[1:100]
} else{
  CpG_selected = rownames(res_methCC_inlimit)
}
rm(res_methCC_inlimit)

#### run the stepwise selection model on the selected CpGs
data_tmp = data.frame(t(methylationdata[CpG_selected, ]),
                      age = samples$age, 
                      sex = samples$sex,
                      smoking = samples$smoking)
colnames(data_tmp) = c(CpG_selected, "age", "sex", "smoking")

# make complete case before we select the subset
print(paste0("Before CC N=", dim(data_tmp)[1]))
cols_forEWAS = c("age", "sex", "smoking") # only listing the covariates used in the EWAS
sumNA <- which(rowSums(is.na(data_tmp[, cols_forEWAS]))>0) # which rows have any missingness - using only covariates used in the EWAS
data_tmp_nonmissing <- data_tmp[-sumNA, ]
print(paste0("After CC N=", dim(data_tmp_nonmissing)[1]))

# multinomial logistic regression because smoking has 3 unordered levels
# NB we use data_tmp_nonmissing for this regression
model.init <- multinom(as.formula("smoking ~ age + sex"), data_tmp_nonmissing)
model.forward.BIC <- step(model.init, direction="forward", k=log(dim(data_tmp_nonmissing)[1]), 
                          scope = list(lower=as.formula("~age + sex"), 
                                       upper=as.formula(paste("~age + sex +", paste(CpG_selected, collapse=" + "), sep=""))), trace=-1, data=data_tmp_nonmissing)
print(model.forward.BIC)

chosen_CpG = model.forward.BIC$coefnames[grep("cg", model.forward.BIC$coefnames)]
print(paste("Chosen model for adjusted: covariates + ", paste(chosen_CpG, collapse=" + "), sep=""))

#### create bins with the remaining CpG sites
CpG_index <- which(!CpGs%in%chosen_CpG)

# calculate remaining size after the chosen_CpGs
remainingbinsize <- binsize - length(chosen_CpG)

# we order the CpG sites randomly, and bin into remainingbinsize groups
bins = matrix(CpG_index[order(runif(length(CpG_index)))], nrow=remainingbinsize) 
# will give warning about multiple of rows, so we correct for that with NAs
lastcol = length(CpG_index)%%remainingbinsize
bins[(lastcol+1):remainingbinsize, dim(bins)[2]] = NA # we don't divide evenly by remainingbinsize, so we correct for this
# add chosenCpGs to every bin
bins <- rbind(bins, matrix(which(CpGs%in%chosen_CpG), nrow=length(chosen_CpG), ncol=ncol(bins)))

saveRDS(bins, file=paste0("WuBins_size", binsize, ".rds"))



