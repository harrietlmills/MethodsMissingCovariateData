################### clean code to release for Wu Bins method
#### Author: Harriet L Mills
## Steps:
# read in random bins (created using Clean_WuBins_MakeBins.R)
# Load data and samples information
# Loop over the bins to impute for each, and use that data for the EWAS for all except the selected CpGs
# just once - run the imputation and EWAS for the selected CpGs only
################### 

rm(list = ls())

#### load packages required
library(survey)
library(mice)

#### source useful functions
source('UsefulFunctions_ForRelease.R')

#### set details
n_imputations = 100
binsize = 95

#### load the random CpG bins as created from file "WuBins_MakeBins.R"
bins = readRDS(paste("WuBins_size", binsize, ".rds", sep=""))

#### WHILE THIS CAN BE RUN ALL AT ONCE (USING bin_index = 1:dim(bins)[2]), IT WILL TAKE A LONG TIME SO IT IS BETTER 
#### TO SET THIS UP AS A TEMPLATE TO DIVIDE THE BINS INTO BATCHES WHICH CAN BE PROCESSED SEPARATELY
## set details for this script (use a .sh script to automate creating copies of this R file with the index "ba=" filled in)
ba =  #determine index set within the iteration
## 
batch_size = 50
bin_index = ((ba-1)*batch_size+1):(ba*batch_size)
if (max(bin_index)>dim(bins)[2]){bin_index = ((ba-1)*batch_size+1):dim(bins)[2] }

#### load the samples data
samples = readRDS(file="samplesdata.rds")
table(samples$smoking, useNA="ifany")
# one row for every sample, with covariates as columns

#### load the methylation data
methylationdata <- readRDS(file="methylationdata_standardised.rds")
# should contain one row for each CpG site, one column for each sample
CpG_fulllist = rownames(methylationdata)

chosen_CpG = CpG_fulllist[intersect(bins[, 1], bins[, 2])]

#### run the imputations and the regressions
to_extract = c("est", "se", "lo 95", "hi 95", "fmi", "lambda") # to extract from the pooled fit regression
for (B in bin_index){ #### impute for each bin using all CpGs and then regress on each CpG sep
  
  # identify the CpGs in this bin
  CpG_bin = CpG_fulllist[bins[, B]]; CpG_bin = CpG_bin[!is.na(CpG_bin)]
  
  # create temporary dataframe with those CpG sites and the sample data
  data_tmp = data.frame(t(methylationdata[CpG_bin, ]),
                        age = samples$age, 
                        sex = samples$sex,
                        smoking = samples$smoking)
  colnames(data_tmp) = c(CpG_bin, "age", "sex", "smoking")
  
  # create a prediction matrix indicating which variables will be used for the imputation (see mice package notes)
  predMat = matrix(1, ncol(data_tmp), ncol(data_tmp))
  diag(predMat) <- 0
  rownames(predMat) <- colnames(data_tmp) 
  predMat[c(CpG_bin, "age", "sex"), ] = 0  #only smoking is missing, so other variables do not need imputing and can set to 0
  
  # run the imputation
  imp <- mice(data = data_tmp[, c(CpG_bin, "age", "sex", "smoking")], 
              method = c(rep("", length(CpG_bin)), "", "", "polyreg"), 
              predictorMatrix = predMat, 
              m=n_imputations, 
              print=FALSE, 
              ridge=0)
  
  # FORCE R TO CORRECT FACTORS
  imp$data$smoking <- as.factor(as.character(imp$data$smoking))
  # NB can check this: table(imp$data$smoking) & table(samples$smoking) should be the same
  
  # remove chosen_CpG from the bin for the regression (otherwise we do this many times for the chosen CpGs)
  CpG_bin = CpG_bin[!CpG_bin%in%chosen_CpG]
  
  # run the EWAS for all remaining CpG sites in the bin
  res_imputed = matrix(NA, nrow=length(CpG_bin), ncol=2*length(to_extract)+1)
  colnames(res_imputed) = c(paste(rep(c("smoking1", "smoking2"), each=length(to_extract)), to_extract, sep="_")   , "smoking_pvalue")
  rownames(res_imputed) = CpG_bin
  for (cp in 1:length(CpG_bin)){
    CpG = CpG_bin[cp]
    
    # the regression
    fit <- with(imp, lm(as.formula(paste0(CpG, " ~ age + sex + smoking"))))
    tab <- summary(pool(fit, "rubin"))
    Fvalues = sapply(1:length(fit$analyses), function(e) regTermTest(fit$analyses[[e]],"smoking")$F) ; 
    df = regTermTest(fit$analyses[[1]],"smoking")$df #degrees of freedom for numerator
    Fvalues_comb = micombine.F(Fvalues, df=df, display=FALSE); # df should be the numerator degrees of freedom (see function code, and http://www.stata.com/support/faqs/statistics/chi-squared-and-f-distributions/)
    pvalue = Fvalues_comb["p"]
    res_imputed[cp, ] = c(tab["smoking1", to_extract], tab["smoking2", to_extract], pvalue)
  }
  
  if (B==1){
    RES_IMPUTED = res_imputed
  } else{
    RES_IMPUTED = rbind(RES_IMPUTED, res_imputed) 
  }

}

#### if ba==1 we perform imputation just for the Wu selected CpG sites (i.e. we only do this once for the chosen CpGs)
if (ba==1){
  
  CpG_bin = chosen_CpG
  
  # create temporary dataframe with those CpG sites and the sample data
  data_tmp = data.frame(t(methylationdata[CpG_bin, ]),
                        age = samples$age, 
                        sex = samples$sex,
                        smoking = samples$smoking)
  colnames(data_tmp) = c(CpG_bin, "age", "sex", "smoking")
  
  # create a prediction matrix indicating which variables will be used for the imputation (see mice package notes)
  predMat = matrix(1, ncol(data_tmp), ncol(data_tmp))
  diag(predMat) <- 0
  rownames(predMat) <- colnames(data_tmp)
  predMat[c(CpG_bin, "age", "sex"), ] = 0  #only smoking is missing, so other variables do not need imputing and can set to 0
  
  # run the imputation
  imp <- mice(data = data_tmp[, c(CpG_bin, "age", "sex", "smoking")], 
              method = c(rep("", length(CpG_bin)), "", "", "polyreg"), 
              predictorMatrix = predMat, 
              m=n_imputations, 
              print=FALSE, 
              ridge=0)
  
  imp$data$smoking <- as.factor(as.character(imp$data$smoking)) # FORCE R TO CORRECT ITS WEIRD FACTOR MESS THAT IT SOMEHOW DID
  # NB can check this: table(imp$data$smoking) & table(samples$smoking) should be the same
  
  res_imputed = matrix(NA, nrow=length(CpG_bin), ncol=2*length(to_extract)+1)
  colnames(res_imputed) = c(paste(rep(c("smoking1", "smoking2"), each=length(to_extract)), to_extract, sep="_")   , "smoking_pvalue")
  rownames(res_imputed) = CpG_bin
  for (cp in 1:length(CpG_bin)){
    CpG = CpG_bin[cp]
    
    ## run the regressions for each CpG in this bin - ADJUSTED
    fit <- with(imp, lm(as.formula(paste0(CpG, " ~ age + sex + smoking"))))
    tab <- summary(pool(fit, "rubin"))
    Fvalues = sapply(1:length(fit$analyses), function(e) regTermTest(fit$analyses[[e]],"smoking")$F) ; 
    df = regTermTest(fit$analyses[[1]],"smoking")$df #degrees of freedom for numerator
    Fvalues_comb = micombine.F(Fvalues, df=df, display=FALSE); # df should be the numerator degrees of freedom (see function code, and http://www.stata.com/support/faqs/statistics/chi-squared-and-f-distributions/)
    pvalue = Fvalues_comb["p"]
    res_imputed[cp, ] = c(tab["smoking1", to_extract], tab["smoking2", to_extract], pvalue)
  }
  RES_IMPUTED = rbind(RES_IMPUTED, res_imputed) 

}

#### save the results
save(RES_IMPUTED, file=paste("Results_WuBins_size", binsize, "_part", ba, ".Rdata", sep=""))
