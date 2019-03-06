################### clean code to release for CC analysis (which is then used in WuBins method)
#### Author: Harriet L Mills
## Steps:
# load the data
# run the CC analysis
################### 

#### load packages required
library(survey)

##### load samples data
samples = readRDS(file="samplesdata.rds")

##### load methylation data
methylationdata <- readRDS(file="methylationdata_standardised.rds")
# should contain one row for each CpG site, one column for each sample
CpGs = rownames(methylationdata)

##### create temporary dataframe with CpG sites and the sample data
data_tmp = data.frame(t(methylationdata),
                      age = samples$age, 
                      sex = samples$sex,
                      smoking = samples$smoking)
colnames(data_tmp) = c(CpGs, "sex", "age", "smoking")
data_tmp = data_tmp[is.na(data_tmp$smoking), ] #removed those rows which do not have smoking

##### run EWAS
# NB this is slow, probably best to run on an HPC, and consider splitting the methylation matrix into evenly sized chunks to load and process (as you never need all loaded together)
LM = lapply(1:length(CpGs), function(i) lm(as.formula(paste0(CpGs[i], " ~ sex + age + smoking")), data = data_tmp))
pval = sapply(LM, function(e) regTermTest(e,"smoking")$p) # returns only the pvalue for smoking, and not for the whole model
se = sapply(LM, function(e) coef(summary(e))[, "Std. Error"])
se_smoking1 = se["smoking1", ]
se_smoking2 = se["smoking2", ]
coefsmoking1 = sapply(LM, function(e) e$coefficients["smoking1"]) 
coefsmoking2 = sapply(LM, function(e) e$coefficients["smoking2"]) 

res_methCC = data.frame(coefsmoking1 = coefsmoking1, coefsmoking2 = coefsmoking2, se_smoking1 = se_smoking1, se_smoking2=se_smoking2, pval = pval)
rownames(res_methCC) = CpGs

##### save results  
save(res_methCC, file="EWAS_results_CC.Rdata")

