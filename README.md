# MethodsMissingCovariateData
Code released with the paper "Methods for dealing with missing covariate data in EWAS studies"

Code is provided for the Random Bins and Wu Bins method only and uses the Tsaprouni dataset downloadable from the NCBI Gene Expression Omnibus (http://www.ncbi.nlm.nih.gov/geo/) accession number GSE50660.

## File description
### SimulationStudyExploringBias.R
Simple code to run a simulation study to explore the effects of including or excluding variables from the imputation model on the coefficients from the regression model.

### Clean_DataPreparation.R
Download data from Tsaprouni paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4623553/ and make some missing.
This data is used to illustrate the methods recommended in the paper

Steps:
1.	download data from Tsaprouni from the NCBI Gene Expression Omnibus (http://www.ncbi.nlm.nih.gov/geo/) accession number GSE50660
2.	read in samples data from Tsaprouni, make some samples have missing data and save modified samples data
3.	read in methylation data, standardise for each CpG and save it

### Clean_CompleteCase_EWAS.R
Run a complete case analysis on the data (for use in the WuBins method)

### Clean_RandomBins_MakeBins.R
Divide CpG sites into the random bins
Bin size needs to be set between 10/1 and 3/1 ratio of samples to variables

### Clean_RandomBins_ImputeAndEWAS.R
Run the EWAS for each bin
1.	Read in random bins (created using Clean_RandomBins_MakeBins.R)
2.	Load data and samples information
3.	Loop over the bins to impute for each, and use that data for the EWAS

### Clean_WuBins_MakeBins.R
Divide CpG sites into the WuBins

Steps:
1.	Load the results of the CC analysis
2.	Pull out the top significant CpG sites
3.	Select model with forward stepwise selection with BIC
4.	Create bins always including these selected sites

### Clean_WuBins_ImputeAndEWAS.R
Run theEWAS on the WuBins

Steps:
1.	Read in random bins (created using Clean_WuBins_MakeBins.R)
2.	Load data and samples information
3.	Loop over the bins to impute for each, and use that data for the EWAS for all except the selected CpGs
4.	Just once - run the imputation and EWAS for the selected CpGs only

