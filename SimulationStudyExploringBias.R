##### simulation of MI to demonstrate bias
### Author: Harriet L Mills harriet.mills@bristol.ac.uk
### Methods and Results described in the supplementary information, section 1

rm(list = ls())

# setwd() # wherever you wish things to be saved

library(mice) # install with install.packages("mice") first - only need to install once

samplesize <- 1000 # set the sample size

set.seed(10) # fix the random seed - can be anything, useful to fix so is reproducible

N_sims <- 1000 # number of simulations to run

N_vars <- 7 # number of variables in our model

n_imputations <- 5 # number of imputations to run

# generate matrices to save the coefficients and imputation indicators in
TrueCoef <- CCCoef <- ImpInd <- ImpCoef <- TrueCoefSE <- CCCoefSE <- ImpCoefSE <- matrix(NA, N_sims, N_vars)

for (sim in 1:N_sims){
  
  # generate the model
  model <- data.frame(matrix(NA, samplesize, 8))
  colnames(model) <- c(paste0("x", 1:N_vars), "smoke")
  model$smoke <- runif(samplesize, 0, 1)<=0.8
  model$x1=rnorm(samplesize,0,1)+0.2*model$smoke 
  model$x2=rnorm(samplesize,0,1)*4+0.2*model$smoke  
  model$x3=rnorm(samplesize,0,1)*2+0.2*model$smoke
  model$x4=rnorm(samplesize,0,1)*2+0.2*model$smoke
  model$x5=rnorm(samplesize,0,1)+0.2*model$smoke
  model$x6=rnorm(samplesize,0,1)+0.2*model$smoke
  model$x7=rnorm(samplesize,0,1)*0.5+0.2*model$smoke
  
  # find the true coefficients
  for (i in 1:N_vars){
    
    lm_mod <- lm(as.formula(paste0("x", i, "~smoke")), data=model)
    sum_lm_mod <- summary(lm_mod)
    TrueCoef[sim, i] <- sum_lm_mod$coefficients["smokeTRUE", "Estimate"]
    TrueCoefSE[sim, i] <- sum_lm_mod$coefficients["smokeTRUE", "Std. Error"]
  }
  
  # add missingness
  missing <- runif(samplesize)>=0.5
  model$smoke[missing==1] <- NA
  
  # complete case analysis, and whether we include variables in the imputation
  CC_model <- model[!is.na(model$smoke), ]
  for (i in 1:N_vars){
    
    lm_mod <- lm(as.formula(paste0("x", i, "~smoke")), data=CC_model)
    sum_lm_mod <- summary(lm_mod)
    CCCoef[sim, i] <- sum_lm_mod$coefficients["smokeTRUE", "Estimate"]
    CCCoefSE[sim, i] <- sum_lm_mod$coefficients["smokeTRUE", "Std. Error"]
    
    if (i!=1){
      ttest <- sum_lm_mod$coefficients["smokeTRUE", "Estimate"]/sum_lm_mod$coefficients["smokeTRUE", "Std. Error"]
      ImpInd[sim, i] <- (ttest>=1.96)|(ttest<=-1.96)
      # variable x1 never appears in the imputation
    }
    
  }
  ImpInd[, 1] <- 0 # variable x1 never appears in the imputation
  
  imp <- mice(data = model, m=n_imputations, 
              method=c(rep("", ncol(model)-1), "logreg"),
              predictorMatrix = matrix(c(rep(0, (ncol(model))-1*ncol(model)),
                                         ImpInd[sim, ], 0), ncol(model), ncol(model), byrow=TRUE), 
              print=FALSE)
  # note that the predictorMatrix uses ImpInd[sim, ] to determine which variables are included in this imputation
  # note that we impute smoke using logreg - logistic regression (suitable for binary)
  
  # run the lm with the imputed data
  for (i in 1:7){
    fit <- with(imp, lm(as.formula(paste0("x", i, "~smoke"))))
    tab <- summary(pool(fit, "rubin"))    
    ImpCoef[sim, i] <- tab["smoke", "est"]
    ImpCoefSE[sim, i] <- tab["smoke", "se"]
  }

}


######## Results:
### assess beta for each variable
beta_table <- matrix(NA, nrow=10, ncol=7)
rownames(beta_table) <- c("True", "TrueSE", "CC", "CCSE", "Mean bias CC", "SD bias CC", "Imp","ImpSE", "Mean bias Imp", "SD bias Imp")
colnames(beta_table) <- paste0("x", 1:7)
beta_table["True", ] <- colMeans(TrueCoef)
beta_table["CC", ] <- colMeans(CCCoef)
beta_table["Imp", ] <- colMeans(ImpCoef)
beta_table["TrueSE", ] <- colMeans(TrueCoefSE)
beta_table["CCSE", ] <- colMeans(CCCoefSE)
beta_table["ImpSE", ] <- colMeans(ImpCoefSE)
beta_table["Mean bias CC", ] <- beta_table["CC",] - beta_table["True",]
beta_table["Mean bias Imp", ] <- beta_table["Imp",] - beta_table["True",]
beta_table["SD bias CC", ] <- sapply(1:7, function(e) sd(CCCoef[,e] - TrueCoef[, e]))
beta_table["SD bias Imp", ] <- sapply(1:7, function(e) sd(ImpCoefSE[,e] - TrueCoef[, e]))
print(beta_table)

## make nice table
beta_table_nice <- rbind("True"=paste0(formatC(beta_table["True", ], digits = 4, format = "f"), " (", formatC(beta_table["TrueSE", ], digits = 4, format = "f"), ")"),
                        "CC"=paste0(formatC(beta_table["CC", ], digits = 4, format = "f"), " (", formatC(beta_table["CCSE", ], digits = 4, format = "f"), ")"),
                        "CC bias"=paste0(formatC(beta_table["Mean bias CC", ], digits = 4, format = "f"), " (", formatC(beta_table["SD bias CC", ], digits = 4, format = "f"), ")"),
                        "Imp"=paste0(formatC(beta_table["Imp", ], digits = 4, format = "f"), " (", formatC(beta_table["ImpSE", ], digits = 4, format = "f"), ")"),
                        "Imp bias"=paste0(formatC(beta_table["Mean bias Imp", ], digits = 4, format = "f"), " (", formatC(beta_table["SD bias Imp", ], digits = 4, format = "f"), ")"))
colnames(beta_table_nice) <- paste0("x", 1:7)
write.csv(beta_table_nice, "beta_table_nice.csv")

png(filename="Simulation_betaplots.png", width=90, height=90, units="mm", res=500)
par(mar=c(2.5, 2.5, 1.5, 1), cex=0.6, tcl=-0.25) #
plot(1:7, beta_table["True", ], pch=15, col="black", ylab="", xlab="", axes=FALSE, ylim=c(-0.3, 0.7)) #ylim=c(0.9*min(beta_table[c("True", "CC", "Imp"), ]), 1.1*max(beta_table[c("True", "CC", "Imp"), ]) 
points(1:7-0.1, beta_table["CC", ], pch=16, col="grey33")
points(1:7+0.1, beta_table["Imp", ], pch=17, col="grey66")
arrows(1:7, beta_table["True", ]-beta_table["TrueSE", ], 1:7, beta_table["True", ]+beta_table["TrueSE", ],  col="black", length=0.02, angle=90, code=3)
arrows(1:7-0.1, beta_table["CC", ]-beta_table["CCSE", ], 1:7-0.1, beta_table["CC", ]+beta_table["CCSE", ],  col="grey33", length=0.02, angle=90, code=3)
arrows(1:7+0.1, beta_table["Imp", ]-beta_table["ImpSE", ], 1:7+0.1, beta_table["Imp", ]+beta_table["ImpSE", ],  col="grey66", length=0.02, angle=90, code=3)
mtext(side=2, "Beta", line=1.3, cex=0.8); mtext(side=1, "Variable", line=1.2, cex=0.8)
axis(1, padj=-1.8, cex.axis=0.9, at=1:7, labels=c("x1*", paste0("x", 2:7))) 
axis(2, las=2, hadj=0.4, cex.axis=0.9); box()
legend("topright", legend=c("True", "CC", "Imp"), pch=c(15, 16,17), col=c("black", "grey33", "grey66"), bty="n")
dev.off()

### assess bias depending on whether they were included or not
bias <- matrix(NA, nrow=8, ncol=7)
rownames(bias) <- c("Imp included", "Imp included SD", "Imp not included", "Imp not included SD", "CC included", "CC included SD", "CC not included", "CC not included SD")
colnames(bias) <- paste0("x", 1:7)
for (i in 1:7){
  bias["Imp included", i] <- mean((ImpCoef[, i] - TrueCoef[, i])[ImpInd[, i]==1])
  bias["Imp included SD", i] <- sd((ImpCoef[, i] - TrueCoef[, i])[ImpInd[, i]==1])
  bias["Imp not included", i] <- mean((ImpCoef[, i] - TrueCoef[, i])[ImpInd[, i]==0])
  bias["Imp not included SD", i] <- sd((ImpCoef[, i] - TrueCoef[, i])[ImpInd[, i]==0])
  bias["CC included", i] <- mean((CCCoef[, i] - TrueCoef[, i])[ImpInd[, i]==1])
  bias["CC included SD", i] <- sd((CCCoef[, i] - TrueCoef[, i])[ImpInd[, i]==1])
  bias["CC not included", i] <- mean((CCCoef[, i] - TrueCoef[, i])[ImpInd[, i]==0])
  bias["CC not included SD", i] <- sd((CCCoef[, i] - TrueCoef[, i])[ImpInd[, i]==0])
}
print(bias)

## make nice table
bias_table_nice <- rbind("Imp included"=paste0(formatC(bias["Imp included", ], digits = 4, format = "f"), " (", formatC(bias["Imp included SD", ], digits = 4, format = "f"), ")"),
                         "Imp not included"=paste0(formatC(bias["Imp not included", ], digits = 4, format = "f"), " (", formatC(bias["Imp not included SD", ], digits = 4, format = "f"), ")"),
                         "CC included"=paste0(formatC(bias["CC included", ], digits = 4, format = "f"), " (", formatC(bias["CC included SD", ], digits = 4, format = "f"), ")"),
                         "CC not included"=paste0(formatC(bias["CC not included", ], digits = 4, format = "f"), " (", formatC(bias["CC not included SD", ], digits = 4, format = "f"), ")"))
write.csv(bias_table_nice, "bias_table_nice.csv")                         

