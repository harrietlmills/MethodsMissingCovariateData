################### taken from package miceadds https://cran.r-project.org/web/packages/miceadds/miceadds.pdf
# and adjusted using http://www.stata.com/support/faqs/statistics/chi-squared-and-f-distributions/ because I think they'd swapped the num/denom df

micombine.chisquare <-function (dk, df, display = TRUE)
{
  M <- length(dk)
  mean.dk <- mean(dk)
  sdk.square <- stats::var(sqrt(dk))
  Dval <- (mean.dk/df - (1 - 1/M) * sdk.square)/(1 + (1 + 1/M) * sdk.square)
  df2 <- (M - 1)/df^(3/M) * (1 + M/(M + 1/M)/sdk.square)^2
  pval <- stats::pf(Dval, df1 = df, df2 = df2, lower.tail = FALSE) # NB: in pf - df1 is numerator and df2 is denominator - 
  chisq.approx <- Dval * df
  p.approx <- 1 - stats::pchisq(chisq.approx, df = df)
  res <- c(D = Dval, p = pval, df = df, df2 = df2)
  if (display) {
    cat("Combination of Chi Square Statistics for Multiply Imputed Data\n")
    cat(paste("Using", M, "Imputed Data Sets\n"))
    cat(paste("F(", df, ",", round(df2, 2), ")", "=", round(Dval,
                                                            3), "     p=", round(pval, 5), sep = ""), "\n")
  }
  invisible(res)
}

micombine.F <- function (Fvalues, df1, display = TRUE)
{
  #M <- length(Fvalues)
  dk <- df1 * Fvalues # actually: chi-squared = (numerator degrees of freedom) * F - i.e. df1 should actually be the NUMERATOR, not the denominator as specified in the package notes
  micombine.chisquare(dk = dk, df = df1, display = display)
}