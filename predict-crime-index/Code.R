# SETTING THE STAGE

# set the environment
graphics.off()
rm(list = ls())

# import packages
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(ks)
library(rjags)
library(runjags)

# import MCMC utilities
source("DBDA2E-utilities.R")

#=============================================================

# FUNCTION TO FIND MCMC SUMMARY

summary_MCMC = function(codaSamples, compVal = NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples, chains = TRUE)
  paramName = colnames(mcmcMat)
  for (pName in paramName) {
    if (pName %in% colnames(compVal)) {
      if (!is.na(compVal[pName])) {
        summaryInfo = rbind(summaryInfo,
                            summarizePost(
                              paramSampleVec = mcmcMat[, pName],
                              compVal = as.numeric(compVal[pName])
                            ))
      }
      else {
        summaryInfo = rbind(summaryInfo,
                            summarizePost(paramSampleVec = mcmcMat[, pName]))
      }
    } else {
      summaryInfo = rbind(summaryInfo,
                          summarizePost(paramSampleVec = mcmcMat[, pName]))
    }
  }
  rownames(summaryInfo) = paramName
  
  return(summaryInfo)
}

#=============================================================

# FUNCTION TO PLOT MCMC RESULTS

plot_MCMC = function(codaSamples,
                     data,
                     xName = "x",
                     yName = "y",
                     showCurve = FALSE,
                     compVal = NULL,
                     saveName = NULL,
                     saveType = "jpg") {
  y = data[, yName]
  x = as.matrix(data[, xName])
  
  mcmcMat = as.matrix(codaSamples, chains = TRUE)
  chainLength = NROW(mcmcMat)
  zbeta0 = mcmcMat[, "zbeta0"]
  zbeta  = mcmcMat[, grep("^zbeta$|^zbeta\\[", colnames(mcmcMat))]
  
  if (ncol(x) == 1) {
    zbeta = matrix(zbeta, ncol = 1)
  }
  
  zVar = mcmcMat[, "zVar"]
  beta0 = mcmcMat[, "beta0"]
  beta  = mcmcMat[, grep("^beta$|^beta\\[", colnames(mcmcMat))]
  
  if (ncol(x) == 1) {
    beta = matrix(beta, ncol = 1)
  }
  
  tau = mcmcMat[, "tau"]
  
  pred  = mcmcMat[, grep("^pred$|^pred\\[", colnames(mcmcMat))]
  
  if (ncol(x) == 1) {
    pred = matrix(pred, ncol = 1)
  }
  
  #-----------------
  
  # Compute R-squared value
  
  Rsq = zbeta %*% matrix(cor(y, x), ncol = 1)
  
  #-----------------
  
  # Marginal histograms
  
  plot_hist = function(panelCount,
                       saveName,
                       finished = FALSE,
                       nRow = 2,
                       nCol = 3) {
    
    # If finishing a set:
    if (finished == TRUE) {
      if (!is.null(saveName)) {
        saveGraph(file = paste0(saveName, ceiling((
          panelCount - 1
        ) / (nRow * nCol))),
        type = saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ((panelCount %% (nRow * nCol)) == 1) {
        # If previous graph was open, save previous one:
        if (panelCount > 1 & !is.null(saveName)) {
          saveGraph(file = paste0(saveName, (panelCount %/% (nRow * nCol))),
                    type = saveType)
        }
        
        # Open new graph
        openGraph(width = nCol * 7.0 / 3, height = nRow * 2.0)
        layout(matrix(1:(nRow * nCol), nrow = nRow, byrow = TRUE))
        par(mar = c(4, 4, 2.5, 0.5),
            mgp = c(2.5, 0.7, 0))
      }
      
      # Increment and return panel count:
      panelCount = panelCount + 1
      return(panelCount)
    }
  }
  
  #-----------------
  
  # Original scale:
  
  panelCount = 1
  if (!is.na(compVal["beta0"])) {
    panelCount = plot_hist(panelCount , saveName = paste0(saveName, "PostMarg_beta0"))
    histInfo = plotPost(
      beta0,
      cex.lab = 1.75,
      showCurve = showCurve,
      xlab = bquote(beta[0]),
      main = "Intercept",
      compVal = as.numeric(compVal["beta0"])
    )
  } else {
    histInfo = plotPost(
      beta0,
      cex.lab = 1.75,
      showCurve = showCurve,
      xlab = bquote(beta[0]),
      main = "Intercept"
    )
  }
  
  for (bIdx in 1:ncol(beta)) {
    panelCount = plot_hist(panelCount, saveName = paste0(saveName, "PostMarg_", bIdx))
    if (!is.na(compVal[paste0("beta[", bIdx, "]")])) {
      histInfo = plotPost(
        beta[, bIdx],
        cex.lab = 1.75,
        showCurve = showCurve,
        xlab = bquote(beta[.(bIdx)]),
        main = xName[bIdx],
        compVal = as.numeric(compVal[paste0("beta[", bIdx, "]")])
      )
    } else {
      histInfo = plotPost(
        beta[, bIdx],
        cex.lab = 1.75,
        showCurve = showCurve,
        xlab = bquote(beta[.(bIdx)]),
        main = xName[bIdx]
      )
    }
  }
  
  panelCount = plot_hist(panelCount, saveName = paste0(saveName, "PostMarg_scale"))
  histInfo = plotPost(
    tau,
    cex.lab = 1.75,
    showCurve = showCurve,
    xlab = bquote(tau),
    main = paste("Scale")
  )
  
  panelCount = plot_hist(panelCount, saveName = paste0(saveName, "PostMarg_rsq"))
  histInfo = plotPost(
    Rsq,
    cex.lab = 1.75,
    showCurve = showCurve,
    xlab = bquote(R ^ 2),
    main = paste("R-squared")
  )
  
  for (predIdx in 1:ncol(pred)) {
    panelCount = plot_hist(panelCount, saveName = paste0(saveName, "PostMarg_", predIdx))
    histInfo = plotPost(
      pred[, predIdx],
      cex.lab = 1.75,
      showCurve = showCurve,
      xlab = bquote(pred[.(predIdx)]),
      main = paste("Prediction", predIdx)
    )
    
  }
  
  #-----------------
  
  # Standardized scale:
  
  panelCount = 1
  panelCount = plot_hist(panelCount, saveName = paste0(saveName, "PostMargZ_intercept"))
  histInfo = plotPost(
    zbeta0,
    cex.lab = 1.75,
    showCurve = showCurve,
    xlab = bquote(z * beta[0]),
    main = "Intercept"
  )
  
  for (bIdx in 1:ncol(beta)) {
    panelCount = plot_hist(panelCount, saveName = paste0(saveName, "PostMargZ_", bIdx))
    histInfo = plotPost(
      zbeta[, bIdx],
      cex.lab = 1.75,
      showCurve = showCurve,
      xlab = bquote(z * beta[.(bIdx)]),
      main = xName[bIdx]
    )
  }
  
  panelCount = plot_hist(panelCount, saveName = paste0(saveName, "PostMargZ_scale"))
  histInfo = plotPost(
    zVar,
    cex.lab = 1.75,
    showCurve = showCurve,
    xlab = bquote(z * tau),
    main = paste("Scale")
  )
  
  panelCount = plot_hist(panelCount, saveName = paste0(saveName, "PostMargZ_rsq"))
  histInfo = plotPost(
    Rsq,
    cex.lab = 1.75,
    showCurve = showCurve,
    xlab = bquote(R ^ 2),
    main = paste("R-squared")
  )
  
  panelCount = plot_hist(panelCount,
                         finished = TRUE,
                         saveName = paste0(saveName, "PostMargZ"))
  
}

#=============================================================

# GETTING DATA INTO R ENVIRONMENT

crime_demo_dat <- read_csv("clinton1.csv")

#=============================================================

# DESCRIPTIVE CHECK

# Kernel density estimation - dependent variable
ggplot(crime_demo_dat, aes(x = crime_index)) +
  geom_density(color = "darkblue", fill = "lightblue") +
  ggtitle("\nKernel Density Plot for Crime Index\n") +
  xlab("Crime Index") +
  ylab("Density")

#-------------------

#Kernel density estimation - independent variables
kde_plot_1 <- ggplot(crime_demo_dat, aes(x = median_age)) +
  geom_density(color = "darkblue", fill = "lightblue") +
  xlab("Median Age") +
  ylab("Density")

kde_plot_2 <- ggplot(crime_demo_dat, aes(x = savings)) +
  geom_density(color = "darkblue", fill = "lightblue") +
  xlab("Savings") +
  ylab("Density")

kde_plot_3 <- ggplot(crime_demo_dat, aes(x = per_capita_income)) +
  geom_density(color = "darkblue", fill = "lightblue") +
  xlab("Per Capita Income") +
  ylab("Density")

kde_plot_4 <- ggplot(crime_demo_dat, aes(x = poverty_prcnt)) +
  geom_density(color = "darkblue", fill = "lightblue") +
  xlab("Poverty Percentage") +
  ylab("Density")

kde_plot_5 <- ggplot(crime_demo_dat, aes(x = veterans_prcnt)) +
  geom_density(color = "darkblue", fill = "lightblue") +
  xlab("Veterans Percentage") +
  ylab("Density")

kde_plot_6 <- ggplot(crime_demo_dat, aes(x = population_density)) +
  geom_density(color = "darkblue", fill = "lightblue") +
  xlab("Population Density") +
  ylab("Density")

fig_kde <- ggarrange(
  kde_plot_1,
  kde_plot_2,
  kde_plot_3,
  kde_plot_4,
  kde_plot_5,
  kde_plot_6,
  nrow = 2,
  ncol = 3
)
fig_kde

#-------------------

# Scatter plots
scatter_plot_1 <- ggplot(crime_demo_dat, aes(x = median_age, y = crime_index)) +
  geom_point() +
  xlab("Median Age") +
  ylab("Crime Index")

scatter_plot_2 <- ggplot(crime_demo_dat, aes(x = savings, y = crime_index)) +
  geom_point() +
  xlab("Savings") +
  ylab("Crime Index")

scatter_plot_3 <- ggplot(crime_demo_dat, aes(x = per_capita_income, y = crime_index)) +
  geom_point() +
  xlab("Per Capita Income") +
  ylab("Crime Index")

scatter_plot_4 <- ggplot(crime_demo_dat, aes(x = poverty_prcnt, y = crime_index)) +
  geom_point() +
  xlab("Poverty Percentage") +
  ylab("Crime Index")

scatter_plot_5 <- ggplot(crime_demo_dat, aes(x = veterans_prcnt, y = crime_index)) +
  geom_point() +
  xlab("Veterans Percentage") +
  ylab("Crime Index")

scatter_plot_6 <- ggplot(crime_demo_dat, aes(x = population_density, y = crime_index)) +
  geom_point() +
  xlab("Population Density") +
  ylab("Crime Index")

fig_scatter <- ggarrange(
  scatter_plot_1,
  scatter_plot_2,
  scatter_plot_3,
  scatter_plot_4,
  scatter_plot_5,
  scatter_plot_6,
  nrow = 2,
  ncol = 3
)
fig_scatter

#-------------------

#box plots
meltData <- melt(crime_demo_dat)

ggplot(meltData, aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue") +
  facet_wrap( ~ variable, scale = "free") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

#-------------------

#correlation between dependent variables
cormat <- round(cor(crime_demo_dat[, -7]), 2)

get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(cormat, na.rm = TRUE)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    size = 12,
    hjust = 1
  )) +
  coord_fixed() + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

#=============================================================

# DATA PREPARATION

#extract sample
set.seed(999)
crime_demo_dat_sample <- sample_n(crime_demo_dat,
                                  round(dim(crime_demo_dat)[1] * 0.1, 0))

#separate dependent and independent variables
y = as.matrix(crime_demo_dat[, "crime_index"])
x = as.matrix(crime_demo_dat[, -7])

#specify input values for dependent variables
xPred = array(NA, dim = c(6, 6))
xPred[1, ] = c(33, 87000, 50000, 7.02, 10, 200)
xPred[2, ] = c(18, 30000, 15000, 17.19, 19.22, 1)
xPred[3, ] = c(60, 690000, 4300, 50.89, 29, 35430)
xPred[4, ] = c(42, 50000, 26700, 3.6, 4.58, 8080.59)
xPred[5, ] = c(50, 400000, 16300, 34.77, 7.77, 3000)
xPred[6, ] = c(27, 10800, 12200, 5.2, 18.59, 131)

#wrap the data in list for JAGS
dataList <- list(
  x = x ,
  y = y ,
  xPred = xPred ,
  Nx = dim(x)[2] ,
  Ntotal = dim(x)[1]
)

# specify initial values for MCMC model
initsList <- list(zbeta0 = 100,
                  zbeta = c(100, 1, 1, 1, 1, 1),
                  Var = 1000)

#=============================================================

# MODEL SPECIFICATION

modelString = "

# Standardize the data:
data {
  y_sd <- sd(y)
  for (i in 1:Ntotal) {
    zy[i] <- y[i,] / y_sd
  }
  for (j in 1:Nx) {
    x_sd[j] <- sd(x[, j])
    for (i in 1:Ntotal) {
      zx[i, j] <- x[i, j] / x_sd[j]
    }
  }
}

#-------------------

# Specify the model for scaled data:
model {

  for (i in 1:Ntotal) {
    zy[i] ~ dgamma((mu[i]^2) / zVar, mu[i] / zVar) 
    mu[i] <- zbeta0 + sum(zbeta[1:Nx] * zx[i, 1:Nx]) 
  }
  
  # Priors on standardized scale:
  zbeta0 ~ dnorm(0, 1 / (20^2))
  zbeta[1] ~ dnorm(-9 / x_sd[1], 1 / (2.7/x_sd[1]^2))
  zbeta[2] ~ dnorm(0 / x_sd[2], 1 / (2.7/x_sd[2]^2))
  zbeta[3] ~ dnorm(0.000001 / x_sd[3], 1 / (2.7/x_sd[3]^2))
  zbeta[4] ~ dnorm(0.8 / x_sd[4], 1 / (0.3/x_sd[4]^2))
  zbeta[5] ~ dnorm(-4 / x_sd[5], 1 / (15/x_sd[5]^2))
  zbeta[6] ~ dnorm(0.0001 / x_sd[6], 1 / (2.7/x_sd[6]^2))
  zVar ~ dgamma(0.01, 0.01)
  
  #------------------
  
  # Transform to original scale:
  beta[1:Nx] <- (zbeta[1:Nx] / x_sd[1:Nx]) * y_sd
  beta0 <- zbeta0 * y_sd
  tau <- zVar * (y_sd)^2

  #------------------

  # Compute predictions at every step of the MCMC
  for (i in 1:6){
    pred[i] <- beta0 + beta[1] * xPred[i, 1] + beta[2] * xPred[i, 2] + beta[3] * xPred[i, 3] + 
               beta[4] * xPred[i, 4] + beta[5] * xPred[i,5] + beta[6] * xPred[i,6]
  }

}
"

#=============================================================

# MCMC SETTINGS AND RUNNING THE MODEL

# Specify MCMC settings
adaptSteps = 1000
burnInSteps = 6000
nChains = 5
thinSteps = 135
numSavedSteps = 2100

# Run the model
start_time <- Sys.time()

runJagsOut <- run.jags(
  method = "parallel",
  model = modelString,
  monitor = c("beta0", "beta", "tau", "zbeta0", "zbeta", "zVar", "pred"),
  data = dataList,
  inits = initsList,
  n.chains = nChains,
  adapt = adaptSteps,
  burnin = burnInSteps,
  sample = numSavedSteps,
  thin = thinSteps,
  summarise = FALSE,
  plots = FALSE
)

codaSamples = as.mcmc.list(runJagsOut)

end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time)

#save image
save.image(file = "run5_p.RData")

#load required image
load(file = "run4_p.RData")

#=============================================================

# DIAGNOSTIC CHECKS

#MCMC check
diagMCMC(codaSamples, parName = "beta0")
diagMCMC(codaSamples, parName = "beta[1]")
diagMCMC(codaSamples, parName = "beta[2]")
diagMCMC(codaSamples, parName = "beta[3]")
diagMCMC(codaSamples, parName = "beta[4]")
diagMCMC(codaSamples, parName = "beta[5]")
diagMCMC(codaSamples, parName = "beta[6]")
diagMCMC(codaSamples, parName = "tau")

diagMCMC(codaSamples, parName = "pred[1]")
diagMCMC(codaSamples, parName = "pred[2]")
diagMCMC(codaSamples, parName = "pred[3]")
diagMCMC(codaSamples, parName = "pred[4]")
diagMCMC(codaSamples, parName = "pred[5]")
diagMCMC(codaSamples, parName = "pred[6]")

#summary
summary <- summary_MCMC(codaSamples = codaSamples)
summary

#=============================================================

# RESULTS

#specify comparison values
compVal <- data.frame(
  "beta0" = NA,
  "beta[1]" = NA,
  "beta[2]" = NA,
  "beta[3]" = NA,
  "beta[4]" =  NA,
  "beta[5]" =  NA,
  "beta[6]" =  NA,
  "tau" = NA,
  check.names = FALSE
)

#------------------

#plot posterior distributions
plot_MCMC(
  codaSamples = codaSamples,
  data = crime_demo_dat,
  xName = c(
    "median_age",
    "savings",
    "per_capita_income",
    "poverty_prcnt",
    "veterans_prcnt",
    "population_density"
  ),
  yName = "crime_index",
  compVal = compVal
)

#------------------

#compute required values
coefficients <- summary[c(2:8), 3]
Variance <- summary[9, 3]

meanGamma <-
  as.matrix(cbind(rep(1, nrow(x)), x)) %*% as.vector(coefficients)
randomData <- rgamma(n = length(y),
                     shape = meanGamma ^ 2 / Variance,
                     rate = meanGamma / Variance)

predicted <- data.frame(crime_index = randomData)
observed <- data.frame(elapsed = y)

predicted$type <- "Predicted"
observed$type <- "Observed"
dataPred <- rbind(predicted, observed)

#------------------

# Display the density plot of observed data and predicted
ggplot(dataPred,
       aes(crime_index, fill = type)) +
  geom_density(alpha = 0.2) +
  xlab("No. of Observations") +
  ylab("Density") +
  ggtitle("Density Estimate: Observed v/s Predicted") +
  theme(plot.title = element_text(hjust = 0.5))
