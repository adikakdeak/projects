Summary:
- Performed descriptive analysis of the data by plotting the time series, scatter plot, Auto Correlation Function(ACF) and Partial ACF plots.
- Decided direction of analysis by first fitting deterministic seasonal model and then the stochastic SARIMA model. Finalized the SARIMA(p, d, q)X(P, D, Q)_s model to model this data.
- Came up with candidate models by finding possible values of p, d, q, P, D and Q by using multiple criterions like ACF and PACF plots, EACF matrix and BIC tables. This was done by handling the seasonal and non-seasonal component of the time series separately.
- Performed a diagnostic check on selected candidate models using methods like residual analysis, Shapiro-Wilk test, Ljung-Box test and parameter estimation using the Conditional Sum of Squares method.
- Finalized the model by referring to the results of diagnostic checks. Using this model, forecasted temperature values for next 10 months. 

Input File:
TempWorld2.csv

Code.Rmd
This file contains all the R codes.

Analysis.html
This file contains detailed description of every step and also the analysis of results.

Report.pdf
This file contains report of all the analysis along with all the results and outputs.
