# R Code for Forecasting High-Frequency Trade Durations: A Regime-Switching Approach with Flexible Hazard Functions

Authors: Tan Yiing Fei; Koh You Beng; Ng Kok Haur  
Year: 2026  

## Description
This repository contains the R code used to produce the empirical results in the manuscript.

FOLDER 1_Simulation:
Note: The numbering in Files, Required Packages, and Instructions corresponds to each other.
## Files
1) 1.1_Sim data EGIG.SCD.R - Simulates EGIG.SCD model data
   2.1_Sim data EGIG.TSCD.R - Simulates EGIG.TSCD model data
2) 1.2_Sim Est EGIG.SCD_SIR.R - Parameter estimation for EGIG.SCD model
   2.2_Sim Est EGIG.TSCD_SIR.R - Parameter estimation for EGIG.TSCD model

## Software Requirements
R version 2025.09.0  
Required packages:
1) GeneralizedHyperbolic - work with heavy-tailed/skewed distributions (density, random, estimation).
2) optimx - advanced numerical optimisation with multiple methods and diagnostics.

## Instructions
- Set working directory
1) To generate a simulated EGIG.SCD data
2) To perform parameter estimation via constrained optimisation

## No empirical data required.


FOLDER 2_Diurnally Adj Dur:
Note: The numbering in Files, Required Packages, and Instructions corresponds to each other.
## Files
1) In-sample Adj Dur.R 
   Out-sample Adj Dur.R 
   - Processes trade-by-trade data to generate diurnally adjusted event-based

## Software Requirements
R version 2025.09.0  
Required packages:
1) MASS - providing functions for applied statistics

## Instructions
- Set working directory
1) To perform robust linear regression

## Data are obtained from Bloomberg (licensed access required).
## Due to license restrictions, raw data are not provided.


FOLDER 3_In-sample:
Note: The numbering in Files, Required Packages, and Instructions corresponds to each other.
## Files
1) 0.2_Combine all run files.R - Combines valid estimation outputs into a single dataset for analysis
2) 0.3_Emp Gamma kernel plots.R - Computes empirical PDF, CDF, Survival, and Hazard functions using Gamma kernel density estimation and generates plots
3) 1.1_Est_Wei.SCD.R
   2.1_Est_GG.SCD.R
   3.1_Est_Burr.SCD.R
   4.1_Est_GB2.SCD.R
   5.1_Est_EGIG.SCD.R
   6.1_Est_EGIG.TSCD.R 
   - Parameter estimation for different models
4) 1.3_Cond vs Emp Plots_Wei.SCD.R
   2.3_Cond vs Emp Plots_GG.SCD.R 
   3.3_Cond vs Emp Plots_Burr.SCD.R
   4.3_Cond vs Emp Plots_GB2.SCD.R
   5.3_Cond vs Emp Plots_EGIG.SCD.R
   6.3_Cond vs Emp Plots_EGIG.TSCD.R
   - Plot conditional functions and empirical Gamma-kernel estimate functions for comparison

## Software Requirements
R version 2025.09.0
Required packages:
1) dplyr - fast, consistent and readable data manipulation
2) (a) ggplot2 - elegant and versatile plotting system
   (b) pracma - practical numerical math functions
   (c) gridExtra - arranges multiple plots on a single page/grid
3) optimx - advanced numerical optimisation with multiple methods and diagnostics
4) (a) ggplot2, (b) pracma, (c) gridExtra - refer to 2)
   (d) MASS - provides functions for applied statistics
   (e) gsl - provides special mathematical functions like gamma, Bessel, and other advanced

## Instructions
- Set working directory
1) To combine multiple estimation output files into a single dataset
2) (a) To generate the corresponding plots
   (b) To compute empirical density, distribution, survival, and hazard functions using Gamma kernel estimation
   (c) To combine multiple plots
3) To perform parameter estimation via constrained optimisation
4) To plot conditional model functions against empirical Gamma-kernel estimates for comparison

## Data are obtained from Bloomberg (licensed access required).
## Due to license restrictions, raw data are not provided.


FOLDER 4_Out-sample:
Note: The numbering in Files, Required Packages, and Instructions corresponds to each other.
## Files
1) 1.1_Roll win est_Wei.SCD.R
   2.1_Roll win est_GG.SCD.R
   3.1_Roll win est_Burr.SCD.R
   4.1_Roll win est_GB2.SCD.R
   5.1_Roll win est_EGIG.SCD.R
   6.1_Roll win est_EGIG.TSCD.R
   - Rolling window parameter estimation for different models
2) 1.3_TaR_Wei.SCD.R
   2.3_TaR_GG.SCD.R
   3.3_TaR_Burr.SCD.R
   4.3_TaR_GB2.SCD.R
   5.3_TaR_EGIG.SCD.R
   6.3_TaR_EGIG.TSCD.R
   - Produces TaR duration forecasts for different models and evaluates their reliability
3) DM test.R - Compares forecast accuracy of multiple models using the DM Test.
4) MCS.R - Evaluates and compares forecast performance using loss functions and the MCS procedure
5) TaR Plot_EGIG.TSCD.R - visualises TaR duration forecasts alongside diurnally adjusted durations

## Software Requirements
R version 2025.09.0
Required packages:
1) optimx - advanced numerical optimisation with multiple methods and diagnostics
2) (a) GAS - used for backtesting Value-at-Risk (VaR)
   (b) rmutil - used for the generalised Gamma quantile function
3) sandwich - used for robust variance estimation in the DM test
4) MCS - used for performing the MCS procedure
5) ggplot2 - elegant and versatile plotting system

## Instructions
- Set working directory
1) To perform parameter estimation via constrained optimisation
2) (a) To compute TaR forecasts, violation statistics, and perform risk backtests
   (b) To compute GG distribution quantiles for TaR forecasts
   (c) To perform statistical backtests of TaR forecasts (Kupiec, Dynamic Quantile, Christoffersen)
3) To perform pairwise forecast comparison using the DM Test
4) To evaluate and compare forecast performance of multiple models using the MCS procedure
5) To generate the corresponding plots

## Data are obtained from Bloomberg (licensed access required).
## Due to license restrictions, raw data are not provided.
