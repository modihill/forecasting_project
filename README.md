# Forecasting_Project

There are Three tasks in this Report.

  1. **Data**: The weekly average Disease Specific Mortality rate in Paris as dependent variable for time series regression model and city's local climate weekly measures including temperature degrees Fahrenheit), size of pollutants and levels of noxious chemical emissions from cars and industry in the air - for the period of 2010 - 2020. Total 508 time points, with mortality as dependent series and temp, chem1, chem2, particle.size as independent series. 
  The main aim of the first task is to predict the `Disease Specific Mortality Rate in Paris` for upcoming four weeks based on the time series analysis. Time Series Regression Models with independent variables, dynlm, exponential smoothing and state-space models are used for this analysis. Then, the appropriate measures of resisual diagnostic checks and model performance measures MASE, AIC and BIC will be considered to find the best optimal model to predict the 4 weeks ahead forecast of `Mortality Rate`.


  2. **Data**: The first day of flowering(FFD) of specific plant and contemporaneous yearly averaged climate variables including Temperature, Rainfall, Radiation Level and Relative Humidity measured from 1984– 2014 (31 years).
  
  The aim of the second task is to analyse whether the FFD of given plant is impacted by such climate factors and give prediction of next 4 years of FFD with optimal univariate model decided by residual checkings and MASE() score. Provide Prediction for optimal univariate Model for each methods (Finite DLM, Poly DLM, koyck DLM, ARDL DLM, Dynamic Model, Exponential Model, State-space Models).
  
  3_a. **Data**: The Rank-based Order similarity metric (RBO) of 81 species of plant and contemporaneous yearly averaged climate variables including Temperature, Rainfall, Radiation Level and Relative Humidity measured from 1984– 2014 (31 years). 
  
  Higher RBO values indicate higher similarity of the order of the first flowering occurrence (based on FFD) of the 81 species from 1983 compared to each of the subsequent years.
  
  The aim of the this task is to analyse whether the decreasing RBO from the base year 1984 is due to such climate factors and give prediction of next 4 years of RBO with optimal univariate model decided by residual checkings and MASE() score. Provide Prediction for univariate Model(Finite DLM, Poly DLM, koyck DLM, ARDL DLM) with lowest `MASE()` score.
  
  3_b.
  
  The aim of the this task is to analyse whether the decreasing RBO from the base year 1984 is due to such climate factors and give prediction of next 4 years of RBO with optimal univariate model decided by residual checkings and MASE() score with consideration of 1996 year as an intervention point due to drought. Provide Prediction for univariate dynamic model from dynlm package with lowest `MASE()` score.
