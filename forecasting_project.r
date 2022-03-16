##### Importing the necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(knitr)
library(TSA)
library(tseries)
library(car)
library(dynlm)
library(Hmisc)
library(forecast)
library(xts)
library(ggplot2)
library(AER)
library(x12)
library(dLagM)
library(kableExtra)


##### Function Definations 

# Descriptive Analysis function
descriptive_analysis <- function(ts, object)
{
  plot(ts,
       ylab = c(paste0(toString(object))),
       main = c(paste0("Time Series Plot of ",toString(object))),
       type="o")
  
  par(mfrow=c(2,1))
  acf(ts,
      lag.max = 48,
      main = c(paste0("ACF plot of ",toString(object))))
  
  pacf(ts,
      lag.max = 48,
      main = c(paste0("PACF plot of ",toString(object))))
  par(mfrow=c(1,1))

  
  print(adf.test(ts))
  print(pp.test(ts))

}


# Function for Decomposition
decom <- function(ts, ts_series)
{
  decomposition <- stl(ts, t.window=15, s.window = "periodic", robust=TRUE)
  plot(decomposition,
       main = c(paste0("Weekly ", toString(ts_series)), " STL Decomposed Series"))
}


# Function for Summary and Residual Analysis
summary_residual_analysis <-  function(m)
{
  summary(m, diagnostics = TRUE)
  checkresiduals(m$model)
  print(bgtest(m$model))
  print(vif(m$model))
  print(shapiro.test(m$model$residuals))
}

#---------------------------------------------
#------------------ Task-1--------------------
#---------------------------------------------

## Data Preparation

mort <- read.csv("mort.csv")
head(mort)
class(mort)
mort <- mort[,-c(1)]
mort_TS <- ts(mort, start = c(2010,1), frequency = 52)


mortality_TS <- ts(mort$mortality, start = c(2010,1), frequency = 52)
temp_TS <- ts(mort$temp, start = c(2010,1), frequency = 52)
chem_1_TS <- ts(mort$chem1, start = c(2010,1), frequency = 52)
chem_2_TS <- ts(mort$chem2, start = c(2010,1), frequency = 52)
particle_size_TS <- ts(mort$particle.size, start = c(2010,1), frequency = 52)

mort_TS %>% head()
class(mort_TS)


##### 1. Mortality
head(mortality_TS)
class(mortality_TS)

##### 2. Temperature
head(temp_TS)
class(temp_TS)

##### 3. Chemical 1
head(chem_1_TS)
class(chem_1_TS)


##### 4. Chemical 2
head(chem_2_TS)
class(chem_2_TS)


##### 5. Partical Size
head(particle_size_TS)
class(particle_size_TS)


## Descriptive Analysis

##### 1. Mortality
descriptive_analysis(mortality_TS, "Averaged weekly disease specific mortality")

##### 2. Temperature
descriptive_analysis(temp_TS, "Averaged weekly temperature")

##### 3. Chemical 1
descriptive_analysis(chem_1_TS, "Averaged weekly Chem1 Emissions")


##### 4. Chemical 2
descriptive_analysis(chem_2_TS, "Averaged weekly Chem2 Emissions")

##### 5. Partical Size
descriptive_analysis(particle_size_TS, "Averaged weekly Pollutants Particle Size")

##### 6. Combined Scaled Time Series Plot
combined <- scale(mort_TS)

plot(combined, 
     plot.type="s", 
     col = c("#05386b","#f01b1d","#5c3c92","#d2601a","#1b6535"), 
     main = "Scaled Time Series Plot of Mortality, Paris's Climate and pollutanats")

legend("topleft", 
       lty=1, 
       col = c("#05386b","#f01b1d","#5c3c92","#d2601a","#1b6535"), 
c("Mortality", "Temperature", "Chemical 1", "Chemical 2" , "Partical Size"))


# Correlation
cor(mort_TS) %>% round(3)


## Decomposition
##### 1. Mortality
decom(mortality_TS, "Average Weekly Disease Specific Mortality")

##### 2. Temperature
decom(temp_TS, "Average Weekly Temperature")

##### 3. Chemical 1
decom(temp_TS, "Average Weekly Chemical 1 Emissions")

##### 4. Chemical 2
decom(temp_TS, "Average Weekly Chemical 2 Emissions")

##### 5. Partical Size
decom(particle_size_TS, "Average Weekly Pollutants Particle Size")


## Time Series Regression Methods
### Finite Distributed Lag Model
for (i in 1:8)
{
  model_dlm <- dlm(formula = mortality ~ temp + chem1 + chem2 + particle.size, data = data.frame(mort), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}

task1_model1_1 <- dlm(formula = mortality ~ temp + chem1 + chem2 + particle.size, 
                      data = data.frame(mort), 
                      q = 8)

summary_residual_analysis(task1_model1_1)


for (i in 1:8)
{
  model_dlm <- dlm(formula = mortality ~ temp + chem1, data = data.frame(mort), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}

task1_model1_2 <- dlm(formula = mortality ~ temp + chem1, 
                      data = data.frame(mort), 
                      q = 8)

summary_residual_analysis(task1_model1_2)


### PolyNomial Distributed Lag Model
task_1_model2 = polyDlm(y = as.vector(mort$mortality), 
                        x = as.vector(mort$chem1) + as.vector(mort$temp), 
                        q = 8, 
                        k = 2, 
                        show.beta = TRUE)

summary_residual_analysis(task_1_model2)


### Koyck Distributed Lag Model
task_1_model3 = koyckDlm(as.vector(mort$chem1) + as.vector(mort$temp), 
                         y = as.vector(mort$mortality))
summary_residual_analysis(task_1_model3)


### Autoregressive Distributed Lag Model
for (i in 1:5)
{
  for(j in 1:5)
  { 
    model_ardlm <- ardlDlm(formula = mortality ~ temp + chem1, 
                           data = data.frame(mort), 
                           p = i,
                           q = j)
    cat("p =", i, 
        "q =" , j,
        "AIC =", AIC(model_ardlm$model), 
        "BIC =", BIC(model_ardlm$model),
        "MASE =", MASE(model_ardlm)$MASE, "\n")
  }
}

#ardlDLM(3,2)
task_1_model4_1 = ardlDlm(formula = mortality ~ temp + chem1,
                          data = data.frame(mort), 
                          p = 3,
                          q = 2)

summary_residual_analysis(task_1_model4_1)

#ardlDLM(4,2)
task_1_model4_2 = ardlDlm(formula = mortality ~ temp + chem1,
                          data = data.frame(mort), 
                          p = 4,
                          q = 2)

summary_residual_analysis(task_1_model4_2)

#ardlDLM(5,2)
task_1_model4_3 = ardlDlm(formula = mortality ~ temp + chem1,
                          data = data.frame(mort), 
                          p = 5,
                          q = 2)

summary_residual_analysis(task_1_model4_3)


### Dynamic Models
#dynlm_1
task_1_model5_1 = dynlm(mortality_TS ~ L(mortality_TS , k = 1 ) + trend(mortality_TS))
summary(task_1_model5_1)
checkresiduals(task_1_model5_1)

#dynlm_2
task_1_model5_2 = dynlm(mortality_TS ~ L(mortality_TS , k = 1) + season(mortality_TS))
summary(task_1_model5_2)
checkresiduals(task_1_model5_2)

#dynlm_3
task_1_model5_3 = dynlm(mortality_TS ~ L(mortality_TS , k = 1 ) + trend(mortality_TS) + season(mortality_TS))
summary(task_1_model5_3)
checkresiduals(task_1_model5_3)

#dynlm_4
task_1_model5_4 = dynlm(mortality_TS ~ L(mortality_TS , k = 1 ) + L(mortality_TS , k = 2 ) + trend(mortality_TS) + season(mortality_TS))
summary(task_1_model5_4)
checkresiduals(task_1_model5_4)

#dynlm_5
task_1_model5_5 = dynlm(mortality_TS ~ L(mortality_TS , k = 1 ) + L(mortality_TS , k = 2 ) + temp_TS + trend(mortality_TS) + season(mortality_TS))
summary(task_1_model5_5)
checkresiduals(task_1_model5_5)

#dynlm_6
task_1_model5_6 = dynlm(mortality_TS ~ L(mortality_TS , k = 1 ) + L(mortality_TS , k = 2 ) + season(mortality_TS) + trend(mortality_TS) + temp_TS + L(temp_TS , k = 1 ))
summary(task_1_model5_6)
checkresiduals(task_1_model5_6)

#dynlm_7
task_1_model5_7 = dynlm(mortality_TS ~ L(mortality_TS , k = 1 ) + L(mortality_TS , k = 2 ) +  L(mortality_TS , k = 3 ) + trend(mortality_TS) + chem_1_TS + temp_TS + L(chem_1_TS, k=1))
summary(task_1_model5_7)
checkresiduals(task_1_model5_7)


#Model Comparison
attr(task_1_model3$model,"class") = "lm"
T1_models <- c("Finite DLM (all predictors)", "Finite DLM (chem1 + temp)", "Poly DLM", "Koyck", "ARDL_3_2", "ARDL_4_2", "ARDL_5_2", "dynlm_1", "dynlm_2", "dynlm_3", "dynlm_4", "dynlm_5" ,"dynlm_6", "dynlm_7")

T1_aic <- AIC(task1_model1_1$model, task1_model1_2$model, task_1_model2$model, task_1_model3$model, task_1_model4_1$model, task_1_model4_2$model, task_1_model4_3$model, task_1_model5_1,task_1_model5_2,task_1_model5_3, task_1_model5_4, task_1_model5_5, task_1_model5_6, task_1_model5_7)$AIC

T1_bic <- BIC(task1_model1_1$model, task1_model1_2$model, task_1_model2$model, task_1_model3$model, task_1_model4_1$model, task_1_model4_2$model, task_1_model4_3$model, task_1_model5_1,task_1_model5_2,task_1_model5_3, task_1_model5_4, task_1_model5_5, task_1_model5_6, task_1_model5_7)$BIC

T1_mase <- MASE(task1_model1_1$model, task1_model1_2$model, task_1_model2$model, task_1_model3$model, task_1_model4_1, task_1_model4_2, task_1_model4_3, lm(task_1_model5_1), lm(task_1_model5_2), lm(task_1_model5_3), lm(task_1_model5_4), lm(task_1_model5_5), lm(task_1_model5_6), lm(task_1_model5_7))$MASE

T1_Model_Comparison <- data.frame(T1_models, T1_mase, T1_aic, T1_bic)
colnames(T1_Model_Comparison) <- c("Model","MASE","AIC", "BIC")

kbl(T1_Model_Comparison) %>% kable_paper()


## Exponential smoothing methods
yearly_mortality_TS = aggregate(zoo(mortality_TS),as.yearmon,sum)
yearly_mortality_TS

yearly_mortality_TS[yearly_mortality_TS == 171.34] <- NA
yearly_mortality_TS = na.omit(yearly_mortality_TS)
yearly_mortality_TS


T1_HW_models = c("Holt_Winter additive method",
                 "Holt_Winter multiplicative method with exponential trend",
                 "Holt_Winter multiplicative method",
                 "Holt_Winter additive method",
                 "Holt_Winter multiplicative method with exponential trend",
                 "Holt_Winter multiplicative method")

T1_exponential = c(TRUE,FALSE)
T1_seasonality = c("additive","multiplicative")
T1_damped = c(TRUE,FALSE)
T1_exponential_models <- expand.grid(T1_exponential, T1_seasonality, T1_damped)
T1_exponential_models <- T1_exponential_models[-c(1,5),]

T1_HW_AIC <- array(NA, 6)
T1_HW_BIC <- array(NA, 6)
T1_HW_MASE <- array(NA, 6)
T1_levels <- array(NA, dim=c(6,3))


for (i in 1:6)
{
  T1_HW_model <- hw(yearly_mortality_TS,
                    exponential = T1_exponential_models[i,1],
                    seasonal = toString(T1_exponential_models[i,2],
                                        damped = T1_exponential_models[i,3]))
  T1_HW_AIC[i] <- T1_HW_model$model$aic
  T1_HW_BIC[i] <- T1_HW_model$model$bic
  T1_HW_MASE[i] <- accuracy(T1_HW_model)[6]
  T1_levels[i,1] <- T1_exponential_models[i,1]
  T1_levels[i,2] <- toString(T1_exponential_models[i,2])
  T1_levels[i,3] <- T1_exponential_models[i,3]
  summary(T1_HW_model)
  checkresiduals(T1_HW_model)
  print(shapiro.test(T1_HW_model$model$residuals))
}

T1_results_HW = data.frame(T1_HW_models, T1_levels, T1_HW_MASE, T1_HW_AIC, T1_HW_BIC)
colnames(T1_results_HW) = c("Model", "Exponential","Seasonality","Damped","MASE","AIC", "BIC")

kbl(T1_results_HW) %>% kable_paper()


#Formating T1_results_HW table 
T1_results_HW$Damped <- factor(T1_results_HW$Damped,
                               levels = c(TRUE, FALSE),
                               labels = c("damped"," "))

T1_results_HW <- unite(T1_results_HW,
                       "Model", c("Model","Damped"), sep = "_")

T1_results_HW <- T1_results_HW[,-c(2,3)]
kbl(T1_results_HW) %>% kable_paper()



## State-space Models
T1_ets_models = c("AAA", "MAA", "MAM", "MMM")
T1_damped = c(TRUE,FALSE)
T1_ETS_models <- expand.grid(T1_ets_models, T1_damped)

T1_ETS_AIC <- array(NA, 8)
T1_ETS_BIC <- array(NA, 8)
T1_ETS_MASE <- array(NA, 8)
T1_levels <- array(NA, dim=c(8,2))

for (i in 1:8)
{
  T1_ETS <- ets(yearly_mortality_TS,
                model = toString(T1_ETS_models[i, 1]), damped = T1_ETS_models[i,2])
  T1_ETS_AIC[i] <- T1_ETS$aic
  T1_ETS_BIC[i] <- T1_ETS$bic
  T1_ETS_MASE[i] <- accuracy(T1_ETS)[6]
  T1_levels[i,1] <- toString(T1_ETS_models[i,1])
  T1_levels[i,2] <- T1_ETS_models[i,2]
  summary(T1_ETS)
  checkresiduals(T1_ETS)
  print(shapiro.test(T1_ETS$residuals))
}


T1_results_ETS = data.frame(T1_levels, T1_ETS_MASE, T1_ETS_AIC, T1_ETS_BIC)
colnames(T1_results_ETS) = c("Model","Damped","MASE","AIC", "BIC")

kbl(T1_results_ETS) %>% kable_paper()

#Formating T1_results_ETS table 
T1_results_ETS$Damped <- factor(T1_results_ETS$Damped,
                                levels = c(TRUE, FALSE),
                                labels = c("damped"," "))

T1_results_ETS <- unite(T1_results_ETS,
                        "Model", c("Model","Damped"), sep = "_")

kbl(T1_results_ETS) %>% kable_paper()


## Model Comparison
T1_Model_Comparison <- rbind(T1_Model_Comparison, T1_results_ETS, T1_results_HW)

T1_sorted_MASE <- T1_Model_Comparison %>% arrange(MASE)
kbl(T1_sorted_MASE) %>%
  kable_paper()


## Forecasting
### Monthly Forcasting from ETS model

#Comparison of Charts
prediction <- ets(yearly_mortality_TS,
                  model="MAM",
                  damped = T)
prediction <- forecast(prediction)

prediction1 <- ets(yearly_mortality_TS,
                   model="MAA",
                   damped = TRUE)
prediction1 <- forecast(prediction1)

plot(prediction,
     main = "Next Two years prediction of mortality using ETS",
     ylab = "Mortality Rate",
     fcol = "#f01b1d")

lines(fitted(prediction1), col = "#05386b")

lines(fitted(prediction), col = "#f01b1d")


legend("topleft", 
       lty = 1, 
       col = c("black", "#f01b1d", "#05386b"), 
       c("Data", "MAM_damped", "MAA_damped"))


#Prediction Chart
plot(prediction1, fcol = "#b20238", 
     main = "Forecasting of Mortality in 2019 and 2020",
     ylab = "Radiation")
lines(fitted(prediction1), col = "#f01b1d")

legend("topleft", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))



#Prediction Points
kbl(prediction) %>% kable_paper()



### weekly Forcasting from dylnm model
#Prediction Chart
q = 4
n = nrow(task_1_model5_6$model)
m.frc = array(NA , (n + q))
m.frc[1:n] = mortality_TS[3:length(mortality_TS)]
trend = array(NA,q)
trend.start = task_1_model5_6$model[n,"trend(mortality_TS)"]
trend = seq(trend.start , trend.start + q/12, 1/12)

for(i in 1:q){
  weeks = array(0,51)
  weeks[(i-12)%%52] = 1
  print(weeks)
  
  data.new =c(1,m.frc[n-1+i],m.frc[n-2+i],weeks,trend[i],1,1)
  m.frc[n+i] = as.vector(task_1_model5_6$coefficients) %*% data.new
}



plot(mortality_TS,
     ylab='Mortality Rate',
     xlab='Year',
     main = "Next 4 weeks (41 - 44) prediction of disease specific Mortality Rate in Paris")

lines(ts(m.frc[(n+1):(n+q)],
         start=c(2019,41), frequency = 52),
      col = "#f01b1d")

legend("topleft", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))


#Prediction Points
t1_Time_point <- c("41, 2019", "42, 2019", "43, 2019", "44, 2019")
t1_points <- m.frc[(n+1):(n+q)]

t1_prediction_points <- data.frame(t1_Time_point, t1_points)
colnames(t1_prediction_points) = c("Time","Point Forecast")


kable(t1_prediction_points) %>% kable_paper()

#------------------------------------------------------#
#----------------------Task - 2------------------------#
#------------------------------------------------------#

## Data Preparation

FFD_df <- read_csv("FFD.csv")
FFD_prediction <- read_csv("Covariate x-values for Task 2.csv")
FFD_prediction = na.omit(FFD_prediction)

head(FFD_df)
class(FFD_df)

FFD_df <- FFD_df[,-c(1)]
FFD_df_TS <- ts(FFD_df, start = 1984)


temp_2_TS <- ts(FFD_df$Temperature, start = 1984)
rain_TS <- ts(FFD_df$Rainfall, start = 1984)
radiation_TS <- ts(FFD_df$Radiation, start = 1984)
humidity_TS <- ts(FFD_df$RelHumidity, start = 1984)
FFD_TS <- ts(FFD_df$FFD, start = 1984)


temp_next <- ts(FFD_prediction$Temperature,start =2015)
rain_next <- ts(FFD_prediction$Rainfall,start =2015)
rad_next <- ts(FFD_prediction$Radiation,start =2015)
hum_next <- ts(FFD_prediction$RelHumidity,start =2015)


FFD_df_TS %>% head()
class(FFD_df_TS)


##### 1. Temperature
temp_2_TS %>% head()
class(temp_2_TS)


##### 2. Rainfall
rain_TS %>% head()
class(rain_TS)


##### 3. Radiation Level
radiation_TS %>% head()
class(radiation_TS)

##### 4. Relative Humidity
humidity_TS %>% head()
class(humidity_TS)

##### 5. Fist Day of Flowering
FFD_TS %>% head()
class(FFD_TS)


## Descriptive Analysis

##### 1. First Day of Flowering
descriptive_analysis(FFD_TS, "day of first flowering")


##### 2. Temperature
descriptive_analysis(temp_2_TS, "Averaged yearly temperature")

##### 3. Rainfall
descriptive_analysis(rain_TS, "Averaged  yearly rainfall")

##### 4. Radiation Level
descriptive_analysis(radiation_TS, "Averaged yearly radiation level")


##### 5. Relative Humidity
descriptive_analysis(humidity_TS, "Averaged yearly relative humidity")


##### 6. Combine Scaled Time Series Plot
combined <- scale(FFD_df_TS)

plot(combined, 
     plot.type="s", 
     col = c("#05386b","#f01b1d","#5c3c92","#d2601a","#1b6535"), 
     main = "Scaled Time Series Plot of plant's day of first flowering and 
     contemporaneous averaged yearly climate variables")

legend("topleft", 
       lty=1, 
       col = c("#05386b","#f01b1d","#5c3c92","#d2601a","#1b6535"), 
       c("Temperature", "Rainfall", "Radiation Level", "Relative Humidity" , "FFD"))

#Correlation
cor(FFD_df_TS) %>% round(3)

## Time Series Regression Methods
### Finite Distributed Lag Model
##### 1.Temperature
for (i in 1:8)
{
  model_dlm <- dlm(formula = FFD ~ Temperature, data = data.frame(FFD_df), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}

task2_temp1 <- dlm(formula = FFD ~ Temperature, 
                   data = data.frame(FFD_df),
                   q = 8)

summary_residual_analysis(task2_temp1)

#Temperature without Intercept
task2_temp1_a <- dlm(formula = FFD ~ Temperature -1, 
                     data = data.frame(FFD_df),
                     q = 8)

summary_residual_analysis(task2_temp1_a)


##### 2.Rainfall
for (i in 1:8)
{
  model_dlm <- dlm(formula = FFD ~ Rainfall, data = data.frame(FFD_df), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}

task2_rain1 <- dlm(formula = FFD ~ Rainfall, 
                   data = data.frame(FFD_df),
                   q = 8)

summary_residual_analysis(task2_rain1)

#Rainfall without Intercept
task2_rain1_a <- dlm(formula = FFD ~ Rainfall -1, 
                     data = data.frame(FFD_df),
                     q = 8)

summary_residual_analysis(task2_rain1_a)


##### 3.Radiation Level
for (i in 1:8)
{
  model_dlm <- dlm(formula = FFD ~ Radiation, data = data.frame(FFD_df), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}

task2_rad1 <- dlm(formula = FFD ~ Radiation, 
                  data = data.frame(FFD_df),
                  q = 8)

summary_residual_analysis(task2_rad1)

#Radiation Level without Intercept
task2_rad1_a <- dlm(formula = FFD ~ Radiation -1, 
                    data = data.frame(FFD_df),
                    q = 8)

summary_residual_analysis(task2_rad1_a)


##### 4.Relative Humidity
for (i in 1:8)
{
  model_dlm <- dlm(formula = FFD ~ RelHumidity, data = data.frame(FFD_df), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}

task2_hum1 <- dlm(formula = FFD ~ RelHumidity, 
                  data = data.frame(FFD_df),
                  q = 8)

summary_residual_analysis(task2_hum1)

#Relative Humidity without Intercept
task2_hum1_a <- dlm(formula = FFD ~ RelHumidity -1, 
                    data = data.frame(FFD_df),
                    q = 8)

summary_residual_analysis(task2_hum1_a)

#### Model Comparison
DLM_models <- c("Temperature", "Temperature without Intercept", 
                "Rainfall", "Rainfall without Intercept",
                "Radiation Level" , "Radiation Level without Intercept",
                "Relative Humidity","Relative Humidity without Intercept")

DLM_aic <- AIC(task2_temp1$model, task2_temp1_a$model, task2_rain1$model, task2_rain1_a$model, task2_rad1$model,
               task2_rad1_a$model, task2_hum1$model, task2_hum1_a$model)$AIC

DLM_bic <- BIC(task2_temp1$model, task2_temp1_a$model, task2_rain1$model, task2_rain1_a$model, task2_rad1$model,
               task2_rad1_a$model, task2_hum1$model, task2_hum1_a$model)$BIC

DLM_mase <- MASE(task2_temp1$model, task2_temp1_a$model, task2_rain1$model, task2_rain1_a$model, task2_rad1$model,
                 task2_rad1_a$model, task2_hum1$model, task2_hum1_a$model)$MASE

DLM_Model_Comparison <- data.frame(DLM_models, DLM_mase, DLM_aic, DLM_bic)
colnames(DLM_Model_Comparison) <- c("DLM Model","MASE","AIC", "BIC")

kbl(DLM_Model_Comparison) %>% kable_paper()



### PolyNomial Distributed Lag Model
##### 1.Temperature
task_2_temp2 = polyDlm(y = as.vector(FFD_df$FFD), 
                       x = as.vector(FFD_df$Temperature), 
                       q = 8, 
                       k = 2, 
                       show.beta = TRUE)

summary_residual_analysis(task_2_temp2)

##### 2.Rainfall
task_2_rain2 = polyDlm(y = as.vector(FFD_df$FFD), 
                       x = as.vector(FFD_df$Rainfall), 
                       q = 8, 
                       k = 2, 
                       show.beta = TRUE)

summary_residual_analysis(task_2_rain2)

##### 3. Radiation Level
task_2_rad2 = polyDlm(y = as.vector(FFD_df$FFD), 
                      x = as.vector(FFD_df$Radiation), 
                      q = 8, 
                      k = 2, 
                      show.beta = TRUE)

summary_residual_analysis(task_2_rad2)

##### 4. Relative Humidity
task_2_hum2 = polyDlm(y = as.vector(FFD_df$FFD), 
                      x = as.vector(FFD_df$RelHumidity), 
                      q = 8, 
                      k = 2, 
                      show.beta = TRUE)

summary_residual_analysis(task_2_hum2)

#### Model Comparison
poly_models <- c("Temperature", "Rainfall", "Radiation Level", "Relative Humidity")

poly_aic <- AIC(task_2_temp2$model, task_2_rain2$model, task_2_rad2$model, task_2_hum2$model)$AIC

poly_bic <- BIC(task_2_temp2$model, task_2_rain2$model, task_2_rad2$model, task_2_hum2$model)$BIC

poly_mase <- MASE(task_2_temp2$model, task_2_rain2$model, task_2_rad2$model, task_2_hum2$model)$MASE

poly_Model_Comparison <- data.frame(poly_models, poly_mase, poly_aic, poly_bic)
colnames(poly_Model_Comparison) <- c("Poly Model","MASE","AIC", "BIC")

kbl(poly_Model_Comparison) %>% kable_paper()


### Koyck Distributed Lag Model
##### 1.Temperature
task_2_temp3 = koyckDlm(x = as.vector(FFD_df$Temperature), 
                        y = as.vector(FFD_df$FFD))
summary_residual_analysis(task_2_temp3)


##### 2.Rainfall
task_2_rain3 = koyckDlm(x = as.vector(FFD_df$Rainfall), 
                        y = as.vector(FFD_df$FFD))
summary_residual_analysis(task_2_rain3)

##### 3.Radiation Level
task_2_rad3 = koyckDlm(x = as.vector(FFD_df$Radiation), 
                       y = as.vector(FFD_df$FFD))
summary_residual_analysis(task_2_rad3)

##### 4.Relative Humidity
task_2_hum3 = koyckDlm(x = as.vector(FFD_df$RelHumidity), 
                       y = as.vector(FFD_df$FFD))
summary_residual_analysis(task_2_hum3)

#### Model Comparison
koyck_models <- c("Temperature", "Rainfall", "Radiation Level", "Relative Humidity")

koyck_aic <- c(AIC(task_2_temp3), AIC(task_2_rain3), AIC(task_2_rad3), AIC(task_2_hum3))

koyck_bic <- c(BIC(task_2_temp3), BIC(task_2_rain3), BIC(task_2_rad3), BIC(task_2_hum3))

koyck_mase <- MASE(lm(task_2_temp3$model), lm(task_2_rain3$model), lm(task_2_rad3$model), lm(task_2_hum3$model))$MASE

koyck_Model_Comparison <- data.frame(koyck_models, koyck_mase, koyck_aic, koyck_bic)
colnames(koyck_Model_Comparison) <- c("Koyck Model","MASE","AIC", "BIC")

kbl(koyck_Model_Comparison) %>% kable_paper()


### Autoregressive Distributed Lag Model
##### 1. Temperature
for (i in 1:5)
{
  for(j in 1:5)
  { 
    model_ardlm <- ardlDlm(formula = FFD ~ Temperature, 
                           data = data.frame(FFD_df), 
                           p = i,
                           q = j)
    cat("p =", i, 
        "q =" , j,
        "AIC =", AIC(model_ardlm$model), 
        "BIC =", BIC(model_ardlm$model),
        "MASE =", MASE(model_ardlm)$MASE, "\n")
  }
}

#ardlDLM(3,1)
task_2_temp4_1 = ardlDlm(formula = FFD ~ Temperature, 
                         data = data.frame(FFD_df), 
                         p = 3,
                         q = 1)

summary_residual_analysis(task_2_temp4_1)


task_2_temp4_1_a = ardlDlm(formula = FFD ~ Temperature -1, 
                           data = data.frame(FFD_df), 
                           p = 3,
                           q = 1)

summary_residual_analysis(task_2_temp4_1_a)


#ardlDLM(4,1)
task_2_temp4_2 = ardlDlm(formula = FFD ~ Temperature, 
                         data = data.frame(FFD_df), 
                         p = 4,
                         q = 1)

summary_residual_analysis(task_2_temp4_2)



task_2_temp4_2_a = ardlDlm(formula = FFD ~ Temperature -1, 
                           data = data.frame(FFD_df), 
                           p = 4,
                           q = 1)

summary_residual_analysis(task_2_temp4_2_a)


#ardlDLM(5,1)
task_2_temp4_3 = ardlDlm(formula = FFD ~ Temperature,
                         data = data.frame(FFD_df),
                         p = 5,
                         q = 1)

summary_residual_analysis(task_2_temp4_3)


task_2_temp4_3_a = ardlDlm(formula = FFD ~ Temperature - 1, 
                           data = data.frame(FFD_df), 
                           p = 5,
                           q = 1)

summary_residual_analysis(task_2_temp4_3_a)


##### 2.Rainfall
for (i in 1:5)
{
  for(j in 1:5)
  { 
    model_ardlm <- ardlDlm(formula = FFD ~ Rainfall, 
                           data = data.frame(FFD_df), 
                           p = i,
                           q = j)
    cat("p =", i, 
        "q =" , j,
        "AIC =", AIC(model_ardlm$model), 
        "BIC =", BIC(model_ardlm$model),
        "MASE =", MASE(model_ardlm)$MASE, "\n")
  }
}

#ardlDLM(3,1)
task_2_rain4_1 = ardlDlm(formula = FFD ~ Rainfall, 
                         data = data.frame(FFD_df), 
                         p = 3,
                         q = 1)

summary_residual_analysis(task_2_rain4_1)



task_2_rain4_1_a = ardlDlm(formula = FFD ~ Rainfall -1, 
                           data = data.frame(FFD_df), 
                           p = 3,
                           q = 1)

summary_residual_analysis(task_2_rain4_1_a)


#ardlDLM(4,1)
task_2_rain4_2 = ardlDlm(formula = FFD ~ Rainfall, 
                         data = data.frame(FFD_df), 
                         p = 4,
                         q = 1)

summary_residual_analysis(task_2_rain4_2)



task_2_rain4_2_a = ardlDlm(formula = FFD ~ Rainfall -1, 
                           data = data.frame(FFD_df), 
                           p = 4,
                           q = 1)

summary_residual_analysis(task_2_rain4_2_a)


#ardlDLM(5,1)
task_2_rain4_3 = ardlDlm(formula = FFD ~ Rainfall, 
                         data = data.frame(FFD_df), 
                         p = 5,
                         q = 1)

summary_residual_analysis(task_2_rain4_3)


task_2_rain4_3_a = ardlDlm(formula = FFD ~ Rainfall -1, 
                           data = data.frame(FFD_df), 
                           p = 5,
                           q = 1)

summary_residual_analysis(task_2_rain4_3_a)


##### 3.Radiation Level
for (i in 1:5)
{
  for(j in 1:5)
  { 
    model_ardlm <- ardlDlm(formula = FFD ~ Radiation, 
                           data = data.frame(FFD_df), 
                           p = i,
                           q = j)
    cat("p =", i, 
        "q =" , j,
        "AIC =", AIC(model_ardlm$model), 
        "BIC =", BIC(model_ardlm$model),
        "MASE =", MASE(model_ardlm)$MASE, "\n")
  }
}
#ardlDLM(3,1)
task_2_rad4_1 = ardlDlm(formula = FFD ~ Radiation, 
                        data = data.frame(FFD_df), 
                        p = 3,
                        q = 1)

summary_residual_analysis(task_2_rad4_1)

task_2_rad4_1_a = ardlDlm(formula = FFD ~ Radiation -1, 
                          data = data.frame(FFD_df), 
                          p = 3,
                          q = 1)

summary_residual_analysis(task_2_rad4_1_a)

#ardlDLM(4,1)
task_2_rad4_2 = ardlDlm(formula = FFD ~ Radiation, 
                        data = data.frame(FFD_df), 
                        p = 4,
                        q = 1)

summary_residual_analysis(task_2_rad4_2)

task_2_rad4_2_a = ardlDlm(formula = FFD ~ Radiation -1, 
                          data = data.frame(FFD_df), 
                          p = 4,
                          q = 1)

summary_residual_analysis(task_2_rad4_2_a)


#ardlDLM(5,1)
task_2_rad4_3 = ardlDlm(formula = FFD ~ Radiation, 
                        data = data.frame(FFD_df), 
                        p = 5,
                        q = 1)

summary_residual_analysis(task_2_rad4_3)

task_2_rad4_3_a = ardlDlm(formula = FFD ~ Radiation, 
                          data = data.frame(FFD_df), 
                          p = 5,
                          q = 1)

summary_residual_analysis(task_2_rad4_3_a)

##### 4.Relative Humidity
for (i in 1:5)
{
  for(j in 1:5)
  { 
    model_ardlm <- ardlDlm(formula = FFD ~ RelHumidity, 
                           data = data.frame(FFD_df), 
                           p = i,
                           q = j)
    cat("p =", i, 
        "q =" , j,
        "AIC =", AIC(model_ardlm$model), 
        "BIC =", BIC(model_ardlm$model),
        "MASE =", MASE(model_ardlm)$MASE, "\n")
  }
}

#ardlDLM(3,1)
task_2_hum4_1 = ardlDlm(formula = FFD ~ RelHumidity, 
                        data = data.frame(FFD_df), 
                        p = 3,
                        q = 1)

summary_residual_analysis(task_2_hum4_1)

task_2_hum4_1_a = ardlDlm(formula = FFD ~ RelHumidity -1, 
                          data = data.frame(FFD_df), 
                          p = 3,
                          q = 1)

summary_residual_analysis(task_2_hum4_1_a)

#ardlDLM(4,1)
task_2_hum4_2 = ardlDlm(formula = FFD ~ RelHumidity, 
                        data = data.frame(FFD_df), 
                        p = 4,
                        q = 1)

summary_residual_analysis(task_2_hum4_2)


task_2_hum4_2_a = ardlDlm(formula = FFD ~ RelHumidity -1, 
                          data = data.frame(FFD_df), 
                          p = 4,
                          q = 1)

summary_residual_analysis(task_2_hum4_2_a)


#ardlDLM(5,1)
task_2_hum4_3 = ardlDlm(formula = FFD ~ RelHumidity, 
                        data = data.frame(FFD_df), 
                        p = 5,
                        q = 1)

summary_residual_analysis(task_2_hum4_3)



task_2_hum4_3_a = ardlDlm(formula = FFD ~ RelHumidity -1, 
                          data = data.frame(FFD_df), 
                          p = 5,
                          q = 1)

summary_residual_analysis(task_2_hum4_3_a)


#### Model Comparison
ARDL_models <- c("Temperature(3,1)", "Temperature(3,1) without Intercept",
                 "Temperature(4,1)", "Temperature(4,1) without Intercept",
                 "Temperature(5,1)", "Temperature(5,1) without Intercept",
                 "Rainfall(3,1)", "Rainfall(3,1) without Intercept",
                 "Rainfall(4,1)", "Rainfall(4,1) without Intercept",
                 "Rainfall(5,1)", "Rainfall(5,1) without Intercept",
                 "Radiation Level(3,1)", "Radiation Level(3,1) without Intercept",
                 "Radiation Level(4,1)", "Radiation Level(4,1) without Intercept",
                 "Radiation Level(5,1)", "Radiation Level(5,1) without Intercept",
                 "Relative Humidity(3,1)", "Relative Humidity(3,1) without Intercept",
                 "Relative Humidity(4,1)", "Relative Humidity(4,1) without Intercept",
                 "Relative Humidity(5,1)", "Relative Humidity(5,1) without Intercept")

ARDL_aic <- AIC(task_2_temp4_1$model, task_2_temp4_1_a$model,
                task_2_temp4_2$model, task_2_temp4_2_a$model,
                task_2_temp4_3$model, task_2_temp4_3_a$model,
                task_2_rain4_1$model, task_2_rain4_1_a$model,
                task_2_rain4_2$model, task_2_rain4_2_a$model,
                task_2_rain4_3$model, task_2_rain4_3_a$model,
                task_2_rad4_1$model, task_2_rad4_1_a$model,
                task_2_rad4_2$model, task_2_rad4_2_a$model,
                task_2_rad4_3$model, task_2_rad4_3_a$model,
                task_2_hum4_1$model, task_2_hum4_1_a$model,
                task_2_hum4_2$model, task_2_hum4_2_a$model,
                task_2_hum4_3$model, task_2_hum4_3_a$model)$AIC

ARDL_bic <- BIC(task_2_temp4_1$model, task_2_temp4_1_a$model,
                task_2_temp4_2$model, task_2_temp4_2_a$model,
                task_2_temp4_3$model, task_2_temp4_3_a$model,
                task_2_rain4_1$model, task_2_rain4_1_a$model,
                task_2_rain4_2$model, task_2_rain4_2_a$model,
                task_2_rain4_3$model, task_2_rain4_3_a$model,
                task_2_rad4_1$model, task_2_rad4_1_a$model,
                task_2_rad4_2$model, task_2_rad4_2_a$model,
                task_2_rad4_3$model, task_2_rad4_3_a$model,
                task_2_hum4_1$model, task_2_hum4_1_a$model,
                task_2_hum4_2$model, task_2_hum4_2_a$model,
                task_2_hum4_3$model, task_2_hum4_3_a$model)$BIC

ARDL_mase <- MASE(task_2_temp4_1, task_2_temp4_1_a, 
                  task_2_temp4_2, task_2_temp4_2_a,
                  task_2_temp4_3,task_2_temp4_3_a,
                  task_2_rain4_1, task_2_rain4_1_a,
                  task_2_rain4_2, task_2_rain4_2_a, 
                  task_2_rain4_3, task_2_rain4_3_a,
                  task_2_rad4_1, task_2_rad4_1_a,
                  task_2_rad4_2, task_2_rad4_2_a,
                  task_2_rad4_3, task_2_rad4_3_a,
                  task_2_hum4_1, task_2_hum4_1_a,
                  task_2_hum4_2, task_2_hum4_2_a,
                  task_2_hum4_3, task_2_hum4_3_a)$MASE


ARDL_Model_Comparison <- data.frame(ARDL_models, ARDL_mase, ARDL_aic, ARDL_bic)
colnames(ARDL_Model_Comparison) <- c("ARDL Model","MASE","AIC", "BIC")

kbl(ARDL_Model_Comparison) %>% kable_paper()



### Dynamic Models

##### 1.FFD
#dynlm_1
task_2_model5_1 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + trend(FFD_TS))
summary(task_2_model5_1)
checkresiduals(task_2_model5_1)

task_2_model5_1_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + trend(FFD_TS) -1)
summary(task_2_model5_1_a)
checkresiduals(task_2_model5_1_a)

#dynlm_2
task_2_model5_2 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + trend(FFD_TS))
summary(task_2_model5_2)
checkresiduals(task_2_model5_2)

task_2_model5_2_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + trend(FFD_TS) -1)
summary(task_2_model5_2_a)
checkresiduals(task_2_model5_2_a)


##### 2.Temperature
#dynlm_3
task_2_temp5_1 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + temp_2_TS + trend(FFD_TS))
summary(task_2_model5_1)
checkresiduals(task_2_model5_1)

task_2_temp5_1_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + temp_2_TS + trend(FFD_TS) -1)
summary(task_2_model5_1_a)
checkresiduals(task_2_model5_1_a)


#dynlm_4
task_2_temp5_2 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + temp_2_TS)
summary(task_2_temp5_2)
checkresiduals(task_2_temp5_2)

task_2_temp5_2_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + temp_2_TS -1)
summary(task_2_temp5_2_a)
checkresiduals(task_2_temp5_2_a)

##### 3.Rainfall
#dynlm_5
task_2_rain5_1 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + rain_TS + trend(FFD_TS))
summary(task_2_rain5_1)
checkresiduals(task_2_rain5_1)

task_2_rain5_1_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + rain_TS + trend(FFD_TS) -1)
summary(task_2_rain5_1_a)
checkresiduals(task_2_rain5_1_a)


#dynlm_6
task_2_rain5_2 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + rain_TS + trend(FFD_TS) )
summary(task_2_rain5_2)
checkresiduals(task_2_rain5_2)

task_2_rain5_2_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + rain_TS + trend(FFD_TS) -1)
summary(task_2_rain5_2_a)
checkresiduals(task_2_rain5_2_a)


#dynlm_7
task_2_rain5_3 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + rain_TS + L(rain_TS, k = 1) + trend(FFD_TS) )
summary(task_2_rain5_3)
checkresiduals(task_2_rain5_3)

task_2_rain5_3_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + rain_TS + L(rain_TS, k = 1) + trend(FFD_TS) -1)
summary(task_2_rain5_3_a)
checkresiduals(task_2_rain5_3_a)


##### 4.Radiation Level
#dynlm_8
task_2_rad5_1 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + radiation_TS + trend(FFD_TS))
summary(task_2_rad5_1)
checkresiduals(task_2_rad5_1)

task_2_rad5_1_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + radiation_TS + trend(FFD_TS) -1)
summary(task_2_rad5_1_a)
checkresiduals(task_2_rad5_1_a)


#dynlm_9
task_2_rad5_2 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + radiation_TS)
summary(task_2_model5_2)
checkresiduals(task_2_rad5_2)

task_2_rad5_2_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + radiation_TS -1)
summary(task_2_model5_2_a)
checkresiduals(task_2_rad5_2_a)


##### 5.Relative Humidity
#dynlm_10
task_2_hum5_1 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + humidity_TS + trend(FFD_TS))
summary(task_2_hum5_1)
checkresiduals(task_2_hum5_1)

task_2_hum5_1_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + humidity_TS + trend(FFD_TS) -1)
summary(task_2_hum5_1_a)
checkresiduals(task_2_hum5_1_a)


#dynlm_11
task_2_hum5_2 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + humidity_TS)
summary(task_2_hum5_2)
checkresiduals(task_2_hum5_2)

task_2_hum5_2_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + humidity_TS -1)
summary(task_2_hum5_2_a)
checkresiduals(task_2_hum5_2_a)


##### 6.Combination of Predictors
#dynlm_12
task_2_comb5_1 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + temp_2_TS + rain_TS + trend(FFD_TS))
summary(task_2_comb5_1)
checkresiduals(task_2_comb5_1)

task_2_comb5_1_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + temp_2_TS + rain_TS + trend(FFD_TS) -1)
summary(task_2_comb5_1_a)
checkresiduals(task_2_comb5_1_a)


#dynlm_13
task_2_comb5_2 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + temp_2_TS + rain_TS)
summary(task_2_comb5_2)
checkresiduals(task_2_comb5_2

task_2_comb5_2_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + temp_2_TS + rain_TS -1)
summary(task_2_comb5_2_a)
checkresiduals(task_2_comb5_2_a)


#dynlm_14
task_2_comb5_3 = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + temp_2_TS + rain_TS + trend(FFD_TS))
summary(task_2_comb5_3)
checkresiduals(task_2_comb5_3)

task_2_comb5_3_a = dynlm(FFD_TS ~ L(FFD_TS , k = 1 ) + L(FFD_TS , k = 2 ) + temp_2_TS + rain_TS + trend(FFD_TS) -1)
summary(task_2_comb5_3_a)
checkresiduals(task_2_comb5_3_a)


#### Model Comparison
dynlm_models <- c("dynlm_1", "dynlm_1 without Intercept",
                  "dynlm_2", "dynlm_2 without Intercept",
                  "dynlm_3", "dynlm_3 without Intercept",
                  "dynlm_4", "dynlm_4 without Intercept",
                  "dynlm_5", "dynlm_5 without Intercept",
                  "dynlm_6", "dynlm_6 without Intercept",
                  "dynlm_7", "dynlm_7 without Intercept",
                  "dynlm_8", "dynlm_8 without Intercept",
                  "dynlm_9", "dynlm_9 without Intercept",
                  "dynlm_10", "dynlm_10 without Intercept",
                  "dynlm_11", "dynlm_11 without Intercept",
                  "dynlm_12", "dynlm_12 without Intercept",
                  "dynlm_13", "dynlm_13 without Intercept",
                  "dynlm_14", "dynlm_14 without Intercept")

dynlm_x <- c(rep("FFD", 4), rep("Temperature", 4), rep("Rainfall", 6), rep("Radiation Level", 4), 
             rep("Relative Humidity", 4), rep("Temperature + Rainfall", 6))

dynlm_aic <- c(AIC(task_2_model5_1), AIC(task_2_model5_1_a), 
               AIC(task_2_model5_2), AIC(task_2_model5_2_a),
               AIC(task_2_temp5_1), AIC(task_2_temp5_1_a),
               AIC(task_2_temp5_2), AIC(task_2_temp5_2_a),
               AIC(task_2_rain5_1), AIC(task_2_rain5_1_a),
               AIC(task_2_rain5_2), AIC(task_2_rain5_2_a), 
               AIC(task_2_rain5_3), AIC(task_2_rain5_3_a),
               AIC(task_2_rad5_1), AIC(task_2_rad5_1_a),
               AIC(task_2_rad5_2), AIC(task_2_rad5_2_a),
               AIC(task_2_hum5_1), AIC(task_2_hum5_1_a),
               AIC(task_2_hum5_2), AIC(task_2_hum5_2_a),
               AIC(task_2_comb5_1), AIC(task_2_comb5_1_a),
               AIC(task_2_comb5_2), AIC(task_2_comb5_2_a),
               AIC(task_2_comb5_3), AIC(task_2_comb5_3_a))

dynlm_bic <- c(BIC(task_2_model5_1), BIC(task_2_model5_1_a), 
               BIC(task_2_model5_2), BIC(task_2_model5_2_a),
               BIC(task_2_temp5_1), BIC(task_2_temp5_1_a),
               BIC(task_2_temp5_2), BIC(task_2_temp5_2_a),
               BIC(task_2_rain5_1), BIC(task_2_rain5_1_a),
               BIC(task_2_rain5_2), BIC(task_2_rain5_2_a), 
               BIC(task_2_rain5_3), BIC(task_2_rain5_3_a),
               BIC(task_2_rad5_1), BIC(task_2_rad5_1_a),
               BIC(task_2_rad5_2), BIC(task_2_rad5_2_a),
               BIC(task_2_hum5_1), BIC(task_2_hum5_1_a),
               BIC(task_2_hum5_2), BIC(task_2_hum5_2_a),
               BIC(task_2_comb5_1), BIC(task_2_comb5_1_a),
               BIC(task_2_comb5_2), BIC(task_2_comb5_2_a),
               BIC(task_2_comb5_3), BIC(task_2_comb5_3_a))

dynlm_mase <- c(MASE(lm(task_2_model5_1))$MASE, MASE(lm(task_2_model5_1_a))$MASE, 
                MASE(lm(task_2_model5_2))$MASE, MASE(lm(task_2_model5_2_a))$MASE,
                MASE(lm(task_2_temp5_1))$MASE, MASE(lm(task_2_temp5_1_a))$MASE,
                MASE(lm(task_2_temp5_2))$MASE, MASE(lm(task_2_temp5_2_a))$MASE,
                MASE(lm(task_2_rain5_1))$MASE, MASE(lm(task_2_rain5_1_a))$MASE,
                MASE(lm(task_2_rain5_2))$MASE, MASE(lm(task_2_rain5_2_a))$MASE, 
                MASE(lm(task_2_rain5_3))$MASE, MASE(lm(task_2_rain5_3_a))$MASE,
                MASE(lm(task_2_rad5_1))$MASE, MASE(lm(task_2_rad5_1_a))$MASE,
                MASE(lm(task_2_rad5_2))$MASE, MASE(lm(task_2_rad5_2_a))$MASE,
                MASE(lm(task_2_hum5_1))$MASE, MASE(lm(task_2_hum5_1_a))$MASE,
                MASE(lm(task_2_hum5_2))$MASE, MASE(lm(task_2_hum5_2_a))$MASE,
                MASE(lm(task_2_comb5_1))$MASE, MASE(lm(task_2_comb5_1_a))$MASE,
                MASE(lm(task_2_comb5_2))$MASE, MASE(lm(task_2_comb5_2_a))$MASE,
                MASE(lm(task_2_comb5_3))$MASE, MASE(lm(task_2_comb5_3_a))$MASE)

dynlm_Model_Comparison <- data.frame(dynlm_models, dynlm_x, dynlm_mase, dynlm_aic, dynlm_bic)
colnames(dynlm_Model_Comparison) <- c("dynlm Model","Predictors","MASE","AIC", "BIC")

kbl(dynlm_Model_Comparison) %>% kable_paper()




## Exponential smoothing methods
T2_holt_models = c("Damped Holt's method with exponential trend",
                   "Damped Holt's method",
                   "Holt's method with exponential trend",
                   "Holt's method")

T2_exponential = c(TRUE,FALSE)
T2_damped = c(TRUE,FALSE)
T2_exponential_models <- expand.grid(T2_exponential, T2_damped)
T2_Holt_AIC <- array(NA, 4)
T2_Holt_BIC <- array(NA, 4)
T2_Holt_MASE <- array(NA, 4)
T2_levels <- array(NA, dim=c(4,2))


for (i in 1:4)
{
  T2_Holt_model <- holt(FFD_TS,
                        exponential = T2_exponential_models[i,1],
                        damped = T2_exponential_models[i,2],
                        h = 4)
  
  T2_Holt_AIC[i] <- T2_Holt_model$model$aic
  T2_Holt_BIC[i] <- T2_Holt_model$model$bic
  T2_Holt_MASE[i] <- accuracy(T2_Holt_model)[6]
  T2_levels[i,1] <- T2_exponential_models[i,1]
  T2_levels[i,2] <- T2_exponential_models[i,2]
  summary(T2_Holt_model)
  checkresiduals(T2_Holt_model)
  print(shapiro.test(T2_Holt_model$model$residuals))
}

T2_results_Holt = data.frame(T2_holt_models, T2_levels, T2_Holt_MASE, T2_Holt_AIC, T2_Holt_BIC)
colnames(T2_results_Holt) = c("Model", "Exponential","Damped","MASE","AIC", "BIC")

kbl(T2_results_Holt) %>% kable_paper()


#Formating T2_results_Holt table
T2_results_Holt$Damped <- factor(T2_results_Holt$Damped,
                                 levels = c(TRUE, FALSE),
                                 labels = c("damped"," "))

T2_results_Holt <- unite(T2_results_Holt,
                         "Model", c("Model","Damped"), sep = "_")

T2_results_Holt <- T2_results_Holt[,-c(2)]
kbl(T2_results_Holt) %>% kable_paper()


## State-space Models

T2_ets_models = c("AAN", "MAN", "MMN")
T2_damped = c(TRUE,FALSE)
T2_ETS_models <- expand.grid(T2_ets_models, T2_damped)

T2_ETS_AIC <- array(NA, 6)
T2_ETS_BIC <- array(NA, 6)
T2_ETS_MASE <- array(NA, 6)
T2_levels <- array(NA, dim=c(6,2))

for (i in 1:6){
  T2_ETS <- ets(FFD_TS, 
                model = toString(T2_ETS_models[i, 1]), 
                damped = T2_ETS_models[i,2])
  T2_ETS_AIC[i] <- T2_ETS$aic
  T2_ETS_BIC[i] <- T2_ETS$bic
  T2_ETS_MASE[i] <- accuracy(T2_ETS)[6]
  T2_levels[i,1] <- toString(T2_ETS_models[i,1])
  T2_levels[i,2] <- T2_ETS_models[i,2]
  print(summary(T2_ETS))
  checkresiduals(T2_ETS)
  print(shapiro.test(T2_ETS$residuals))
}

T2_results_ETS = data.frame(T2_levels, T2_ETS_MASE, T2_ETS_AIC, T2_ETS_BIC)
colnames(T2_results_ETS) = c("Model","Damped","MASE","AIC", "BIC")

kbl(T2_results_ETS) %>% kable_paper()


#Formating T2_results_ETS table 
T2_results_ETS$Damped <- factor(T2_results_ETS$Damped,
                                levels = c(TRUE, FALSE),
                                labels = c("damped"," "))

T2_results_ETS <- unite(T2_results_ETS,
                        "Model", c("Model","Damped"), sep = "_")

kbl(T2_results_ETS) %>% kable_paper()


## Forecasting
### Finite Distributed Lag Model
#Following table display the Finite Distributed Lag modesl with respective `MASE, AIC, BIC`  in ascending order of `MASE`
DLM_sorted_MASE <- DLM_Model_Comparison %>% arrange(MASE)

kbl(DLM_sorted_MASE) %>% kable_paper()

Finite DLM model with Intercept and `Rainfall` as predictor gives the lowest `MASE(), AIC(), BIC()` 

Prediction Chart
  
DLM_forecasting <-dLagM::forecast(task2_rain1,x = as.vector(rain_next),h=4)


plot(ts(c(as.vector(FFD_TS), as.vector(DLM_forecasting$forecasts)),start =1984),
     col ="#f01b1d",
     xlim = c(1984,2019.5),type = 'l',
     ylab ="FFD", 
     main ="Forecast of Fisrt Flowreing Day(DLM), x=rainfall")

lines(FFD_TS,type = 'l')

legend("topleft", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))


#Prediction Points
t2_Time_point <- c("2015", "2016", "2017", "2018")
t2_DLM_points <- DLM_forecasting$forecasts

t2_DLM_prediction_points <- data.frame(t2_Time_point, t2_DLM_points)
colnames(t2_DLM_prediction_points) = c("Time","Point Forecast")

kable(t2_DLM_prediction_points) %>% kable_paper()


### PolyNomial Distributed Lag Model
#Following table display the Polynomial Distributed Lag modesl with respective `MASE, AIC, BIC`  in ascending order of `MASE`
poly_sorted_MASE <- poly_Model_Comparison %>% arrange(MASE)

kbl(poly_sorted_MASE) %>% kable_paper()

#Prediction Chart
poly_forecasting <- dLagM::forecast(task_2_rain2,x = as.vector(rain_next),h=4)

plot(ts(c(as.vector(FFD_TS), as.vector(poly_forecasting$forecasts)),start =1984),
     col ="#f01b1d",
     xlim = c(1984,2019.5),type = 'l',
     ylab ="FFD", 
     main ="Forecast of Fisrt Flowreing Day (polyDLM), x=rainfall")

lines(FFD_TS,type = 'l')

legend("topleft", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))


#Prediction Points
t2_interval_poly = dLagM::forecast(interval = TRUE, task_2_rain2, x = as.vector(rain_next),h=4)
t2_interval_poly$forecasts


### Koyck Distributed Lag Model
#Following table display the Koyck Distributed Lag modesl with respective `MASE, AIC, BIC`  in ascending order of `MASE`
koyck_sorted_MASE <- koyck_Model_Comparison %>% arrange(MASE)

kbl(koyck_sorted_MASE) %>% kable_paper()

#Prediction Chart
koyck_forecasting <- dLagM::forecast(task_2_temp3,x = as.vector(temp_next),h=4)


plot(ts(c(as.vector(FFD_TS), as.vector(koyck_forecasting$forecasts)),start =1984),
     col ="#f01b1d",
     xlim = c(1984,2019.5),type = 'l',
     ylab ="FFD", 
     main ="Forecast of Fisrt Flowreing Day (Koyck DLM), x=temp")

lines(FFD_TS,type = 'l')

legend("topleft", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))

#Prediction Points
  
t2_interval_koyck = dLagM::forecast(interval = TRUE, task_2_temp3, x = as.vector(temp_next),h=4)
t2_interval_koyck$forecasts


### Autoregressive Distributed Lag Model
#Following table display the Autoregressive Distributed Lag modesl with respective `MASE, AIC, BIC`  in ascending order of `MASE`
ARDL_sorted_MASE <- ARDL_Model_Comparison %>% arrange(MASE)

kbl(ARDL_sorted_MASE) %>% kable_paper()

ardlDLM(3,1) with Intercept and `x = temperature` gives the `lowest MASE()`, `2nd lowest AIC()` and `3rd lowest BIC()` 
Prediction Chart
  
  
ardl_forecasting <- dLagM::forecast(task_2_temp4_1,x = as.vector(temp_next),h=4)


plot(ts(c(as.vector(FFD_TS), as.vector(ardl_forecasting$forecasts)),start =1984),
     col ="#f01b1d",
     xlim = c(1984,2019.5),type = 'l',
     ylab ="FFD", 
     main ="Forecast of Fisrt Flowreing Day (ARDL DLM(3,1), x=temp)")

lines(FFD_TS,type = 'l')

legend("topleft", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))

#Prediction Points
  
t2_interval_ardl = dLagM::forecast(interval = TRUE, task_2_temp4_1, x = as.vector(temp_next),h=4)
t2_interval_ardl$forecasts

### Dynamic Models
#Following table display the Dynamic Lag modesl with respective `MASE, AIC, BIC`  in ascending order of `MASE`
dynlm_sorted_MASE <- dynlm_Model_Comparison %>% arrange(MASE)

kbl(dynlm_sorted_MASE) %>% kable_paper()

# Prediction Chart
q = 4
n = nrow(task_2_rain5_3$model)
ffd.frc = array(NA , (n + q))
ffd.frc[1:n] = FFD_TS[3:length(FFD_TS)]
trend = array(NA,q)
trend.start = task_2_rain5_3$model[n,"trend(FFD_TS)"]
trend = seq(trend.start , trend.start + q/12, 1/12)


for(i in 1:q){
  data.new = c(1,ffd.frc[n-1+i],ffd.frc[n-2+i],1,1,trend[i])
  ffd.frc[n+i] = as.vector(task_2_rain5_3$coefficients) %*% data.new
}

plot(FFD_TS,xlim=c(1984,2018),
     ylab='FFD',xlab='Year',
     main = "Next 4 years prediction of First Flowering Day in specific plant")

lines(ts(ffd.frc[(n+1):(n+q)],start=c(2015)),col="#f01b1d")

legend("topleft", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))


#Prediction Points
t2_points_dynlm <- ffd.frc[(n+1):(n+q)]

t2_dynlm_prediction_points <- data.frame(t2_Time_point, t2_points_dynlm)
colnames(t2_dynlm_prediction_points) = c("Time","Point Forecast")


kable(t2_dynlm_prediction_points) %>% kable_paper()

t2_points_dynlm


### Exponential Smoothing Models
#Following table display the Exponential Smoothing modesl with respective `MASE, AIC, BIC`  in ascending order of `MASE`
Holt_sorted_MASE <- T2_results_Holt %>% arrange(MASE)

kbl(Holt_sorted_MASE) %>% kable_paper()

#Prediction Chart
prediction_holt <- holt(FFD_TS,
                        damped = FALSE,
                   h = 4)

plot(prediction_holt,
     main = "Next four years prediction of FFD using Exponential Smoothing Model",
     ylab = "FFD",
     fcol = "#f01b1d")
       
legend("topleft", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))

#Prediction Points
kbl(prediction_holt) %>% kable_paper()


### State-Space Models
#Following table display the State space modesl with respective `MASE, AIC, BIC`  in ascending order of `MASE`
ETS_sorted_MASE <- T2_results_ETS %>% arrange(MASE)

kbl(ETS_sorted_MASE) %>% kable_paper()

#Prediction Chart
prediction_ETS_2 <- ets(FFD_TS,
                        model="AAN",
                        damped = F)
prediction_ETS_2 <- forecast(prediction_ETS_2)

plot(prediction_ETS_2,
     main = "Next four years prediction of FFD using ETS(A,A,N) model",
     ylab = "FFD",
     fcol = "#f01b1d")

legend("topleft", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))


#Prediction Points
kbl(prediction_ETS_2) %>% kable_paper()

#-----------------------------------------------------------#
#------------------------Task - 3(a)------------------------#
#-----------------------------------------------------------#

## Data Preparation
RBO_df <- read_csv("D:/Study_Material/SEM_4/Forecasting/Assignment3/RBO.csv")
RBO_prediction <- read_csv("D:/Study_Material/SEM_4/Forecasting/Assignment3/Covariate x-values for Task 3.csv")

RBO_prediction = na.omit(RBO_prediction)

head(RBO_df)
class(RBO_df)

RBO_df <- RBO_df[,-c(1)]
RBO_df_TS <- ts(RBO_df, start = 1984)

RBO_TS <- ts(RBO_df$RBO, start = 1984)
temp_3_TS <- ts(RBO_df$Temperature, start = 1984)
rain_2_TS <- ts(RBO_df$Rainfall, start = 1984)
radiation_2_TS <- ts(RBO_df$Radiation, start = 1984)
humidity_2_TS <- ts(RBO_df$RelHumidity, start = 1984)


temp_2_next <- ts(RBO_prediction$Temperature,start =2015)
rain_2_next <- ts(RBO_prediction$Rainfall,start =2015)
rad_2_next <- ts(RBO_prediction$Radiation,start =2015)
hum_2_next <- ts(RBO_prediction$RelHumidity,start =2015)


RBO_df_TS %>% head()
class(RBO_df_TS)


##### 1.Temperature
temp_3_TS %>% head()
class(temp_3_TS)

##### 2.Rainfall
rain_2_TS %>% head()
class(rain_2_TS)

##### 3.Radiation Level
radiation_2_TS %>% head()
class(radiation_2_TS)

##### 4.Relative Humidity
humidity_2_TS %>% head()
class(humidity_2_TS)

##### 5.Rank-based Flowering-Order similarity metric (RBO)
RBO_TS %>% head()
class(RBO_TS)


## Descriptive Analysis
##### 1.Rank-based Flowering-Order (RBO)
descriptive_analysis(RBO_TS, "Rank-based Flowering-Order similarity metric (RBO)")


##### 2. Temperature
descriptive_analysis(temp_3_TS, "Averaged yearly temperature")

##### 3. Rainfall
descriptive_analysis(rain_2_TS, "Averaged  yearly rainfall")

##### 4. Radiation Level
descriptive_analysis(radiation_2_TS, "Averaged yearly radiation level")

##### 5. Relative Humidity
descriptive_analysis(humidity_2_TS, "Averaged yearly relative humidity")

##### 6. Combine Scaled Time Series Plot
combined <- scale(RBO_df_TS)

plot(combined, 
     plot.type="s", 
     col = c("#05386b","#f01b1d","#5c3c92","#d2601a","#1b6535"), 
     main = "Scaled Time Series Plot of RBO similarity values and 
     contemporaneous averaged yearly climate variables")

legend("topleft", 
       lty=1, 
       col = c("#05386b","#f01b1d","#5c3c92","#d2601a","#1b6535"), 
       c("RBO","Temperature", "Rainfall", "Radiation Level", "Relative Humidity"))

#Correlation
cor(RBO_df_TS) %>% round(3)

## Time Series Regression Methods
### Finite Distributed Lag Model
##### 1.Temperature
for (i in 1:8)
{
  model_dlm <- dlm(formula = RBO ~ Temperature, data = data.frame(RBO_df), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}

task3_temp1 <- dlm(formula = RBO ~ Temperature, 
                   data = data.frame(RBO_df),
                   q = 8)

summary_residual_analysis(task3_temp1)


##### 2.Rainfall
for (i in 1:8)
{
  model_dlm <- dlm(formula = RBO ~ Rainfall, data = data.frame(RBO_df), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}

task3_rain1 <- dlm(formula = RBO ~ Rainfall, 
                   data = data.frame(RBO_df),
                   q = 8)

summary_residual_analysis(task3_rain1)


##### 3.Radiation Level
for (i in 1:8)
{
  model_dlm <- dlm(formula = RBO ~ Radiation, data = data.frame(RBO_df), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}

task3_rad1 <- dlm(formula = RBO ~ Radiation, 
                  data = data.frame(RBO_df),
                  q = 8)

summary_residual_analysis(task3_rad1)

##### 4.Relative Humidity
for (i in 1:8)
{
  model_dlm <- dlm(formula = RBO ~ RelHumidity, data = data.frame(RBO_df), q = i)
  cat("q = ", i, 
      "AIC = ", AIC(model_dlm$model), 
      "BIC = ", BIC(model_dlm$model), 
      "MASE = ", MASE(model_dlm)$MASE, "\n")
}

task3_hum1 <- dlm(formula = RBO ~ RelHumidity, 
                  data = data.frame(RBO_df),
                  q = 8)

summary_residual_analysis(task3_hum1)


### PolyNomial Distributed Lag Model
##### 1.Temperature
task_3_temp2 = polyDlm(y = as.vector(RBO_df$RBO), 
                       x = as.vector(RBO_df$Temperature), 
                       q = 8, 
                       k = 2, 
                       show.beta = TRUE)

summary_residual_analysis(task_3_temp2)

##### 2.Rainfall
task_3_rain2 = polyDlm(y = as.vector(RBO_df$RBO), 
                       x = as.vector(RBO_df$Rainfall), 
                       q = 8, 
                       k = 2, 
                       show.beta = TRUE)

summary_residual_analysis(task_3_rain2)

##### 3. Radiation Level
task_3_rad2 = polyDlm(y = as.vector(RBO_df$RBO), 
                      x = as.vector(RBO_df$Radiation), 
                      q = 8, 
                      k = 2, 
                      show.beta = TRUE)

summary_residual_analysis(task_3_rad2)

##### 4. Relative Humidity
task_3_hum2 = polyDlm(y = as.vector(RBO_df$RBO), 
                      x = as.vector(RBO_df$RelHumidity), 
                      q = 8, 
                      k = 2, 
                      show.beta = TRUE)

summary_residual_analysis(task_3_hum2)


### Koyck Distributed Lag Model
##### 1.Temperature
task_3_temp3 = koyckDlm(x = as.vector(RBO_df$Temperature), 
                        y = as.vector(RBO_df$RBO))
summary_residual_analysis(task_3_temp3)

##### 2.Rainfall
task_3_rain3 = koyckDlm(x = as.vector(RBO_df$Rainfall), 
                        y = as.vector(RBO_df$RBO))
summary_residual_analysis(task_3_rain3)

##### 3.Radiation Level
task_3_rad3 = koyckDlm(x = as.vector(RBO_df$Radiation), 
                       y = as.vector(RBO_df$RBO))
summary_residual_analysis(task_3_rad3)


##### 4.Relative Humidity
task_3_hum3 = koyckDlm(x = as.vector(RBO_df$RelHumidity), 
                       y = as.vector(RBO_df$RBO))
summary_residual_analysis(task_3_hum3)

### Autoregressive Distributed Lag Model
##### 1. Temperature
for (i in 1:5)
{
  for(j in 1:5)
  { 
    model_ardlm <- ardlDlm(formula = RBO ~ Temperature, 
                           data = data.frame(RBO_df), 
                           p = i,
                           q = j)
    cat("p =", i, 
        "q =" , j,
        "AIC =", AIC(model_ardlm$model), 
        "BIC =", BIC(model_ardlm$model),
        "MASE =", MASE(model_ardlm)$MASE, "\n")
  }
}

#ardlDLM(1,3)
task_3_temp4_1 = ardlDlm(formula = RBO ~ Temperature, 
                         data = data.frame(RBO_df), 
                         p = 1,
                         q = 3)

summary_residual_analysis(task_3_temp4_1)

#ardlDLM(2,3)
task_3_temp4_2 = ardlDlm(formula = RBO ~ Temperature, 
                         data = data.frame(RBO_df), 
                         p = 2,
                         q = 3)

summary_residual_analysis(task_3_temp4_2)

#ardlDLM(3,3)
task_3_temp4_3 = ardlDlm(formula = RBO ~ Temperature, 
                         data = data.frame(RBO_df), 
                         p = 3,
                         q = 3)

summary_residual_analysis(task_3_temp4_3)


##### 2.Rainfall
for (i in 1:5)
{
  for(j in 1:5)
  { 
    model_ardlm <- ardlDlm(formula = RBO ~ Rainfall, 
                           data = data.frame(RBO_df), 
                           p = i,
                           q = j)
    cat("p =", i, 
        "q =" , j,
        "AIC =", AIC(model_ardlm$model), 
        "BIC =", BIC(model_ardlm$model),
        "MASE =", MASE(model_ardlm)$MASE, "\n")
  }
}

#ardlDLM(1,3)
task_3_rain4_1 = ardlDlm(formula = RBO ~ Rainfall, 
                         data = data.frame(RBO_df), 
                         p = 1,
                         q = 3)

summary_residual_analysis(task_3_rain4_1)

#ardlDLM(2,3)
task_3_rain4_2 = ardlDlm(formula = RBO ~ Rainfall, 
                         data = data.frame(RBO_df), 
                         p = 2,
                         q = 3)

summary_residual_analysis(task_3_rain4_2)

#ardlDLM(3,3)
task_3_rain4_3 = ardlDlm(formula = RBO ~ Rainfall, 
                         data = data.frame(RBO_df), 
                         p = 3,
                         q = 3)

summary_residual_analysis(task_3_rain4_3)



##### 3.Radiation Level
for (i in 1:5)
{
  for(j in 1:5)
  { 
    model_ardlm <- ardlDlm(formula = RBO ~ Radiation, 
                           data = data.frame(RBO_df), 
                           p = i,
                           q = j)
    cat("p =", i, 
        "q =" , j,
        "AIC =", AIC(model_ardlm$model), 
        "BIC =", BIC(model_ardlm$model),
        "MASE =", MASE(model_ardlm)$MASE, "\n")
  }
}

#ardlDLM(1,3)
task_3_rad4_1 = ardlDlm(formula = RBO ~ Radiation, 
                        data = data.frame(RBO_df), 
                        p = 1,
                        q = 3)

summary_residual_analysis(task_3_rad4_1)

#ardlDLM(2,3)
task_3_rad4_2 = ardlDlm(formula = RBO ~ Radiation, 
                        data = data.frame(RBO_df), 
                        p = 2,
                        q = 3)

summary_residual_analysis(task_3_rad4_2)

#ardlDLM(3,2)
task_3_rad4_3 = ardlDlm(formula = RBO ~ Radiation, 
                        data = data.frame(RBO_df), 
                        p = 3,
                        q = 2)

summary_residual_analysis(task_3_rad4_3)

##### 4.Relative Humidity
for (i in 1:5)
{
  for(j in 1:5)
  { 
    model_ardlm <- ardlDlm(formula = RBO ~ RelHumidity, 
                           data = data.frame(RBO_df), 
                           p = i,
                           q = j)
    cat("p =", i, 
        "q =" , j,
        "AIC =", AIC(model_ardlm$model), 
        "BIC =", BIC(model_ardlm$model),
        "MASE =", MASE(model_ardlm)$MASE, "\n")
  }
}

#ardlDLM(1,2)
task_3_hum4_1 = ardlDlm(formula = RBO ~ RelHumidity, 
                        data = data.frame(RBO_df), 
                        p = 1,
                        q = 2)

summary_residual_analysis(task_3_hum4_1)

#ardlDLM(2,3)
task_3_hum4_2 = ardlDlm(formula = RBO ~ RelHumidity, 
                        data = data.frame(RBO_df), 
                        p = 2,
                        q = 3)

summary_residual_analysis(task_3_hum4_2)

#ardlDLM(3,2)
task_3_hum4_3 = ardlDlm(formula = RBO ~ RelHumidity, 
                        data = data.frame(RBO_df), 
                        p = 3,
                        q = 2)

summary_residual_analysis(task_3_hum4_3)


## Model Comparison
attr(task_3_temp3$model,"class") = "lm"
attr(task_3_rain3$model,"class") = "lm"
attr(task_3_rad3$model,"class") = "lm"
attr(task_3_hum3$model,"class") = "lm"

T3_models <- c(rep("Finite DLM",4), rep("Poly DLM",4), rep("Koyck", 4), rep("ARDL",12))

T3_x <- c(rep(c("Temperature", "Rainfall", "Radiation Level", "Relative Humidity"),3), 
          "temp(1,3)", "temp(2,3)", "temp(3,3)",
          "rain(1,3)", "rain(2,3)", "rain(3,3)",
          "rad(1,3)", "rad(2,3)", "rad(3,2)",
          "hum(1,2)", "temp(2,3)", "hum(3,2)")

T3_aic <- AIC(task3_temp1$model, task3_rain1$model, task3_rad1$model, task3_hum1$model, 
              task_3_temp2$model, task_3_rain2$model, task_3_rad2$model, task_3_hum2$model,
              task_3_temp3$model, task_3_rain3$model, task_3_rad3$model, task_3_hum3$model,
              task_3_temp4_1$model, task_3_temp4_2$model, task_3_temp4_3$model,
              task_3_rain4_1$model, task_3_rain4_2$model, task_3_rain4_3$model,
              task_3_rad4_1$model, task_3_rad4_2$model, task_3_rad4_3$model,
              task_3_hum4_1$model, task_3_hum4_2$model, task_3_hum4_3$model)$AIC


T3_bic <- BIC(task3_temp1$model, task3_rain1$model, task3_rad1$model, task3_hum1$model, 
              task_3_temp2$model, task_3_rain2$model, task_3_rad2$model, task_3_hum2$model,
              task_3_temp3$model, task_3_rain3$model, task_3_rad3$model, task_3_hum3$model,
              task_3_temp4_1$model, task_3_temp4_2$model, task_3_temp4_3$model,
              task_3_rain4_1$model, task_3_rain4_2$model, task_3_rain4_3$model,
              task_3_rad4_1$model, task_3_rad4_2$model, task_3_rad4_3$model,
              task_3_hum4_1$model, task_3_hum4_2$model, task_3_hum4_3$model)$BIC

T3_mase <- MASE(task3_temp1$model, task3_rain1$model, task3_rad1$model, task3_hum1$model, 
                task_3_temp2$model, task_3_rain2$model, task_3_rad2$model, task_3_hum2$model,
                task_3_temp3$model, task_3_rain3$model, task_3_rad3$model, task_3_hum3$model,
                task_3_temp4_1, task_3_temp4_2, task_3_temp4_3,
                task_3_rain4_1, task_3_rain4_2, task_3_rain4_3,
                task_3_rad4_1, task_3_rad4_2, task_3_rad4_3,
                task_3_hum4_1, task_3_hum4_2, task_3_hum4_3)$MASE

T3_Model_Comparison <- data.frame(T3_models, T3_x, T3_mase, T3_aic, T3_bic)
colnames(T3_Model_Comparison) <- c("Model","Predictors", "MASE","AIC", "BIC")

kbl(T3_Model_Comparison) %>% kable_paper()


## Forecasting
### MASE as optimal Measure
T3_sorted_MASE <- T3_Model_Comparison %>% arrange(MASE)
kbl(T3_sorted_MASE) %>% kable_paper()

#Prediction Chart
DLM_forecasting_2 <-dLagM::forecast(task3_rad1,x = as.vector(rad_2_next),h=3)

plot(ts(c(as.vector(RBO_TS), as.vector(DLM_forecasting_2$forecasts)),start =1984),
     col ="#f01b1d",
     xlim = c(1984,2019.5),type = 'l',
     ylab = "RBO", 
     main = "Next 3 years Forecasting of Flowering-Order (RBO) using Finite DLM (Radiation Level)")

lines(RBO_TS,type = 'l')

legend("topright", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))

#Prediction Points
t3_DLM_points <- DLM_forecasting_2$forecasts

t3_DLM_prediction_points <- data.frame(t2_Time_point[1:3], t3_DLM_points)
colnames(t3_DLM_prediction_points) = c("Time","Point Forecast")

kable(t3_DLM_prediction_points) %>% kable_paper()


### AIC  as optimal Measure
T3_sorted_AIC <- T3_Model_Comparison %>% arrange(AIC)
kbl(T3_sorted_AIC) %>% kable_paper()

# Predication Chart
DLM_forecasting_3 <-dLagM::forecast(task_3_rad4_1,x = as.vector(rad_2_next),h=3)

plot(ts(c(as.vector(RBO_TS), as.vector(DLM_forecasting_3$forecasts)),start =1984),
     col = "#f01b1d",
     xlim = c(1984,2019.5),type = 'l',
     ylab = "RBO", 
     main = "3 years Forecasting of Flowering-Order (RBO) using ARDL (Radiation Level)")

lines(RBO_TS,type = 'l')

legend("topright", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))

#Prediction Points
t3_interval_AIC = dLagM::forecast(interval = TRUE,task_3_rad4_1, x = as.vector(rad_2_next),h=4)
t3_interval_AIC$forecasts


#---------------------------------------------------------#
#----------------------- Task - 3(b)----------------------#
#---------------------------------------------------------#

## Dynamic Models
T = 12
S.T = (1*seq(RBO_TS)>=T)
S.T.1 = Lag(S.T, +1)

#dynlm_1
task_3_model5_1 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + trend(RBO_TS) + S.T )
summary(task_2_model5_1)
checkresiduals(task_2_model5_1)

#dynlm_2
task_3_model5_2 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + trend(RBO_TS) + S.T + S.T.1)
summary(task_3_model5_2)
checkresiduals(task_3_model5_2)

##### 1.Temperature
#dynlm_3
task_3_temp5_1 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + temp_3_TS + trend(RBO_TS) + S.T)
summary(task_3_temp5_1)
checkresiduals(task_3_temp5_1)

#dynlm_4
task_3_temp5_2 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + temp_3_TS + S.T)
summary(task_3_temp5_2)
checkresiduals(task_3_temp5_2)

#dynlm_5
task_3_temp5_3 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + temp_3_TS + L(temp_3_TS, k = 1) + trend(RBO_TS) + S.T)
summary(task_3_temp5_3)
checkresiduals(task_3_temp5_3)

#dynlm_6
task_3_temp5_4 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + temp_3_TS + L(temp_3_TS, k = 1) + S.T + S.T.1 + L(temp_3_TS, k = 2) + trend(RBO_TS))
summary(task_3_temp5_4)
checkresiduals(task_3_temp5_4)

#dynlm_7
task_3_temp5_5 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + L(RBO_TS , k = 3 ) + temp_3_TS + L(temp_3_TS, k = 1) + L(temp_3_TS, k = 2) + S.T + S.T.1)
summary(task_3_temp5_5)
checkresiduals(task_3_temp5_5)

##### 2.Rainfall
#dynlm_8
task_3_rain5_1 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + rain_2_TS + trend(RBO_TS) + S.T)
summary(task_3_rain5_1)
checkresiduals(task_3_rain5_1)

#dynlm_9
task_3_rain5_2 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + rain_2_TS + trend(RBO_TS) + S.T + S.T.1)
summary(task_3_rain5_2)
checkresiduals(task_3_rain5_2)

#dynlm_10
task_3_rain5_3 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + rain_2_TS + L(rain_2_TS, k = 1) + trend(RBO_TS) + S.T)
summary(task_3_rain5_3)
checkresiduals(task_3_rain5_3)

#dynlm_11
task_3_rain5_4 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + L(RBO_TS , k = 3 ) + rain_2_TS + trend(RBO_TS) + S.T)
summary(task_3_rain5_4)
checkresiduals(task_3_rain5_4)

#dynlm_12
task_3_rain5_5 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + L(RBO_TS , k = 3 ) + rain_2_TS + L(rain_2_TS, k = 1) + trend(RBO_TS) + S.T + S.T.1)
summary(task_3_rain5_5)
checkresiduals(task_3_rain5_5)


##### 3.Radiation Level
#dynlm_13
task_3_rad5_1 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + radiation_2_TS + trend(RBO_TS) + S.T)
summary(task_3_rad5_1)
checkresiduals(task_3_rad5_1)

#dynlm_14
task_3_rad5_2 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + radiation_2_TS + trend(RBO_TS) + S.T + S.T.1)
summary(task_3_rad5_2)
checkresiduals(task_3_rad5_2)

#dynlm_15
task_3_rad5_3 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + radiation_2_TS + L(radiation_2_TS, k = 1) + trend(RBO_TS) + S.T + S.T.1)
summary(task_3_rad5_3)
checkresiduals(task_3_rad5_3)

#dynlm_16
task_3_rad5_4 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + L(RBO_TS , k = 3 ) + radiation_2_TS + trend(RBO_TS)+ S.T)
summary(task_3_rad5_4)
checkresiduals(task_3_rad5_4)

#dynlm_17
task_3_rad5_5 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + L(RBO_TS , k = 3 ) + radiation_2_TS + L(radiation_2_TS, k = 1) + S.T + S.T.1 + trend(RBO_TS))
summary(task_3_rad5_5)
checkresiduals(task_3_rad5_5)

#dynlm_18
task_3_rad5_6 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + L(RBO_TS , k = 3 ) + L(RBO_TS , k = 4 ) + radiation_2_TS + L(radiation_2_TS, k = 1) + trend(RBO_TS) + S.T + S.T.1)
summary(task_3_rad5_6)
checkresiduals(task_3_rad5_6)

#dynlm_19
task_3_rad5_7 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + L(RBO_TS , k = 3 ) + L(RBO_TS , k = 4 ) + radiation_2_TS + L(radiation_2_TS, k = 1) + L(radiation_2_TS, k = 2) + trend(RBO_TS) + S.T + S.T.1)
summary(task_3_rad5_7)
checkresiduals(task_3_rad5_7)


##### 4.Relative Humidity
#dynlm_20
task_3_hum5_1 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + humidity_2_TS + trend(RBO_TS) + S.T)
summary(task_3_hum5_1)
checkresiduals(task_3_hum5_1)

#dynlm_21
task_3_hum5_2 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + humidity_2_TS + trend(RBO_TS) + S.T)
summary(task_3_hum5_2)
checkresiduals(task_3_hum5_2)

#dynlm_22
task_3_hum5_3 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + humidity_2_TS + L(humidity_2_TS, k = 1) + trend(RBO_TS) + S.T + S.T.1)
summary(task_3_hum5_3)
checkresiduals(task_3_hum5_3)

#dynlm_23
task_3_hum5_4 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + L(RBO_TS , k = 3 ) + humidity_2_TS + S.T + S.T.1)
summary(task_3_hum5_4)
checkresiduals(task_3_hum5_4)

#dynlm_24
task_3_hum5_5 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + L(RBO_TS , k = 3 ) + humidity_2_TS + L(humidity_2_TS, k = 1) + S.T + S.T.1)
summary(task_3_hum5_5)
checkresiduals(task_3_hum5_5)

#dynlm_25
task_3_hum5_6 = dynlm(RBO_TS ~ L(RBO_TS , k = 1 ) + L(RBO_TS , k = 2 ) + L(RBO_TS , k = 3 ) + humidity_2_TS + L(humidity_2_TS, k = 1) + trend(RBO_TS) + S.T + S.T.1)
summary(task_3_hum5_6)
checkresiduals(task_3_hum5_6)


## Model Comparison
T3_b_models <- c("dynlm_1", "dynlm_2", "dynlm_3", "dynlm_4", "dynlm_5",
                 "dynlm_6", "dynlm_7", "dynlm_8", "dynlm_9", "dynlm_10", 
                 "dynlm_11", "dynlm_12", "dynlm_13", "dynlm_14", "dynlm_15", 
                 "dynlm_16", "dynlm_17", "dynlm_18", "dynlm_19", "dynlm_20",
                 "dynlm_21", "dynlm_22", "dynlm_23", "dynlm_24", "dynlm_25")

T3_b_x <- c(rep("FFD", 2), rep("Temperature", 5), rep("Rainfall", 5), rep("Radiation Level", 7), rep("Relative Humidity", 6)) 

T3_b_aic <- AIC(task_3_model5_1, task_3_model5_2, 
                task_3_temp5_1, task_3_temp5_2, task_3_temp5_3, task_3_temp5_4, task_3_temp5_5, 
                task_3_rain5_1, task_3_rain5_2, task_3_rain5_3, task_3_rain5_4, task_3_rain5_5,
                task_3_rad5_1, task_3_rad5_2, task_3_rad5_3, task_3_rad5_4, task_3_rad5_5, task_3_rad5_6, task_3_rad5_7,
                task_3_hum5_1, task_3_hum5_2, task_3_hum5_3, task_3_hum5_4, task_3_hum5_5, task_3_hum5_6)$AIC

T3_b_bic <- BIC(task_3_model5_1, task_3_model5_2, task_3_temp5_1,
                task_3_temp5_2, task_3_temp5_3, task_3_temp5_4, task_3_temp5_5, 
                task_3_rain5_1, task_3_rain5_2, task_3_rain5_3, task_3_rain5_4, task_3_rain5_5,
                task_3_rad5_1, task_3_rad5_2, task_3_rad5_3, task_3_rad5_4, task_3_rad5_5, task_3_rad5_6, task_3_rad5_7,
                task_3_hum5_1, task_3_hum5_2, task_3_hum5_3, task_3_hum5_4, task_3_hum5_5, task_3_hum5_6)$BIC

T3_b_mase <- MASE(lm(task_3_model5_1), lm(task_3_model5_2), lm(task_3_temp5_1),
                  lm(task_3_temp5_2), lm(task_3_temp5_3), lm(task_3_temp5_4), lm(task_3_temp5_5), 
                  lm(task_3_rain5_1), lm(task_3_rain5_2), lm(task_3_rain5_3), lm(task_3_rain5_4), lm(task_3_rain5_5),
                  lm(task_3_rad5_1), lm(task_3_rad5_2), lm(task_3_rad5_3), lm(task_3_rad5_4), lm(task_3_rad5_5), lm(task_3_rad5_6), lm(task_3_rad5_7),
                  lm(task_3_hum5_1), lm(task_3_hum5_2), lm(task_3_hum5_3), lm(task_3_hum5_4), lm(task_3_hum5_5), lm(task_3_hum5_6))$MASE

T3_b_Model_Comparison <- data.frame(T3_b_models, T3_b_x, T3_b_mase, T3_b_aic, T3_b_bic)
colnames(T3_b_Model_Comparison) <- c("Model","Predictors","MASE","AIC", "BIC")
kbl(T3_b_Model_Comparison) %>% kable_paper()


## Forecasting
T3_b_sorted_MASE <- T3_b_Model_Comparison %>% arrange(MASE)
kbl(T3_b_sorted_MASE) %>% kable_paper()

#Prediction Chart
q = 3
n = nrow(task_3_rad5_7$model)
rbo.frc = array(NA , (n + q))
rbo.frc[1:n] = RBO_TS[5:length(RBO_TS)]
trend = array(NA,q)
trend.start = task_3_rad5_7$model[n,"trend(RBO_TS)"]
trend = seq(trend.start , trend.start + q/12, 1/12)

for(i in 1:q){
  data.new =c(1,rbo.frc[n-1+i],rbo.frc[n-2+i],rbo.frc[n-3+i],rbo.frc[n-4+i],1,1,1,trend[i],1,1)
  rbo.frc[n+i] = as.vector(task_3_rad5_7$coefficients) %*% data.new
}

plot(RBO_TS,xlim=c(1984,2018), 
     ylab='RBO',xlab='Year', 
     main = "Next 3 years forecasting of RBO using dynlm_19 model")

lines(ts(rbo.frc[(n+1):(n+q)],start=c(2015)),col="#f01b1d")

legend("topright", 
       lty = 1, 
       col = c("black", "#f01b1d"), 
       c("Data", "Prediction"))

#Prediction Points
t3_points_dynlm <- rbo.frc[(n+1):(n+q)]

t3_dynlm_prediction_points <- data.frame(t2_Time_point[1:3], t3_points_dynlm)
colnames(t3_dynlm_prediction_points) = c("Time","Point Forecast")


kable(t3_dynlm_prediction_points) %>% kable_paper()
