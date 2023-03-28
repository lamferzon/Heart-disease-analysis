library(MASS)
library(ggplot2)
library(tseries)
library(e1071)
library(dplyr)
library(tree)
library(caret)
library(randomForest)
library(boot)
library(ROSE)
library(car)
library(randomForest)

rm(list = ls())

# ******************************************************************************
### Importation data from "Hearth disease dataset" ----
DS <- read.csv("C:/Users/loren/Desktop/Statistical Learning/Scientific paper/Heart disease dataset.csv",
               header=TRUE, sep=",", dec=".", stringsAsFactors = TRUE)
N <- nrow(DS)

# ******************************************************************************
### Functions ----
create_results_table <- function(data_table, row_names, N, flag = TRUE){
  if (flag == TRUE){
    results_table <- matrix(nrow = length(row_names), ncol = 12)
    colnames(results_table) <- c("Min", "Mean", "Max", "Std", "1st qu.", 
                                 "Median", "2nd qu.", "Kurtosis", "Skewness",
                                 "JB p-value", "CI_low", "CI_up")
  } else{
    results_table <- matrix(nrow = length(row_names), ncol = 4)
    colnames(results_table) <- c("Min", "Mean", "Max", "Std")
  }
  rownames(results_table) <- row_names
  cont <- dim(results_table)[1]
  for (i in 1:cont){
    results_table[i,1] <- round(min(data_table[,i]),2)
    results_table[i,2] <- round(mean(data_table[,i]),2)
    results_table[i,3] <- round(max(data_table[,i]),2)
    results_table[i,4] <- round(sd(data_table[,i]),2)
    if (flag == TRUE){
      results_table[i,5] <- round(quantile(data_table[,i],.25),2)
      results_table[i,6] <- round(median(data_table[,i]),2)
      results_table[i,7] <- round(quantile(data_table[,i],.75),2)
      results_table[i,8] <- round(kurtosis(data_table[,i]),2)
      results_table[i,9] <- round(skewness(data_table[,i]),2)
      results_table[i,10] <- round(jarque.bera.test(data_table[,i])$p.value,2)
      if (results_table[i,10] > 0.05){
        mean <- mean(data_table[,i])
        sigma <- sqrt(var(data_table[,i]))
        N <- nrow(data_table)
        t_percentile <- qt(1-(0.05/2),df = N - 1)
        results_table[i,11] <- round(mean - t_percentile*(sigma/sqrt(N)),2)
        results_table[i,12] <- round(mean + t_percentile*(sigma/sqrt(N)),2)
      }
    }
  }
  return(results_table)
}

remove_outliers <- function(data, alpha){
  BMI_Qs <- quantile(DS$BMI, c(.25, .75))
  BMI_IQR <- IQR(DS$BMI)
  BMI_cond <- DS$BMI > (BMI_Qs[1] - alpha*BMI_IQR) & 
    DS$BMI < (BMI_Qs[2] + alpha*BMI_IQR)
  PhysicalH_Qs <- quantile(DS$PhysicalHealth, c(.25, .75))
  PhysicalH_IQR <- IQR(DS$PhysicalHealth)
  PhysicalH_cond <- DS$PhysicalHealth > (PhysicalH_Qs[1] - alpha*PhysicalH_IQR) & 
    DS$PhysicalHealth < (PhysicalH_Qs[2] + alpha*PhysicalH_IQR)
  MentalH_Qs <- quantile(DS$MentalHealth, c(.25, .75))
  MentalH_IQR <- IQR(DS$MentalHealth)
  MentalH_cond <- DS$MentalHealth > (MentalH_Qs[1] - alpha*MentalH_IQR) & 
    DS$MentalHealth < (MentalH_Qs[2] + alpha*MentalH_IQR)
  SleepT_Qs <- quantile(DS$SleepTime, c(.25, .75))
  SleepT_IQR <- IQR(DS$SleepTime)
  SleepT_cond <- DS$SleepTime > (SleepT_Qs[1] - alpha*SleepT_IQR) & 
    DS$SleepTime < (SleepT_Qs[2] + alpha*SleepT_IQR)
  no_outliers_DS <- subset(DS, BMI_cond & PhysicalH_cond & MentalH_cond & SleepT_cond)
  return(no_outliers_DS)
}
# ******************************************************************************
### Outliers individuation ----
no_outliers_DS <- remove_outliers(DS,15)
nrow(DS) - nrow(no_outliers_DS)

# Statistics both with and without outliers
num_vars_names <- c("BMI", "PhysicalHealth", "MentalHealth", "SleepTime")
num_vars_DS_yes_outliers <- DS[, num_vars_names]
num_vars_statistics_yes_outliers <- create_results_table(num_vars_DS_yes_outliers,
                                                         num_vars_names, flag = TRUE)

num_vars_DS_no_outliers <- no_outliers_DS[, num_vars_names]
num_vars_statistics_no_outliers <- create_results_table(num_vars_DS_no_outliers,
                                                        num_vars_names, flag = TRUE)

# ******************************************************************************
### LOGISTIC REGRESSION
## LR: model selection, oversampling + undersampling, with outliers (TO REVIEW) -----
# First step: most complex model also with non-linear terms
N_rep <- 500
mean_pvalues_1 <- rep(0,48)
mean_aic_yes_outliers_1 <- 0
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "both")$data
  log_model_yes_outliers <- glm(HeartDisease ~. + I(BMI^2) + I(PhysicalHealth^2) + 
                          I(MentalHealth^2) + I(SleepTime^2) + 
                          BMI:PhysicalHealth + BMI:MentalHealth + BMI:SleepTime + 
                          PhysicalHealth:MentalHealth + PhysicalHealth:SleepTime + 
                          MentalHealth:SleepTime, data = balanced_DS, family = binomial)
  mean_pvalues_1 <- mean_pvalues_1 + summary(log_model_yes_outliers)$coefficients[,4]
  mean_aic_yes_outliers_1 <- mean_aic_yes_outliers_1 + summary(log_model_yes_outliers)$aic
}
mean_pvalues_1 <- mean_pvalues_1/N_rep
mean_aic_yes_outliers_1 <- round(mean_aic_yes_outliers_1/N_rep)
for (i in 1:48){
  if(mean_pvalues_1[i] > 0.05){
    print(mean_pvalues_1[i])
  }
}
# Second step: model complexity reduction 
dummy_age <- contrasts(DS$AgeCategory)
dummy_age[2,] <- rep(0,12)
dummy_diabetic <- contrasts(DS$Diabetic)
dummy_diabetic[4,] <- rep(0,3)
dummy_race <- contrasts(DS$Race)
dummy_race[5,] <- rep(0,5)
dummy_race[6,] <- rep(0,5)
mean_aic_yes_outliers_2 <- 0
mean_pvalues_2 <- rep(0,38)
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "both")$data
  log_model_yes_outliers <- glm(HeartDisease ~. - PhysicalActivity + I(BMI^2) + 
                          I(PhysicalHealth^2) + I(SleepTime^2) + 
                          PhysicalHealth:MentalHealth + MentalHealth:SleepTime, 
                        data = balanced_DS, family = binomial, 
                        contrasts = list("AgeCategory" = dummy_age,
                                         "Race" = dummy_race,
                                         "Diabetic" = dummy_diabetic))
  mean_pvalues_2 <- mean_pvalues_2 + summary(log_model_yes_outliers)$coefficients[,4]
  mean_aic_yes_outliers_2 <- mean_aic_yes_outliers_2 + summary(log_model_yes_outliers)$aic
}
mean_pvalues_2 <- mean_pvalues_2/N_rep
mean_aic_yes_outliers_2 <- round(mean_aic_yes_outliers_2/N_rep)
for (i in 1:38){
  if(mean_pvalues_2[i] > 0.05){
    print(mean_pvalues_2[i])
  }
}

# ******************************************************************************
## LR: model selection, oversampling + undersampling, without outliers (TO REVIEW) ----
# First step: most complex model also with non-linear terms
N_rep <- 100
mean_pvalues_1 <- rep(0,48)
mean_aic_no_outliers_1 <- 0
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS, 
                             method = "both")$data
  log_model_no_outliers <- glm(HeartDisease ~. + I(BMI^2) + I(PhysicalHealth^2) + 
                                  I(MentalHealth^2) + I(SleepTime^2) + 
                                  BMI:PhysicalHealth + BMI:MentalHealth + BMI:SleepTime + 
                                  PhysicalHealth:MentalHealth + PhysicalHealth:SleepTime + 
                                  MentalHealth:SleepTime, data = balanced_DS, family = binomial)
  mean_pvalues_1 <- mean_pvalues_1 + summary(log_model_no_outliers)$coefficients[,4]
  mean_aic_no_outliers_1 <- mean_aic_no_outliers_1 + summary(log_model_no_outliers)$aic
}
mean_pvalues_1 <- mean_pvalues_1/N_rep
mean_aic_no_outliers_1 <- round(mean_aic_no_outliers_1/N_rep)
for (i in 1:48){
  if(mean_pvalues_1[i] > 0.05){
    print(mean_pvalues_1[i])
  }
}

# Second step: model complexity reduction 1
dummy_age <- contrasts(DS$AgeCategory)
dummy_age[2,] <- rep(0,12)
mean_aic_no_outliers_2 <- 0
mean_pvalues_2 <- rep(0,39)
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS, 
                             method = "both")$data
  log_model_no_outliers <- glm(HeartDisease ~. - MentalHealth - PhysicalActivity + 
                                  I(PhysicalHealth^2) + I(MentalHealth^2) + 
                                  I(SleepTime^2) + BMI:SleepTime, data = balanced_DS, 
                                family = binomial, 
                                contrasts = list("AgeCategory" = dummy_age))
  mean_pvalues_2 <- mean_pvalues_2 + summary(log_model_no_outliers)$coefficients[,4]
  mean_aic_no_outliers_2 <- mean_aic_no_outliers_2 + summary(log_model_no_outliers)$aic
}
mean_pvalues_2 <- mean_pvalues_2/N_rep
mean_aic_no_outliers_2 <- round(mean_aic_no_outliers_2/N_rep)
for (i in 1:38){
  if(mean_pvalues_2[i] > 0.05){
    print(mean_pvalues_2[i])
  }
}

# Third step: model complexity reduction 2
mean_aic_no_outliers_3 <- 0
mean_pvalues_3 <- rep(0,38)
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS, 
                             method = "both")$data
  log_model_no_outliers <- glm(HeartDisease ~. - MentalHealth - PhysicalActivity + 
                                 I(PhysicalHealth^2) + I(SleepTime^2) + 
                                 BMI:SleepTime, data = balanced_DS, 
                               family = binomial, 
                               contrasts = list("AgeCategory" = dummy_age))
  mean_pvalues_3 <- mean_pvalues_3 + summary(log_model_no_outliers)$coefficients[,4]
  mean_aic_no_outliers_3 <- mean_aic_no_outliers_3 + summary(log_model_no_outliers)$aic
}
mean_pvalues_3 <- mean_pvalues_3/N_rep
mean_aic_no_outliers_3 <- round(mean_aic_no_outliers_3/N_rep)
for (i in 1:38){
  if(mean_pvalues_3[i] > 0.05){
    print(mean_pvalues_3[i])
  }
}

# ******************************************************************************
## LR: model selection, undersampling, with outliers ----
N_rep <- 500
US_log_mod_yes_out <- NULL
US_avg_pvalues_yes_out <- NULL
US_aic_yes_out <- NULL

# 1st step: all variables in the model
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,38)
temp_pvalues <- rep(0,38)
pvalues_count <- rep(0,38) 

temp_aic <- 0
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  temp_log_mod <- glm(HeartDisease ~., data = balanced_DS, family = binomial)
  temp_pvalues <- summary(temp_log_mod)$coefficients[,4]
  for (i in 1:length(temp_pvalues)){
    if (temp_pvalues[i] >= 0.01) {
      pvalues_count[i] = pvalues_count[i] + 1;
    }
  }
  temp_avg_pvalues <- temp_avg_pvalues + temp_pvalues
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
names(pvalues_count) <- names(temp_pvalues)
pvalues_count <- sort(pvalues_count, decreasing = T)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_yes_out[1] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_yes_out[1] <- list(temp_aic)
for (i in 1:38){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}
for (i in 1:38){
  print(pvalues_count[i])
}

# 2nd step: -Race (Other)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,37)
temp_aic <- 0
dummy_race <- contrasts(DS$Race)
dummy_race <- dummy_race[,-4]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  temp_log_mod <- glm(HeartDisease ~., data = balanced_DS, family = binomial,
                      contrasts = list("Race" = dummy_race))
  temp_avg_pvalues <- temp_avg_pvalues + summary(temp_log_mod)$coefficients[,4]
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
US_log_mod_yes_out[2] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_yes_out[2] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_yes_out[2] <- list(temp_aic)
for (i in 1:37){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}

# 3rd step: -PhysicalActivity 
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,36)
temp_aic <- 0
dummy_race <- contrasts(DS$Race)
dummy_race <- dummy_race[,-4]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  temp_log_mod <- glm(HeartDisease ~. - PhysicalActivity, data = balanced_DS, 
                      family = binomial,
                      contrasts = list("Race" = dummy_race))
  temp_avg_pvalues <- temp_avg_pvalues + summary(temp_log_mod)$coefficients[,4]
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
US_log_mod_yes_out[3] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_yes_out[3] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_yes_out[3] <- list(temp_aic)
for (i in 1:36){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}

# 4th step: -AgeCategory (25-29)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,35)
temp_aic <- 0
dummy_race <- contrasts(DS$Race)
dummy_race <- dummy_race[,-4]
dummy_age <- contrasts(DS$AgeCategory)
dummy_age <- dummy_age[,-1]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  temp_log_mod <- glm(HeartDisease ~. - PhysicalActivity, data = balanced_DS, 
                      family = binomial,
                      contrasts = list("Race" = dummy_race,
                                       "AgeCategory" = dummy_age))
  temp_avg_pvalues <- temp_avg_pvalues + summary(temp_log_mod)$coefficients[,4]
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
US_log_mod_yes_out[4] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_yes_out[4] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_yes_out[4] <- list(temp_aic)
for (i in 1:35){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}

# 5th step: -Diabetic (Yes, during pregnancy)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,34)
temp_aic <- 0
dummy_race <- contrasts(DS$Race)
dummy_race <- dummy_race[,-4]
dummy_age <- contrasts(DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(DS$Diabetic)
dummy_diabetic <- dummy_diabetic[,-3]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  temp_log_mod <- glm(HeartDisease ~. - PhysicalActivity, data = balanced_DS, 
                      family = binomial,
                      contrasts = list("Race" = dummy_race,
                                       "AgeCategory" = dummy_age,
                                       "Diabetic" = dummy_diabetic))
  temp_avg_pvalues <- temp_avg_pvalues + summary(temp_log_mod)$coefficients[,4]
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
US_log_mod_yes_out[5] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_yes_out[5] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_yes_out[5] <- list(temp_aic)
for (i in 1:34){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}

# 6th step: -Race (White)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,33)
temp_aic <- 0
dummy_race <- contrasts(DS$Race)
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-4]
dummy_age <- contrasts(DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(DS$Diabetic)
dummy_diabetic <- dummy_diabetic[,-3]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  temp_log_mod <- glm(HeartDisease ~. - PhysicalActivity, data = balanced_DS, 
                      family = binomial,
                      contrasts = list("Race" = dummy_race,
                                       "AgeCategory" = dummy_age,
                                       "Diabetic" = dummy_diabetic))
  temp_avg_pvalues <- temp_avg_pvalues + summary(temp_log_mod)$coefficients[,4]
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
US_log_mod_yes_out[6] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_yes_out[6] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_yes_out[6] <- list(temp_aic)
for (i in 1:33){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}

# 7th step: -Race (Hispanic)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,32)
temp_aic <- 0
dummy_race <- contrasts(DS$Race)
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-3]
dummy_age <- contrasts(DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(DS$Diabetic)
dummy_diabetic <- dummy_diabetic[,-3]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  temp_log_mod <- glm(HeartDisease ~. - PhysicalActivity, data = balanced_DS, 
                      family = binomial,
                      contrasts = list("Race" = dummy_race,
                                       "AgeCategory" = dummy_age,
                                       "Diabetic" = dummy_diabetic))
  temp_avg_pvalues <- temp_avg_pvalues + summary(temp_log_mod)$coefficients[,4]
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
US_log_mod_yes_out[7] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_yes_out[7] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_yes_out[7] <- list(temp_aic)
for (i in 1:32){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}

# 8th step: -Diabetic (No, borderline diabetes)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,31)
temp_aic <- 0
dummy_race <- contrasts(DS$Race)
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-3]
dummy_age <- contrasts(DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(DS$Diabetic)
dummy_diabetic <- dummy_diabetic[,-3]
dummy_diabetic <- dummy_diabetic[,-1]
dummy_diabetic <- t(dummy_diabetic)
dummy_diabetic <- t(dummy_diabetic)
colnames(dummy_diabetic) <- "Yes"
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  temp_log_mod <- glm(HeartDisease ~. - PhysicalActivity, data = balanced_DS, 
                      family = binomial,
                      contrasts = list("Race" = dummy_race,
                                       "AgeCategory" = dummy_age,
                                       "Diabetic" = dummy_diabetic))
  temp_avg_pvalues <- temp_avg_pvalues + summary(temp_log_mod)$coefficients[,4]
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
US_log_mod_yes_out[8] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_yes_out[8] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_yes_out[8] <- list(temp_aic)
for (i in 1:31){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}

# 9th step: previous (best) model with quadratic terms
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,35)
temp_aic <- 0
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  temp_log_mod <- glm(HeartDisease ~. - PhysicalActivity + I(BMI^2) 
                      + I(PhysicalHealth^2) + I(MentalHealth^2) + I(SleepTime^2), 
                      data = balanced_DS, 
                      family = binomial,
                      contrasts = list("Race" = dummy_race,
                                       "AgeCategory" = dummy_age,
                                       "Diabetic" = dummy_diabetic))
  temp_avg_pvalues <- temp_avg_pvalues + summary(temp_log_mod)$coefficients[,4]
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
US_log_mod_yes_out[9] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_yes_out[9] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_yes_out[9] <- list(temp_aic)
for (i in 1:35){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}
# NB: there isn't a marked improvement in terms of AIC, therefore it's not 
# necessary to consider also the quadratic terms.

# ******************************************************************************
## LR: model selection, undersampling, without outliers ----
N_rep <- 500
US_log_mod_no_out <- NULL
US_avg_pvalues_no_out <- NULL
US_not_significant_count_no_out <- NULL
US_aic_no_out <- NULL

# 1st step: all variables in the model
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,38)
temp_not_significant_count <- rep(0,38) 
temp_aic <- 0
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS,
                             method = "under")$data
  temp_log_mod <- glm(HeartDisease ~., data = balanced_DS, family = binomial)
  temp_pvalues <- summary(temp_log_mod)$coefficients[,4]
  for (i in 1:length(temp_pvalues)){
    if (temp_pvalues[i] >= 0.01) {
      temp_not_significant_count[i] = temp_not_significant_count[i] + 1;
    }
  }
  temp_avg_pvalues <- temp_avg_pvalues + temp_pvalues
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
names(temp_not_significant_count) <- names(temp_pvalues)
US_not_significant_count_no_out[1] <- list(temp_not_significant_count)
temp_not_significant_count <- sort(temp_not_significant_count, decreasing = T)
US_log_mod_no_out[1] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_no_out[1] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_no_out[1] <- list(temp_aic)
for (i in 1:38){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}
for (i in 1:38){
  print(temp_not_significant_count[i])
}

# 2nd step: -AgeCategory (25-29)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,37)
temp_not_significant_count <- rep(0,37)
temp_aic <- 0
dummy_age <- contrasts(no_outliers_DS$AgeCategory)
dummy_age <- dummy_age[,-1]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS,
                             method = "under")$data
  temp_log_mod <- glm(HeartDisease ~., data = balanced_DS, family = binomial,
                      contrasts = list("AgeCategory" = dummy_age))
  temp_pvalues <- summary(temp_log_mod)$coefficients[,4]
  for (i in 1:length(temp_pvalues)){
    if (temp_pvalues[i] >= 0.01) {
      temp_not_significant_count[i] = temp_not_significant_count[i] + 1;
    }
  }
  temp_avg_pvalues <- temp_avg_pvalues + temp_pvalues
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
names(temp_not_significant_count) <- names(temp_pvalues)
US_not_significant_count_no_out[2] <- list(temp_not_significant_count)
temp_not_significant_count <- sort(temp_not_significant_count, decreasing = T)
US_log_mod_no_out[2] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_no_out[2] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_no_out[2] <- list(temp_aic)
for (i in 1:37){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}
for (i in 1:37){
  print(temp_not_significant_count[i])
}

# 3rd step: -PhysicalActivity
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,36)
temp_not_significant_count <- rep(0,36)
temp_aic <- 0
dummy_age <- contrasts(no_outliers_DS$AgeCategory)
dummy_age <- dummy_age[,-1]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS,
                             method = "under")$data
  temp_log_mod <- glm(HeartDisease ~.- PhysicalActivity, data = balanced_DS, family = binomial,
                      contrasts = list("AgeCategory" = dummy_age))
  temp_pvalues <- summary(temp_log_mod)$coefficients[,4]
  for (i in 1:length(temp_pvalues)){
    if (temp_pvalues[i] >= 0.01) {
      temp_not_significant_count[i] = temp_not_significant_count[i] + 1;
    }
  }
  temp_avg_pvalues <- temp_avg_pvalues + temp_pvalues
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
names(temp_not_significant_count) <- names(temp_pvalues)
US_not_significant_count_no_out[3] <- list(temp_not_significant_count)
temp_not_significant_count <- sort(temp_not_significant_count, decreasing = T)
US_log_mod_no_out[3] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_no_out[3] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_no_out[3] <- list(temp_aic)
for (i in 1:36){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}
for (i in 1:36){
  print(temp_not_significant_count[i])
}

# 4th step: -Diabetic (Yes, during pregnancy)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,35)
temp_not_significant_count <- rep(0,35)
temp_aic <- 0
dummy_age <- contrasts(no_outliers_DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(no_outliers_DS$Diabetic) 
dummy_diabetic <- dummy_diabetic[,-3]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS,
                             method = "under")$data
  temp_log_mod <- glm(HeartDisease ~.- PhysicalActivity, data = balanced_DS, family = binomial,
                      contrasts = list("AgeCategory" = dummy_age,
                                       "Diabetic" = dummy_diabetic))
  temp_pvalues <- summary(temp_log_mod)$coefficients[,4]
  for (i in 1:length(temp_pvalues)){
    if (temp_pvalues[i] >= 0.01) {
      temp_not_significant_count[i] = temp_not_significant_count[i] + 1;
    }
  }
  temp_avg_pvalues <- temp_avg_pvalues + temp_pvalues
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
names(temp_not_significant_count) <- names(temp_pvalues)
US_not_significant_count_no_out[4] <- list(temp_not_significant_count)
temp_not_significant_count <- sort(temp_not_significant_count, decreasing = T)
US_log_mod_no_out[4] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_no_out[4] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_no_out[4] <- list(temp_aic)
for (i in 1:35){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}
for (i in 1:35){
  print(temp_not_significant_count[i])
}

# 5th step: -Race (White)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,34)
temp_not_significant_count <- rep(0,34)
temp_aic <- 0
dummy_age <- contrasts(no_outliers_DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(no_outliers_DS$Diabetic) 
dummy_diabetic <- dummy_diabetic[,-3]
dummy_race <- contrasts(no_outliers_DS$Race)
dummy_race <- dummy_race[,-5]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS,
                             method = "under")$data
  temp_log_mod <- glm(HeartDisease ~.- PhysicalActivity, data = balanced_DS, family = binomial,
                      contrasts = list("AgeCategory" = dummy_age,
                                       "Diabetic" = dummy_diabetic,
                                       "Race" = dummy_race))
  temp_pvalues <- summary(temp_log_mod)$coefficients[,4]
  for (i in 1:length(temp_pvalues)){
    if (temp_pvalues[i] >= 0.01) {
      temp_not_significant_count[i] = temp_not_significant_count[i] + 1;
    }
  }
  temp_avg_pvalues <- temp_avg_pvalues + temp_pvalues
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
names(temp_not_significant_count) <- names(temp_pvalues)
US_not_significant_count_no_out[5] <- list(temp_not_significant_count)
temp_not_significant_count <- sort(temp_not_significant_count, decreasing = T)
US_log_mod_no_out[5] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_no_out[5] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_no_out[5] <- list(temp_aic)
for (i in 1:34){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}
for (i in 1:34){
  print(temp_not_significant_count[i])
}

# 6th step: -Race (Other)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,33)
temp_not_significant_count <- rep(0,33)
temp_aic <- 0
dummy_age <- contrasts(no_outliers_DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(no_outliers_DS$Diabetic) 
dummy_diabetic <- dummy_diabetic[,-3]
dummy_race <- contrasts(no_outliers_DS$Race)
dummy_race <- dummy_race[,-5]
dummy_race <- dummy_race[,-4]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS,
                             method = "under")$data
  temp_log_mod <- glm(HeartDisease ~.- PhysicalActivity, data = balanced_DS, family = binomial,
                      contrasts = list("AgeCategory" = dummy_age,
                                       "Diabetic" = dummy_diabetic,
                                       "Race" = dummy_race))
  temp_pvalues <- summary(temp_log_mod)$coefficients[,4]
  for (i in 1:length(temp_pvalues)){
    if (temp_pvalues[i] >= 0.01) {
      temp_not_significant_count[i] = temp_not_significant_count[i] + 1;
    }
  }
  temp_avg_pvalues <- temp_avg_pvalues + temp_pvalues
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
names(temp_not_significant_count) <- names(temp_pvalues)
US_not_significant_count_no_out[6] <- list(temp_not_significant_count)
temp_not_significant_count <- sort(temp_not_significant_count, decreasing = T)
US_log_mod_no_out[6] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_no_out[6] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_no_out[6] <- list(temp_aic)
for (i in 1:33){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}
for (i in 1:33){
  print(temp_not_significant_count[i])
}

# 7th step: -Race (Hispanic)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,32)
temp_not_significant_count <- rep(0,32)
temp_aic <- 0
dummy_age <- contrasts(no_outliers_DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(no_outliers_DS$Diabetic) 
dummy_diabetic <- dummy_diabetic[,-3]
dummy_race <- contrasts(no_outliers_DS$Race)
dummy_race <- dummy_race[,-5]
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-3]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS,
                             method = "under")$data
  temp_log_mod <- glm(HeartDisease ~.- PhysicalActivity, data = balanced_DS, family = binomial,
                      contrasts = list("AgeCategory" = dummy_age,
                                       "Diabetic" = dummy_diabetic,
                                       "Race" = dummy_race))
  temp_pvalues <- summary(temp_log_mod)$coefficients[,4]
  for (i in 1:length(temp_pvalues)){
    if (temp_pvalues[i] >= 0.01) {
      temp_not_significant_count[i] = temp_not_significant_count[i] + 1;
    }
  }
  temp_avg_pvalues <- temp_avg_pvalues + temp_pvalues
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
names(temp_not_significant_count) <- names(temp_pvalues)
US_not_significant_count_no_out[7] <- list(temp_not_significant_count)
temp_not_significant_count <- sort(temp_not_significant_count, decreasing = T)
US_log_mod_no_out[7] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_no_out[7] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_no_out[7] <- list(temp_aic)
for (i in 1:32){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}
for (i in 1:32){
  print(temp_not_significant_count[i])
}

# 8th step: -Diabetic (No, borderline diabetes)
temp_log_mod <- NULL
temp_avg_pvalues <- rep(0,31)
temp_not_significant_count <- rep(0,31)
temp_aic <- 0
dummy_age <- contrasts(no_outliers_DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(DS$Diabetic)
dummy_diabetic <- dummy_diabetic[,-3]
dummy_diabetic <- dummy_diabetic[,-1]
dummy_diabetic <- t(dummy_diabetic)
dummy_diabetic <- t(dummy_diabetic)
colnames(dummy_diabetic) <- "Yes"
dummy_race <- contrasts(no_outliers_DS$Race)
dummy_race <- dummy_race[,-5]
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-3]
for (i in 1:N_rep){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS,
                             method = "under")$data
  temp_log_mod <- glm(HeartDisease ~.- PhysicalActivity, data = balanced_DS, family = binomial,
                      contrasts = list("AgeCategory" = dummy_age,
                                       "Diabetic" = dummy_diabetic,
                                       "Race" = dummy_race))
  temp_pvalues <- summary(temp_log_mod)$coefficients[,4]
  for (i in 1:length(temp_pvalues)){
    if (temp_pvalues[i] >= 0.01) {
      temp_not_significant_count[i] = temp_not_significant_count[i] + 1;
    }
  }
  temp_avg_pvalues <- temp_avg_pvalues + temp_pvalues
  temp_aic <- temp_aic + summary(temp_log_mod)$aic
}
names(temp_not_significant_count) <- names(temp_pvalues)
US_not_significant_count_no_out[8] <- list(temp_not_significant_count)
temp_not_significant_count <- sort(temp_not_significant_count, decreasing = T)
US_log_mod_no_out[8] <- list(temp_log_mod)
temp_avg_pvalues <- temp_avg_pvalues/N_rep
US_avg_pvalues_no_out[8] <- list(temp_avg_pvalues)
temp_aic <- round(temp_aic/N_rep)
US_aic_no_out[8] <- list(temp_aic)
for (i in 1:31){
  if(temp_avg_pvalues[i] > 0.01){
    print(temp_avg_pvalues[i])
  }
}
for (i in 1:31){
  print(temp_not_significant_count[i])
}

# ******************************************************************************
## LR: validation 70/30, search for the best threshold, with outliers ----
thresholds <- seq(from = .35, to = .65, by = .005)
N_iter <- 500

MIS_rate_train_yes_out <- matrix(nrow = N_iter, ncol = length(thresholds))
colnames(MIS_rate_train_yes_out) <- thresholds
FP_rate_train_yes_out <- matrix(nrow = N_iter, ncol = length(thresholds))
colnames(FP_rate_train_yes_out) <- thresholds
FN_rate_train_yes_out <- matrix(nrow = N_iter, ncol = length(thresholds))
colnames(FN_rate_train_yes_out) <- thresholds
MIS_rate_test_yes_out <- matrix(nrow = N_iter, ncol = length(thresholds))
colnames(MIS_rate_test_yes_out) <- thresholds
FP_rate_test_yes_out <- matrix(nrow = N_iter, ncol = length(thresholds))
colnames(FP_rate_test_yes_out) <- thresholds
FN_rate_test_yes_out <- matrix(nrow = N_iter, ncol = length(thresholds))
colnames(FN_rate_test_yes_out) <- thresholds

dummy_race <- contrasts(DS$Race)
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-3]
dummy_age <- contrasts(DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(DS$Diabetic)
dummy_diabetic <- dummy_diabetic[,-3]
dummy_diabetic <- dummy_diabetic[,-1]
dummy_diabetic <- t(dummy_diabetic)
dummy_diabetic <- t(dummy_diabetic)
colnames(dummy_diabetic) <- "Yes"

for (i in 1:N_iter){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = no_outliers_DS, method = "under")$data
  train <- sample(1:nrow(balanced_DS), floor(nrow(balanced_DS)*.70))
  best_US_log_mod_yes_out <- glm(HeartDisease ~. - PhysicalActivity,
                                 data = balanced_DS, 
                                 family = binomial,
                                 contrasts = list("Race" = dummy_race,
                                                  "AgeCategory" = dummy_age,
                                                  "Diabetic" = dummy_diabetic),
                                 subset = train)
  pred_p_train <- predict(best_US_log_mod_yes_out, type = "response")
  for (j in 1:length(thresholds)){
    clas_train <- rep("No", length(pred_p_train))
    clas_train[pred_p_train > thresholds[j]] = "Yes";
    MIS_train_table <- table(clas_train, balanced_DS[train,1])
    MIS_rate_train_yes_out[i,j] <- round((MIS_train_table[1,2] + 
                                      MIS_train_table[2,1])/sum(MIS_train_table)*100,2)
    FP_rate_train_yes_out[i,j] <- round(MIS_train_table[2,1]/sum(MIS_train_table[,1])*100,2)
    FN_rate_train_yes_out[i,j] <- round(MIS_train_table[1,2]/sum(MIS_train_table[,2])*100,2)
  }
  pred_p_test <- predict(best_US_log_mod_yes_out, 
                         newdata = balanced_DS[-train,], type = "response")
  for (j in 1:length(thresholds)){
    clas_test <- rep("No", length(pred_p_test))
    clas_test[pred_p_test > thresholds[j]] = "Yes";
    MIS_test_table <- table(clas_test, balanced_DS[-train,1])
    MIS_rate_test_yes_out[i,j] <- round((MIS_test_table[1,2] + 
                                            MIS_test_table[2,1])/sum(MIS_test_table)*100,2)
    FP_rate_test_yes_out[i,j] <- round(MIS_test_table[2,1]/sum(MIS_test_table[,1])*100,2)
    FN_rate_test_yes_out[i,j] <- round(MIS_test_table[1,2]/sum(MIS_test_table[,2])*100,2)
  }
}

# Results
results_MIS_rate_train_yes_out <- create_results_table(MIS_rate_train_yes_out, 
                                                  thresholds, flag = TRUE)
results_FP_rate_train_yes_out <- create_results_table(FP_rate_train_yes_out, 
                                                  thresholds, flag = TRUE)
results_FN_rate_train_yes_out <- create_results_table(FN_rate_train_yes_out, 
                                                      thresholds, flag = TRUE)
results_MIS_rate_test_yes_out <- create_results_table(MIS_rate_test_yes_out, 
                                                       thresholds, flag = TRUE)
results_FP_rate_test_yes_out <- create_results_table(FP_rate_test_yes_out, 
                                                      thresholds, flag = TRUE)
results_FN_rate_test_yes_out <- create_results_table(FN_rate_test_yes_out, 
                                                      thresholds, flag = TRUE)
# Graphs
# Train misclassification error distribution, th = 50%
h <- hist(MIS_rate_train_yes_out[,31], breaks=10, col="lightblue", 
          main = "Train misclassification error distribution",
          include.lowest = TRUE,
          xlab = "MISE (%)",
          ylim = c(0,250),
          cex.lab = 1.5,
          cex.axis = 1.4,
          cex.main = 1.5)
xfit <- seq(min(MIS_rate_train_yes_out[,31]),max(MIS_rate_train_yes_out[,31]),length=40)
yfit <- dnorm(xfit,mean=mean(MIS_rate_train_yes_out[,31]),sd=sd(MIS_rate_train_yes_out[,41]))
yfit <- yfit*diff(h$mids[1:2])*length(MIS_rate_train_yes_out[,31])
lines(xfit, yfit, col="green", lwd=2)
legend("topright", legend = c("Histogram","Normal"),
       fill = c("lightblue","green"),
       cex = 1.5)

# Test misclassification error distribution, th = 50%
h <- hist(MIS_rate_test_yes_out[,31], breaks=10, col="lightblue", 
          main = "Test misclassification error distribution",
          include.lowest = TRUE,
          xlab = "MISE (%)",
          cex.lab = 1.5,
          cex.axis = 1.4,
          cex.main = 1.5)
xfit <- seq(min(MIS_rate_test_yes_out[,31]),max(MIS_rate_test_yes_out[,31]),length=40)
yfit <- dnorm(xfit,mean=mean(MIS_rate_test_yes_out[,31]),sd=sd(MIS_rate_test_yes_out[,41]))
yfit <- yfit*diff(h$mids[1:2])*length(MIS_rate_test_yes_out[,31])
lines(xfit, yfit, col="green", lwd=2)
legend("topright", legend = c("Histogram","Normal"),
       fill = c("lightblue","green"),
       cex = 1.5)

# Train false-negative error distribution, th = 50%
h <- hist(FN_rate_train_yes_out[,31], breaks=10, col="lightblue", 
          main = "Train false-negative error rate distribution",
          include.lowest = TRUE,
          xlab = "FNE (%)",
          ylim = c(0,300),
          cex.lab = 1.5,
          cex.axis = 1.4,
          cex.main = 1.5)
xfit <- seq(min(FN_rate_train_yes_out[,31]),max(FN_rate_train_yes_out[,31]),length=40)
yfit <- dnorm(xfit,mean=mean(FN_rate_train_yes_out[,31]),sd=sd(FN_rate_train_yes_out[,31]))
yfit <- yfit*diff(h$mids[1:2])*length(FN_rate_train_yes_out[,31])
lines(xfit, yfit, col="green", lwd=2)
legend("topright", legend = c("Histogram","Normal"),
       fill = c("lightblue","green"),
       cex = 1.5)

# Test false-negative error distribution, th = 50%
h <- hist(FN_rate_test_yes_out[,31], breaks=10, col="lightblue", 
          main = "Test false-negative error rate distribution",
          include.lowest = TRUE,
          xlab = "FNE (%)",
          ylim = c(0,400),
          cex.lab = 1.5,
          cex.axis = 1.4,
          cex.main = 1.5)
xfit <- seq(min(FN_rate_test_yes_out[,31]),max(FN_rate_test_yes_out[,31]),length=40)
yfit <- dnorm(xfit,mean=mean(FN_rate_test_yes_out[,31]),sd=sd(FN_rate_test_yes_out[,31]))
yfit <- yfit*diff(h$mids[1:2])*length(FN_rate_test_yes_out[,31])
lines(xfit, yfit, col="green", lwd=2)
legend("topright", legend = c("Histogram","Normal"),
       fill = c("lightblue","green"),
       cex = 1.5)

# Train false-positive error distribution, th = 50%
h <- hist(FP_rate_train_yes_out[,31], breaks=10, col="lightblue", 
          main = "Train false-positive error rate distribution",
          include.lowest = TRUE,
          xlab = "FPE (%)",
          ylim = c(0,300),
          cex.lab = 1.5,
          cex.axis = 1.4,
          cex.main = 1.5)
xfit <- seq(min(FP_rate_train_yes_out[,31]),max(FP_rate_train_yes_out[,31]),length=40)
yfit <- dnorm(xfit,mean=mean(FP_rate_train_yes_out[,31]),sd=sd(FP_rate_train_yes_out[,31]))
yfit <- yfit*diff(h$mids[1:2])*length(FP_rate_train_yes_out[,31])
lines(xfit, yfit, col="green", lwd=2)
legend("topright", legend = c("Histogram","Normal"),
       fill = c("lightblue","green"),
       cex = 1.5)

# Test false-positive error distribution, th = 50%
h <- hist(FP_rate_test_yes_out[,41], breaks=10, col="lightblue", 
          main = "Test false-positive error rate distribution",
          include.lowest = TRUE,
          xlab = "FPE (%)",
          ylim = c(0,400),
          cex.lab = 1.5,
          cex.axis = 1.4,
          cex.main = 1.5)
xfit <- seq(min(FP_rate_test_yes_out[,31]),max(FP_rate_test_yes_out[,31]),length=40)
yfit <- dnorm(xfit,mean=mean(FP_rate_test_yes_out[,31]),sd=sd(FP_rate_test_yes_out[,31]))
yfit <- yfit*diff(h$mids[1:2])*length(FP_rate_test_yes_out[,31])
lines(xfit, yfit, col="green", lwd=2)
legend("topright", legend = c("Histogram","Normal"),
       fill = c("lightblue","green"),
       cex = 1.5)

# Threshold's analysis
plot(thresholds*100, results_MIS_rate_train_yes_out[,2], 
     ylim = c(22,30),
     col = "red",
     main = "Train MISE rate and classification threshold",
     ylab = "MISE (%)",
     xlab = "Threshold (%)",
     cex.lab = 1.5,
     cex.axis = 1.4,
     cex.main = 1.5)
lines(thresholds*100, results_MIS_rate_train_yes_out[,1], col = "blue")
lines(thresholds*100, results_MIS_rate_train_yes_out[,3], col = "green")
legend("bottomright", legend = c("Max","Avg", "Min"),
       fill = c("green","red", "blue"),
       cex = 1.5)

plot(thresholds*100, results_MIS_rate_test_yes_out[,2], 
     ylim = c(22,30),
     col = "blue",
     main = "Test MISE and classification threshold",
     ylab = "MISE (%)",
     xlab = "Threshold (%)",
     cex.lab = 1.5,
     cex.axis = 1.4,
     cex.main = 1.5)
lines(thresholds*100, results_MIS_rate_test_yes_out[,1], col = "gold")
lines(thresholds*100, results_MIS_rate_test_yes_out[,3], col = "red")
abline(v=47, col = "green")
legend("topright", legend = c("Max","Avg", "Min"),
       fill = c("red","blue", "gold"),
       cex = 1.5)

plot(thresholds*100, results_FP_rate_test_yes_out[,2],
     col = "blue3",
     ylim = c(10, 40),
     main = "Test FPE - FNE and classification threshold",
     ylab = "FPE - FNE (%)",
     xlab = "Threshold (%)",
     cex.lab = 1.5,
     cex.axis = 1.4,
     cex.main = 1.5,)
lines(thresholds*100, results_FN_rate_test_yes_out[,2], col = "cornflowerblue", 
      type = "p")
abline(v=47, col = "green")
legend("left", legend = c("FPE","FNE"),
       fill = c("blue3","cornflowerblue"),
       cex = 1.5)

# PDF files' dimensions:
# - 233.2 mm --> 9.18 inch
# - 145 mm --> 5.71 inch

# ******************************************************************************
## LR: validation 70/30, search for the best alpha, without outliers ----
alphas <- seq(from = 1, to = 15, by = 1)
N_iter <- 1000

MIS_rate_train_no_out <- matrix(nrow = N_iter, ncol = length(alphas))
colnames(MIS_rate_train_no_out) <- alphas
FP_rate_train_no_out <- matrix(nrow = N_iter, ncol = length(alphas))
colnames(FP_rate_train_no_out) <- alphas
FN_rate_train_no_out <- matrix(nrow = N_iter, ncol = length(alphas))
colnames(FN_rate_train_no_out) <- alphas
MIS_rate_test_no_out <- matrix(nrow = N_iter, ncol = length(alphas))
colnames(MIS_rate_test_no_out) <- alphas
FP_rate_test_no_out <- matrix(nrow = N_iter, ncol = length(alphas))
colnames(FP_rate_test_no_out) <- alphas
FN_rate_test_no_out <- matrix(nrow = N_iter, ncol = length(alphas))
colnames(FN_rate_test_no_out) <- alphas

dummy_race <- contrasts(DS$Race)
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-4]
dummy_race <- dummy_race[,-3]
dummy_age <- contrasts(DS$AgeCategory)
dummy_age <- dummy_age[,-1]
dummy_diabetic <- contrasts(DS$Diabetic)
dummy_diabetic <- dummy_diabetic[,-3]
dummy_diabetic <- dummy_diabetic[,-1]
dummy_diabetic <- t(dummy_diabetic)
dummy_diabetic <- t(dummy_diabetic)
colnames(dummy_diabetic) <- "Yes"

for (i in 1:N_iter){
  print("New iteration")
  print(i)
  for (j in 1:length(alphas)){
    new_DS <- remove_outliers(DS, alphas[j])
    balanced_DS <- ovun.sample(HeartDisease ~., data = new_DS, method = "under")$data
    table(balanced_DS$HeartDisease)
    train <- sample(1:nrow(balanced_DS), floor(nrow(balanced_DS)*.70))
    best_US_log_mod_no_out <- glm(HeartDisease ~. - PhysicalActivity,
                                   data = balanced_DS, 
                                   family = binomial,
                                   contrasts = list("Race" = dummy_race,
                                                    "AgeCategory" = dummy_age,
                                                    "Diabetic" = dummy_diabetic),
                                   subset = train)
    print(alphas[j])
    pred_p_train <- predict(best_US_log_mod_no_out, type = "response")
    clas_train <- rep("No", length(pred_p_train))
    clas_train[pred_p_train > 0.5] = "Yes";
    MIS_train_table <- table(clas_train, balanced_DS[train,1])
    MIS_rate_train_no_out[i,j] <- round((MIS_train_table[1,2] + 
                                            MIS_train_table[2,1])/sum(MIS_train_table)*100,2)
    FP_rate_train_no_out[i,j] <- round(MIS_train_table[2,1]/sum(MIS_train_table[,1])*100,2)
    FN_rate_train_no_out[i,j] <- round(MIS_train_table[1,2]/sum(MIS_train_table[,2])*100,2)
    
    pred_p_test <- predict(best_US_log_mod_no_out, 
                           newdata = balanced_DS[-train,], type = "response")
    clas_test <- rep("No", length(pred_p_test))
    clas_test[pred_p_test > 0.5] = "Yes";
    MIS_test_table <- table(clas_test, balanced_DS[-train,1])
    MIS_rate_test_no_out[i,j] <- round((MIS_test_table[1,2] + 
                                           MIS_test_table[2,1])/sum(MIS_test_table)*100,2)
    FP_rate_test_no_out[i,j] <- round(MIS_test_table[2,1]/sum(MIS_test_table[,1])*100,2)
    FN_rate_test_no_out[i,j] <- round(MIS_test_table[1,2]/sum(MIS_test_table[,2])*100,2)
  }
}

# Results
results_MIS_rate_train_no_out <- create_results_table(MIS_rate_train_no_out, 
                                                       alphas, flag = TRUE)
results_FP_rate_train_no_out <- create_results_table(FP_rate_train_no_out, 
                                                      alphas, flag = TRUE)
results_FN_rate_train_no_out <- create_results_table(FN_rate_train_no_out, 
                                                      alphas, flag = TRUE)
results_MIS_rate_test_no_out <- create_results_table(MIS_rate_test_no_out, 
                                                      alphas, flag = TRUE)
results_FP_rate_test_no_out <- create_results_table(FP_rate_test_no_out, 
                                                     alphas, flag = TRUE)
results_FN_rate_test_no_out <- create_results_table(FN_rate_test_no_out, 
                                                     alphas, flag = TRUE)

# ******************************************************************************
### Random forest (1st validation approach) ----
random_forest_list <- NULL
num_tree <- 200;
num_iter <- 75;
for (i in 1:num_iter){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  rf_model = randomForest(HeartDisease ~. ,
                          data = balanced_DS,
                          mtry = floor(sqrt(ncol(balanced_DS)-1)),
                          importance = TRUE,
                          ntree = num_tree,
                          replace = TRUE)
  random_forest_list[i] <- list(rf_model)
}

clas_matrix <- data.frame(row.names = 1:nrow(DS))
for (i in 1:num_iter){
  print(i)
  clas_matrix[as.character(i)] <- predict(random_forest_list[[i]], 
                             newdata = DS, 
                             type = "class")
}
clas <- rep(0,nrow(DS))
for (i in 1:nrow(DS)){
  print(i)
  count_yes <- 0
  count_no <- 0
  for (j in 1:num_iter){
    if(clas_matrix[i,j] == "Yes")
      count_yes <- count_yes + 1
    else
      count_no <- count_no + 1
  }
  if (count_yes > count_no)
    clas[i] <- "Yes"
  else
    clas[i] <- "No"
}
confusion_matrix_rf <- table(clas, DS$HeartDisease)
MIS_rate_rf <- round((confusion_matrix_rf[1,2] + 
                        confusion_matrix_rf[2,1])/sum(confusion_matrix_rf)*100,2)
FP_rate_rf <- round(confusion_matrix_rf[2,1]/sum(confusion_matrix_rf[,1])*100,2)
FN_rate_rf <- round(confusion_matrix_rf[1,2]/sum(confusion_matrix_rf[,2])*100,2)
table(DS$HeartDisease)
confusion_matrix_rf

# ******************************************************************************
### Random forest (2nd validation approach) ----
N_iter <- 100
N_tree <- 300

var_names <- c("MISE_train", "MISE_test", "FPR_train", "FPR_test", 
               "FNR_train", "FNR_test")
rf_val_table <- matrix(nrow = N_iter, ncol = length(var_names))
colnames(rf_val_table) <- var_names

for (i in 1:N_iter){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  train <- sample(1:nrow(balanced_DS), floor(nrow(balanced_DS)*.70))
  rf_model = randomForest(HeartDisease ~. ,
                          data = balanced_DS,
                          subset = train,
                          mtry = floor(sqrt(ncol(balanced_DS)-1)),
                          importance = TRUE,
                          ntree = N_tree,
                          replace = TRUE)
  pred_train <- predict(rf_model, type = "class")
  train_confusion_matrix <- table(pred_train, balanced_DS[train,1])
  rf_val_table[i,1] <- round((train_confusion_matrix[1,2] + 
                                train_confusion_matrix[2,1])/sum(train_confusion_matrix)*100,2)
  rf_val_table[i,3] <- round(train_confusion_matrix[2,1]/sum(train_confusion_matrix[,1])*100,2)
  rf_val_table[i,5] <- round(train_confusion_matrix[1,2]/sum(train_confusion_matrix[,2])*100,2)
  
  pred_test <- predict(rf_model, newdata = balanced_DS[-train,], type = "class")
  test_confusion_matrix <- table(pred_test, balanced_DS[-train,1])
  rf_val_table[i,2] <- round((test_confusion_matrix[1,2] + 
                                test_confusion_matrix[2,1])/sum(test_confusion_matrix)*100,2)
  rf_val_table[i,4] <- round(test_confusion_matrix[2,1]/sum(test_confusion_matrix[,1])*100,2)
  rf_val_table[i,6] <- round(test_confusion_matrix[1,2]/sum(test_confusion_matrix[,2])*100,2)
}

# Results
results_rf_validation <- create_results_table(rf_val_table, 
                                              var_names, flag = TRUE)

# ******************************************************************************
### QDA ----
N_iter <- 1000

var_names <- c("MISE_train", "MISE_test", "FPR_train", "FPR_test", 
               "FNR_train", "FNR_test")
qda_val_table <- matrix(nrow = N_iter, ncol = length(var_names))
colnames(qda_val_table) <- var_names

for (i in 1:N_iter){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  train <- sample(1:nrow(balanced_DS), floor(nrow(balanced_DS)*.70))
  qda_model = qda(HeartDisease ~. ,
                 data = balanced_DS,
                 subset = train)
  pred_train <- predict(qda_model, type = "response")
  train_confusion_matrix <- table(pred_train$class, balanced_DS[train,1])
  qda_val_table[i,1] <- round((train_confusion_matrix[1,2] + 
                                train_confusion_matrix[2,1])/sum(train_confusion_matrix)*100,2)
  qda_val_table[i,3] <- round(train_confusion_matrix[2,1]/sum(train_confusion_matrix[,1])*100,2)
  qda_val_table[i,5] <- round(train_confusion_matrix[1,2]/sum(train_confusion_matrix[,2])*100,2)
  
  pred_test <- predict(qda_model, newdata = balanced_DS[-train,], type = "response")
  test_confusion_matrix <- table(pred_test$class, balanced_DS[-train,1])
  qda_val_table[i,2] <- round((test_confusion_matrix[1,2] + 
                                test_confusion_matrix[2,1])/sum(test_confusion_matrix)*100,2)
  qda_val_table[i,4] <- round(test_confusion_matrix[2,1]/sum(test_confusion_matrix[,1])*100,2)
  qda_val_table[i,6] <- round(test_confusion_matrix[1,2]/sum(test_confusion_matrix[,2])*100,2)
}

# Results
results_qda_validation <- create_results_table(qda_val_table, 
                                              var_names, flag = TRUE)

# ******************************************************************************
### FDA ----
N_iter <- 1000

var_names <- c("MISE_train", "MISE_test", "FPR_train", "FPR_test", 
               "FNR_train", "FNR_test")
lda_val_table <- matrix(nrow = N_iter, ncol = length(var_names))
colnames(lda_val_table) <- var_names

for (i in 1:N_iter){
  print(i)
  balanced_DS <- ovun.sample(HeartDisease ~., data = DS, method = "under")$data
  train <- sample(1:nrow(balanced_DS), floor(nrow(balanced_DS)*.70))
  lda_model = lda(HeartDisease ~. ,
                  data = balanced_DS,
                  subset = train)
  pred_train <- predict(lda_model, type = "response")
  train_confusion_matrix <- table(pred_train$class, balanced_DS[train,1])
  lda_val_table[i,1] <- round((train_confusion_matrix[1,2] + 
                                 train_confusion_matrix[2,1])/sum(train_confusion_matrix)*100,2)
  lda_val_table[i,3] <- round(train_confusion_matrix[2,1]/sum(train_confusion_matrix[,1])*100,2)
  lda_val_table[i,5] <- round(train_confusion_matrix[1,2]/sum(train_confusion_matrix[,2])*100,2)
  
  pred_test <- predict(lda_model, newdata = balanced_DS[-train,], type = "response")
  test_confusion_matrix <- table(pred_test$class, balanced_DS[-train,1])
  lda_val_table[i,2] <- round((test_confusion_matrix[1,2] + 
                                 test_confusion_matrix[2,1])/sum(test_confusion_matrix)*100,2)
  lda_val_table[i,4] <- round(test_confusion_matrix[2,1]/sum(test_confusion_matrix[,1])*100,2)
  lda_val_table[i,6] <- round(test_confusion_matrix[1,2]/sum(test_confusion_matrix[,2])*100,2)
}

# Results
results_lda_validation <- create_results_table(lda_val_table, 
                                               var_names, flag = TRUE)
