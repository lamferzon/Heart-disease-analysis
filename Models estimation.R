
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("caret")
#install.packages("randomForest")
library(MASS)
library(ggplot2)
library(tseries)
library(e1071)
library(dplyr)
library(tree)
library(caret)
library(randomForest)
library(boot)

rm(list = ls())

#Importation data from "Hearth desease dataset"
DS <- read.csv("C:/Users/loren/Desktop/Statistical Learning/Scientific paper/Heart desease dataset.csv",
              header=TRUE, sep=",", dec=".")
DS <- na.omit(DS)
N <- nrow(DS)
cat_var = colnames(DS)
cat_var = cat_var[-c(2,6,7,15)]
for (i in cat_var){
  DS[,i] <- factor(DS[,i])
}

## STATISTICS CALCULATION
#Binary variables statistics calculation
bin_var <- colnames(DS[,c(1,3,4,5,8,13,16,17,18)])
bin_var_stat <- data.frame(Yes_per=1:9, No_per=1:9, row.names = bin_var)
for(i in 1:9){
  temp <- DS %>% group_by(eval(parse(text = bin_var[i]))) %>% summarise(count = n())
  bin_var_stat[i,1] <- round((temp[2,2]/N)*100,2)
  bin_var_stat[i,2] <- round((temp[1,2]/N)*100,2)
}

#Numeric variables statistics calculation
num_var <- colnames(DS[,c(2,6,7,15)])
num_var_stat <- data.frame(Mean=1:4, Median=1:4, Var=1:4, Std_dev=1:4,
                           Kurtosis=1:4, Skewness=1:4,
                           Pvalue_JB=1:4, H0_rejection=1:4,
                           row.names = num_var)
for (i in 1:4){
  num_var_stat[i,1] <- round(mean(DS[,num_var[i]]),2)
  num_var_stat[i,2] <- round(median(DS[,num_var[i]]),2)
  num_var_stat[i,3] <- round(var(DS[,num_var[i]]),2)
  num_var_stat[i,4] <- round(sd(DS[,num_var[i]]),2)
  num_var_stat[i,5] <- round(kurtosis(DS[,num_var[i]]),2)
  num_var_stat[i,6] <- round(skewness(DS[,num_var[i]]),2)
  num_var_stat[i,7] <- jarque.bera.test(DS[,num_var[i]])$p.value
  if (jarque.bera.test(DS[,num_var[i]])$p.value > 0.05){
    num_var_stat[i,8] <- 0
  }else{
    num_var_stat[i,8] <- 1
  }
}

#Diabetic variable statistics calculation
diabetic_stat <- DS %>% group_by(Diabetic) %>% summarise(count = n())
colnames(diabetic_stat)[2] <- c("Percentage")
diabetic_stat[,2] <- round((diabetic_stat[,2]/N)*100,2)

#Genre variable statistics calculation
genre_stat <- DS %>% group_by(Sex) %>% summarise(count = n())
colnames(genre_stat)[2] <- c("Percentage")
genre_stat[,2] <- round((genre_stat[,2]/N)*100,2)

#Race variable statistics calculation
race_stat <- DS %>% group_by(Race) %>% summarise(count = n())
colnames(race_stat)[2] <- c("Percentage")
race_stat[,2] <- round((race_stat[,2]/N)*100,2)

#General health statistics calculation
gen_health_stat <- DS %>% group_by(GenHealth) %>% summarise(count = n())
colnames(gen_health_stat)[2] <- c("Percentage")
gen_health_stat[,2] <- round((gen_health_stat[,2]/N)*100,2)

#Age statistics calculation
age_stat <- DS %>% group_by(AgeCategory) %>% summarise(count = n())
colnames(age_stat)[2] <- c("Percentage")
age_stat[,2] <- round((age_stat[,2]/N)*100,2)
rm(temp)

SIGMA = round(cov(DS[,num_var]),2)
R = matrix(nrow = 4, ncol = 4)
for (i in 1:4){
  for (j in 1:4){
    R[i,j] = round(SIGMA[i,j]/sqrt(SIGMA[i,i]*SIGMA[j,j]),2)
  }
}

## BARPLOTS
# Binary variables' barplots
bin_labels = c("Heart disease", "Smoking", "Alcohol drinking", "Stroke", 
          "Difficulty walking", "Physical activity", "Asthma",
          "Kidney disease", "Skin cancer")
for (i in 1:9){
  mybar <- barplot(table(DS[,bin_var[i]]),
                  width = 2,
                  main = bin_labels[i],
                  xlab = "Number of people",
                  xlim = c(0,350000),
                  names.arg = c("No", "Yes"),
                  col=topo.colors(2),
                  horiz = TRUE,
                  cex.lab = 1.5,
                  cex.axis = 1.4,
                  cex.names = 1.5,
                  cex.main = 1.5)
  if(i == 6){
    pos <- "right"
  }else{
    pos <- "topright"
  }
  legend(pos, paste(c(bin_var_stat[i,1],bin_var_stat[i,2]), "%"),
         fill = rev(topo.colors(2)), title = "Percentage", cex = 1.5)
  text(table(DS[,bin_var[i]])+20000, mybar, labels = table(DS[,bin_var[i]]), cex = 1.5)
}

# Diabetic variable's barplot
mybar <- barplot(table(DS[,"Diabetic"]),
                width = 2,
                main = "Diabetic",
                xlab = "Number of people",
                xlim = c(0,350000),
                names.arg = c(" ", " ", " ", " "),
                col=topo.colors(4),
                horiz = TRUE,
                cex.lab = 1.5,
                cex.axis = 1.4,
                cex.names = 1.5,
                cex.main = 1.5)
legend("topright", paste(c(diabetic_stat[4,2], diabetic_stat[3,2], 
                           diabetic_stat[2,2], diabetic_stat[1,2]), "%"),
       fill = rev(topo.colors(4)), title = "Percentage", cex = 1.5)
text(table(DS[,"Diabetic"])+c(35000, 80000, 35000, 80000), mybar, 
     labels = paste( c("No", "No, bordeline diabetes", "Yes", 
                       "Yes, during pregnancy"), "-", table(DS[,"Diabetic"])), 
     cex = 1.5)

# Genre variable's statistic
mybar <- barplot(table(DS[,"Sex"]),
                width = 2,
                main = "Genre",
                xlab = "Number of people",
                xlim = c(0,350000),
                names.arg = c("Female", "Male"),
                col=topo.colors(2),
                horiz = TRUE,
                cex.lab = 1.5,
                cex.axis = 1.4,
                cex.names = 1.5,
                cex.main = 1.5)
legend("topright", paste(c(genre_stat[2,2], genre_stat[1,2]), "%"),
       fill = rev(topo.colors(2)), title = "Percentage", cex = 1.5)
text(table(DS[,"Sex"])+20000, mybar, labels = table(DS[,"Sex"]), cex = 1.5)

# Race variable's barplot
mybar <- barplot(table(DS[,"Race"]),
                width = 2,
                main = "Race",
                xlab = "Number of people",
                xlim = c(0,350000),
                names.arg = c(" ", " ", " ", " ", " ", " "),
                col=topo.colors(6),
                horiz = TRUE,
                cex.lab = 1.5,
                cex.axis = 1.4,
                cex.names = 1.5,
                cex.main = 1.5)
legend("right", paste(c(race_stat[6,2], race_stat[5,2], 
                           race_stat[4,2], race_stat[3,2],
                           race_stat[2,2], race_stat[1,2]), "%"),
       fill = rev(topo.colors(6)), title = "Percentage", cex = 1.5)
text(table(DS[,"Race"])+c(98000, 35000, 38000, 45000, 38000, 42000), mybar, 
     labels = paste( c("American indian/Alaskan native", "Asian", "Black", "Hispanic", "Other",
                       "White"), "-", table(DS[,"Race"])), 
     cex = 1.5)

# General health's barplot
mybar <- barplot(table(DS[,"GenHealth"]),
                width = 2,
                main = "General health",
                xlab = "Number of people",
                xlim = c(0,350000),
                names.arg = c(" ", " ", " ", " ", " "),
                col=topo.colors(5),
                horiz = TRUE,
                cex.lab = 1.5,
                cex.axis = 1.4,
                cex.names = 1.5,
                cex.main = 1.5)
legend("topright", paste(c(gen_health_stat[5,2], gen_health_stat[4,2], 
                           gen_health_stat[3,2], gen_health_stat[2,2],
                           gen_health_stat[1,2]), "%"),
       fill = rev(topo.colors(5)), title = "Percentage", cex = 1.5)
text(table(DS[,"GenHealth"])+c(47000, 33000, 37000, 36000, 52000), mybar, 
     labels = paste( c("Excellent", "Fair", "Good", "Poor", "Very good"), 
                     "-", table(DS[,"GenHealth"])), 
     cex = 1.5)

a <- table(DS[,"AgeCategory"])

# Age's barplot
mybar <- barplot(table(DS[,"AgeCategory"]),
                width = 2,
                main = "Age category",
                xlab = "Number of people",
                xlim = c(0,60000),
                names.arg = c("18-24","25-29","30-34","35-39","40-44","45-49",
                              "50-54","55-59","60-64","65-69","70-74","75-79",
                              ">= 80"),
                col=topo.colors(13),
                horiz = TRUE,
                cex.lab = 1.5,
                cex.axis = 1.4,
                cex.names = 1.4,
                cex.main = 1.5,
                las = 1)
legend("topright", paste(c(age_stat[13,2], age_stat[12,2], age_stat[11,2], age_stat[10,2],
                           age_stat[9,2], age_stat[8,2], age_stat[7,2], age_stat[6,2],
                           age_stat[5,2], age_stat[4,2], age_stat[3,2],
                           age_stat[2,2], age_stat[1,2]), "%"),
       fill = rev(topo.colors(13)), title = "Percentage", cex = 1.25)
text(table(DS[,"AgeCategory"])+4000, mybar, 
     labels = paste(table(DS[,"AgeCategory"])), 
     cex = 1.4)

## HISTOGRAMS
# BMI histogram
h <- hist(DS$BMI, breaks=80, col="lightblue", 
        main = "Body Mass Index",
        xlab = "BMI [kg/m^2]",
        xlim = c(0,60),
        ylim = c(0,30000),
        include.lowest = TRUE,
        cex.lab = 1.5,
        cex.axis = 1.4,
        cex.main = 1.5)
xfit <- seq(min(DS$BMI),max(DS$BMI),length=40)
yfit <- dnorm(xfit,mean=mean(DS$BMI),sd=sd(DS$BMI))
yfit <- yfit*diff(h$mids[1:2])*length(DS$BMI)
lines(xfit, yfit, col="green", lwd=2)
legend("topright", legend = c("Histogram","Normal"),
                              fill = c("lightblue","green"),
       cex = 1.5)
summary(DS$BMI)

# Physical health histogram
h <- hist(DS$PhysicalHealth, breaks=30, col="lightblue", 
          main = "Physical health",
          xlab = "Days",
          xlim = c(0,30),
          ylim = c(0,20000),
          include.lowest = TRUE,
          cex.lab = 1.5,
          cex.axis = 1.4,
          cex.main = 1.5)
xfit <- seq(min(DS$PhysicalHealth),max(DS$PhysicalHealth),length=50)
yfit <- dnorm(xfit,mean=mean(DS$PhysicalHealth),sd=sd(DS$PhysicalHealth))
yfit <- yfit*diff(h$mids[1:2])*length(DS$PhysicalHealth)
lines(xfit, yfit, col="green", lwd=2)
legend("top", legend = c("Histogram","Normal"),
       fill = c("lightblue","green"),
       cex = 1.5)
text(3.5, 20000, 
     labels = "<-- 237078", 
     cex = 1.4)
summary(DS$PhysicalHealth)

# Mental health histogram
h <- hist(DS$MentalHealth, breaks=30, col="lightblue", 
          main = "Mental health",
          xlab = "Days",
          xlim = c(0,30),
          ylim = c(0,20000),
          include.lowest = TRUE,
          cex.lab = 1.5,
          cex.axis = 1.4,
          cex.main = 1.5)
xfit <- seq(min(DS$MentalHealth),max(DS$MentalHealth),length=50)
yfit <- dnorm(xfit,mean=mean(DS$MentalHealth),sd=sd(DS$MentalHealth))
yfit <- yfit*diff(h$mids[1:2])*length(DS$MentalHealth)
lines(xfit, yfit, col="green", lwd=2)
legend("top", legend = c("Histogram","Normal"),
       fill = c("lightblue","green"),
       cex = 1.5)
text(3.5, 20000, 
     labels = "<-- 214692", 
     cex = 1.4)
summary(DS$MentalHealth)

# Sleep time histogram
h <- hist(DS$SleepTime, breaks=24, col="lightblue", 
          main = "Sleep time",
          xlab = "Hours",
          xlim = c(0,25),
          ylim = c(0,85000),
          include.lowest = TRUE,
          cex.lab = 1.5,
          cex.axis = 1.4,
          cex.main = 1.5)
xfit <- seq(min(DS$SleepTime),max(DS$SleepTime),length=40)
yfit <- dnorm(xfit,mean=mean(DS$SleepTime),sd=sd(DS$SleepTime))
yfit <- yfit*diff(h$mids[1:2])*length(DS$SleepTime)
lines(xfit, yfit, col="green", lwd=2)
legend("topright", legend = c("Histogram","Normal"),
       fill = c("lightblue","green"),
       cex = 1.5)
text(c(4,10), c(80000,80000), 
     labels = c("97751 -->", "<-- 97602"), 
     cex = 1.4)
summary(DS$SleepTime)

## CLASSIFICATION TREE
counters <- table(DS$HeartDisease)
index_yes <- array(dim = counters[2])
index_no <- array(dim = counters[1])
idx_yes <- 1
idx_no <- 1
for (i in 1:N){
  if(DS[i,1] == "No"){
    index_no[idx_no] <- i
    idx_no <- idx_no+1
  }else{
    index_yes[idx_yes] <- i
    idx_yes <- idx_yes+1
  }
}
train1 <- sample(index_yes, floor(length(index_yes)*0.5))
train2 <- sample(index_no, floor(length(index_yes)*0.5))
train <- sort(append(train1, train2))
table(DS[train, "HeartDisease"])

random_forest_model <- randomForest(HeartDisease ~ . ,data = DS, subset = train,
                           mtry = floor(sqrt(ncol(DS)-1)), importance = TRUE, ntree = 100, replace = TRUE)
random_forest_model
plot(random_forest_model)

#test1 <- sample(setdiff(index_yes, train1), floor(length(index_no)*0.50), replace = TRUE)
#test2 <- setdiff(index_no, train2)
#test <- union(test1, test2)

yhat <- predict(random_forest_model, newdata = DS[-train,], type = "class")
table(yhat,DS$HeartDisease[-train])

## LOGISTIC REGRESSION
# Creation of a balanced training subset
counters <- table(DS$HeartDisease)
index_yes <- array(dim = counters[2])
index_no <- array(dim = counters[1])
idx_yes <- 1
idx_no <- 1
for (i in 1:N){
  if(DS[i,1] == "No"){
    index_no[idx_no] <- i
    idx_no <- idx_no+1
  }else{
    index_yes[idx_yes] <- i
    idx_yes <- idx_yes+1
  }
}
train1 <- sample(index_yes, floor(length(index_yes)*0.5))
train2 <- sample(index_no, floor(length(index_yes)*0.5))
train <- sort(append(train1, train2))

logistic_model <- glm(HeartDisease ~., 
                      data = DS, family = binomial, 
                      subset = train)
summary(logistic_model)
fit <- predict(logistic_model, type = "response")
clas <- rep( "No" , 27372)
clas[fit > .50] = "Yes";

test1 <- setdiff(index_yes, train1)
test2 <- sample(setdiff(index_no, train2), floor(length(index_yes)*0.5))
test <- sort(append(test1, test2))

yhat <- predict(logistic_model, newdata = DS[test,], type = "response")
clas2 <- rep( "No" , 27373)
clas2[yhat > .50] = "Yes";
table(clas, DS[train,1])
table(clas2, DS[test,1])

step_model <- stepAIC(logistic_model, direction = "backward", trace = FALSE)
summary(step_model)

fit <- predict(step_model, type = "response")
clas <- rep( "No" , 27372)
clas[fit > .50] = "Yes";

yhat <- predict(step_model, newdata = DS[test,], type = "response")
clas2 <- rep( "No" , 27373)
clas2[yhat > .50] = "Yes";
table(clas, DS[train,1])
table(clas2, DS[test,1])

complex_model <- glm(HeartDisease ~. -PhysicalActivity
                     +poly(BMI,2) +poly(PhysicalHealth,2)
                     +poly(SleepTime,2), 
                     data = DS, family = binomial, subset = train)
summary(complex_model)
fit <- predict(complex_model, type = "response")
clas <- rep( "No" , 41059)
clas[fit > .50] = "Yes";

yhat <- predict(complex_model, newdata = DS[test,], type = "response")
clas2 <- rep( "No" , 41060)
clas2[yhat > .50] = "Yes";
table(clas, DS[train,1])
table(clas2, DS[test,1])

# LDA
lda_fit <- lda(HeartDisease ~.,
                data = DS, subset = train); 
summary(lda_fit)
fit <- predict(lda_fit, type = "response")
yhat <- predict(lda_fit, newdata = DS[test,], type = "response")

table(fit$class, DS[train,1])
table(yhat$class, DS[test,1])


# QDA
qda_fit <- qda(HeartDisease ~.,
               data = DS, subset = train); 
summary(qda_fit)
fit <- predict(qda_fit, type = "response")
yhat <- predict(qda_fit, newdata = DS[test,], type = "response")

table(fit$class, DS[train,1])
table(yhat$class, DS[test,1])


# Naive Bayes
library(e1071)
nb_fit <- naiveBayes(HeartDisease ~.,
               data = DS, subset = train); 
summary(nb_fit)
fit <- predict(nb_fit, newdata = DS[train,])
yhat <- predict(nb_fit, newdata = DS[test,])

table(fit, DS[train,1])
table(yhat, DS[test,1]) 

#KNN
library(class)
matrix <- cbind(DS$BMI, DS$Smoking, DS$AlcoholDrinking, DS$Stroke, DS$PhysicalHealth,
                DS$MentalHealth, DS$DiffWalking, DS$Sex, DS$AgeCategory, DS$Race,
                DS$Diabetic, DS$PhysicalActivity, DS$GenHealth, DS$SleepTime,
                DS$Asthma, DS$KidneyDisease, DS$SkinCancer); 
train <- sample(1:nrow(DS), floor(nrow(DS)*0.5)) 
data_train <- matrix[train,];
data_test <- matrix[-train,];
label_train <- DS$HeartDisease[train]; 
knn_fit <- knn(data_train,data_test,label_train, k=3)

table(knn_fit, DS$HeartDisease[-train])
