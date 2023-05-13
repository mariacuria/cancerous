#Mariia KIM   r0831824
#Soniya LAMA   r0715145
#Statistical Methods for Bioinformatics

library("Hmisc")
library("pastecs")
library("doBy")
library("stats")
#library("raster")
library("XLConnect")

library(reshape2)
library(corrplot)
library(leaps)

#load packages
library(glmnet)
library(tidyverse)
library(boot)
library(corrplot)
library(onewaytests)
require(methods)
library(pROC)
library(doBy)
library(leaps)
library(splines)
library(gam)
library(janitor)
library(pls)

#Load the data
load('prostate2.Rdata')
head(prostate)

######################################

#1. Study and describe the predictor variables. Do you see any issues that are relevant for making predictions? 

#############Question 1###############
#Overview
dim(prostate)
summary(prostate)
#Descriptive statistics
descrip.prostate<-stat.desc(prostate[,names(prostate)],basic=TRUE, desc=TRUE)
descrip.prostate
#Calculating the covariance
covar.prostate <- cov(prostate)
covar.prostate
# Cscore (response variable) has the highest correlation with lpsa (0.7215889).
# The highest correlation is between lpsa and lcavol (0.7344603).

#Box plots for each variable
par(mfrow = c(1, 2))
boxplot(prostate$Cscore,main="Cscore",ylab="Cscore", col = "gray")
#boxplot(prostate$lcavol,main="lcavol",ylab="lcavol", col = "gray")
#par(mfrow = c(1, 2))
#boxplot(prostate$lweight,main="lweight",ylab="lweight", col = "gray")
boxplot(prostate$age,main="age",ylab="age", col = "gray")
par(mfrow = c(1, 2))
boxplot(prostate$lbph,main="lbph",ylab="lbph", col = "gray")
boxplot(prostate$lcp,main="lcp",ylab="lcp", col = "gray")
par(mfrow = c(1, 2))
boxplot(prostate$lpsa,main="lpsa",ylab="lpsa", col = "gray")
boxplot(prostate$svi,main="svi",ylab="svi", col = "gray")
# Calculate the correlation matrix
cormat = cor(prostate)
cormat
#correlation plot
ggplot(data = melt(cormat), aes(x=Var1, y=Var2, fill=value)) + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation")+ geom_tile()

# The boxplot for Cscore indicates a high variability of cancer progression between patients.
# There are many high outliers in Cscore (patients who have a higher progression of cancer).
# The median is positioned lower than the middle of the box, suggesting that the distribution is skewed to the left
# (negatively skewed). In other words, the data has a tail on the left side, with more values on the higher end.
# The third IQR (the difference between the third quartile and the median) is twice as big as the first IQR (the
# difference between the median and the first quartile), it indicates that the upper half of the data has more variability
# or dispersion compared to the lower half. This suggests that there is a wider range of values in the upper part of the
# distribution.

#determine what predictors have a correlation greater than 0.5
m1 <- apply(cormat>0.5,2,sum)
m1
max(m1)
length(which(m1>1))/length(m1) 

#about 63% of the variables have a correlation greater than 0.5 with another variable.

#plot the numeric predictors against one another
pairs(prostate)
#there is a pattern of linear association between lcavol and lcp

#finally let's look at a summary of the response variable
summary(prostate[,"Cscore"])

######################################

#2. Generate your best linear regression model using only linear effects. Are there any indications that
# assumptions underlying inferences with the model are violated? Evaluate the effect of any influential
# point, or outlier.

#############Question 2###############

#Best Subset Selection
#The regsubsets() function (part of the leaps library) performs best subset selection by identifying the
#best model that contains a given number of predictors, where best is quantified using RSS. The summary()
#command outputs the best set of variables for each model size.
regfit.full = regsubsets(Cscore~., prostate)
summary(regfit.full)
regfit.full = regsubsets(Cscore~., data = prostate)
reg.summary = summary(regfit.full) #returns R2, RSS, adjusted R2, Cp, and BIC
names(reg.summary)
reg.summary$rsq
#R2 statistic increases from 52%, when only one variable is included in the model, to 59%, when all
#variables are included. As expected, the R2 statistic increases monotonically as more variables are
#included.

#Plotting RSS, adjusted R2, Cp, and BIC for all of the models at once will help us decide which model to
#select.
par(mfrow=c(2,2))
plot(reg.summary$rss ,xlab="Number of Variables ",ylab="RSS", type="l")
plot(reg.summary$adjr2 ,xlab="Number of Variables ", ylab="Adjusted RSq",type="l")
#What is the model with the largest adjusted R2 statistic?
which.max(reg.summary$adjr2)
points(5,reg.summary$adjr2[5], col="red",cex=2,pch=20)
#plot the Cp and BIC statistics, and indicate the models with the smallest statistic using which.min()
plot(reg.summary$cp, xlab="Number of Variables ", ylab="Cp", type='l')
which.min(reg.summary$cp)
points(5,reg.summary$cp [5],col="red",cex=2,pch=20)
which.min(reg.summary$bic)
plot(reg.summary$bic ,xlab="Number of Variables ",ylab="BIC", type='l')
points(2,reg.summary$bic [2],col="red",cex=2,pch=20)
#display the selected variables for the best model with a given number of predictors, ranked according to
#the BIC, Cp, adjusted R2, or AIC
plot(regfit.full,scale="r2")
plot(regfit.full,scale="adjr2")
plot(regfit.full,scale="Cp")
plot(regfit.full,scale="bic")
#the model with the highest adjusted R2 is the five-variable model
coef(regfit.full,5)
#the model with the lowest BIC is the two-variable model
coef(regfit.full,2)

#We got ambiguous results: 5 vs 2 predictor variables. Moreover, it estimates the training error better than
#the test error whereas we are more interested in the test error.
#Let's use cross-validation (page 213 of the book) which provides a direct estimate of the test error.

### Choosing Among Models Using Cross-Validation ###

#split the observations into a training set and a test set
set.seed(1)
train = sample(c(TRUE,FALSE), nrow(prostate), rep=TRUE)
test = (!train)
#apply regsubsets() to the training set in order to perform best subset selection
regfit.best = regsubsets(Cscore~.,data=prostate[train,])
#compute the validation set error for the best model of each model size
#1) make a model matrix from the test data
test.mat = model.matrix(Cscore~., data = prostate[test,]) #build an “X” matrix from data
#2) for each size i, extract the coefficients from regfit.best for the best model of that size, multiply
#   them into the appropriate columns of the test model matrix to form the predictions, and compute the test
#   MSE.
val.errors = rep(NA,7)
for(i in 1:7) {
  coefi = coef(regfit.best, id = i)
  pred = test.mat[,names(coefi)]%*%coefi
  val.errors[i] = mean((prostate$Cscore[test]-pred)^2)
}
#the best model is the one that contains ___ variables.
val.errors
which.min(val.errors)












#2022

train = sample(1:nrow(prostate), nrow(prostate)*2/3)
train
test=(!train)
test

#preform a LOOCV and k-fold of 10 CV for 
#different degree polynomials
loocv.error=rep(0,3)
cv10.error=rep(0,3)
for (i in 1:3){
  glm.fit=glm(caroBlood~poly(alcohol,i),data=prostate[train,])
  loocv.error[i]=cv.glm(prostate[train,],glm.fit)$delta[1]
  cv10.error[i]=cv.glm(prostate[train,],glm.fit,K=10)$delta[1]
}

#plot the results
plot(loocv.error,type="b",ylim=c(min(loocv.error,cv10.error),max(loocv.error,cv10.error)))
lines(cv10.error, type="b",col=2)

#Seemingly the model does not improve by using a higher order
#polynomial, since delta (the MSE) does not reduce after
#a linear model. Thus, we will stick with the linear model.


#make a general model to test data on
glm.mod <- glm(caroBlood~alcohol, data = prostate[train,])
glm.test <- predict(glm.mod, newdata = prostate[test,])
glm.MSE <- mean((prostate[test,'caroBlood'] - glm.test)^2)

#now we can see if the outlier is influencing the model
caroAlc.glm <- glm(caroBlood~alcohol ,data=prostate)
summary(caroAlc.glm)

#will test the leverage of the points
plot(hatvalues(caroAlc.glm))
which.max (hatvalues(caroAlc.glm))
prostate[62,]

#Just for fun let's rerun our analysis without the outlier
#to see if the model becomes significant
rmOut <- prostate[-62,]
boxplot(rmOut[,'alcohol'])

#let's test the leverage of the outlier
#an observation has significant leverage if it
#has a hatvalue larger than 2*p/ n
threshold <- 2*(dim(prostate)[2]) / (dim(prostate)[1])
threshold

hatvalues(caroAlc.glm)[62]
#the outlier appears to not have significant leverage 

#re-do testing wth the outlier removed just ot be safe
set.seed(1)
trainRM = sample(1:nrow(rmOut), nrow(rmOut)*2/3)
trainRM
testRM=(-trainRM)

#CV again
loocv.error=rep(0,3)
cv10.error=rep(0,3)
for (i in 1:3){
  glm.fit=glm(caroBlood~poly(alcohol,i),data=rmOut[trainRM,])
  loocv.error[i]=cv.glm(rmOut[trainRM,],glm.fit)$delta[1]
  cv10.error[i]=cv.glm(rmOut[trainRM,],glm.fit,K=10)$delta[1]
}

#plot the results
plot(loocv.error,type="b",ylim=c(min(loocv.error,cv10.error),max(loocv.error,cv10.error)))
lines(cv10.error, type="b",col=2)
#even with the outlier removed there linear model still seems to have
#the lowest MSE

#make glm model with all of the data
caroAlc.glm.rm <- glm(caroBlood~alcohol, data = rmOut)
summary(caroAlc.glm.rm)
#the p-value for the slope actually increased, so the
#outlier is likely not directly influencing the impact
#of the relationship between caroBlood and alcohol alone


#now the final model w/ the outlier is
caroAlc.glm

plot(y=prostate[,'caroBlood'], x=prostate[,'alcohol'], ylab = "caroBlood",
    xlab = "alcohol")
abline(caroAlc.glm, col = 2)

#now we can try with with bootstrapping since
#the alcohol variable does not follow a normal distribution

####create a boot function for a linear model
boot.fn=function (d ,index)
  return (coef(glm(caroBlood~alcohol ,data=d,subset =index)))


#obtain boot strap results
boot.mod.lin <- boot(prostate[train,], boot.fn, 1000)
boot.mod.lin
#produes the same results as the glm.mod


####create a boot function for a quadratic model
boot.fn=function (d ,index )
  coefficients(glm(caroBlood~alcohol +I( alcohol ^2) ,data=d ,subset =index))

#obtain boot strap results
boot.mod.quad <- boot(prostate[train,] ,boot.fn ,1000)
boot.mod.quad
#the standard error of the second degree model is
#seemingly narrow range around 0, suggesting
#this term is not statistically significant for the model
#thus we will not consider it when comparing models

#seemingly the linear model's coefficients have a smaller
#standard error sum than with the quadratic. We again
#find that the linear model is the best approximation.
#However, from the model generated we don't have 
#significant evidence to think there is an association between
#alcohol and caroBlood alone.
summary(caroAlc.glm)


####test w/ smoothing and natrual spline
#begin with smoothing
smooth.gam <- gam(caroBlood~s(alcohol,4), data = prostate[train,])

#obtain predictions
smooth.pred <- predict(smooth.gam, newdata = prostate[test,])
smooth.MSE <- mean((prostate[test,'caroBlood'] - smooth.pred)^2)
  
#now test natural spline
natural.gam <- gam(caroBlood~ns(alcohol,4), data= prostate[train,])

#obtain predictions 
natural.pred <- predict(natural.gam, newdata = prostate[test,])
natural.MSE <- mean((prostate[test,'caroBlood'] - natural.pred)^2)


#compare the models
glm.MSE
smooth.MSE
natural.MSE
anova( glm.mod,smooth.gam, natural.gam, test = "Chisq")
#seemingly the spline models do not differ
#so we go with the generalized linear model as
#this has the same capability of expliaing the observed
#variance and the easiest interpretability 
######################################







#2023

#3. Make an appropriate LASSO model, with the appropriate link and error function, and evaluate the
#   prediction performance. Do you see evidence that over-learning is an issue?














#2022
  
#############Question 3###############
set.seed(1)
#As seen during the exploratory analysis there is a drastic
#difference in the number of males and females.

#begin by comparing box plots
boxplot(prostate[,'caroBlood'] ~ prostate[,'gender'], ylab = 'caroBlood',
        xlab = 'gender', col = c("#FFE0B2", "#FFA726"))

#now let's compare the variance 
#Welch's test for uneuqal variance does not assume
#the two groups to have equal size or variance 
welch.test(caroBlood~gender, prostate[])

#we find there is likely a difference in the mean value of
#carotene in blood for males and females 

#let's build a logistic regression model
#to see if caroBlood is a good predictor for sex
sex.glm <- glm(gender~caroBlood, data=prostate[train,], family = 'binomial')
summary(sex.glm)

#make predictions
sex.glm.test <- predict(sex.glm, prostate[train,])
sex.glm.mse <- mean((sex.glm.test - prostate[train,'caroBlood'])^2)

#let's also build a model to see if gender can
#predict caroblood levels (as the goal of the study
#is to see which predictors are good at forecasting 
#a pateitns caroBlood levels)
cb.glm<- glm(caroBlood~gender,data=prostate[train,])
summary(cb.glm)

#make predictions
cb.glm.test <- predict(cb.glm, prostate[train,])
cb.glm.mse <- mean((cb.glm.test - prostate[train,'caroBlood'])^2)

#compare the models
sex.glm.mse
cb.glm.mse
#It seems like there is predictive power in both models
#However, we can see that prediciting caroBlood levels
#from gender alone is better

#So the means, and thus variance, of the males
#and females are statistically different according to 
#Welch's T-test. Additionaly, there is seemingly
#evidence of statistical significance in there being a
#predictive relationship between caroBlood and gender.
#however, since the size of males samples was rather low, there
#could in fact be a bias since maybe only certain types
#of men were recruited or tempted by the study....

######################################  

#4. Make your best model for predicting carotene levels in the blood. 
#Make and compare a ridge regression model, 
#a linear model and a generalized additive model (GAM).

#############Question 4###############
set.seed(1)

#convert categorical variables into dummy variables
data_dummy <- model.matrix( ~ .-1, prostate[, c(2,3,5)])
colnames(data_dummy)
dim(data_dummy)
dim(numericData)

newData <- cbind(numericData,data_dummy)
newData <- newData %>% clean_names()
names(newData)
#smoking: former (smokingFormer = 1, smokingCurrent = 0)
#smoking: never (smokingFormer = 0, smokingCurrent = 0)
#smoking: current (smokingFormer = 0, smokingCurrent = 1)
#gender: female (genderfemale = 1, gendermale = 0)
#gender: male (genderfemale = 0, gendermale = 1)
#vitalSuppl: yes,fairly often (not often = 0, No =0)
#vitalSuppl: yes, not often (not often = 1, No = 0)
#vitalSuppl: no (not often = 0, No =0)

###ridge regression
#generate train and test sets
set.seed(1)
train = sample(1:nrow(newData), nrow(newData)*2/3)
train
test=(-train)

#split response and other predictors
cb <- newData[,10]
preds <- as.matrix(newData[,-10])

#let's run a CV of 10-fold
ridge.mod <- glmnet(y=cb[train], x = preds[train,], alpha = 0)
plot(ridge.mod)
cv.ridge <- cv.glmnet(x = preds[train,], y = cb[train], alpha = 0)
plot(cv.ridge)
cv.ridge$lambda.min
#make a prediction
ridge.pred <- predict(ridge.mod,s=cv.ridge$lambda.min,newx=preds[test,],type="response")
plot(ridge.pred~cb[test])

#RSS
ridge.rss <- sum((ridge.pred - cb[test])^2)
ridge.rss

###linear model 
#let's do a best subset selection
#again we will do a K-fold CV on the training data
#to determine which size model we should select for
#our linear model
set.seed(1)
k = 10
trainData <- prostate[train,]
folds <- sample(1:k, nrow(trainData), replace = T)
cv.errors <- matrix(NA,k,15, dimnames= list(NULL< paste(1:15)))

#perform cross validation on best subset selection
for(j in 1:k)
{
  best.fit <- regsubsets(caroBlood~.,data = trainData[folds!=j,],
                         nvmax=15)
  for(i in 1:15)
  {
    pred = predict(best.fit, trainData[folds==j,], id = i) #changed from predict to predict.regsubsets
    cv.errors[j,i]=sum((trainData$caroBlood[folds==j]-pred)^2)
  }
}

#calculate the RSS
RSS.cv.errors = apply(cv.errors,2,sum)
RSS.cv.errors
#plot and determine which model is the best
plot(RSS.cv.errors, type = 'b')
numCoef <- which.min(RSS.cv.errors)

#seems like a model with 10 variables has the best RSS
full.best <- regsubsets(caroBlood~., data = trainData, nvmax = 15)
full.coefs <- coef(full.best, numCoef)
#determine the test errror
full.RSS <- sum((predict(full.best, prostate[test,], id = numCoef) - (prostate[test,'caroBlood']))^2)

#we can also do a PCR to construct a linear model
pcr.fit <- pcr(caroBlood~., data = prostate[train,],
               scale = T, validation = 'CV')

validationplot(pcr.fit,val.type="MSEP")
#we see choosing 5 components is the best

#make predictions and obtain RSS for PCR 
pcr.pred <- predict(pcr.fit, prostate[test,], ncomp = 5)
pcr.RSS <- sum((prostate[test,'caroBlood'] - pcr.pred)^2)

##GAM
set.seed(1)
#maybe we should only test the variables found to be
#strongest in predicting for the linear model
nD <- newData[,c("age", "genderfemale", "smoking_former", "smoking_current_smoker",
                 "weight_over_height_sq","vita_suppl_yes_not_often", "vita_suppl_no",
                 "fiber", "caro_diet", "retinol_blood", "caro_blood")]
names(nD)

#now plot a gam with the features found in best subset selection
plot(gam(caro_blood~s(age,4) + genderfemale + smoking_former + smoking_current_smoker +
           s(weight_over_height_sq, 4) + vita_suppl_yes_not_often + vita_suppl_no +
           s(fiber,4) + s(caro_diet,4) + s(retinol_blood,4),
    k =10,data=nD[train,]), se = T, col = 'blue')

#from the plots alone it seems like age, weight_over_height_sq, fiber, caro_diet,
#and retinol_blood may have non linear relations

gam.mod <- gam(caro_blood~s(age,4) + genderfemale + smoking_former + smoking_current_smoker +
                 s(weight_over_height_sq, 4) + vita_suppl_yes_not_often + vita_suppl_no +
                 s(fiber,4) + s(caro_diet,4) + s(retinol_blood,4),
               k =10,data=nD[train,])

#run a summary to confirm
summary(gam.mod)
#from the nonparametric effects it seems like only fiber
#and caro_diet follow a non linear model


gam.mod <- gam(caro_blood~age + genderfemale + smoking_former + smoking_current_smoker +
                 weight_over_height_sq + vita_suppl_yes_not_often + vita_suppl_no +
                 s(fiber,4) + s(caro_diet,4) + retinol_blood,
               k =10,data=nD[train,])

summary(gam.mod)
#from the two parametric summaries we see that age, smoking_former and 
#vita_suppl_yes_not_often were no significant in either approach
#so we can remove these to reduce the complexity of the model

gam.mod <- gam(caro_blood~genderfemale +  smoking_current_smoker +
                 weight_over_height_sq + vita_suppl_no +
                 s(fiber,4) + s(caro_diet,4) + retinol_blood,
               k =10,data=nD[train,])
#verify that removing the variable did not change the 
#significance of the others 
summary(gam.mod)
#retional blood may not be significant but we leave it in there
#just in case

#now test to so if a natural spline may be a better chouse
gam.ns.mod <- gam(caro_blood~genderfemale +  smoking_current_smoker +
                 weight_over_height_sq + vita_suppl_no +
                 ns(fiber,4) + ns(caro_diet,4) + retinol_blood,
               k =10,data=nD[train,])

summary(gam.ns.mod)

#now predit
gam.pred <- predict(gam.mod, newData[test,])
gam.RSS <- sum((gam.pred - newData[test,'caro_blood'])^2)

gam.ns.pred <- predict(gam.ns.mod, newData[test,])
gam.ns.RSS <- sum((gam.ns.pred - newData[test, 'caro_blood'])^2)

#compare the two splines
gam.ns.RSS
gam.RSS

###comparison of the ridge, lm and gam models from part 4
ridge.rss
full.RSS
pcr.RSS
gam.RSS
gam.ns.RSS

#the GAM model has the highest RSS and it's interpretablility
#is low so therefore we exclude it when comparing the linear
#models obtained from ridge, pcr and best subset selection

#However, ridge accounts for interacting predictors better
#than a normal linear model as produced with the best subset selection.
#But PCR also account for associations in the variables. But
#we choose to select ridge regression as the best model,
# because although we lose some interpretability compared to 
#subset seletion we lose more interpretability with PCR. 
coef(ridge.mod)[,which.min(cv.ridge$lambda)]

######################################

#5. Summarize your analyses. 
#Discuss which is the most important dietary contribution, 
#and the effect of correlations between variables.

#############Question 5###############
#we can use a summary on our all models
#but it seems like there is no clear winner
#of which dierary contribution is most important
summary(ridge.mod)
summary(full.best)
summary(gam.mod)

#however, we can test for the most influential dietary contribution
#with a random forest approach
library(randomForest)

set.seed(1)
#perform randomforrest
rf.d <- randomForest(caroBlood~., data = prostate[test,], 
                      mtry = sqrt(dim(prostate)[2]-1), ntree = 5000, importantce = T)

#results
rf.d

#let's see how well the model did 
yhat.bag <- predict(rf.d, newddata=prostate[test,])
sum((yhat.bag - prostate[test,'caroBlood'])^2)
#This RSS value is the worst for any model
#we have constructed yet, but bagging does
#account for interactions 

#look at the importance of each variable 
importance(rf.d)
varImpPlot(rf.d)

#when considering effects of correlations between variables
#we can look more into the interactions via fitting a linear
#model that accounts for all possible interactions and then calling anova
#to compare the effects of the variables to one another 
interactions <- lm(caroBlood ~ (.)^2,data=prostate)
anova(interactions)
#seemingly there may be interactions between age and caroDiet, smoking
#and cal, vitaSuppl and fat. However, this a rough approach towards assessing
#interactions between variables.

#PCR analysis could also give more insight
#we will use the full data set as we are just
#trying to understand the correlations and interactions
#between vairbale:
pcr.mod <- pcr(caroBlood~., data = prostate, scale = T,
               validation = "CV")
validationplot(pcr.mod, val.type = "MSEP")
summary(pcr.mod)
#seemingly there is not too much collinearity
#This is because we have to use 4 components
#to explain about 50% of the variation over
#15 predictors 
#additionally CV selects a PCR model with 7 components
######################################