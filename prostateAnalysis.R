#Mariia KIM   r0831824
#Soniya LAMA   r0715145
#Statistical Methods for Bioinformatics

#load packages
library("Hmisc")
library("pastecs")
library("doBy")
library("stats")
#library("raster")
library("XLConnect")

library("ggplot2")
library(reshape2)
library(corrplot)
library(leaps)

library(glmnet)
library(tidyverse)
library(boot)
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

############# Question 1 ###############
#Overview
dim(prostate)
summary(prostate)
#make a subset of the data w/ only quantitative values
prostateNumeric <- subset(prostate, select = -c(svi))
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
boxplot(prostate$lcavol,main="lcavol",ylab="lcavol", col = "gray")

boxplot(prostate$lweight,main="lweight",ylab="lweight", col = "gray")
boxplot(prostate$age,main="age",ylab="age", col = "gray")

boxplot(prostate$lbph,main="lbph",ylab="lbph", col = "gray")
boxplot(prostate$lcp,main="lcp",ylab="lcp", col = "gray")

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

################### question 2 ################### 

# Create boxplots for all variables in your data frame
par(mfrow=c(1,1))
boxplot(prostateNumeric)
# Identify the outliers based on the boxplot
outliers <- boxplot(prostateNumeric)$out
outliersSvi <- boxplot(prostate)$out
# Remove the outliers from the data frame
prostateCleanAgeSvi <- prostate[!(prostate$age %in% outliers), ]
prostateCleanAge <- prostateNumeric[!(prostateNumeric$age %in% outliers), ]
boxplot(prostateCleanAge)
prostateCleanCscoreSvi <- prostate[!(prostate$Cscore %in% outliers), ]
prostateCleanCscore <- prostateNumeric[!(prostateNumeric$Cscore %in% outliers), ]
boxplot(prostateCleanCscore)
prostateCleanLpsaSvi <- prostate[!(prostate$lpsa %in% outliers), ]
prostateCleanLpsa <- prostateNumeric[!(prostateNumeric$lpsa %in% outliers), ]
boxplot(prostateCleanLpsa)
cleanProstateSvi <- prostate[!(c(prostate$age, prostate$Cscore, prostate$lpsa) %in% outliers), ]
cleanProstate <- prostateNumeric[!(c(prostateNumeric$age, prostateNumeric$Cscore, prostateNumeric$lpsa) %in% outliers), ]
boxplot(cleanProstate)



#Best Subset Selection
#The regsubsets() function (part of the leaps library) performs best subset selection by identifying the
#best model that contains a given number of predictors, where best is quantified using RSS. The summary()
#command outputs the best set of variables for each model size.
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
plot(reg.summary$rss ,xlab="Number of Variables",ylab="RSS", type="l")
plot(reg.summary$adjr2 ,xlab="Number of Variables", ylab="Adjusted RSq",type="l")
#What is the model with the largest adjusted R2 statistic?
which.max(reg.summary$adjr2)
points(5,reg.summary$adjr2[5], col="red",cex=2,pch=20)
#plot the Cp and BIC statistics, and indicate the models with the smallest statistic using which.min()
plot(reg.summary$cp, xlab="Number of Variables", ylab="Cp", type='l')
which.min(reg.summary$cp)
points(5,reg.summary$cp [5],col="red",cex=2,pch=20)
which.min(reg.summary$bic)
plot(reg.summary$bic,xlab="Number of Variables",ylab="BIC", type='l')
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

#Best Subset Selection with only numeric variables (no svi)
regfit.full.numeric = regsubsets(Cscore~., prostateNumeric)
summary(regfit.full.numeric)
regfit.full.numeric = regsubsets(Cscore~., data = prostateNumeric)
reg.summary.numeric = summary(regfit.full.numeric) #returns R2, RSS, adjusted R2, Cp, and BIC
names(reg.summary.numeric)
reg.summary.numeric$rsq
par(mfrow=c(2,2))
plot(reg.summary.numeric$rss, xlab="Number of Variables", ylab="RSS", type="l")
points(4, reg.summary.numeric$rss[4], col="red", cex=2, pch=20)
plot(reg.summary.numeric$adjr2, xlab="Number of Variables", ylab="Adjusted RSq", type="l")
#What is the model with the largest adjusted R2 statistic?
which.max(reg.summary.numeric$adjr2)
points(4, reg.summary.numeric$adjr2[4], col="red", cex=2, pch=20)
#plot the Cp and BIC statistics, and indicate the models with the smallest statistic using which.min()
plot(reg.summary.numeric$cp, xlab="Number of Variables", ylab="Cp", type='l')
which.min(reg.summary.numeric$cp)
points(4, reg.summary.numeric$cp[4], col="red", cex=2, pch=20)
plot(reg.summary.numeric$bic, xlab="Number of Variables", ylab="BIC", type='l')
which.min(reg.summary.numeric$bic)
points(2, reg.summary.numeric$bic[2], col="red", cex=2, pch=20)

#Best Subset Selection with only numeric variables (no svi), no age outliers
regfit.full.ageout = regsubsets(Cscore~., prostateCleanAge)
summary(regfit.full.ageout)
regfit.full.ageout = regsubsets(Cscore~., data = prostateCleanAge)
reg.summary.ageout = summary(regfit.full.ageout) #returns R2, RSS, adjusted R2, Cp, and BIC
reg.summary.ageout$rsq
par(mfrow=c(2,2))
plot(reg.summary.ageout$rss, xlab="Number of Variables", ylab="RSS", type="l")
points(4, reg.summary.ageout$rss[4], col="red", cex=2, pch=20)
plot(reg.summary.ageout$adjr2, xlab="Number of Variables", ylab="Adjusted RSq", type="l")
#What is the model with the largest adjusted R2 statistic?
which.max(reg.summary.ageout$adjr2)
points(4, reg.summary.ageout$adjr2[4], col="red", cex=2, pch=20)
#plot the Cp and BIC statistics, and indicate the models with the smallest statistic using which.min()
plot(reg.summary.ageout$cp, xlab="Number of Variables", ylab="Cp", type='l')
which.min(reg.summary.ageout$cp)
points(4, reg.summary.ageout$cp[4], col="red", cex=2, pch=20)
plot(reg.summary.ageout$bic, xlab="Number of Variables", ylab="BIC", type='l')
which.min(reg.summary.ageout$bic)
points(2, reg.summary.ageout$bic[2], col="red", cex=2, pch=20)

#Best Subset Selection with only numeric variables (no svi), no Cscore outliers
regfit.full.cscoreout = regsubsets(Cscore~., prostateCleanCscore)
summary(regfit.full.cscoreout)
regfit.full.cscoreout = regsubsets(Cscore~., data = prostateCleanCscore)
reg.summary.cscoreout = summary(regfit.full.cscoreout) #returns R2, RSS, adjusted R2, Cp, and BIC
reg.summary.cscoreout$rsq
par(mfrow=c(2,2))
plot(reg.summary.cscoreout$rss, xlab="Number of Variables", ylab="RSS", type="l")
points(4, reg.summary.cscoreout$rss[4], col="red", cex=2, pch=20)
plot(reg.summary.cscoreout$adjr2, xlab="Number of Variables", ylab="Adjusted RSq", type="l")
#What is the model with the largest adjusted R2 statistic?
which.max(reg.summary.cscoreout$adjr2)
points(5, reg.summary.cscoreout$adjr2[5], col="red", cex=2, pch=20)
#plot the Cp and BIC statistics, and indicate the models with the smallest statistic using which.min()
plot(reg.summary.cscoreout$cp, xlab="Number of Variables", ylab="Cp", type='l')
which.min(reg.summary.cscoreout$cp)
points(2, reg.summary.cscoreout$cp[2], col="red", cex=2, pch=20)
plot(reg.summary.cscoreout$bic, xlab="Number of Variables", ylab="BIC", type='l')
which.min(reg.summary.cscoreout$bic)
points(2, reg.summary.cscoreout$bic[2], col="red", cex=2, pch=20)

#Best Subset Selection with only numeric variables (no svi), no lpsa outliers
regfit.full.lpsaout = regsubsets(Cscore~., prostateCleanLpsa)
summary(regfit.full.lpsaout)
regfit.full.lpsaout = regsubsets(Cscore~., data = prostateCleanLpsa)
reg.summary.lpsaout = summary(regfit.full.lpsaout) #returns R2, RSS, adjusted R2, Cp, and BIC
reg.summary.lpsaout$rsq
par(mfrow=c(2,2))
plot(reg.summary.lpsaout$rss, xlab="Number of Variables", ylab="RSS", type="l")
points(4, reg.summary.lpsaout$rss[4], col="red", cex=2, pch=20)
plot(reg.summary.lpsaout$adjr2, xlab="Number of Variables", ylab="Adjusted RSq", type="l")
#What is the model with the largest adjusted R2 statistic?
which.max(reg.summary.lpsaout$adjr2)
points(4, reg.summary.lpsaout$adjr2[4], col="red", cex=2, pch=20)
#plot the Cp and BIC statistics, and indicate the models with the smallest statistic using which.min()
plot(reg.summary.lpsaout$cp, xlab="Number of Variables", ylab="Cp", type='l')
which.min(reg.summary.lpsaout$cp)
points(4, reg.summary.lpsaout$cp[4], col="red", cex=2, pch=20)
plot(reg.summary.lpsaout$bic, xlab="Number of Variables", ylab="BIC", type='l')
which.min(reg.summary.lpsaout$bic)
points(2, reg.summary.lpsaout$bic[2], col="red", cex=2, pch=20)

#Best Subset Selection with only numeric variables (no svi), excluding age, Cscore and lpsa outliers
regfit.full.clean = regsubsets(Cscore~., cleanProstate)
summary(regfit.full.clean)
regfit.full.clean = regsubsets(Cscore~., data = cleanProstate)
reg.summary.clean = summary(regfit.full.clean) #returns R2, RSS, adjusted R2, Cp, and BIC
names(reg.summary.clean)
reg.summary.clean$rsq
par(mfrow=c(2,2))
plot(reg.summary.clean$rss, xlab="Number of Variables", ylab="RSS", type="l")
points(4, reg.summary.clean$rss[4], col="red", cex=2, pch=20)
plot(reg.summary.clean$adjr2, xlab="Number of Variables", ylab="Adjusted RSq", type="l")
#What is the model with the largest adjusted R2 statistic?
which.max(reg.summary.clean$adjr2)
points(4, reg.summary.clean$adjr2[4], col="red", cex=2, pch=20)
#plot the Cp and BIC statistics, and indicate the models with the smallest statistic using which.min()
plot(reg.summary.clean$cp, xlab="Number of Variables", ylab="Cp", type='l')
which.min(reg.summary.clean$cp)
points(4, reg.summary.clean$cp[4], col="red", cex=2, pch=20)
plot(reg.summary.clean$bic, xlab="Number of Variables", ylab="BIC", type='l')
which.min(reg.summary.clean$bic)
points(2, reg.summary.clean$bic[2], col="red", cex=2, pch=20)

### Choosing Among Models Using Cross-Validation ###

#split the observations into a training set and a test set
set.seed(1)
train = sample(c(TRUE,FALSE), nrow(prostate), rep=TRUE)
test = (!train)
#apply regsubsets() to the training set in order to perform best subset selection
regfit.best = regsubsets(Cscore~., data = prostate[train,])
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
#the best model is the one that contains 1 variable.
val.errors
which.min(val.errors)

#LOOCV
library(boot)
glm.fit = glm(Cscore ~ ., data = prostate)
cv.err = cv.glm(prostate, glm.fit)
cv.err$delta
cv.error = rep(0,5)
for (i in 1:5) {
  glm.fit = glm(Cscore ~ poly(lpsa,i), data = prostate)
  cv.error[i] = cv.glm(prostate, glm.fit)$delta[1]
}
cv.error
#1481.7016  777.7890  803.7298  891.4756 1008.7859
#sharp drop in the estimated test MSE between the linear and quadratic fits, but then no clear
#improvement from using higher-order polynomials

#k-Fold Cross-Validation with k = 10
set.seed(17)
cv.error.10 = rep(0,10)
for (i in 1:10){
  glm.fit = glm(Cscore ~ poly(lpsa,i), data = prostate)
  cv.error.10[i] = cv.glm(prostate, glm.fit, K = 10)$delta[1]
}
cv.error.10
#1529.8534   769.2539   795.8788   893.9604  1021.9147  1257.1746  1613.5471  2925.3239 24967.8074 71946.0660


#### Multiple Linear Regression ####
lm.fit=lm(Cscore ~ ., data = prostate)
summary(lm.fit) #lpsa very significant
sum(resid(lm.fit)^2) #RSS
lm.fit.numeric = lm(Cscore ~ ., data = prostateNumeric)
summary(lm.fit.numeric) #lcavol somewhat significant, lcp significant, lpsa very significant
lm.fit.ageout = lm(Cscore ~ ., data = prostateCleanAge)
summary(lm.fit.ageout) #lcavol somewhat significant, lcp significant, lpsa very significant
lm.fit.ageout.svi = lm(Cscore ~ ., data = prostateCleanAgeSvi)
summary(lm.fit.ageout.svi) #lcavol somewhat significant, lcp significant, lpsa very significant
lm.fit.cscoreout = lm(Cscore ~ ., data = prostateCleanCscore)
summary(lm.fit.cscoreout) #lcp and lpsa very significant
lm.fit.cscoreout.svi = lm(Cscore ~ ., data = prostateCleanCscoreSvi)
summary(lm.fit.cscoreout.svi) #lcp somewhat significant, lpsa very significant
lm.fit.lpsaout = lm(Cscore ~ ., data = prostateCleanLpsa)
summary(lm.fit.lpsaout) #lcp and lpsa very significant
lm.fit.lpsaout.svi = lm(Cscore ~ ., data = prostateCleanLpsaSvi)
summary(lm.fit.lpsaout.svi) #lcp significant, lpsa very significant
lm.fit.clean = lm(Cscore ~ ., data = cleanProstate)
summary(lm.fit.clean) #lcavol somewhat significant, lcp significant, lpsa very significant
lm.fit.clean.svi = lm(Cscore ~ ., data = cleanProstateSvi)
summary(lm.fit.clean.svi) #lcavol somewhat significant, lcp significant, lpsa very significant
#compute variance inflation factors
library(car)
vif(lm.fit)
vif(lm.fit.numeric)
vif(lm.fit.ageout)
vif(lm.fit.cscoreout)
vif(lm.fit.lpsaout)
vif(lm.fit.clean)



################### question 3 ################### 

#Make an appropriate LASSO model with the appropriate link and error function, 
#and evaluate the prediction performance. Do you see any evidence that over-learning is an issue? 

x = model.matrix(Cscore~.,prostateNumeric)[,-1]
y = prostateNumeric$Cscore

set.seed(1)
train=sample(1:nrow(x), nrow(x)/2)
test =(-train)
y.test = y[test]
library(glmnet)
grid = 10^seq(10,-2,length=100) 

lasso.mod = glmnet(x[train ,], y[train], alpha=1, lambda=grid)
plot(lasso.mod, label = TRUE)

set.seed (2)
cv.out = cv.glmnet (x[train ,],y[train],alpha =1) #10-fold cross validation
plot(cv.out)
bestlam =cv.out$lambda.min 
bestlam
#[1] 1.552688
lasso.pred=predict (lasso.mod ,s=bestlam ,newx=x[test ,])
mean(( lasso.pred -y.test)^2)
#[1] 709.4205

out=glmnet (x,y,alpha =1, lambda =grid)
lasso.coef=predict(out ,type = "coefficients",s=bestlam )[1:7 ,]
lasso.coef
#(Intercept)      lcavol     lweight         age        lbph         svi         lcp        lpsa 
#-8.320718   -2.369483   -6.331806    0.000000    0.000000   18.835766    4.882504   27.274831 



################### question 4 ###################
#Look at the coefficient for “lcavol” in your LASSO model. Does this coefficient correspond
#to how well it can predict Cscore? Explain your observation. 

#LASSO shrinks the coefficient estimate for "lcavol" to zero. Hence, "lcavol" is not part
#of the variable selection by LASSO. Let's determine if it corresponds to how well lcavol
#predicts the Cscore.

#1. Visualize the relationship between lcavol and Cscore
#plot the data: scatter plot of lcavol (predictor var) against Cscore (response var)
plot(prostate$lcavol, prostate$Cscore, xlab="Log cancer volume", ylab="Cscore", pch=16, cex=1.5, main="Cancer volume and cancer progression")
abline(lm(prostate$Cscore ~ prostate$lcavol), lwd=2, col="red")
lines(lowess(prostate$lcavol, prostate$Cscore, f=0.9), lwd=2, col="blue")
abline(lm(prostate$Cscore ~ prostate$lcavol + log(prostate$lcavol)), lwd=2, col="blue")
#The relationship does not seem to be linear. However, there seems to be some sort of an
#exponential relationship. This suggests that the predictor variable may be a good predictor
#of the response variable, just not in the linear model.

#2. Calculate the correlation coefficient between lcavol and Cscore
cov.lcavol.cscore <- cov(prostate$lcavol, prostate$Cscore)
cov.lcavol.cscore
corr.lcavol.cscore <- cor(prostate$lcavol, prostate$Cscore)
corr.lcavol.cscore #0.511
#The correlation coefficient between lcavol and Cscore is around 0.51 - not a strong indicator
#of a linear relationship.

#3. Build a regression model
#See section "Multiple Linear Regression". "lcavol" has been shown to be somewhat
#significant (0.024 < p-value < 0.043 numeric) in some of the models excluding some
#outliers, in others not significant at all.
res <- lm(Cscore ~ lcavol, data = prostate)
prostate.anova <- anova(res)
prostate.anova
prostate.summary <- summary(res)
prostate.summary
plot(prostate$lcavol, prostate$Cscore, xlab="Cancer volume", ylab="Cscore", main="Cscore vs cancer volume")
abline(res, col ="blue", lwd=2)
#lcavol alone seems to be a significant predictor of the Cscore, and the relationship seems
#to be linear. However, this regression model does not perform well for this data. The
#adjusted R-squared statistic is not very high - around 0.25.
#calculate mean squared error
mse <- mean(res$residuals^2)
rmse <- sqrt(mse)
range(prostate$Cscore)
median(prostate$Cscore)
rmse
#The root mean squared error is rather high, around 45 - higher than the median and the
#mean Cscores.
#This information combined indicates that the present model is not well-fitted for our data.

#In conclusion, the coefficient for lcavol estimated by LASSO does not entirely correspond
#to how well this variable can predict Cscore. This predictor has some significance for
#the prediction of Cscore, and their relationship seems to be non-linear.



################ Question 5 ####################
library(gam)
gam.m5 <- gam(Cscore ~ s(lcavol, 4) + s(lweight, 4) + s(lcp, 4) + s(lpsa, 4) + s(age, 4) + svi + s(lbph, 4), data = prostate)
head(prostate)

head(prostateNumeric)
##GAM
head(prostate)
set.seed(1)
#maybe we should only test the variables found to be
#strongest in predicting for the linear model
nD <- prostate[,c("Cscore", "lcavol", "lweight", "age", "lbph","svi",
                 "lcp","lpsa")]
names(nD)

library(gam)
#now plot a gam with the features found in best subset selection
plot(gam(Cscore ~ s(lcavol, 4) + s(lweight, 4) + s(lcp, 4) + s(lpsa, 4) + s(age, 4) + svi + s(lbph, 4),
         k = 10,data=prostate[train,]), se = T, col = 'blue')


#from the plots alone it seems like age, lbph, lcavol, lweight may have non-linear
#relationships
gam.mod <- gam(Cscore ~ s(lcavol, 4) + s(lweight, 4) + s(lcp, 4) + s(lpsa, 4) + s(age, 4) + svi + s(lbph, 4),
               k = 10, data = prostate[train,])

#run a summary to confirm
summary(gam.mod)
#from the nonparametric effects it seems like only lpsa follows a non-linear model

gam.mod <- gam(Cscore ~ lcavol + lweight + lcp  + age + svi + lbph + s(lpsa, 4),
               k = 10, data = prostate[train,])
summary(gam.mod)
#from the two parametric summaries we see that svi and lbph were not significant in either
#approach, so we can remove these two to reduce the complexity of the model
gam.mod <- gam(Cscore ~ lcavol + lweight + lcp  + age + s(lpsa, 4),
               k = 10, data = prostate[train,])
#verify that removing the variable did not change the 
#significance of the others 
summary(gam.mod)

#now test if a natural spline may be a better choice
gam.ns.mod <- gam(Cscore ~ lcavol + lweight + lcp  + age + ns(lpsa, 4),
               k = 10, data = prostate[train,])
summary(gam.ns.mod)

#predict
gam.pred <- predict(gam.mod, prostate[test,])
gam.RSS <- sum((gam.pred - prostate[test,'Cscore'])^2)

gam.ns.pred <- predict(gam.ns.mod, prostate[test,])
gam.ns.RSS <- sum((gam.ns.pred - prostate[test, 'Cscore'])^2)

#compare the splines
gam.ns.RSS
gam.RSS
