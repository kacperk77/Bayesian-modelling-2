#Autor: Kacper Kalinowski 76975

library(tidyverse)
library(Metrics)
if(!require("mvtnorm")) {install.packages("mvtnorm"); library(mvtnorm)}
if(!require("coda")) {install.packages("coda"); library(coda)}
library(bridgesampling)
remove(list=ls())


setwd('C:\\Users\\lenovo\\Desktop\\ekonometria bayesowska\\pd2')


data <- read.csv('test_scores.csv')

#przekodowanie zmiennych
data$Urban<- ifelse(data$school_setting %in% c('Urban'), 1, 0)
data$Suburban<- ifelse(data$school_setting %in% c('Suburban'), 1, 0)

data$Public<- ifelse(data$school_type %in% c('Public'), 1, 0)

#model dla wszystkich obserwacji
lm <- lm(data$posttest~ +data$n_student+data$Public+data$Urban+data$Suburban)

summary(lm)



#losowanei zbioru
set.seed(1234)
sample_data <- data %>%  sample_n(100)


sample_data$Urban<- ifelse(sample_data$school_setting %in% c('Urban'), 1, 0)
sample_data$Suburban<- ifelse(sample_data$school_setting %in% c('Suburban'), 1, 0)


ols.results <- lm(sample_data$posttest~ +sample_data$n_student+sample_data$school_type+sample_data$Urban+sample_data$Suburban)

summary(ols.results)




#Ustalamy parametry rozk³adu a posteriori ze skrajnie nieinformacyjnym rozk³adem a priori N-G
ols.beta <- ols.results$coefficients
ols.sigma <- vcov(ols.results)
ols.sum.sq.residuals <- sum(ols.results$residuals ^ 2)

v.posterior <- nrow(sample_data) - length(ols.beta)
beta.posterior <- ols.beta
U.posterior <- ols.sigma / (ols.sum.sq.residuals / v.posterior)
s2.posterior <- ols.sum.sq.residuals / v.posterior


#Próbkowanie z funkcji q
library(mvtnorm)
S <- 10 ^ 5

losowanie.beta <- rmvt(S, sigma = s2.posterior * U.posterior,
                       df = v.posterior, delta = beta.posterior,
                       type = "shifted")
colnames(losowanie.beta) <- rownames(summary(ols.results)$coefficients)
head(losowanie.beta) #Ten sam wspó³czynnik wystêpuje z ró¿nymi znakami

beta0 <- losowanie.beta[, 1] # sta³a
beta1 <- losowanie.beta[, 2] # liczba studentow
beta2 <- losowanie.beta[, 3] # szkola publiczna
beta3 <- losowanie.beta[, 4] # szkola zlokalizowana w miescie
beta4 <- losowanie.beta[, 5] # szkola zlokalizowana na przedmiesciach


t(t(apply(losowanie.beta,2,mean)))

#Importance sampling
important.beta <- losowanie.beta[beta1 < 0 & beta2 < 0 & beta3 > 0 & beta4 > 0 , ]

t(t(apply(important.beta,2,mean)))[1,]

#Prawdopodobieñstwo a posteriori spe³nienia restrykcji

#################################################################

#Metoda 1
nrow(important.beta)/nrow(losowanie.beta)

#Du¿a liczba obserwacji spe³nia restrykcje



#Ilustracja gêstoœci parametru przy zmiennej n_student
par(mfrow = c(2, 2))

important.beta1 <- important.beta[, 2]


unrestricted.h <- hist(beta1,  col = "gray", 
                       xlab = "liczba uczniow w klasie: parametr", 
                       main = "Rozk³ad a posteriori przy skrajnie nieinformacyjny rozk³adzie a priori oraz
                       przy uciêtym rozk³adzie a priori")

restricted.h <- hist(important.beta1, 
                     col = "red", 
                     xlab = "liczba uczniow w klasie: parametr", 
                     main = "Rozk³ad a posteriori przy uciêtym rozk³adzie a priori", add = T)

     


#Ilustracja gêstoœci parametru przy zmiennej Public
important.beta2 <- important.beta[, 3]


unrestricted.h <- hist(beta2,  col = "gray", 
                       xlab = "szko³a publiczna: parametr", 
                       main = "Rozk³ad a posteriori przy skrajnie nieinformacyjny rozk³adzie a priori oraz
                       przy uciêtym rozk³adzie a priori")

restricted.h <- hist(important.beta2, 
                     col = "red", 
                     xlab = "szko³a publiczna: parametr", 
                     main = "Rozk³ad a posteriori przy uciêtym rozk³adzie a priori oraz
                       przy uciêtym rozk³adzie a priori", add = T,
                     breaks=unrestricted.h$breaks)



#Ilustracja gêstoœci parametru przy zmiennej Urban
important.beta3 <- important.beta[, 4]


unrestricted.h <- hist(beta3,  col = "gray", 
                       xlab = "szko³a miejska: parametr", 
                       main = "Rozk³ad a posteriori przy skrajnie nieinformacyjny rozk³adzie a priori oraz
                       przy uciêtym rozk³adzie a priori")

restricted.h <- hist(important.beta3, 
                     col = "red", 
                     xlab = "szko³a miejska: parametr", 
                     main = "Rozk³ad a posteriori przy uciêtym rozk³adzie a priori oraz
                       przy uciêtym rozk³adzie a priori", add = T,
                     breaks=unrestricted.h$breaks)


#Ilustracja gêstoœci parametru przy zmiennej Suburban

important.beta4 <- important.beta[, 5]


unrestricted.h <- hist(beta4,  col = "gray", 
                       xlab = "szko³a podmiejska: parametr", 
                       main = "Rozk³ad a posteriori przy skrajnie nieinformacyjny rozk³adzie a priori oraz
                       przy uciêtym rozk³adzie a priori")

restricted.h <- hist(important.beta4, 
                     col = "red", 
                     xlab = "szko³a podmiejska: parametr", 
                     main = "Rozk³ad a posteriori przy uciêtym rozk³adzie a priori", add = T,
                     breaks=unrestricted.h$breaks)


library(coda)  #Pakiet s³u¿¹cy do wyliczenia przedzia³ów HPD (Plummer et al., 2006)

beta.restricted <- apply(important.beta, 2, mean)
hpd.restricted <- HPDinterval(mcmc(important.beta)) #funkcja licz¹ca bayesowskie przedzia³y HPD (95%)
restricted.results <- cbind(hpd.restricted[,1], beta.restricted, hpd.restricted[,2])
rownames(restricted.results) <- names(ols.beta)
colnames(restricted.results) <- c("Dolna granica HPD","Oszacowanie z ograniczeniami","Górna granica HPD")
restricted.results <- round(restricted.results, digits = 3)


frequentist.results <- cbind(confint(ols.results)[,1],ols.results$coefficients,confint(ols.results)[,2])
colnames(frequentist.results) <- c("Dolna granica 95% p.u.","Oszacowanie bez ograniczeñ", "Górna granica 95% p.u.");
frequentist.results <- round(frequentist.results, digits=3)

summary(ols.results)


#ocena predykcji model OLS vs zaproponowany model

'%nin%' <- Negate('%in%')

data_pred <- data %>% filter(student_id %nin% sample_data$student_id )

data_pred2 <- data_pred %>%  sample_n(100)



  
data_pred2$prediction_OLS <- predict(ols.results, newdata = data_pred2)
data_pred2$prediction_proposed_model <- restricted.results[1,2]+data_pred2$n_student*restricted.results[2,2]+
                                        data_pred2$Public*restricted.results[3,2]+data_pred2$Urban*restricted.results[4,2]+
                                        data_pred2$Suburban*restricted.results[5,2]     




RMSE_OLS <- rmse(data_pred2$posttest, data_pred2$prediction_OLS)

RMSE_proposed_model <- rmse(data_pred2$posttest, data_pred2$prediction_proposed_model)

mean(data_pred2$posttest)

