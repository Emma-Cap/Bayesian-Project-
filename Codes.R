library(readxl)
library(ggplot2)
library(cowplot)
library(leaps)
library(boot)
library(gridExtra)
library(nortest)
library(corrplot) 
library(copula)
library(copBasic)
library(brms)
library(mvtnorm)
library(truncnorm)
library(coda)
library(tibble)
library(knitr)
library(RColorBrewer)

data<-"BayesianDataset.xlsx"
data<- read_excel(data, sheet=1)


### Histogram ###
# of the number of companies based on their rating level. 
ggplot(data, aes(x = factor(Binary_Response))) +
  geom_bar(fill = c("Non-Investment Grade" = "steelblue", "Investment Grade" = "darkorange"), color = "black") +
  labs(x = "Rating", y = "Count", title = "Binary Response Variable") +
  scale_x_discrete(labels = c("Non-Investment Grade", "Investment Grade")) +
  theme_minimal() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))



data_filtered <- data[!is.na(data$"Debt/Equity"), ]


### Boxplot ###
# with the Binary Response and the covariate variable Debt/Equity
ggplot(data_filtered, aes(x = factor(Binary_Response), y = Debt_Equit, fill = factor(Binary_Response))) +
  geom_boxplot() +
  labs(x = "Rating", y = "Debt to Equity Ratio",
       title = "Debt to Equity Ratio by Rating", fill = "Legend") +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "darkorange"),
                    labels = c("0" = "Non-Investment Grade", "1" = "Investment Grade")) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1500)) +
  theme(plot.title = element_text(hjust = 0.5))


### Correlation Matrix ###

selected_data = data[, 4:ncol(data)]
complete_data = selected_data[complete.cases(selected_data),]
correlation_matrix <- cor(complete_data)
color_palette <- colorRampPalette(c("#7F0000", "white", "#0033CC"))(200)

corrplot(correlation_matrix, method="number", col=color_palette, 
         number.cex=0.65, number.digits=2, 
         tl.col="black", tl.srt=45, # Aggiunge etichette inclinate per evitare sovrapposizioni
         # Assegna le etichette alle righe e alle colonne
         tl.cex=0.6, cl.cex=0.6) # Controlla la dimensione del testo delle etichette




######################################
### LOGISTIC REGRESSION WITH JAGS  ###
######################################

data = data[complete.cases(data),]
y=data$Binary_Response  # Response Variable
X = as.data.frame(data[,c(4:14)])
X.x= model.matrix(y ~., X) 
p = ncol(X.x)

X = X.x

# Standardize the predictors:

for(j in 2:p){
  
  X[,j] = scale(X.x[,j])
}

# Specification of the data:

library(R2jags)
set.seed(1)

jags_data = with(data, list(y = y, X = X, n = length(y), p = ncol(X)))
logistic_regr_jags = function(){
  
  # Likelihood:
  
  for(i in 1:n){
    
    y[i] ~ dbern(pi[i])
    
    logit(pi[i]) = beta%*%X[i,]
    
  }
  
  # Priors:
  
  for(j in 1:p){
    
    beta[j] ~ dnorm(0, 0.01) # where 0.01 is the precision, so the variance is 100
    
  }
  
}

# Initial values for beta:


init_values = function(){
  
  list(beta = rep(0, p))
  
}


# We want to do posterior inference on beta:

params = c("beta")

# Jags 

jags_posterior = jags(data = jags_data,
                      inits = init_values,
                      parameters.to.save = params,
                      model.file = logistic_regr_jags,
                      n.chains = 1,
                      n.iter = 5000,
                      n.burnin = 500,
                      n.thin = 1)

post_theta=jags_posterior$BUGSoutput



###############################
### CONVERGENCE DIAGNOSTICS ###
###############################
post_beta = post_theta$sims.list$beta

plot(post_beta[,7], ylab = expression(paste(beta[i])),ylim=c(-2.5,0.5),xab = "Iterazione", type = "l")
plot(post_beta[,8], ylab = expression(paste(beta[i])),ylim=c(-1,8),  xlab = "Iterazione", type = "l")

par(oma = c(2, 2, 2, 2)) # Riduci i margini esterni
par(mfrow = c(1, 2))

# Grafico della funzione di autocorrelazione
for(i in 7:8){
  acf(post_beta[,i], xlab = "lag", main = "", mar = c(4, 4, 2, 2)) # Riduci i margini interni
}


# Geweke test
geweke_result <- geweke.diag(post_beta, frac1 = 0.3, frac2 = 0.4)
geweke_values <- geweke_result$z
geweke_matrix <- matrix(geweke_values, nrow = 1)
colnames(geweke_matrix) <- paste("Beta", 1:length(geweke_values))
kable(geweke_matrix, caption = "Geweke Test", align = 'c')

# Effective Sample Size (ESS) 
ess_values <- numeric(p)
for(i in 1:p){
  ess_values[i] <- round(coda::effectiveSize(mcmc(post_beta[,i])), 1)
}
ess_matrix <- matrix(ess_values, nrow = 1)
colnames(ess_matrix) <- paste("Beta", 1:p)
kable(ess_matrix, caption = "Effective Sample Size", align = 'c')

# Parameter Analysis
boxplot(jags_posterior$BUGSoutput$sims.array[,1,1:11], ylim=c(-6,23), col="orange")
abline(h=0, col="blue")

# Quantiles
names(post_theta$sims.list$beta) = gsub("\\s+", "", tolower(names(post_theta$sims.list$beta)))
res = t(apply(post_theta$sims.list$beta, 2, quantile, prob=c(.025,.975)))
rownames(res) = colnames(X)
colnames(res) = c("2.5 % Quantile","97.5 % Quantile")
knitr::kable(res)

# Posterior Distribution of the parameters
resu = colMeans(post_theta$sims.list$beta)
res = t(apply(post_theta$sims.list$beta, 2, quantile, prob=c(.025,.975)))
rownames(res) = colnames(X)
par(mfrow=c(1,2)) 
num_parametri_beta <- ncol(post_theta$sims.list$beta)
i <- 6
hist(post_theta$sims.list$beta[,i],
     main=paste("Posterior distribution of Net Income"),
     xlab = rownames(res)[i], cex=0.8, prob = TRUE, cex.main = .9, col="steelblue")

abline(v=res[i,1], col="red", lty=2, lwd=3)
abline(v=res[i,2], col="red", lty=2, lwd=3)
abline(v=resu[i], col="orange", lwd=3)
i <- 8
hist(post_theta$sims.list$beta[,i],
     main=paste("Posterior distribution of Accounts Receivable"),
     xlab = rownames(res)[i], cex=0.8, prob = TRUE, cex.main = .9, col="steelblue")

abline(v=res[i,1], col="red", lty=2, lwd=3)
abline(v=res[i,2], col="red", lty=2, lwd=3)
abline(v=resu[i], col="orange", lwd=3)





##########################################################
### Bayesian Model Averaging with Spike and Slab prior ###
##########################################################

data = data[complete.cases(data),]
y=data$Binary_Response  # Response Variable
X = as.data.frame(data[,c(4:14)])
X.x= model.matrix(y ~., X) 
p = ncol(X.x)

X = X.x

# Standardize the predictors:

for(j in 2:p){
  
  X[,j] = scale(X.x[,j])
}


# Specification of the data:

library(R2jags)

jags_data = with(data, list(y = y, X = X, n = length(y), p = ncol(X)))
logistic_regr_jags = function(){
  
  # Likelihood:
  
  for(i in 1:n){
    
    y[i] ~ dbern(pi[i])
    
    logit(pi[i]) = (gamma*beta)%*%X[i,]
    
  }
  
  # Priors:
  
  for(j in 1:p){
    
    beta[j] ~ dnorm(0, 0.01) # where 0.01 is the precision, so the variance is 100
    
  }
  
  for(j in 1:p){
    
    gamma[j] ~ dbern(w)
    
  }
  
  w ~ dbeta(1, 1)
  
}

# Initial values for beta:

init_values = function(){
  
  list(beta = rep(0, p), gamma = rep(1, p))
  
}


# We want to do posterior inference on beta:

params = c("beta", "gamma")

# Jags 

jags_posterior = jags(data = jags_data,
                      inits = init_values,
                      parameters.to.save = params,
                      model.file = logistic_regr_jags,
                      n.chains = 1,
                      n.iter = 5000,
                      n.burnin = 500,
                      n.thin = 1)

out = jags_posterior$BUGSoutput

## Extract samples from the posterior of beta and gamma

beta_post  = out$sims.list$beta
gamma_post = out$sims.list$gamma
S = nrow(gamma_post)
prob_inclusion = colMeans(gamma_post)

names(prob_inclusion) = c("Intercept",
                          "Debt/Equity",
                          "Curr_Ratio",
                          "EPS",
                          "REV/BAS_sh",
                          "Ni/Profit",
                          "Cf/Ni",
                          "A/R_Days",
                          "Px_to_Bas_EPS",
                          "Dvd12Myld_Gross",
                          "ROA",
                          "ROIC")


par(mfrow = c(1,1), mar = c(5,5,2,2)) 
barplot(prob_inclusion, col = "steelblue", horiz = TRUE, xlab = expression(hat(p)[j]), las = 1) 


datatest<-"testDataset.xlsx"
datatest<- read_excel(datatest, sheet=1)
datatest = datatest[complete.cases(datatest),]
y=datatest$Binary_Response  # Response Variable
X = as.data.frame(datatest[,c(4:14)])
X.x= model.matrix(y ~., X) 
p = ncol(X.x)
X = X.x

# Standardize the predictors:

for(j in 2:p){
  
  X[,j] = scale(X.x[,j])
}

n = nrow(X)

S = dim(beta_post)[1]
eta = matrix(0, nrow=n, ncol=S)

# Linear Predictor for each row. So we will have 110 rows and 4500 columns 
# Beta_post sarebbe l'output della SPIKE AND SLAP
for(i in 1:n){
  for(s in 1:S){
    eta[i,s] = (gamma_post[s,]*beta_post[s,])%*%X[i,]
  }
}

pi.star = list()
pi.star = exp(eta)/(1+exp(eta))

y_hat = round(rowMeans(pi.star))

library(reshape2)
library(caret)
library(pheatmap)
library(pROC)

y_actual <- factor(y)  
y_pred <- factor(y_hat)  
cm <- confusionMatrix(factor(y_pred), factor(y_actual), dnn = c("Prediction", "Reference"))

plt <- as.data.frame(cm$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
        geom_tile() + geom_text(aes(label=Freq)) +
        scale_fill_gradient(low="white", high="steelblue") +
        labs(x = "Reference",y = "Prediction") +
        scale_x_discrete(labels=c("1","0","1","0")) +
        scale_y_discrete(labels=c("0","1","1","0"))

y_actual <- factor(y) 
y_pred <- factor(y_hat)  

# Ensure that y_pred is a factor with the correct levels
y_pred <- factor(y_pred, levels = levels(y_actual))

# Calculate the ROC curve
roc_obj <- roc(y_actual, as.numeric(y_pred))

# Create a dataframe for ggplot
roc_data <- data.frame(
  cutpoints = roc_obj$thresholds,
  specificity = 1 - roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

# Plot the ROC 
ggplot(roc_data, aes(x = specificity, y = sensitivity)) +
  geom_line() +
  geom_abline(linetype = "dashed") +
  labs(
    title = "ROC Curve",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_minimal()

# Calculate the area under the ROC curve (AUC)
auc_value <- auc(roc_obj)
print(paste("The Area Under the Curve (AUC) is:", auc_value))





#########################
### Probit Regression ###
#########################

jags_data = with(data, list(y = y, X = X, n = length(y), p = ncol(X)))

probit_regr_jags = function(){
  
  # Likelihood:
  for(i in 1:n){
    y[i] ~ dbern(pi[i])
    probit(pi[i]) = (gamma*beta)%*%X[i,]
  }
  
  # Priors:
  for(j in 1:p){
    beta[j] ~ dnorm(0, 0.01)
  }
  
  for(j in 1:p){
    gamma[j] ~ dbern(w)
  }
  
  w ~ dbeta(1, 1)
  
}

# Set initial values for beta and tau
init_values = function(){
  
  list(beta = rep(0, p), gamma = rep(1, p))
  
}

params = c("beta","gamma") # parameters on which we are interested

# jags function
jags_posterior = jags(data = jags_data, inits = init_values, 
                      parameters.to.save = params, n.chains = 1, n.iter = 5000,
                      n.burnin = 1000, n.thin = 1, model.file = probit_regr_jags)


out = jags_posterior$BUGSoutput

## Extract samples from the posterior of beta and gamma

beta_post  = out$sims.list$beta
gamma_post = out$sims.list$gamma

S = nrow(gamma_post)

## Estimate the posterior probability of inclusion of each predictor Xj
## i.e. proportion of times gammaj = 1

prob_inclusion = colMeans(gamma_post)

names(prob_inclusion) = c("Intercept",
                          "Debt/Equity",
                          "Curr_Ratio",
                          "EPS",
                          "REV/BAS_sh",
                          "Ni/Profit",
                          "Cf/Ni",
                          "A/R_Days",
                          "Px_to_Bas_EPS",
                          "Dvd12Myld_Gross",
                          "ROA",
                          "ROIC")


par(mfrow = c(1,1), mar = c(5,5,2,2)) 
barplot(prob_inclusion, col = "steelblue", horiz = TRUE, xlab = expression(hat(p)[j]), las = 1)

# GEWEKE TEST 
# Per verificare che la chain abbia raggiunto stationarity. 
geweke_result <- geweke.diag(beta_post, frac1 = 0.3, frac2 = 0.4)
geweke_values <- geweke_result$z
geweke_matrix <- matrix(geweke_values, nrow = 1)
colnames(geweke_matrix) <- paste("Beta", 1:length(geweke_values))
kable(geweke_matrix, caption = "Geweke Test", align = 'c')


datatest = datatest[complete.cases(datatest),]
y=datatest$Binary_Response  # Response Variable
X = as.data.frame(datatest[,c(4:14)])
X.x= model.matrix(y ~., X) 
p = ncol(X.x)
X = X.x

# Standardize the predictors:

for(j in 2:p){
  
  X[,j] = scale(X.x[,j])
}

# Numero di righe nel  dataset
n = nrow(X)

# Numero di simulazioni
S = dim(beta_post)[1]

# Inizializza una matrice per i predittori lineari
eta = matrix(0, nrow=n, ncol=S)

for(i in 1:n){
  for(s in 1:S){
    eta[i,s] = (gamma_post[s,]*beta_post[s,])%*%X[i,]
  }
}

#  funzione media e la distribuzione predittiva
pi.star = list()
pi.star = exp(eta)/(1+exp(eta))

y_hat = round(rowMeans(pi.star))

cm <- confusionMatrix(factor(y_pred), factor(y_actual), dnn = c("Prediction", "Reference"))

plt <- as.data.frame(cm$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="steelblue") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=c("1","0","1","0")) +
  scale_y_discrete(labels=c("0","1","1","0"))

y_pred <- factor(y_pred, levels = levels(y_actual))

roc_obj <- roc(y_actual, as.numeric(y_pred))

roc_data <- data.frame(
  cutpoints = roc_obj$thresholds,
  specificity = 1 - roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

ggplot(roc_data, aes(x = specificity, y = sensitivity)) +
  geom_line() +
  geom_abline(linetype = "dashed") +
  labs(
    title = "ROC Curve",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_minimal()

# Calcolo dell'area sotto la ROC curve (AUC)
auc_value <- auc(roc_obj)
print(paste("L'Area Under the Curve (AUC) è:", auc_value))
