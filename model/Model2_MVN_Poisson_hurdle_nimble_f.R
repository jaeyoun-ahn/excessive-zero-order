##########################################################
####### Model2(MVN Poisson hurdle model) #################
##########################################################

source("Data_for_monotonicity_f.R")

# 0. set constants and some functions
nn = dim(Y_mat)[1] # number of policyholders (train data) = 409
tt = dim(Y_mat)[2] # number of t (train data) = 5
kk = dim(XX_train)[2] # number of exponential variables = 8
zero_mean=rep(0, kk) # for prior(mean)
diag_cov= diag(rep(1.0, kk)) # for prior(covariance)
sigmoid <-function(t){
  1/(1+exp(-t))
}

# 1. nimble code for Model2
library(nimble)
data_list <- list(
  I  = I_mat, #I=Z ~ Ber(sigmoid(theta[i,1]))
  N  = N_mat  #N ~ Pois(exp(theta[i,2]))
)
constants <- list(
  nn = nn, tt = tt,
  kk = kk, 
  zero_mean=zero_mean, 
  diag_cov=diag_cov, 
  X=XX_train, 
  wI_mat=wI_mat, wN_mat=wN_mat # exposure
)

set.seed(101)
code <- nimbleCode({
  
  for (i in 1:nn) {
    theta[i,1:2] ~ dmnorm(mean=mu[i,1:2], cov = Sigma[1:2,1:2]) # theta ~ MVN((mu1,mu2)^T,Sigma)
    mu[i,1] <- inprod(beta1[1:kk], X[i, 1:kk])  # mu1_i=<beta1,X>
    mu[i,2] <- inprod(beta2[1:kk], X[i, 1:kk])  # mu2_i=<beta2,X>
  }
  
  beta1[1:kk] ~ dmnorm(mean = zero_mean[], cov = diag_cov[,])
  beta2[1:kk] ~ dmnorm(mean = zero_mean[], cov = diag_cov[,])
  
  sigma1 ~ dgamma(1,1)
  sigma2 ~ dgamma(1,1)
  
  rho ~ dunif(-0.99, 0.99)
  Sigma[1,1] <- sigma1^2
  Sigma[2,2] <- sigma2^2
  Sigma[1,2] <- rho* sigma1*sigma2
  Sigma[2,1] <- rho* sigma1*sigma2
  
  for (i in 1:nn) {
    for (t in 1:tt) {
      I[i, t] ~ dbern( wI_mat[i,t] * (1/(1+exp(-theta[i,1]))) )
    }
    
    for (t in 1:tt) {
      N[i, t] ~ dpois( wN_mat[i,t] * exp(theta[i,2]) )
    }
  } 
})


# 2. MCMC
samples_mod2 <- nimbleMCMC(
  code      = code,
  constants = constants,
  data      = data_list,
  monitors  = c("beta1","beta2","sigma1","sigma2","rho","theta"),
  niter     = 30000, nburnin = 10000, thin = 10, nchains = 3
)

# 3. Results of MCMC
library(MCMCvis)
summary_mcmc_mod2 = MCMCsummary(object = samples_mod2, round = 2, params = c("beta1","beta2","sigma1","sigma2","rho"))
MCMCtrace(object = samples_mod2,
          pdf = FALSE, # no export to PDF
          ind = TRUE,
          params = c("beta1","beta2","sigma1","sigma2","rho"))


# 4. estimators of coefficients, thetas
beta1_hat = summary_mcmc_mod2[1:kk,1]
beta2_hat = summary_mcmc_mod2[(kk+1):(kk+kk),1]
sigma1_hat = summary_mcmc_mod2[(kk+kk+1),1]
sigma2_hat = summary_mcmc_mod2[(kk+kk+2),1]
rho_hat = summary_mcmc_mod2[(kk+kk+3),1]
kkk = dim(summary_mcmc_mod2)[1]
theta1_sim_train <-rbind(samples_mod2$chain1, samples_mod2$chain2, samples_mod2$chain3)[,(kkk+1):(kkk+nn)]
theta2_sim_train <-rbind(samples_mod2$chain1, samples_mod2$chain2, samples_mod2$chain3)[,(kkk+nn+1):(kkk+nn+nn)]

# 5. Test validation (test MSE, MAE)
theta1_sim_test= theta1_sim_train[, ID_train %in% ID_test] # theta1 for test set
theta2_sim_test= theta2_sim_train[, ID_train %in% ID_test] # theta2 for test set
# XX_test = XX_train[ID_train %in% ID_test, ] # X for test set
# mu1hat_test = XX_test %*% beta1_hat
# mu2hat_test = XX_test %*% beta2_hat

Y_test_hat = rep(NA, dim(XX_test)[1])
for(i in 1:dim(XX_test)[1]){
  Y_test_hat[i] = mean(sigmoid(theta1_sim_test[,i])*(1 + exp(theta2_sim_test[,i])))
}

MSE <- mean((Y_test - Y_test_hat)^2)
cat("MSE =", round(MSE, 4), "\n") # MSE = 1.0489 

MAE <- mean(abs(Y_test - Y_test_hat))
cat("MAE =", round(MAE, 4), "\n") # MAE = 0.599

# 5.1. Test posterior predictive log-likelihood (LPD)

log_mean_exp <- function(logv){
  m <- max(logv)
  m + log(mean(exp(logv - m))) # log( mean(exp(logv)) )
}

Sdraw <- nrow(theta1_sim_test)     # posterior draws
ntest <- ncol(theta1_sim_test)     # test

## I_test_mat, N_test_mat
I_test <- as.integer(I_test_mat[,1])
N_test <- N_test_mat[,1]           

stopifnot(length(I_test) == ntest)

eps <- .Machine$double.eps
lpd_i <- numeric(ntest)

for(i in 1:ntest){
  ## (p, lambda)
  p_s <- sigmoid(theta1_sim_test[, i])
  p_s <- pmin(pmax(p_s, eps), 1 - eps)         # 0/1
  lambda_s <- exp(theta2_sim_test[, i])
  
  ## Bernoulli part: I_i | theta1
  loglik_I_s <- dbinom(I_test[i], size = 1, prob = p_s, log = TRUE)
  
  ## Poisson part: N_i | theta2 
  if(I_test[i] == 1L){
    loglik_N_s <- dpois(N_test[i], lambda = lambda_s, log = TRUE)
  } else {
    loglik_N_s <- rep(0, Sdraw) 
  }
  
  ## p(test_i | train) ≈ mean_s p(test_i | theta^(s))
  lpd_i[i] <- log_mean_exp(loglik_I_s + loglik_N_s)
}

test_log_pred <- sum(lpd_i)

cat("Test log predictive likelihood (sum LPD) =", round(test_log_pred, 4), "\n")
cat("Mean test LPD per policy =", round(mean(lpd_i), 4), "\n")

#Test log predictive likelihood (sum LPD) = -208.9472 
#Mean test LPD per policy = -0.8259 

####################################################################################################
# 6. Comparison of (E[Y6 | Y5=0, y_{1:4}] vs E[Y6 | Y5=1, y_{1:4}])
####################################################################################################
# 6.1. simulate E[Y6 | Y5=0, y_{1:4}]
# 6.1.1. nimble code for comparison Model2
data_list_c_0 <- list(
  I  = I_mat_c_0, 
  N  = N_mat_c_0 
)
constants_c_0 <- list(
  nn = nn, tt = tt,
  kk = kk, 
  X=XX_train, 
  wI_mat=wI_mat_c_0, wN_mat=wN_mat_c_0, # exposure
  beta1=beta1_hat, beta2=beta2_hat,
  sigma1 = sigma1_hat, sigma2 = sigma2_hat, rho = rho_hat
)

set.seed(101)
code_comparison <- nimbleCode({
  
  for (i in 1:nn) {
    theta[i,1:2] ~ dmnorm(mean=mu[i,1:2], cov = Sigma[1:2,1:2]) # theta ~ MVN((mu1,mu2)^T,Sigma)
    mu[i,1] <- inprod(beta1[1:kk], X[i, 1:kk])  # mu1_i=<beta1,X>
    mu[i,2] <- inprod(beta2[1:kk], X[i, 1:kk])  # mu2_i=<beta2,X>
  }
  Sigma[1,1] <- sigma1^2
  Sigma[2,2] <- sigma2^2
  Sigma[1,2] <- rho* sigma1*sigma2
  Sigma[2,1] <- rho* sigma1*sigma2
  
  for (i in 1:nn) {
    for (t in 1:tt) {
      I[i, t] ~ dbern( wI_mat[i,t] * (1/(1+exp(-theta[i,1]))) )
    }
    
    for (t in 1:tt) {
      N[i, t] ~ dpois( wN_mat[i,t] * exp(theta[i,2]) )
    }
  } 
})

# 6.1.2. MCMC
samples_mod2_c_0 <- nimbleMCMC(
  code      = code_comparison,
  constants = constants_c_0,
  data      = data_list_c_0,
  monitors  = c("theta"),
  niter     = 50000, nburnin = 10000, thin = 10, nchains = 3
)
theta1_sim_train_c_0 <-rbind(samples_mod2_c_0$chain1, samples_mod2_c_0$chain2, samples_mod2_c_0$chain3)[,1:nn]
theta2_sim_train_c_0 <-rbind(samples_mod2_c_0$chain1, samples_mod2_c_0$chain2, samples_mod2_c_0$chain3)[,(nn+1):(nn+nn)]

# 6.1.3. estimation of E[Y6 | Y5=0, y_{1:4}]
Y_train_hat_c_0 = rep(NA, dim(XX_train)[1])
for(i in 1:dim(XX_train)[1]){
  Y_train_hat_c_0[i] = mean(sigmoid(theta1_sim_train_c_0[,i])*(1 + exp(theta2_sim_train_c_0[,i])))
}


# 6.2. simulate E[Y6 | Y5=1, y_{1:4}]
# 6.2.1. nimble code for comparison Model2
data_list_c_1 <- list(
  I  = I_mat_c_1, 
  N  = N_mat_c_1  
)
constants_c_1 <- list(
  nn = nn, tt = tt,
  kk = kk, 
  X=XX_train, 
  wI_mat=wI_mat_c_1, wN_mat=wN_mat_c_1, # exposure
  beta1=beta1_hat, beta2=beta2_hat,
  sigma1 = sigma1_hat, sigma2 = sigma2_hat, rho = rho_hat
)
set.seed(101)

# 6.2.2. MCMC
samples_mod2_c_1 <- nimbleMCMC(
  code      = code_comparison,
  constants = constants_c_1,
  data      = data_list_c_1,
  monitors  = c("theta"),
  niter     = 50000, nburnin = 10000, thin = 10, nchains = 3
)
theta1_sim_train_c_1 <-rbind(samples_mod2_c_1$chain1, samples_mod2_c_1$chain2, samples_mod2_c_1$chain3)[,1:nn]
theta2_sim_train_c_1 <-rbind(samples_mod2_c_1$chain1, samples_mod2_c_1$chain2, samples_mod2_c_1$chain3)[,(nn+1):(nn+nn)]

# 6.2.3. estimation of E[Y_6 | Y_5=0, y_{1:4}]
Y_train_hat_c_1 = rep(NA, dim(XX_train)[1])
for(i in 1:dim(XX_train)[1]){
  Y_train_hat_c_1[i] = mean(sigmoid(theta1_sim_train_c_1[,i])*(1 + exp(theta2_sim_train_c_1[,i])))
}

# 6.3. comparison of E[Y6 | Y5=0, y_{1:4}] vs E[Y6 | Y5=1, y_{1:4}]
dff_Y_train_hat_c = Y_train_hat_c_1 - Y_train_hat_c_0 
neg_ratio <- mean(dff_Y_train_hat_c < 0)
neg_ratio

summary(dff_Y_train_hat_c)
min(dff_Y_train_hat_c)
sort(dff_Y_train_hat_c)[1:10]

# > neg_ratio
# [1] 0.09290954
# > summary(dff_Y_train_hat_c)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -2.56286  0.09618  0.12822  0.07346  0.14351  0.26776 
# > min(dff_Y_train_hat_c)
# [1] -2.562855
# > sort(dff_Y_train_hat_c)[1:10]
# [1] -2.5628554 -1.9770556 -1.4781614 -1.3604918 -1.3007370 -1.0180066 -0.8984751 -0.7245999 -0.6790919
# [10] -0.5801572



##########################################################################################
# 7. Comparison of (E[(Y6-d)^+ | Y5=0, y_{1:4}] vs E[(Y6-d)^+ | Y5=1, y_{1:4}]), d=1, 2
##########################################################################################
sigmoid <- function(x) 1 / (1 + exp(-x))

Ds    <- c(1, 2)  # d
n_sim <- 10          # simulation number

compute_EY_excess_mc <- function(theta1_mat, theta2_mat, d, n_sim) {
  M      <- nrow(theta1_mat)
  n_obs  <- ncol(theta1_mat)
  EY_vec <- numeric(n_obs)
  
  for (i in 1:n_obs) {
    temp <- numeric(M)
    for (j in 1:M) {
      th1 <- theta1_mat[j, i]
      th2 <- theta2_mat[j, i]
      
      p   <- sigmoid(th1)
      lam <- exp(th2)
      
      Z6 <- rbinom(n_sim, size = 1, prob = p)
      N6 <- rpois(n_sim, lambda = lam)
      Y6 <- Z6 * (1 + N6)
      
      temp[j] <- mean(pmax(Y6 - d, 0))
    }
    EY_vec[i] <- mean(temp)  
  }
  
  EY_vec
}

## 7.1. E[(Y6 - d)^+|Y5=0]
EY_excess_c_0 <- sapply(Ds, function(d) {
  compute_EY_excess_mc(
    theta1_mat = theta1_sim_train_c_0,
    theta2_mat = theta2_sim_train_c_0,
    d          = d,
    n_sim      = n_sim
  )
})
colnames(EY_excess_c_0) <- paste0("d", Ds)

## 7.2. E[(Y6 - d)^+|Y5=1]
EY_excess_c_1 <- sapply(Ds, function(d) {
  compute_EY_excess_mc(
    theta1_mat = theta1_sim_train_c_1,
    theta2_mat = theta2_sim_train_c_1,
    d          = d,
    n_sim      = n_sim
  )
})
colnames(EY_excess_c_1) <- paste0("d", Ds)

## 7.3. E[(Y6-d)^+ | Y5=1] vs E[(Y6-d)^+ | Y5=0]
for (k in seq_along(Ds)) {
  d <- Ds[k]
  diff_d <- EY_excess_c_1[, k] - EY_excess_c_0[, k]
  neg_ratio_d <- mean(diff_d < 0)
  
  cat("\n===== d =", d, "=====\n")
  cat("neg_ratio =", neg_ratio_d, "\n")
  print(summary(diff_d))
  cat("min diff =", min(diff_d), "\n")
  cat("10 smallest diffs:\n")
  print(sort(diff_d)[1:10])
}


# ===== d = 1 =====
#   neg_ratio = 0.1075795 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -2.56456  0.01962  0.03040 -0.01687  0.04292  0.11130 
# min diff = -2.564558 
# 10 smallest diffs:
#   [1] -2.5645583 -1.9888000 -1.4800667 -1.3678083 -1.3139000 -1.0414083 -0.9298500 -0.7330750 -0.6974417
# [10] -0.5962083
# 
# ===== d = 2 =====
#   neg_ratio = 0.1271394 
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -2.546050  0.002600  0.005042 -0.042137  0.009008  0.026942 
# min diff = -2.54605 
# 10 smallest diffs:
#   [1] -2.5460500 -1.9860333 -1.4752500 -1.3929417 -1.3331750 -0.9769333 -0.9345083 -0.7489333 -0.6575333
# [10] -0.5944167


####################################################################################################
# 8. Comparison of (E[min(y6,d) | Y5=0, y_{1:4}] vs E[min(y6,d) | Y5=1, y_{1:4}]), d=1, 2
####################################################################################################
sigmoid <- function(x) 1 / (1 + exp(-x))

M       <- nrow(theta1_sim_train_c_0)
n_train <- ncol(theta1_sim_train_c_0)   

Ds    <- c(1, 2)   # d
n_sim <- 10           # simulation number

compute_Emin_mc <- function(theta1_mat, theta2_mat, d, n_sim = 10) {
  M      <- nrow(theta1_mat)
  n_obs  <- ncol(theta1_mat)
  EY_vec <- numeric(n_obs)
  
  for (i in 1:n_obs) {
    temp <- numeric(M)
    for (j in 1:M) {
      th1 <- theta1_mat[j, i]
      th2 <- theta2_mat[j, i]
      
      p   <- sigmoid(th1)
      lam <- exp(th2)
      
      Z6 <- rbinom(n_sim, size = 1, prob = p)
      N6 <- rpois(n_sim, lambda = lam)
      Y6 <- Z6 * (1 + N6)
      
      temp[j] <- mean(pmin(Y6, d))
    }
    EY_vec[i] <- mean(temp)  
  }
  
  EY_vec
}

## 8.1. E[min(Y6, d)|Y5=0]
EY_min_c_0 <- sapply(Ds, function(d) {
  compute_Emin_mc(
    theta1_mat = theta1_sim_train_c_0,
    theta2_mat = theta2_sim_train_c_0,
    d          = d,
    n_sim      = n_sim
  )
})
colnames(EY_min_c_0) <- paste0("d", Ds)

## 8.2. E[min(Y6, d)|Y5=1]
EY_min_c_1 <- sapply(Ds, function(d) {
  compute_Emin_mc(
    theta1_mat = theta1_sim_train_c_1,
    theta2_mat = theta2_sim_train_c_1,
    d          = d,
    n_sim      = n_sim
  )
})
colnames(EY_min_c_1) <- paste0("d", Ds)

## 8.3. E[min(Y6,d) | Y5=1] vs E[min(Y6,d) | Y5=0] 
for (k in seq_along(Ds)) {
  d <- Ds[k]
  diff_d <- EY_min_c_1[, k] - EY_min_c_0[, k]
  neg_ratio_d <- mean(diff_d < 0)
  
  cat("\n===== d =", d, "=====\n")
  cat("neg_ratio =", neg_ratio_d, "\n")
  print(summary(diff_d))
  cat("min diff =", min(diff_d), "\n")
  cat("10 smallest diffs:\n")
  print(sort(diff_d)[1:10])
}

# ===== d = 1 =====
#   neg_ratio = 0 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000925 0.074917 0.092425 0.090221 0.105783 0.177500 
# min diff = 0.000925 
# 10 smallest diffs:
#   [1] 0.000925000 0.007366667 0.010425000 0.010508333 0.011425000 0.012725000 0.014058333 0.015433333
# [9] 0.015783333 0.015841667
# 
# ===== d = 2 =====
#   neg_ratio = 0.009779951 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.03231  0.09304  0.12104  0.11577  0.13578  0.24818 
# min diff = -0.03230833 
# 10 smallest diffs:
#   [1] -0.032308333 -0.021616667 -0.021575000 -0.003333333  0.015508333  0.017683333  0.018391667  0.018858333
# [9]  0.019033333  0.020100000


