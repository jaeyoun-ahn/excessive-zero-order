##########################################################
####### Model4(Comonotonic Poisson hurdle model) #########
##########################################################

source("Data_for_monotonicity_f.R")

####### Data for classical zero-inflated model ###########
# I_mat  <- I_train
# N_mat  <- N_train
# wI_mat <- !is.na(I_mat) * 1L     # 1 = 관측, 0 = NA
# wN_mat <- !is.na(N_mat) * 1L
# I_mat[is.na(I_mat)] <- 0L        # NA → 0 (더미)
# N_mat[is.na(N_mat)] <- 0L
# T_I <- ncol(I_mat);  T_N <- ncol(N_mat)
# Y_mat
# ID_train
# XX_train
#########################################################

# 0. set constants and some functions
nn = dim(Y_mat)[1] # number of policyholders (train data) = 409
tt = dim(Y_mat)[2] # number of t (train data) = 5
kk = dim(XX_train)[2] # number of exponential variables = 8
zero_mean=rep(0, kk) # for prior of coefficients (mean)
diag_cov= diag(rep(3.0, kk)) # for prior of coefficients (covariance)
sigmoid <-function(t){
  1/(1+exp(-t))
}
softplus <- function(t){
  log(1+exp(t))
}


# 1. nimble code for Model4
library(nimble)
data_list <- list(
  I  = I_mat, #I=Z ~ Ber(sigmoid(c_i+theta_i))
  N  = N_mat  #N ~ Pois(lambda(d_i+theta_i)) where lambda(d_i+theta_i)=ln(1+exp(d_i+theta_i))
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
    theta[i] ~ dnorm(0, sd = sigma_theta) # theta ~ N(0,sigma_theta^2)
  }
  sigma_theta ~ dgamma(shape=0.3, rate=0.3) # prior of sigma_theta
  
  for (i in 1:nn) {
    c_i[i] <- inprod(beta_c[1:kk], X[i, 1:kk]) # c_i=<beta_c,X>
    d_i[i] <- inprod(beta_d[1:kk], X[i, 1:kk]) # d_i=<beta_d,X>
    
    for (t in 1:tt) {
      I[i, t] ~ dbern( wI_mat[i,t] * (1/(1+exp(-(c_i[i]+theta[i])))) )
    }
    
    for (t in 1:tt) {
      N[i, t] ~ dpois( wN_mat[i,t] * log(1+exp(d_i[i] + theta[i])) )
    }
  } 
  beta_c[1:kk] ~ dmnorm(mean = zero_mean[1:kk], cov = diag_cov[1:kk,1:kk]) # prior of coefficients(c_i)
  beta_d[1:kk] ~ dmnorm(mean = zero_mean[1:kk], cov = diag_cov[1:kk,1:kk]) # prior of coefficients(d_i)
})


# 2. MCMC
samples_mod4 <- nimbleMCMC(
  code      = code,
  constants = constants,
  data      = data_list,
  monitors  = c("beta_c","beta_d","theta","sigma_theta"),
  niter     = 30000, nburnin = 10000, thin = 10, nchains = 3
)


# 3. Results of MCMC
library(MCMCvis)
summary_mcmc_mod4 = MCMCsummary(object = samples_mod4, round = 2, params = c("beta_c","beta_d","sigma_theta"))
MCMCtrace(object = samples_mod4,
          pdf = FALSE, # no export to PDF
          ind = TRUE,
          params = c("beta_c","beta_d","sigma_theta"))


# 4. estimators of coefficients, thetas
beta_c_hat = summary_mcmc_mod4[1:kk,1]
beta_d_hat = summary_mcmc_mod4[(kk+1):(kk+kk),1]
sigma_theta_hat = summary_mcmc_mod4[(kk+kk+1),1]
kkk = dim(summary_mcmc_mod4)[1]
theta_sim_train <-rbind(samples_mod4$chain1, samples_mod4$chain2, samples_mod4$chain3)[,(kkk+1):(kkk+nn)]


# 5. Test validation (test MSE, MAE)
theta_sim_test= theta_sim_train[, ID_train %in% ID_test] # thetas for test set
XX_test = XX_train[ID_train %in% ID_test, ] # X for test set
chat_test = XX_test %*% beta_c_hat
dhat_test = XX_test %*% beta_d_hat

# mean(ID_train %in% ID_test)
# mean(ID_test %in% ID_train)
# Check whether the data is the set
# id_train<-c(1,2,3,4,5)
# id_test<-c(   2,  4,5)
# theta_train<-c(6,7,8,9,10)
# id_train %in% id_test # [1] FALSE  TRUE FALSE  TRUE  TRUE
# theta_test = theta_train[id_train %in% id_test] #[1]  7  9 10

Y_test_hat = rep(NA, dim(XX_test)[1])
for(i in 1:dim(XX_test)[1]){
  Y_test_hat[i] = mean(sigmoid(chat_test[i] + theta_sim_test[,i])*(1 + softplus(dhat_test[i] + theta_sim_test[,i])))
}

MSE <- mean((Y_test - Y_test_hat)^2)
cat("MSE =", round(MSE, 4), "\n") # MSE = 0.9382

MAE <- mean(abs(Y_test - Y_test_hat))
cat("MAE =", round(MAE, 4), "\n") # MAE = 0.5994

# 6. Test posterior predictive log-likelihood (Model4)

log_mean_exp <- function(logv){
  m <- max(logv)
  m + log(mean(exp(logv - m)))
}

Sdraw <- nrow(theta_sim_test)
ntest <- ncol(theta_sim_test)

# test hurdle (ntest x 1)
I_test <- as.integer(I_test_mat[,1])
N_test <- N_test_mat[,1]   

stopifnot(length(I_test) == ntest)

eps <- .Machine$double.eps
lpd_i <- numeric(ntest)

for(i in 1:ntest){
  # p, lambda
  p_s <- sigmoid(as.numeric(chat_test[i]) + theta_sim_test[, i])
  p_s <- pmin(pmax(p_s, eps), 1 - eps)
  
  lambda_s <- softplus(as.numeric(dhat_test[i]) + theta_sim_test[, i])
  lambda_s <- pmax(lambda_s, eps)
  
  # Bernoulli part: I
  loglik_I_s <- dbinom(I_test[i], size = 1, prob = p_s, log = TRUE)
  
  # Poisson part: N 
  if(I_test[i] == 1L){
    loglik_N_s <- dpois(N_test[i], lambda = lambda_s, log = TRUE)
  } else {
    loglik_N_s <- rep(0, Sdraw)
  }
  
  lpd_i[i] <- log_mean_exp(loglik_I_s + loglik_N_s)
}

test_log_pred <- sum(lpd_i)
cat("Test log predictive likelihood (sum LPD) =", round(test_log_pred, 4), "\n")
cat("Mean test LPD per policy =", round(mean(lpd_i), 4), "\n")

# > cat("Test log predictive likelihood (sum LPD) =", round(test_log_pred, 4), "\n")
# Test log predictive likelihood (sum LPD) = -208.5055 
# > cat("Mean test LPD per policy =", round(mean(lpd_i), 4), "\n")
# Mean test LPD per policy = -0.8241 




