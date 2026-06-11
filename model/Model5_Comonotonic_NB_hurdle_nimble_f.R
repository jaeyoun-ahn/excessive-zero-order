##########################################################
#### Model5(Comonotonic Negative Binomial hurdle model) ###
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
kk = dim(XX_train)[2] # number of explanatory variables = 8
zero_mean=rep(0, kk) # for prior of coefficients (mean)
diag_cov= diag(rep(3.0, kk)) # for prior of coefficients (covariance)
eps_nb <- 1.0E-12 # small number for numerical stability of NB probability

sigmoid <-function(t){
  1/(1+exp(-t))
}
softplus <- function(t){
  log(1+exp(t))
}


# 1. nimble code for Model5
library(nimble)
data_list <- list(
  I  = I_mat,
  N  = N_mat 
)
constants <- list(
  nn = nn, tt = tt,
  kk = kk,
  zero_mean=zero_mean,
  diag_cov=diag_cov,
  X=XX_train,
  wI_mat=wI_mat, wN_mat=wN_mat,
  eps_nb=eps_nb
)
set.seed(101)
code <- nimbleCode({

  for (i in 1:nn) {
    theta[i] ~ dnorm(0, sd = sigma_theta) 
  }
  sigma_theta ~ dgamma(shape=0.3, rate=0.3) 
  r_N ~ dgamma(shape=0.3, rate=0.3)         

  for (i in 1:nn) {
    c_i[i] <- inprod(beta_c[1:kk], X[i, 1:kk]) # c_i=<beta_c,X>
    d_i[i] <- inprod(beta_d[1:kk], X[i, 1:kk]) # d_i=<beta_d,X>

    for (t in 1:tt) {
      I[i, t] ~ dbern( wI_mat[i,t] * (1/(1+exp(-(c_i[i]+theta[i])))) )
    }

    for (t in 1:tt) {
      eta_N[i, t] <- d_i[i] + theta[i]/r_N
      q_N_raw[i, t] <- 1/(1+exp(eta_N[i, t]))
      q_N_obs[i, t] <- eps_nb + (1 - 2*eps_nb) * q_N_raw[i, t]

      prob_N[i, t] <- (1-wN_mat[i,t]) + wN_mat[i,t] * q_N_obs[i, t]
      N[i, t] ~ dnegbin(prob = prob_N[i, t], size = r_N)
    }
  }
  beta_c[1:kk] ~ dmnorm(mean = zero_mean[1:kk], cov = diag_cov[1:kk,1:kk]) 
  beta_d[1:kk] ~ dmnorm(mean = zero_mean[1:kk], cov = diag_cov[1:kk,1:kk])
})


# 2. MCMC
inits_mod5 <- function(){
  list(
    beta_c = rep(0, kk),
    beta_d = rep(0, kk),
    theta = rep(0, nn),
    sigma_theta = 1,
    r_N = 1
  )
}

samples_mod5 <- nimbleMCMC(
  code      = code,
  constants = constants,
  data      = data_list,
  inits     = inits_mod5,
  monitors  = c("beta_c","beta_d","theta","sigma_theta","r_N"),
  niter     = 30000, nburnin = 10000, thin = 10, nchains = 3
)


# 3. Results of MCMC
library(MCMCvis)
summary_mcmc_mod5 = MCMCsummary(object = samples_mod5, round = 2, params = c("beta_c","beta_d","sigma_theta","r_N"))
MCMCtrace(object = samples_mod5,
          pdf = FALSE, # no export to PDF
          ind = TRUE,
          params = c("beta_c","beta_d","sigma_theta","r_N"))


# 4. estimators of coefficients, thetas, and r
samples_mat_mod5 <- rbind(samples_mod5$chain1, samples_mod5$chain2, samples_mod5$chain3)

beta_c_hat = summary_mcmc_mod5[paste0("beta_c[", 1:kk, "]"), 1]
beta_d_hat = summary_mcmc_mod5[paste0("beta_d[", 1:kk, "]"), 1]
sigma_theta_hat = summary_mcmc_mod5["sigma_theta", 1]
r_N_hat = summary_mcmc_mod5["r_N", 1]

theta_cols <- paste0("theta[", 1:nn, "]")
theta_sim_train <- samples_mat_mod5[, theta_cols, drop = FALSE]
r_N_sim <- samples_mat_mod5[, "r_N"]


# 5. Test validation (test MSE, MAE)
theta_sim_test= theta_sim_train[, ID_train %in% ID_test, drop = FALSE] # thetas for test set
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
  p_I_s <- sigmoid(as.numeric(chat_test[i]) + theta_sim_test[,i])

  eta_N_s <- as.numeric(dhat_test[i]) + theta_sim_test[,i]/r_N_sim
  q_N_s <- sigmoid(-eta_N_s) # R/NIMBLE prob = 1 - p_N
  q_N_s <- pmin(pmax(q_N_s, eps_nb), 1 - eps_nb)
  mean_N_s <- r_N_sim * (1 - q_N_s) / q_N_s

  Y_test_hat[i] = mean(p_I_s * (1 + mean_N_s))
}

MSE <- mean((Y_test - Y_test_hat)^2)
cat("MSE =", round(MSE, 4), "\n")
# MSE = 0.9624 
MAE <- mean(abs(Y_test - Y_test_hat))
cat("MAE =", round(MAE, 4), "\n")
# MAE = 0.606

# 6. Test posterior predictive log-likelihood (Model5)

log_mean_exp <- function(logv){
  m <- max(logv)
  m + log(mean(exp(logv - m)))
}

Sdraw <- nrow(theta_sim_test)
ntest <- ncol(theta_sim_test)

# test hurdle (ntest x 1)
I_test <- as.integer(I_test_mat[,1])
N_test <- as.integer(N_test_mat[,1])

stopifnot(length(I_test) == ntest)
stopifnot(length(r_N_sim) == Sdraw)

eps <- eps_nb
lpd_i <- numeric(ntest)

for(i in 1:ntest){
 
  p_s <- sigmoid(as.numeric(chat_test[i]) + theta_sim_test[, i])
  p_s <- pmin(pmax(p_s, eps), 1 - eps)

  r_s <- pmax(r_N_sim, eps)
  eta_N_s <- as.numeric(dhat_test[i]) + theta_sim_test[, i]/r_s
  prob_N_s <- sigmoid(-eta_N_s)
  prob_N_s <- pmin(pmax(prob_N_s, eps), 1 - eps)

  loglik_I_s <- dbinom(I_test[i], size = 1, prob = p_s, log = TRUE)

  if(I_test[i] == 1L){
    loglik_N_s <- dnbinom(N_test[i], size = r_s, prob = prob_N_s, log = TRUE)
  } else {
    loglik_N_s <- rep(0, Sdraw)
  }

  lpd_i[i] <- log_mean_exp(loglik_I_s + loglik_N_s)
}

test_log_pred <- sum(lpd_i)
cat("Test log predictive likelihood (sum LPD) =", round(test_log_pred, 4), "\n")
# Test log predictive likelihood (sum LPD) = -214.7405 
cat("Mean test LPD per policy =", round(mean(lpd_i), 4), "\n")
# Mean test LPD per policy = -0.8488 
