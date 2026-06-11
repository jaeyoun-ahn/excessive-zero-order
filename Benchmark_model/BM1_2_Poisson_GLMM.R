########################################################
####### BM1_2(Poisson GLMM) ############################
########################################################

source("Data_for_monotonicity_f.R")

# 0. set constants and some functions
nn = dim(Y_mat)[1] # number of policyholders (train data) = 409
tt = dim(Y_mat)[2] # number of t (train data) = 5
kk = dim(XX_train)[2] # number of exponential variables = 8
zero_mean=rep(0, kk) # for prior(mean)
diag_cov= diag(rep(1.0, kk)) # for prior(covariance)

# 1. nimble code for Model6
library(nimble)
data_list <- list(
  Y = Y_mat #Y ~ Pois(lambda[i]*exp(R[i])
)
constants <- list(
  nn = nn, tt = tt,
  kk = kk,
  zero_mean=zero_mean,
  diag_cov=diag_cov,
  X=XX_train,
  wY_mat=wY_mat # exposure
)

set.seed(101)
code <- nimbleCode({

  for (i in 1:nn) {
    R[i] ~ dnorm(mean=(-1/2)*(d^2), sd=d) # R ~ N((-1/2)*(d^2), d^2)
    lambda[i] <- exp(inprod(beta[1:kk], X[i, 1:kk]))  # lambda_i=exp(<beta,X>)
  }
  beta[1:kk]~ dmnorm(mean = zero_mean[], cov = diag_cov[,])

  log_d ~ dnorm(0, sd = 1.5)
  d     <- exp(log_d)

  for (i in 1:nn) {
    for (t in 1:tt) {
      Y[i, t] ~ dpois( wY_mat[i,t] * (lambda[i]*(exp(R[i])) ) )
    }
  }
})

# 2. MCMC
samples_bm1_2 <- nimbleMCMC(
  code      = code,
  constants = constants,
  data      = data_list,
  monitors  = c("beta","d","R"),
  niter     = 30000, nburnin = 10000, thin = 10, nchains = 3
)

# 3. Results of MCMC
library(MCMCvis)
summary_mcmc_bm1_2 = MCMCsummary(object = samples_bm1_2, round = 2, params = c("beta","d"))
MCMCtrace(object = samples_bm1_2,
          pdf = FALSE, # no export to PDF
          ind = TRUE,
          params = c("beta","d"))

# 4. estimators of coefficients, thetas
beta_hat = summary_mcmc_bm1_2[1:kk,1]
d_hat = summary_mcmc_bm1_2[(kk+1),1]
R_sim_train <- rbind(samples_bm1_2$chain1, samples_bm1_2$chain2, samples_bm1_2$chain3)[,1:nn]

# 5. Test validation (test MSE, MAE)
R_sim_test= R_sim_train[, ID_train %in% ID_test] # R for test set
XX_test = XX_train[ID_train %in% ID_test, ] # X for test set
log_lamda_hat_test = XX_test %*% beta_hat

Y_test_hat = rep(NA, dim(XX_test)[1])
for(i in 1:dim(XX_test)[1]){
  Y_test_hat[i] = mean(exp(log_lamda_hat_test[i])*exp(R_sim_test[,i]))
}

MSE <- mean((Y_test - Y_test_hat)^2)
cat("MSE =", round(MSE, 4), "\n") #MSE = 1.1171

MAE <- mean(abs(Y_test - Y_test_hat))
cat("MAE =", round(MAE, 4), "\n") #MAE = 0.6063


# 6. Test posterior predictive log-likelihood (BM1_2)
mcmc_bm1_2 <- rbind(samples_bm1_2$chain1, samples_bm1_2$chain2, samples_bm1_2$chain3)

beta_cols <- grep("^beta\\[", colnames(mcmc_bm1_2))
R_cols    <- grep("^R\\[",    colnames(mcmc_bm1_2))

stopifnot(length(beta_cols) == kk)
stopifnot(length(R_cols) == nn)

beta_sim <- mcmc_bm1_2[, beta_cols, drop = FALSE]      # Sdraw x kk
R_sim_train <- mcmc_bm1_2[, R_cols, drop = FALSE]      # Sdraw x nn
R_sim_test  <- R_sim_train[, ID_train %in% ID_test, drop = FALSE]

log_mean_exp <- function(logv){
  m <- max(logv)
  m + log(mean(exp(logv - m)))
}

lambda_draws <- exp(beta_sim %*% t(XX_test))      # Sdraw x ntest

eps <- .Machine$double.eps
ntest <- length(Y_test)
lpd_i <- numeric(ntest)

for(i in 1:ntest){
  mu_s <- lambda_draws[, i] * exp(R_sim_test[, i])
  mu_s <- pmax(mu_s, eps)
  lpd_i[i] <- log_mean_exp(dpois(Y_test[i], lambda = mu_s, log = TRUE))
}

test_log_pred <- sum(lpd_i)
cat("Test log predictive likelihood (sum LPD) =", round(test_log_pred, 4), "\n")
cat("Mean test LPD per policy =", round(mean(lpd_i), 4), "\n")

# > cat("Test log predictive likelihood (sum LPD) =", round(test_log_pred, 4), "\n")
# Test log predictive likelihood (sum LPD) = -208.6458
# > cat("Mean test LPD per policy =", round(mean(lpd_i), 4), "\n")
# Mean test LPD per policy = -0.8221
