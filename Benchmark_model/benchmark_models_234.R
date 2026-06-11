########################## DATA setting ########################################
source("repeated_data_management.R")

#Inspect the data
head(mydata, 1)              #training data
head(data.valid,1)           #test data
table(mydata$col.Cov.Idx)
# n:frequency, s:aggregate loss, Type: one-hot encoding, col.Cov.Idx: 0,1,2
# id: PolicyNum
mean(mydata$PolicyNum %in% data.valid$PolicyNum) # not all training data is in test data
mean(data.valid$PolicyNum %in% mydata$PolicyNum) # all test data is in training data

mean(mydata$TypeCity)
length( unique(mydata$PolicyNum) )
length( unique(data$PolicyNum) )

install.packages("writexl")
library(writexl)
write_xlsx(mydata, path = "mydata.xlsx")
write_xlsx(data.valid, path = "testdata.xlsx")


######################################################################
####### Benchmark model 2. Poisson GLM ###############################
######################################################################
packages <- c("readxl","dplyr","yardstick","tibble","forcats","broom","sandwich","lmtest")
to_install <- setdiff(packages, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install)
lapply(packages, library, character.only = TRUE)

train_path <- "mydata.xlsx"
test_path  <- "testdata.xlsx"

df_train <- readxl::read_excel(train_path)
df_test  <- readxl::read_excel(test_path)

vars_keep <-
c("PolicyNum","n","TypeCity","TypeCounty","TypeMisc","TypeSchool","TypeTown","col.Cov.Idx")
df_train <- df_train |> dplyr::select(all_of(vars_keep)) |> dplyr::filter(!is.na(n))
df_test  <- df_test  |> dplyr::select(all_of(vars_keep)) |> dplyr::filter(!is.na(n))

df_train$.__split__ <- "train"
df_test$.__split__  <- "test"
df_all <- dplyr::bind_rows(df_train, df_test)

factor_cols <-
c("TypeCity","TypeCounty","TypeMisc","TypeSchool","TypeTown","col.Cov.Idx","PolicyNum")
df_all <- df_all |>
    mutate(across(all_of(factor_cols), ~ as.factor(.x)))

df_train <- df_all |> filter(.__split__ == "train") |> select(-.__split__)
df_test  <- df_all |> filter(.__split__ == "test")  |> select(-.__split__)

form_glm <- n ~ TypeCity + TypeCounty + TypeMisc + TypeSchool + TypeTown + factor(col.Cov.Idx)

fit_pois <- glm(
    formula = form_glm,
    data    = df_train,
    family  = poisson(link = "log")
)

summary(fit_pois)

phi_hat <- sum(residuals(fit_pois, type = "pearson")^2) / df.residual(fit_pois)
vcov_cl <- sandwich::vcovCL(fit_pois, cluster = df_train$PolicyNum)

pred_train <- predict(fit_pois, type = "response")
pred_test  <- predict(fit_pois, newdata = df_test, type = "response")

eval_train <- tibble(truth = df_train$n, estimate = pred_train)
eval_test  <- tibble(truth = df_test$n,  estimate = pred_test)

train_rmse <- yardstick::rmse(eval_train, truth = truth, estimate = estimate)$.estimate
train_mse  <- train_rmse^2
train_mae  <- yardstick::mae (eval_train, truth = truth, estimate = estimate)$.estimate

test_rmse  <- yardstick::rmse(eval_test,  truth = truth, estimate = estimate)$.estimate
test_mse   <- test_rmse^2
test_mae   <- yardstick::mae (eval_test,  truth = truth, estimate = estimate)$.estimate

cat(sprintf("[TRAIN] MSE=%.4f, MAE=%.4f\n", train_mse, train_mae))
cat(sprintf("[ TEST] MSE=%.4f, MAE=%.4f\n",  test_mse,  test_mae))

# > cat(sprintf("[TRAIN] MSE=%.4f, MAE=%.4f\n", train_mse, train_mae))
# [TRAIN] MSE=2.9777, MAE=0.8690
# > cat(sprintf("[ TEST] MSE=%.4f, MAE=%.4f\n", test_mse, test_mae))
# [ TEST] MSE=2.0241, MAE=0.8575

## ===== Test log-likelihood (Poisson GLM plug-in) =====
eps <- .Machine$double.eps
mu_test <- pmax(pred_test, eps)

test_loglik <- sum(dpois(df_test$n, lambda = mu_test, log = TRUE))
mean_test_loglik <- test_loglik / nrow(df_test)

cat(sprintf("[ TEST] logLik=%.4f, mean logLik=%.4f\n", test_loglik, mean_test_loglik))

# > cat(sprintf("[ TEST] logLik=%.4f, mean logLik=%.4f\n", test_loglik, mean_test_loglik))
# [ TEST] logLik=-266.5104, mean logLik=-1.0534


######################################################################
####### Benchmark model 3. Poisson Hurdle ############################
######################################################################
packages <- c("readxl","dplyr","yardstick","tibble","glmmTMB")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

train_path <- "mydata.xlsx"
test_path  <- "testdata.xlsx"

df_train <- readxl::read_excel(train_path)
df_test  <- readxl::read_excel(test_path)

vars_keep <- c("PolicyNum","n",
               "TypeCity","TypeCounty","TypeMisc","TypeSchool","TypeTown",
               "col.Cov.Idx")

df_train <- df_train |>
  dplyr::select(all_of(vars_keep)) |>
  dplyr::filter(!is.na(n))

df_test <- df_test |>
  dplyr::select(all_of(vars_keep)) |>
  dplyr::filter(!is.na(n))

to01 <- function(x) {
  if (is.numeric(x)) return(as.numeric(x > 0))
  if (is.logical(x)) return(as.numeric(x))
  z <- trimws(as.character(x))
  z <- tolower(z)
  out <- rep(NA_real_, length(z))
  out[z %in% c("1","y","yes","t","true")] <- 1
  out[z %in% c("0","n","no","f","false")] <- 0
  is_num <- suppressWarnings(!is.na(as.numeric(z)))
  out[is_num & is.na(out)] <- as.numeric(as.numeric(z[is_num]) > 0)
  out[is.na(out)] <- 0
  out
}

bin_cols <- c("TypeCity","TypeCounty","TypeMisc","TypeSchool","TypeTown")
for (nm in bin_cols) {
  if (nm %in% names(df_train)) df_train[[nm]] <- to01(df_train[[nm]])
  if (nm %in% names(df_test))  df_test[[nm]]  <- to01(df_test[[nm]])
}

df_train$col.Cov.Idx <- suppressWarnings(as.numeric(df_train$col.Cov.Idx))
df_test$col.Cov.Idx  <- suppressWarnings(as.numeric(df_test$col.Cov.Idx))

need_cols <- c("n", bin_cols, "col.Cov.Idx")
df_train <- tidyr::drop_na(df_train, dplyr::any_of(need_cols))
df_test  <- tidyr::drop_na(df_test,  dplyr::any_of(need_cols))

form_count <- n ~ TypeCity + TypeCounty + TypeMisc + TypeSchool + TypeTown + col.Cov.Idx

form_zi <- ~ 1

fit_hurdle <- glmmTMB::glmmTMB(
  formula   = form_count,
  ziformula = form_zi,
  family    = truncated_poisson(),
  data      = df_train
  # , offset = log(exposure)
)

cat("\n=== Hurdle(Poisson) model summary ===\n")
print(summary(fit_hurdle))

pred_train <- predict(fit_hurdle, type = "response")
pred_test  <- predict(fit_hurdle, newdata = df_test, type = "response")

eval_train <- tibble::tibble(truth = df_train$n, estimate = pred_train)
eval_test  <- tibble::tibble(truth = df_test$n, estimate = pred_test)

train_rmse <- yardstick::rmse(eval_train, truth = truth, estimate = estimate)$.estimate
train_mse  <- train_rmse^2
train_mae  <- yardstick::mae (eval_train, truth = truth, estimate = estimate)$.estimate

test_rmse  <- yardstick::rmse(eval_test,  truth = truth, estimate = estimate)$.estimate
test_mse   <- test_rmse^2
test_mae   <- yardstick::mae (eval_test,  truth = truth, estimate = estimate)$.estimate

cat(sprintf("\n[TRAIN] MSE=%.4f, MAE=%.4f\n", train_mse, train_mae))
cat(sprintf("[ TEST] MSE=%.4f, MAE=%.4f\n",  test_mse,  test_mae))

# > cat(sprintf("\n[TRAIN] MSE=%.4f, MAE=%.4f\n", train_mse, train_mae))
#
# [TRAIN] MSE=3.7118, MAE=0.9641
# > cat(sprintf("[ TEST] MSE=%.4f, MAE=%.4f\n", test_mse, test_mae))
# [ TEST] MSE=2.4239, MAE=0.9158


## ===== Test log-likelihood for Poisson hurdle (glmmTMB) =====
eps <- .Machine$double.eps

y <- df_test$n

## 1) eta -> mu = exp(eta)
eta <- predict(fit_hurdle, newdata = df_test, type = "link")
mu <- exp(eta)

## 2) hurdle : pi0 = P(Y=0)  (ziformula)
pi0 <- predict(fit_hurdle, newdata = df_test, type = "zprob")
pi0 <- pmin(pmax(pi0, eps), 1 - eps)

## 3) truncated Poisson: 1 - P(Pois(mu)=0) = 1 - exp(-mu)
denom <- -expm1(-mu)                 # = 1 - exp(-mu)
denom <- pmax(denom, eps)

logp <- numeric(length(y))

is0 <- (y == 0)
logp[is0] <- log(pi0[is0])

## y>0:  P(Y=y) = (1-pi0) * Pois(y;mu) / (1-exp(-mu))
logp[!is0] <- log1p(-pi0[!is0]) +
  dpois(y[!is0], lambda = mu[!is0], log = TRUE) -
  log(denom[!is0])

test_loglik <- sum(logp)
mean_test_loglik <- mean(logp)

cat(sprintf("[ TEST] logLik=%.4f, mean logLik=%.4f\n", test_loglik, mean_test_loglik))

# > cat(sprintf("[ TEST] logLik=%.4f, mean logLik=%.4f\n", test_loglik, mean_test_loglik))
# [ TEST] logLik=-301.3224, mean logLik=-1.1910


######################################################################
####### Benchmark model 4. Poisson Zero-inflated #####################
######################################################################
packages <- c("readxl","dplyr","yardstick","tibble","forcats","broom.mixed","glmmTMB")
to_install <- setdiff(packages, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install)
lapply(packages, library, character.only = TRUE)

train_path <- "mydata.xlsx"
test_path  <- "testdata.xlsx"

df_train <- readxl::read_excel(train_path)
df_test  <- readxl::read_excel(test_path)

vars_keep <-
c("PolicyNum","n","TypeCity","TypeCounty","TypeMisc","TypeSchool","TypeTown","col.Cov.Idx")
df_train <- df_train |> dplyr::select(all_of(vars_keep)) |> dplyr::filter(!is.na(n))
df_test  <- df_test  |> dplyr::select(all_of(vars_keep)) |> dplyr::filter(!is.na(n))

df_train$.__split__ <- "train"
df_test$.__split__  <- "test"
df_all <- dplyr::bind_rows(df_train, df_test)

df_train <- df_all |>
  dplyr::filter(.__split__ == "train") |>
  dplyr::select(-all_of(".__split__"))

df_test <- df_all |>
  dplyr::filter(.__split__ == "test") |>
  dplyr::select(-all_of(".__split__"))

factor_cols <-
c("TypeCity","TypeCounty","TypeMisc","TypeSchool","TypeTown","col.Cov.Idx","PolicyNum")
df_all <- df_all |> mutate(across(all_of(factor_cols), ~ as.factor(.x)))

df_train <- df_all |>
  dplyr::filter(.__split__ == "train") |>
  dplyr::select(-dplyr::all_of(".__split__"))

df_test <- df_all |>
  dplyr::filter(.__split__ == "test") |>
  dplyr::select(-dplyr::all_of(".__split__"))

form_count <- n ~ TypeCity + TypeCounty + TypeMisc + TypeSchool + TypeTown + factor(col.Cov.Idx)
form_zi    <- ~ 1

fit_zip <- glmmTMB::glmmTMB(
  formula   = form_count,      # count part
  ziformula = form_zi,         # zero-inflation part
  family    = poisson(link = "log"),
  data      = df_train
  # , offset = log(exposure)
)

summary(fit_zip)

pred_train <- predict(fit_zip, type = "response")
pred_test  <- predict(fit_zip, newdata = df_test, type = "response")

eval_train <- tibble(truth = df_train$n, estimate = pred_train)
eval_test  <- tibble(truth = df_test$n,  estimate = pred_test)

train_rmse <- yardstick::rmse(eval_train, truth = truth, estimate = estimate)$.estimate
train_mse  <- train_rmse^2
train_mae  <- yardstick::mae (eval_train, truth = truth, estimate = estimate)$.estimate

test_rmse  <- yardstick::rmse(eval_test,  truth = truth, estimate = estimate)$.estimate
test_mse   <- test_rmse^2
test_mae   <- yardstick::mae (eval_test,  truth = truth, estimate = estimate)$.estimate

cat(sprintf("[TRAIN] MSE=%.4f, MAE=%.4f\n", train_mse, train_mae))
cat(sprintf("[ TEST] MSE=%.4f, MAE=%.4f\n",  test_mse,  test_mae))

# > cat(sprintf("[TRAIN] MSE=%.4f, MAE=%.4f\n", train_mse, train_mae))
# [TRAIN] MSE=3.0002, MAE=0.8659
# > cat(sprintf("[ TEST] MSE=%.4f, MAE=%.4f\n", test_mse, test_mae))
# [ TEST] MSE=1.9762, MAE=0.8396

## ===== Test log-likelihood for Poisson ZIP (glmmTMB) =====
eps <- .Machine$double.eps
y <- df_test$n

## 1) count part: eta -> mu
eta <- predict(fit_zip, newdata = df_test, type = "link")
mu <- exp(eta)

## 2) zero-inflation probability pi = P(structural zero)
pi <- predict(fit_zip, newdata = df_test, type = "zprob")
pi <- pmin(pmax(pi, eps), 1 - eps)

logp <- numeric(length(y))
is0 <- (y == 0)

## y=0: log( pi + (1-pi)*exp(-mu) )
a <- log(pi[is0])
b <- log1p(-pi[is0]) - mu[is0]            # log((1-pi)*exp(-mu))
m <- pmax(a, b)
logp[is0] <- m + log(exp(a - m) + exp(b - m))

## y>0: log(1-pi) + log Pois(y;mu)
logp[!is0] <- log1p(-pi[!is0]) + dpois(y[!is0], lambda = mu[!is0], log = TRUE)

test_loglik <- sum(logp)
mean_test_loglik <- mean(logp)

cat(sprintf("[ TEST] logLik=%.4f, mean logLik=%.4f\n", test_loglik, mean_test_loglik))

# > cat(sprintf("[ TEST] logLik=%.4f, mean logLik=%.4f\n", test_loglik, mean_test_loglik))
# [ TEST] logLik=-256.4331, mean logLik=-1.0136
