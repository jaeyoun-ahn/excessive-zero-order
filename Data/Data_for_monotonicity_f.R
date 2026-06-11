load("data.RData")
head(data)

dim(data)
length( unique(data$PolicyNum) )
data$col.freq<-data$FreqCN+data$FreqCO         #Collision old and new
data$col.Cov<-data$CoverageCN+data$CoverageCO  #Coverage old and new

index<-(data$CoverageCN>0 & data$CoverageCO>0)
data.train<-data[index,]
dim(data.train)
length(unique(data.train$PolicyNum))

mycov<-sort(data.train$col.Cov)
length(mycov)

#data.train$col.Cov.Idx<- as.factor((data.train$col.Cov>0.635)  +1)
data.train$col.Cov.Idx<- as.factor((data.train$col.Cov>mycov[489])+(data.train$col.Cov>mycov[978])  )
unique(data.train$col.Cov.Idx)

#Define frequency and severity of collision of old and new
data.train$n <-data.train$col.freq
data.train$s <- data.train$yAvgCN * data.train$FreqCN   +  data.train$yAvgCO * data.train$FreqCO

#Define frequency and severity of collision of old and new
index2<-(data.train$col.freq>0)
data.train.sev<-data.train[index2,]
data.train.sev$m<- (data.train.sev$yAvgCN * data.train.sev$FreqCN   +  data.train.sev$yAvgCO * data.train.sev$FreqCO)/(data.train.sev$FreqCN+data.train.sev$FreqCO)

load("dataout.RData")
head(dataout)

# Check whether the data is the set
# x<-c(1,2,3,4,5)
# y<-c(2,4,6,2)
# # match(y, x)
# y %in% x

# Out of test data, use index which is from the training data & coverage>0
out.idx1 <- dataout$PolicyNum %in% data.train$PolicyNum
out.idx2 <- dataout$CoverageCO>0 & dataout$CoverageCN>0
out.idx <- out.idx1 & out.idx2
length( unique(data.train$PolicyNum) )
sum(out.idx)
data.valid<-dataout[out.idx,] #test data

# Now, define explanatory variable in test data. 
# (I am not sure the explanatory variable will be ever used or not.)
out.idx.sev<- data.valid$PolicyNum %in% data.train$PolicyNum
data.valid$col.freq<-data.valid$FreqCN+data.valid$FreqCO
data.valid$col.Cov<-data.valid$CoverageCN+data.valid$CoverageCO
data.valid$col.Cov.Idx<-as.factor((data.valid$col.Cov>mycov[489])+(data.valid$col.Cov>mycov[978])  )
#data.valid$col.Cov.Idx<- as.factor((data.valid$col.Cov>0.635)  +1)
data.valid$s <- data.valid$yAvgCN * data.valid$FreqCN   +  data.valid$yAvgCO * data.valid$FreqCO
data.valid$n <- data.valid$FreqCN   +   data.valid$FreqCO

###For the Research
typeof(data.train$s)
head(data.train, 1)
idx<-data.train$Year==2010
sum(data.train$s[idx] )/sum(data.train[idx,]$n )
sum(data.valid$s)/sum(data.valid$n)

#Inspect the data
head(data.train, 1)           #training data 
head(data.valid,1)        #test data
table(data.train$col.Cov.Idx) 
# n:frequency, s:aggregate loss, Type: one-hot encoding, col.Cov.Idx: 0,1,2
# id: PolicyNum
mean(data.train$PolicyNum %in%  data.valid$PolicyNum) # not all training data is in test data
mean(data.valid$PolicyNum %in% data.train$PolicyNum) # all test data is in training data

unique_policy = unique(data.train$PolicyNum)
n_pol = length(unique_policy)
Y_train = matrix(NA, nrow=n_pol, ncol=5)
ID_train = rep(NA, n_pol)
X_train = matrix(NA, nrow=n_pol, ncol=7)
Length_observed_train = rep(NA, n_pol)

for(i in 1:n_pol){
  temp_idx = which(unique_policy[i]== data.train$PolicyNum  )
  year_idx = data.train$Year[temp_idx]-2005
  len_temp = length(temp_idx)
  for(t in 1:5){
    Y_train[i, year_idx] = c(data.train$n[temp_idx])
    #Y_train[i, 1:5] = c(data.train$n[temp_idx], rep(NA, 5-len_temp))
  }
  X_train[i, 1:6] = with(c(TypeCity[temp_idx[1]], TypeCounty[temp_idx[1]], TypeSchool[temp_idx[1]], TypeTown[temp_idx[1]], 
                           TypeVillage[temp_idx[1]], col.Cov.Idx[temp_idx[1]]), data=data.train)
  ID_train[i] = unique_policy[i]
  Length_observed_train[i] = len_temp
}

X_train <- as.data.frame(X_train, stringsAsFactors = FALSE)
colnames(X_train) <- c(
  "CityType",
  "CountyType",
  "SchoolType",
  "TownType",
  "VillageType",
  "OtherCov"
)
XX_train = model.matrix(~ CityType+CountyType+SchoolType+TownType+VillageType+factor(OtherCov), data=X_train)
head(XX_train)


n_test = length(data.valid$PolicyNum)
X_test = matrix(NA, nrow=n_pol, ncol=6)
Y_test = rep(NA, nrow=n_test)
ID_test = rep(NA, nrow=n_test)
i=100
which(data.valid$PolicyNum[i]== unique_policy)

for(i in 1:n_test){
  temp_idx = which(data.valid$PolicyNum[i]== unique_policy  )
  Y_test[i] = data.valid$n[i]
  ID_test[i] = unique_policy[temp_idx]
  X_test[i, 1:6] = with(c(TypeCity[i], TypeCounty[i], TypeSchool[i], TypeTown[i], 
                          TypeVillage[i], col.Cov.Idx[i]), data=data.valid)
}

X_test <- as.data.frame(X_test, stringsAsFactors = FALSE)
colnames(X_test) <- c(
  "CityType",
  "CountyType",
  "SchoolType",
  "TownType",
  "VillageType",
  "OtherCov"
)
XX_test = model.matrix(~ CityType+CountyType+SchoolType+TownType+VillageType+factor(OtherCov), data=X_test)
head(XX_test, 1)

# Input DATA
head(Y_train)
head(XX_train)
IDX = !is.na(Y_train)
num_obs= rowSums(!is.na(Y_train))
n = dim(Y_train)[1]
p = dim(XX_train)[2]
k = dim(XX_train)[2]

betaMean  <- rep(0, p)        
betaCov   <- diag(1, p)
dim(XX_train)
IDXint = matrix(NA, nrow=n, ncol=5)
for(i in 1:n){
  obsCols    <- which(IDX[i,])
  IDXint[i, 1:num_obs[i]] = obsCols
}
length(betaMean)
dim(betaCov)

#Data for hurdle model(Y:=Z(1+N))
num_obs1 = num_obs
IDXint1 = IDXint

I_train = (Y_train>0)
N_train = Y_train-1
N_train[N_train==-1] <- NA

IDX2 = !is.na(N_train) 
IDXint2 = matrix(NA, nrow=n, ncol=5)

IDX3 = which(rowSums(!is.na(N_train))>0)
length(IDX3)

num_obs2= rowSums(!is.na(N_train))
for(i in IDX3){
  # print(i)
  obsCols2    <- which(IDX2[i,])
  IDXint2[i, 1:num_obs2[i]] = obsCols2
}



###############################
#### Final data to be used ####
###############################

## Data for classical zero-inflated model
I_mat  <- I_train
N_mat  <- N_train
wI_mat <- !is.na(I_mat) * 1L     # 1 = 관측, 0 = NA
wN_mat <- !is.na(N_mat) * 1L
I_mat[is.na(I_mat)] <- 0L        # NA → 0 (더미)
N_mat[is.na(N_mat)] <- 0L
T_I <- ncol(I_mat);  T_N <- ncol(N_mat)

## Data for state-space copula model
# n   <- nrow(Y_train)
# k   <- ncol(XX_train)



#############
# Train data: 
#############
Y_mat <- Y_train
wY_mat <- !is.na(Y_mat) * 1L
Y_mat[!wY_mat] <- 0L
ID_train

# Y_mat
# n \times 5 matrix of counting observations

# XX_train
# n \times k matrix of explanatory variable including the intercept.
# Assume there are k-1 explanatory variables and one intercept so that it has k columns
# Assume that explanatory variable does not change with time so that it is not a 3 dim array.

# wY_mat
# n \times 5 matrix of 0,1 exposures

# ID_train
# length n vector of ID to be linked with test data


#############
# Test data:  
#############
Y_test
XX_test
ID_test

# Y_test
# n_test vector of counting observations

# XX_test
# n_test \times k matrix of explanatory variable including the intercept.
# Assume there are k-1 explanatory variables and one intercept so that it has k columns
# Assume that explanatory variable does not change with time so that it is not a 3 dim array.

# ID_test
# length n_test vector of ID to be linked with train data

# Y_test: n_test vector of counts (0,1,2,...)
I_test <- as.integer(Y_test > 0)
N_test <- Y_test - 1L
N_test[Y_test == 0] <- NA_integer_   # if 0, then NA

# matrix form (n_test x 1)
I_test_mat <- matrix(I_test, ncol = 1)
N_test_mat <- matrix(N_test, ncol = 1)

# exposure matrix
wI_test_mat <- matrix(1L, nrow = nrow(I_test_mat), ncol = ncol(I_test_mat))
wN_test_mat <- matrix(as.integer(!is.na(N_test_mat)), nrow = nrow(N_test_mat), ncol = ncol(N_test_mat))

# Nimble data(NA -> 0)
I_test_mat_nim <- I_test_mat
N_test_mat_nim <- N_test_mat
I_test_mat_nim[is.na(I_test_mat_nim)] <- 0L
N_test_mat_nim[is.na(N_test_mat_nim)] <- 0L


########################################################################################################


#Data for comparison : (y_{1:4}, Y_5=0), (y_{1:4}, Y_5=1)
num_obs1 = num_obs
IDXint1 = IDXint

Y_train_c_0 = Y_train
Y_train_c_0[,5] <- 0
Y_train_c_0

Y_train_c_1 = Y_train
Y_train_c_1[,5] <- 1 
Y_train_c_1

I_train_c_0 = (Y_train_c_0>0)
N_train_c_0 = Y_train_c_0-1
N_train_c_0[N_train_c_0==-1] <- NA

I_train_c_1 = (Y_train_c_1>0)
N_train_c_1 = Y_train_c_1-1
N_train_c_1[N_train_c_1==-1] <- NA

IDX2_c_0 = !is.na(N_train_c_0) 
IDXint2_c_0 = matrix(NA, nrow=n, ncol=5)
IDX3_c_0 = which(rowSums(!is.na(N_train_c_0))>0)
length(IDX3_c_0)
num_obs2_c_0= rowSums(!is.na(N_train_c_0))
for(i in IDX3_c_0){
  # print(i)
  obsCols2_c_0    <- which(IDX2_c_0[i,])
  IDXint2_c_0[i, 1:num_obs2_c_0[i]] = obsCols2_c_0
}

IDX2_c_1 = !is.na(N_train_c_1) 
IDXint2_c_1 = matrix(NA, nrow=n, ncol=5)
IDX3_c_1 = which(rowSums(!is.na(N_train_c_1))>0)
length(IDX3_c_1)
num_obs2_c_1= rowSums(!is.na(N_train_c_1))
for(i in IDX3_c_1){
  # print(i)
  obsCols2_c_1    <- which(IDX2_c_1[i,])
  IDXint2_c_1[i, 1:num_obs2_c_1[i]] = obsCols2_c_1
}

## Final Data for comparison of (E[Y_6 | Y_5=0, y_{1:4}] vs E[Y_6 | Y_5=1, y_{1:4}])
I_mat_c_0  <- I_train_c_0
N_mat_c_0  <- N_train_c_0
wI_mat_c_0 <- !is.na(I_mat_c_0) * 1L     # 1 = 관측, 0 = NA
wN_mat_c_0 <- !is.na(N_mat_c_0) * 1L
I_mat_c_0[is.na(I_mat_c_0)] <- 0L        # NA → 0 (더미)
N_mat_c_0[is.na(N_mat_c_0)] <- 0L
T_I_c_0 <- ncol(I_mat_c_0);  T_N_c_0 <- ncol(N_mat_c_0)

I_mat_c_1  <- I_train_c_1
N_mat_c_1  <- N_train_c_1
wI_mat_c_1 <- !is.na(I_mat_c_1) * 1L     # 1 = 관측, 0 = NA
wN_mat_c_1 <- !is.na(N_mat_c_1) * 1L
I_mat_c_1[is.na(I_mat_c_1)] <- 0L        # NA → 0 (더미)
N_mat_c_1[is.na(N_mat_c_1)] <- 0L
T_I_c_1 <- ncol(I_mat_c_1);  T_N_c_1 <- ncol(N_mat_c_1)

Y_mat_c_0 <- Y_train_c_0
wY_mat_c_0 <- !is.na(Y_mat_c_0) * 1L
Y_mat_c_0[!wY_mat_c_0] <- 0L

Y_mat_c_1 <- Y_train_c_1
wY_mat_c_1 <- !is.na(Y_mat_c_1) * 1L
Y_mat_c_1[!wY_mat_c_1] <- 0L



















