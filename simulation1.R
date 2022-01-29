# Ina Klasner
# 1/26/2022
# Simulation Study on PCR with ME Version 1

library(matlib)
library(pls)

#function for reading in filenames
file <- function(gamma, rep) {
  #Data_rho0.75_uGAMMA_e0.5_repREP.csv
  name <- paste0('Data_rho0.75_u', gamma, '_e0.5_rep', rep, '.csv')
  return(name)
}

directory <- '~/Desktop/michigan_tech/research/biostats/plsda/pcr/Simulation/Data/Rho075/Sige05/'

##### Test on one #####
data<- read.csv('Data_rho0.75_u1_e0.5_rep1.csv')

# initialize matrices 
Y <- data$Y
W <- as.matrix(data[,2:6]) #variable w/ME
X <- as.matrix(data[,7:11]) #true variable
cov_W <- t(W)%*%W/(nrow(data)-1)
B_true <- data$true.beta

# Side Test
# scale_x <- apply(X,2,scale)
# cov_X <- t(scale_x)%*%scale_x/(nrow(data)-1)
# svd(cov_X)

# Run PCR with Y and W
# get singular value decomp of W
W_svd <- svd(cov_W)
V <- W_svd$v
U <- W_svd$u

# create the diagonal matrix of singular values
D <- matrix(diag(W_svd$d), ncol=length(W_svd$d))

#eigenvalue decomp
W_ed <- eigen(cov_W)

# B = V.(D^-1).(U^T).Y
B <- V%*%solve(D, tol = 1e-17)%*%t(U)%*%Y
# t <- U%*%D # orthogonal scores
# P <- V
# Z <- U%*%D
# A <- solve(D, tol = 1e-17)%*%t(U)%*%Y
# fitted <- Z%*%A

# Check MSE of estimated parameters with B_true
mean(t(B_true[1:length(B)]-B)%*%(B_true[1:length(B)]-B))

# compare to pcr coefficients from package
pcr_mod <-pcr(Y ~ (W.1+W.2+W.3+W.4+W.5), data=data, scale=TRUE, validation='CV')
summary(pcr_mod)
pcr_mod$fitted.values

##### repeat on the rest #####
i <- c(1,5)
perf <- array(0, c(100, 2))
for(g in 1:2) { #gamma = 0.1, 0.5
  for(r in 1:100) { #rep
    dir <- paste0(directory, 'ME0', i[g]) #add folder
    setwd(dir)
    if (g==1) {
      data <- read.csv(file(i[g], r)) #read data file
    } else {
      data <- read.csv(file(0.5, r)) #read data file
    }
    
    # initialize matrices 
    Y <- matrix(data$Y,ncol=1)
    W <- as.matrix(data[,2:6]) #variable w/ME
    X <- as.matrix(data[,7:11]) #true variable
    B_true <- data$true.beta
    cov_W <- t(W)%*%W/(nrow(data)-1)
    
    # Run PCR with Y and W
    W_svd <- svd(W) # get singular value decomp of W
    V <- W_svd$v
    U <- W_svd$u
    D <- matrix(diag(W_svd$d), ncol=length(W_svd$d)) # create the diagonal matrix of singular values
    B <- V%*%solve(D)%*%t(U)%*%Y # B = V.(D^-1).(U^T).Y
    
    # Check MSE of estimated parameters with B_true
    perf[r,g] <- t(B_true[1:length(B)]-B)%*%(B_true[1:length(B)]-B)
  }
}

##### Review of Performance #####
mean(perf[,1]) #ME 1
mean(perf[,2]) #ME 0.5
