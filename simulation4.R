# Ina Klasner
# 2/10/2022
# Simulation Study on PCR with ME Version 4
# with estimator of beta

##### Setup #####
#function for reading in filenames
file <- function(gamma, rep, rho, sige) {
  #Data_rho0.RHO_uGAMMA_eSIGE_repREP.csv
  name <- paste0('Data_rho0.', rho, '_u', gamma, '_e', sige,'_rep', rep, '.csv')
  return(name)
}

directory <- '~/Desktop/michigan_tech/research/biostats/plsda/pcr/Simulation/Data/'
i <- c('05','1')
i_file <- c('0.5','1')
j <- c(75, 85, 99)
perf <- array(0, c(100, 2, 2, 3))
newperf<- array(0, c(100, 2, 2, 3))

##### Test on one #####
data<- read.csv('Data_rho0.75_u1_e0.5_rep1.csv')

# initialize matrices 
Y <- data$Y
W <- as.matrix(data[,2:6]) #variable w/ME
X <- as.matrix(data[,7:11]) #true variable
cov_W <- t(W)%*%W/(nrow(data)-1) #same as cov(W)
cov_WY <- t(W)%*%Y/(nrow(data)-1) #xy covariance
B_true <- data$true.beta

# Run PCR with Y and W
# get singular value decomp of W
W_svd <- svd(cov_W)
V <- W_svd$v
U <- W_svd$u

# create the diagonal matrix of singular values
S <- matrix(diag(W_svd$d), ncol=length(W_svd$d))

K <- solve(cov_W)%*%(cov_W-nrow(data)*S)# find K to calibrate the estimation of B

# Bhat = (Sxx-nSIGMA)^-1 Sxy
B_hat <- solve(cov_W-nrow(data)*S)%*%cov_WY

# B0 = y-B'X
B_hat0 <- Y-t(t(B_hat)%*%t(W)) #B_hat or conjugate transpose of B_hat ?
T_scores <- W%*%V

# B*=T(T'SxxT)^-1 T'Sxy
B_hat22 <- T_scores%*%solve(T_scores%*%cov_W%*%t(T_scores))%*%t(T_scores)%*%cov_WY

# B = V.(D^-1).(U^T).Y
B <- matrix(V[,1],ncol=1)%*%solve(S[1,1])%*%t(matrix(U[,1],ncol=1))%*%Y

# Check MSE of estimated parameters with B_true
mean(t(B_true[1:length(B)]-B)%*%(B_true[1:length(B)]-B))


##### Iterate through files and Compute B #####
for(p in 1:3) { #rho = 0.75, 0.85, 0.99
  for(s in 1:2) { #sige = 1, 0.5
    for(g in 1:2) { #gamma = 1, 0.5
      for(r in 1:100) { #rep
        if(r==31 && p==2 && s==2 && g==2) {
          break
        }
        dir <- paste0(directory, 'Rho0', j[p], '/Sige', i[s],'/ME', i[g]) #add folder
        setwd(dir)
        data <- read.csv(file(i_file[g], r, j[p], i_file[s])) 
        
        # initialize matrices 
        Y <- matrix(data$Y,ncol=1)
        W <- as.matrix(data[,2:6]) #variable w/ME
        X <- as.matrix(data[,7:11]) #true variable
        cov_W <- t(W)%*%W/(nrow(data)-1) #same as cov(W)
        cov_WY <- t(W)%*%Y/(nrow(data)-1) #xy covariance
        B_true <- data$true.beta
        
        # Run PCR with Y and W
        W_svd <- svd(W) # get singular value decomp of W
        V <- W_svd$v
        U <- W_svd$u
        S <- matrix(diag(W_svd$d), ncol=length(W_svd$d)) # create the diagonal matrix of singular values
        
        B_hat <- solve(cov_W-nrow(data)*S)%*%cov_WY # Bhat = (Sxx-nSIGMA)^-1 Sxy
        B_hat0 <- Y-t(t(B_hat)%*%t(W)) # B0 = y-B'X
        B <- matrix(V[,1],ncol=1)%*%solve(D[1,1])%*%t(matrix(U[,1],ncol=1))%*%Y # B = V.(D^-1).(U^T).Y
        
        # y=b0+b'x+e

        # Check MSE of estimated parameters with B_true
        perf[r,g,s,p] <- t(B_true[1:length(B)]-B)%*%(B_true[1:length(B)]-B)
        newperf[r,g,s,p] <- t(B_true[1:length(B_hat)]-B_hat)%*%(B_true[1:length(B_hat)]-B_hat)
      }
    }
  }
}

##### Analyze Performance #####
tab1 <- matrix(c(mean(perf[,1,1,1]), mean(perf[,1,1,2]), mean(perf[,1,1,3]),
              mean(perf[,2,1,1]), mean(perf[,2,1,2]), mean(perf[,2,1,3])), ncol=6)
tab <- rbind(tab1, c(mean(perf[,1,2,1]), mean(perf[,1,2,2]), mean(perf[,1,2,3]),
                     mean(perf[,2,2,1]), mean(perf[,2,2,2]), mean(perf[,2,2,3])))
colnames(tab) <- c('p=0.75', 'p=0.85', 'p=0.99', 'p=0.75', 'p=0.85', 'p=0.99')
rownames(tab) <- c('Sige 0.5', 'Sige 1')
tab2 <- matrix(c(mean(newperf[,1,1,1]), mean(newperf[,1,1,2]), mean(newperf[,1,1,3]),
                 mean(newperf[,2,1,1]), mean(newperf[,2,1,2]), mean(newperf[,2,1,3])), ncol=6)
tab_hat <- rbind(tab2, c(mean(newperf[,1,2,1]), mean(newperf[,1,2,2]), mean(newperf[,1,2,3]),
                     mean(newperf[,2,2,1]), mean(newperf[,2,2,2]), mean(newperf[,2,2,3])))
colnames(tab) <- c('p=0.75', 'p=0.85', 'p=0.99', 'p=0.75', 'p=0.85', 'p=0.99')
rownames(tab) <- c('Sige 0.5', 'Sige 1')

par(mfrow=c(1,2))
x <- c(0.75, 0.85, 0.99)
mse_05_05 <- tab[1,1:3]
mse_05_1 <- tab[1,4:6]
mse_1_05 <- tab[2,1:3]
mse_1_1 <- tab[2,4:6]

plot(x, mse_05_05,
     main='Sige 0.5 MSE of Beta',
     xlab='Rho', ylab='MSE',
     xlim=c(0.75, 1), ylim=c(0,0.01),
     col='red')
points(x, mse_05_1, col='blue')
legend("topleft",
       inset=0.02,
       legend = c("ME = 0.5", "ME = 1"),
       fill = c('red', 'blue'),
       cex=0.7)

plot(x, mse_1_05,
     main='Sige 1 MSE of Beta',
     xlab='Rho', ylab='MSE',
     xlim=c(0.75, 1), ylim=c(0,0.02),
     col='red')
points(x, mse_1_1, col='blue')
