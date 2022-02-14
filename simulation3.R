# Ina Klasner
# 2/4/2022
# Simulation Study on PCR with ME Version 3
# New rho, variance of random error, and measurement error scale
# update: using 1 principle component

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
        B_true <- data$true.beta
        
        # Run PCR with Y and W
        W_svd <- svd(W) # get singular value decomp of W
        V <- W_svd$v
        U <- W_svd$u
        D <- matrix(diag(W_svd$d), ncol=length(W_svd$d)) # create the diagonal matrix of singular values
        B <- matrix(V[,1],ncol=1)%*%solve(D[1,1])%*%t(matrix(U[,1],ncol=1))%*%Y # B = V.(D^-1).(U^T).Y
        
        # Check MSE of estimated parameters with B_true
        perf[r,g,s,p] <- t(B_true[1:length(B)]-B)%*%(B_true[1:length(B)]-B)
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
