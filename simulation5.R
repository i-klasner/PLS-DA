# Ina Klasner
# 2/16/2022
# Simulation Study on PCR with ME Version 5
# Beta estimated with Equation 9 and 22, Applying Kaiser rule

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
perf_n  <- array(0, c(100, 2, 2, 3))
perf_nrest <- array(0, c(100, 2, 2, 3))
perf_22 <- array(0, c(100, 2, 2, 3))
perf_n_alt <- array(0, c(100, 2, 2, 3))
perf_nrest_pc <- array(0, c(100, 2, 2, 3))

##### Iterate through files and Compute Betas #####
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
        
        ##### initialize matrices #####
        Y <- matrix(data$Y,ncol=1)
        W <- as.matrix(data[,2:6]) #variable w/ME
        X <- as.matrix(data[,7:11]) #true variable
        cov_W <- t(W)%*%W/(nrow(data)-1) #same as cov(W)
        cov_WY <- t(W)%*%Y/(nrow(data)-1) #xy covariance
        cov_dd <- t(W-X)%*%(W-X)/(nrow(data)-1) #covariance of measurement error, delta
        B_true <- data$true.beta
        
        ##### Singular Value Decomposition ####
        W_svd <- svd(W) # get singular value decomp of W
        V <- W_svd$v
        U <- W_svd$u
        S <- matrix(diag(W_svd$d), ncol=length(W_svd$d)) # create the diagonal matrix of singular values
        
        ###### Eigenvalues and Choosing Number of PC #####
        lambda <-eigen(cov_W)$values
        num_pc <- 1
        for(l in 1:length(lambda)) {
          if (lambda[l] < 1) {
            num_pc <- l-1
            break
          }
        }
        
        ##### Apply Equations 9, 22 to find Betas#####
        B_n <- solve(cov_W-nrow(data)*cov_dd)%*%cov_WY # Equation 9
        sigma_data <- c(0.254, 0.003, 0.002, -0.005, -0.001, 
                        0.003, 0.252, 0.005, 0.004, 0.001, 
                        0.002, 0.005, 0.252, 0.004, -0.001, 
                        -0.005, 0.004, 0.004, 0.247, 0,
                        -0.001, 0.001, -0.001, 0, 0.257) # variance-covariance matrix on page 13
        sigma <- matrix(sigma_data,nrow=5,ncol=5,byrow=TRUE)
        B_n_alt <- solve(cov_W-nrow(data)*sigma)%*%cov_WY
        #### WITH ALL PC 
        # B_22 <- V%*%solve(V%*%cov_W%*%t(V))%*%t(V)%*%cov_WY # Equation 22
        #### WITH 1 PC V(V'SxxV)^-1 V'Sxy
        V1 <- matrix(V[,1],ncol=1)
        B_22 <- V1%*%solve(t(V1)%*%cov_W%*%V1)%*%t(V1)%*%cov_WY # Equation 22
        B_nrest <- V1%*%solve(S[1,1])%*%t(matrix(U[,1],ncol=1))%*%Y # B = V.(S^-1).(U^T).Y
        B_nrest_pc <- matrix(V[,num_pc],ncol=num_pc)%*%solve(S[num_pc,num_pc])%*%t(matrix(U[,num_pc],ncol=num_pc))%*%Y 

        ##### Find MSE between estimated Betas and true Beta, add to matrices ####
        perf_n[r,g,s,p] <- t(B_true[1:length(B_n)]-B_n)%*%(B_true[1:length(B_n)]-B_n)
        perf_nrest[r,g,s,p] <- t(B_true[1:length(B_nrest)]-B_nrest)%*%(B_true[1:length(B_nrest)]-B_nrest)
        perf_22[r,g,s,p] <- t(B_true[1:length(B_22)]-B_22)%*%(B_true[1:length(B_22)]-B_22)
        perf_n_alt[r,g,s,p] <- t(B_true[1:length(B_n_alt)]-B_n_alt)%*%(B_true[1:length(B_n_alt)]-B_n_alt)
        perf_nrest_pc[r,g,s,p] <- t(B_true[1:length(B_nrest_pc)]-B_nrest_pc)%*%(B_true[1:length(B_nrest_pc)]-B_nrest_pc)
      }
    }
  }
}

##### Analyze Performance #####
perf_table <- function(perf) {
  tab1 <- matrix(c(mean(perf[,1,1,1]), mean(perf[,1,1,2]), mean(perf[,1,1,3]),
                   mean(perf[,2,1,1]), mean(perf[,2,1,2]), mean(perf[,2,1,3])), ncol=6)
  tab <- rbind(tab1, c(mean(perf[,1,2,1]), mean(perf[,1,2,2]), mean(perf[,1,2,3]),
                       mean(perf[,2,2,1]), mean(perf[,2,2,2]), mean(perf[,2,2,3])))
  colnames(tab) <- c('p=0.75', 'p=0.85', 'p=0.99', 'p=0.75', 'p=0.85', 'p=0.99')
  rownames(tab) <- c('Sige 0.5', 'Sige 1')
  return(tab)
}
tab_n <- perf_table(perf_n) #beta_n
tab_nrest <- perf_table(perf_nrest) #beta_nrest
tab_22 <- perf_table(perf_22) #beta from equation 22
tab_n_alt <- perf_table(perf_n_alt) #beta_n using sigma from page 13
tab_nrest_pc <- perf_table(perf_nrest_pc) #beta_nrest using Kaiser rule to select number of components

par(mfrow=c(1,2))
x <- c(0.75, 0.85, 0.99)
mse_05_05 <- tab[1,1:3]
mse_05_1 <- tab[1,4:6]
mse_1_05 <- tab[2,1:3]
mse_1_1 <- tab[2,4:6]
