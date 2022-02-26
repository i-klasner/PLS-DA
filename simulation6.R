# Ina Klasner
# 2/22/2022
# Simulation Study on PCR with ME Version 6
# Include more estimated betas to fill Table 1

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
perf <- array(0, c(100, 2, 2, 3, 8)) #matrix of performances for all B

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
        cov_X <- t(X)%*%X/nrow(data) #maximum likelihood estimator
        B_true <- data$true.beta
        
        # create sigma matrices
        if (g==1) { # ME=0.5
          cov_dd <- diag(0.5,5,5)
        } else {
          cov_dd <- diag(1,5,5)
        }
        
        # define reliability matrix K
        # K <- solve(cov_W + cov_dd)%*%cov_W
        K_hat <- solve(cov_W)%*%(cov_W-nrow(data)*cov_dd)
        
        # define H and h
        H <- matrix(c(1,1,1,1,1), nrow=1, ncol=5)
        h <- H%*%B_true[1:length(B_n)] + nrow(data)^(-1/2)
        
        # define P
        W_c <- W - mean(W) # centered
        l <- eigen(t(W_c)%*%W_c)$values #eigenvalues of Wc'W
        L <- matrix(diag(l), ncol=length(l)) # diagonal matrix of eigenvalues
        P <- solve(W_c[1:5,])%*%t(W)
        
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
        
        ##### Find Beta n, nr, nrest, nr,rest#####
        # Equation 9
        B_n <- solve(cov_W-nrow(data)*cov_dd)%*%cov_WY 
        
        # restricted least squares estimator Equation 13
        B_nrest <- B_n - solve(t(K_hat)%*%cov_W%*%K_hat)%*%t(H)%*%solve(H%*%solve(t(K_hat)%*%cov_W%*%K_hat)%*%t(H))%*%(H%*%B_n - h) # Equation 13
        
        # principal components estimator Equation 24
        # select number of PC => r = num_pc
        P <- V
        P_r <- matrix(P[,num_pc],ncol=num_pc)
        B_nr <- solve(cov_W-nrow(data)*cov_dd)%*%cov_W%*%P_r%*%solve(t(P_r)%*%cov_W%*%P_r)%*%t(P_r)%*%cov_WY # Equation 24 (21)
        # Equation 23
        B_r <- solve(K_hat)%*%P_r%*%solve(t(P_r)%*%cov_W%*%P_r)%*%t(P_r)%*%cov_WY
        
        # principle components estimator with restriction Equation 39
        B_nrrest <- B_nr - solve(K_hat)%*%P_r%*%solve(t(P_r)%*%cov_W%*%P_r)%*%t(P_r)%*%solve(t(K_hat))%*%t(H)%*%
          solve(H%*%solve(K_hat)%*%P_r%*%solve(t(P_r)%*%cov_W%*%P_r)%*%t(P_r)%*%solve(t(K_hat))%*%t(H))%*%(H%*%B_nr - h)
        # Equation 37
        B_rrest <- B_r - solve(K_hat)%*%P_r%*%solve(t(P_r)%*%cov_W%*%P_r)%*%t(P_r)%*%t(H)%*%solve(H%*%P_r%*%solve(t(P_r)%*%cov_W%*%P_r)%*%t(P_r)%*%t(H))%*%(H%*%K%*%B_r - h)
        
        new_V <- matrix(V[,num_pc],ncol=num_pc)
        B_pcr <- V%*%solve(S)%*%t(U)%*%Y # B = V.(S^-1).(U^T).Y with all PCs
        B_pcr_pc <- new_V%*%solve(S[num_pc,num_pc])%*%t(matrix(U[,num_pc],ncol=num_pc))%*%Y #selecting num of PCs

        ##### Find MSE between estimated Betas and true Beta, add to performance matrix ####
        betas <- list(B_n, B_nr, B_nrest, B_nrrest, B_r, B_rrest, B_pcr, B_pcr_pc)
        for (b in 1:length(betas)) {
          B <- betas[b]
          perf[r,g,s,p,b] <- t(B_true[1:length(B)]-B[[1]])%*%(B_true[1:length(B)]-B[[1]])
        }
      }
    }
  }
}

##### Analyze Performance #####
perf_table <- function(s) {
  tab <- matrix(c(mean(perf[,1,s,1,1]), mean(perf[,1,s,2,1]), mean(perf[,1,s,3,1]),
                   mean(perf[,2,s,1,1]), mean(perf[,2,s,2,1]), mean(perf[,2,s,3,1])), ncol=6)
  for (e in 2:8) {
    tab <- rbind(tab, c(mean(perf[,1,s,1,e]), mean(perf[,1,s,2,e]), mean(perf[,1,s,3,e]),
                     mean(perf[,2,s,1,e]), mean(perf[,2,s,2,e]), mean(perf[,2,s,3,e])))
  }
  colnames(tab) <- c('p=0.75', 'p=0.85', 'p=0.99', 'p=0.75', 'p=0.85', 'p=0.99')
  rownames(tab) <- c('B_n', 'B_nr', 'B_nrest', 'B_nr,rest', 'B_r', 'B_rrest', 'B_pcr', 'B_pcr_pc')
  return(tab)
}
tab_0.5 <- perf_table(1)
tab_1 <- perf_table(2)
