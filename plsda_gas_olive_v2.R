# Ina Klasner
# 1/13/2022
# Simulation Study on PLSDA with Olive Oil & Gasoline Data

##### Add Packages and Data #####
library('pls')
library(caret)
library(Lock5Data)
data('oliveoil')
data('gasoline')


##### Generate Indices with Measurement Error #####
# gasoline data
set.seed(123)
n <- 60
gas_me <- sample(c(1:n), round(n*0.3), replace=FALSE)

# oliveoil data
olive_me <- c(1,2,6,7,12,13)


##### Read in Data Files #####
# create a function for reading in data files
file_names <- function(name, gamma, rep) {
  #iris_sl_me_0.GAMMA_rep_REP.csv
  name <- paste0(name, '_me_0.', gamma, '_rep_', rep, '.csv')
  return(name)
}

# dataframes for analysis
olive_perf <- array(0, c(100, 2*3, 3))
gas_perf <- array(0, c(100, 2*3, 3))
num_components <- data.frame()


##### Perform Analysis #####
for (x in 1:2) { #iterate through olive and gas folders
  # create string for file directory
  if(x == 1) { #go through Olive ME files first
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/plsda/Olive'
    name <- 'olive'
    index <- olive_me #the indices with measurement error
  } else {
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/plsda/Gas'
    name <- 'gas'
    index <- gas_me # the indices with measurement error
  }
  
  i <- seq(1,5,2)
  # iterate through the folder for each gamma and rep number
  for (b in 1:3) { #through gamma values 0.1, 0.3, 0.5
    rmsep <- array(0, c(100,3))
    rsq <- array(0, c(100,3))
    for (j in 1:100) { #through rep values 1, 2, ..., 100
      for (case in 1:3) { #through each measurement error case
        # read in data
        directory <- paste0(string, '/ME0', i[b], '/')
        data <- read.csv(paste0(directory, file_names(name, i[b], j))) #read in data file
        
        ##### Create Training and Testing Sets #####
        # create training and testing sets based on case number
        # for olive, testing is 37.5% of data
        # for gas, testing is 30% of data
        if (case == 1) { # ME all in training
          if (name == 'gas') {
            rest_train <- sample(as.integer(row.names(data[-index,])), nrow(data) - (2*length(index)), replace=FALSE)
          } else {
            rest_i <- sample(3:5, 1, replace=FALSE)
            rest_i2 <- sample(c(11, 14, 15, 16), 2, replace=FALSE)
            rest_train <- c(rest_i, rest_i+5, rest_i2)
          }
          train <- data[c(index,rest_train),]
          test <- data[-as.integer(row.names(train)),]
        } else if (case == 2) { # ME all in testing
          test <- data[index,] #put all ME indices into testing
          train <- data[-index,] #rest into training
        } else { # ME in both
          if (name == 'gas') {
            half_i <- sample(index, length(index)/2, replace=FALSE)
            rest_test <- sample(as.integer(row.names(data[-index,])), length(index)/2, replace=FALSE)
          } else {
            half_i <- c(1,6,12)
            rest_i <- sample(3:5, 1, replace=FALSE)
            rest_i2 <- sample(c(11, 14, 15, 16), 1, replace=FALSE)
            rest_test <- c(rest_i, rest_i+5, rest_i2)
          }
          test <- data[c(half_i,rest_test),]
          train <- data[-as.integer(row.names(test)),]
        }
        
        # fit plsr model with leave-one-out cross-validated predictions
        if (name == 'olive') {
          model <- plsr(cbind(train$sensory.yellow, train$sensory.green, train$sensory.brown, 
                              train$sensory.glossy, train$sensory.transp, train$sensory.syrup) ~ ., data=train, validation='LOO')
        } else {
          model <- plsr(octane ~ ., data=train, validation='LOO')
        }
        
        ##### Component Selection #####
        perc_comp <- cumsum(explvar(model))
        
        # method 1: percentage rule, with cutoff at 90%
        num_comps1 <- 1 #number components by method 1
        for (comp in 1:length(perc_comp)) { #iterate through the principle components
          if (perc_comp[comp] > 90) { #stop when reach 90% mark
            num_comps1 <- comp #set number of components to use
            break
          }
        }
        
        # method 2: slope cutoff
        num_comps2 <- 1
        # plot(RMSEP(model), legendpos='topright')
        if (name=='olive') {
          cutoff <- 3
        } else {
          cutoff <- 0.01
        }
        for (comp in 1:(length(perc_comp)-1)) {
          diff <- abs(perc_comp[comp]-perc_comp[comp+1])
          if (diff < cutoff) {
            num_comps2 <- comp
            break
          }
        }
        
        # method 3: rmsep
        rmsep_subset <- RMSEP(model)$val[1,,] #root mean squared error of prediction
        num_comps3 <- 1
        rmsep_diffs <- abs(diff(rmsep_subset))
        if (name=='olive') {
          cutoff <- 4.5
        } else {
          cutoff <- 1e-03
        }
        for (comp in 1:(length(rmsep_diffs))) {
          if (rmsep_diffs[comp] < cutoff) {
            num_comps3 <- comp-1 
            break
          }
        }
        
        num_components <- rbind(num_components, c(num_comps1, num_comps2, num_comps3))
        num_pc <- round(mean(c(num_comps1, num_comps2, num_comps3)))
        
        
        ##### Inspect Fit #####
        # predict response variables of new observations
        plsda.pred <- predict(model, ncomp=num_pc, newdata=test)
        
        # get R squared and RMESP https://rpubs.com/omicsdata/pls
        eval <- data.frame(obs=test[,1], pred=plsda.pred[,1,1])
        rsq[j,case] <- defaultSummary(eval)[2]
        rmsep[j,case] <- defaultSummary(eval)[1]
        
        if(j == 100) {
          if (name=='olive') {
            olive_perf[,,b] <- cbind(rmsep, rsq)
          }
          else {
            gas_perf[,,b] <- cbind(rmsep, rsq)
          }
        }
      }
    }
  }
}


##### Format Output #####
compare <- data.frame()

for (y in 1:3) {
  s <- 0
  for (col in 1:6) {
    s <- s+1
    compare[y,s] <- mean(olive_perf[,col,y])
    s <- s+1
    compare[y,s] <- mean(gas_perf[,col,y])
  }
}
row.names(compare) <- c('0.1', '0.3', '0.5')
colnames(compare) <- c('Olive RMSEP 1', 'Gas RMSEP 1', 'Olive RMSEP 2', 'Gas RMSEP 2', 
                       'Olive RMSEP 3', 'Gas RMSEP 3', 'Olive R21', 'Gas R2 1', 
                       'Olive R2 2', 'Gas R2 2', 'Olive R2 3', 'Gas R2 3')
compare


##### Models without Measurement Error #####
# make gas plsda model without ME
gas_train <- gasoline[-gas_me,] #training set
gas_test <- gasoline[gas_me,] # testing set
gas1 <- plsr(octane ~ NIR, data=gas_train, validation='LOO')
gas_pred <- predict(gas1, ncomp=2, newdata=gas_test)
gas_eval <- data.frame(obs=gas_test[,1], pred=gas_pred[,1,1])
gas_rsq <- defaultSummary(gas_eval)[2]
gas_rmse <- defaultSummary(gas_eval)[1]

# make olive oil plsda model without ME
olive_train <- oliveoil[-olive_me,] #training set
olive_test <- oliveoil[olive_me,] # testing set
olive1 <- plsr(sensory~chemical, scale=TRUE, data=oliveoil, validation='LOO')
olive_pred <- predict(olive1, ncomp=3, newdata=olive_test)
olive_eval <- data.frame(obs=olive_test[,1], pred=olive_pred[,1,1])
olive_rsq <- defaultSummary(olive_eval)[2] ##### NA
olive_rmse <- defaultSummary(olive_eval)[1]


##### Visualize Results #####
par(mfrow=c(1,2))
x <- c(0.1, 0.3, 0.5)
olive1_rmsep <- unlist(compare[,1], use.names = FALSE)
olive2_rmsep <- unlist(compare[,3], use.names = FALSE)
olive3_rmsep <- unlist(compare[,5], use.names = FALSE)
olive1_r2 <- unlist(compare[,7], use.names = FALSE)
olive2_r2 <- unlist(compare[,9], use.names = FALSE)
olive3_r2 <- unlist(compare[,11], use.names = FALSE)
gas1_rmsep <- unlist(compare[,2], use.names = FALSE)
gas2_rmsep <- unlist(compare[,4], use.names = FALSE)
gas3_rmsep <- unlist(compare[,6], use.names = FALSE)
gas1_r2 <- unlist(compare[,8], use.names = FALSE)
gas2_r2 <- unlist(compare[,10], use.names = FALSE)
gas3_r2 <- unlist(compare[,12], use.names = FALSE)

plot(x, olive1_rmsep,
     main='Olive Oil',
     xlab='Gamma Value', ylab='RMSEP',
     col='red', ylim=c(16,22))
points(x, olive2_rmsep, col='blue')
points(x, olive3_rmsep, col='forest green')
legend("topleft",
       inset=0.03,
       legend = c("Case 1", "Case 2", 'Case 3'),
       fill = c('red', 'blue', 'forest green'),
       cex=0.7)

plot(x, olive1_r2,
     main='Olive Oil',
     xlab='Gamma Value', ylab='R Squared',
     col='red', ylim=c(0.1,0.3))
points(x, olive2_r2, col='blue')
points(x, olive3_r2, col='forest green')

plot(x, gas1_rmsep,
     main='Gasoline',
     xlab='Gamma Value', ylab='RMSEP',
     col='red', ylim=c(0.3,0.9))
points(x, gas2_rmsep, col='blue')
points(x, gas3_rmsep, col='forest green')
legend("topleft", 
       inset=0.03,
       legend = c("Case 1", "Case 2", 'Case 3'),
       fill = c('red', 'blue', 'forest green'),
       cex=0.7)

plot(x, gas1_r2,
     main='Gasoline',
     xlab='Gamma Value', ylab='R Squared',
     col='red', ylim=c(0.5,1))
points(x, gas2_r2, col='blue')
points(x, gas3_r2, col='forest green')
