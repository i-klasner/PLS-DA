# Ina Klasner
# 12/1/2021
# Simulation Study on PLSDA with Olive Oil & Gasoline Data - All ME


##### Add Packages and Data #####
library('pls')
library(caret)
library(Lock5Data)
data('gasoline')


##### Read in Data Files #####
# create a function for reading in data files
file_names <- function(name, gamma, rep) {
  #iris_sl_me_0.GAMMA_rep_REP.csv
  name <- paste0(name, '_me_0.', gamma, '_rep_', rep, '.csv')
  return(name)
}

# dataframes for analysis
gas.rmsep <- data.frame()
num_components <- data.frame()


##### Perform Analysis #####

string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/plsda/Gas/ME_All'
name <- 'gas'
# iterate through the folder for each gamma and rep number
for (i in seq(from=1, to=5, by=2)) { #through gamma values 0.1, 0.3, 0.5
  for (j in 1:100) { #through rep values 1, 2, ..., 100
    # read in data
    directory <- paste0(string, '/ME0', i, '/')
    data <- read.csv(paste0(directory, file_names(name, i, j))) #read in data file
    #data <- read.csv('/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/plsda/Gas/ME_All/ME01/gas_me_0.1_rep_1.csv')
    
    ##### Create Training and Testing Sets #####
    set.seed(100)
    index <- createDataPartition(data[,1], p=0.75, list=FALSE)
    train <- data[index,] #create a training data set of 75% of the original data
    test <- data[-index,] #create a testing set of all other data
    
    # fit plsr model with leave-one-out cross-validated predictions
    y <- train$octane
    model <- plsr(y ~ ., data=train, validation='LOO')
    
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
    
    # method 3: slope cutoff
    num_comps3 <- 1
    cutoff <- 1
    for (comp in 1:(length(perc_comp)-1)) {
      diff <- abs(perc_comp[comp]-perc_comp[comp+1])
      if (diff < cutoff) {
        num_comps3 <- comp
        break
      }
    }
    
    # method 4: rmsep
    rmsep_subset <- RMSEP(model)$val[1,,] #root mean squared error of prediction
    num_comps4 <- 1
    rmsep_diffs <- abs(diff(rmsep_subset))
    cutoff <- 1e-03
    for (comp in 1:(length(rmsep_diffs))) {
      if (rmsep_diffs[comp] < cutoff) {
        num_comps4 <- comp-1 
        break
      }
    }
    
    num_components <- rbind(num_components, c(num_comps1, num_comps3, num_comps4))
    num_pc <- round(mean(c(num_comps1, num_comps3, num_comps4)))
    
    
    ##### Inspect Fit #####
    # get RMSEPs
    train.rmsep <- RMSEP(model, ncomp=num_pc)$val[2,,]
    
    # predict response variables of new observations
    plsda.pred <- predict(model, ncomp=num_pc, newdata=test)
    
    # get R squared
    # https://rpubs.com/omicsdata/pls
    eval <- data.frame(obs=test[,1], pred=plsda.pred[,1,1])
    rsq <- defaultSummary(eval)[2]
    
    # manually find RMSEP for test
    test.rmsep <- sqrt(mean((plsda.pred-test$octane)^2))
    gas.rmsep <- rbind(gas.rmsep, test.rmsep, rsq)
    
  }
}


##### Format Output #####
compare <- data.frame()
for (y in 1:3) {
  compare[y,1] <- mean(gas.rmsep[seq(from=(y-1)*600+1, to=y*600, by=6),1]) #gas rmsep case 1
  compare[y,2] <- mean(gas.rmsep[seq(from=(y-1)*600+2, to=y*600, by=6),1]) #gas r2 case 1
}

row.names(compare) <- c('0.1', '0.3', '0.5')
colnames(compare) <- c('Gas_RMSEP', 'Gas_R2')
compare


##### Models without Measurement Error #####
# make gas plsda model without ME
set.seed(100)
index <- createDataPartition(data[,1], p=0.75, list=FALSE)
gas_train <- gasoline[index,] #training set
gas_test <- gasoline[-index,] # testing set
gas1 <- plsr(octane ~ NIR, data=gas_train, validation='LOO')
#plot(RMSEP(gas1), legendpos='topright')
gas_pred <- predict(gas1, ncomp=2, newdata=gas_test)
gas_eval <- data.frame(obs=gas_test[,1], pred=gas_pred[,1,1])
gas_rsq <- defaultSummary(gas_eval)[2]
gas_rmse <- defaultSummary(gas_eval)[1]


##### Visualize Results #####
par(mfrow=c(1,2))
x <- c(0.1, 0.3, 0.5)
x1 <- c(0, 0.1, 0.3, 0.5)
gas_rmsep <- unlist(compare[,1], use.names = FALSE)
gas_rmsep1 <- c(gas_rmse, gas_rmsep)
gas_r2 <- unlist(compare[,2], use.names = FALSE)
gas_r21 <- c(gas_rsq, gas_r2)

# plot the measurement error cases
plot(x, gas_rmsep,
     main='Gasoline RMSE',
     xlab='Gamma Value', ylab='RMSEP',
     col='red', ylim=c(0.023,0.025))

plot(x, gas_r2,
     main='Gasoline R Squared',
     xlab='Gamma Value', ylab='R Squared',
     col='blue', ylim=c(0.99978,0.99979))

# plot measurement error cases with no measurement error
plot(x1, gas_rmsep1,
     main='Gasoline RMSE',
     xlab='Gamma Value', ylab='RMSEP',
     col='red', ylim=c(0.023,0.695))

plot(x1, gas_r21,
     main='Gasoline R Squared',
     xlab='Gamma Value', ylab='R Squared',
     col='blue', ylim=c(0.778,0.99979))
