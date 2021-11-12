# Ina Klasner
# 11/09/2021
# PLS-DA with Measurement Error

##### Add Packages and Data #####
library('pls')
library(caret)
library(Lock5Data)
data('oliveoil')
data('gasoline')


##### Read in Data Files #####
#create a function for reading in data files
file_names <- function(name, gamma, rep) {
  #iris_sl_me_0.GAMMA_rep_REP.csv
  name <- paste0(name, '_me_0.', gamma, '_rep_', rep, '.csv')
  return(name)
}

# dataframes for analysis
olive.rmsep <- data.frame()
gas.rmsep <- data.frame()
num_components <- data.frame()


##### Perform Analysis #####
for (x in 1:2) { #iterate through olive and gas folders
  # create string for file directory
  if(x == 1) { #go through Olive ME files first
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/plsda/Olive'
    name <- 'olive'
  } else {
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/plsda/Gas'
    name <- 'gas'
  }
  
  # iterate through the folder for each gamma and rep number
  for (i in seq(from=1, to=5, by=2)) { #through gamma values 0.1, 0.3, 0.5
    for (j in 1:100) { #through rep values 1, 2, ..., 100
      
      # read in data
      directory <- paste0(string, '/ME0', i, '/')
      data <- read.csv(paste0(directory, file_names(name, i, j))) #read in data file
      #data <- read.csv('/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/plsda/Olive/ME01/olive_me_0.1_rep_1.csv')
      #data <- read.csv('/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/plsda/Gas/ME01/gas_me_0.1_rep_1.csv')
      #name <- 'olive'
      
      # create training and testing data sets
      set.seed(100)
      index <- createDataPartition(data[,1], p=0.75, list=FALSE)
      train <- data[index,] #create a training data set of 75% of the original data
      test <- data[-index,] #create a testing set of all other data

      # fit plsr model with leave-one-out cross-validated predictions
      if (name == 'olive') {
        # xnam <- paste0('chemical.', c('Acidity', 'Peroxide', 'K232', 'K270', 'DK'))
        # ynam <- paste0('sensory.', c('yellow', 'green', 'brown', 'glossy', 'transp', 'syrup'))
        # (fmla <- as.formula(paste(paste(ynam, collapse= '+'), '~' , paste(xnam, collapse= '+'))))
        y <- cbind(train$sensory.yellow, train$sensory.green, train$sensory.brown, 
                   train$sensory.glossy, train$sensory.transp, train$sensory.syrup)
      } else {
        # NIR <- train[,2:402]
        # xnam <- paste0('NIR.', seq(from=900, to=1700, by=2), '.nm')
        # (fmla <- as.formula(paste('octane ~ ', paste(xnam, collapse= '+'))))
        y <- train$octane
        
      }
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
      
      # # method 2: kaiser rule, adapt to introduced bias
      # pca <- princomp(~., data=train)
      # pca$loadings
      # eigen <- eigenvals(pca)
      # num_comps2 <- 1
      # for (comp in 1:(length(eigen))) {
      #   if(eigen[comp] < mean(eigen)*0.7) {
      #     num_comps2 <- comp
      #     break
      #   }
      # }
      
      # method 3: slope cutoff
      num_comps3 <- 1
      # plot(RMSEP(model), legendpos='topright')
      if (name=='olive') {
        cutoff <- 3
      } else {
        cutoff <- 1
      }
      for (comp in 1:(length(perc_comp)-1)) {
        diff <- abs(perc_comp[comp]-perc_comp[comp+1])
        #print(diff)
        if (diff < cutoff) {
          num_comps3 <- comp
          break
        }
      }
      
      # method 4: rmsep
      rmsep_subset <- RMSEP(model)$val[1,,] #root mean squared error of prediction
      num_comps4 <- 1
      rmsep_diffs <- abs(diff(rmsep_subset))
      if (name=='olive') {
        cutoff <- 1
      } else {
        cutoff <- 1e-03
      }
      for (comp in 1:(length(rmsep_diffs))) {
        #print(rmsep_diffs[comp])
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
      # test_rmsep <- RMSEP(model, newdata=test)$val[2,,]
      
      # predict response variables of new observations
      plsda.pred <- predict(model, ncomp=num_pc, newdata=test)
      
      # manually find RMSEP for test
      if (name=='olive') {
        #test.rmsep <- sqrt(mean((as.data.frame(plsda.pred)-test[,1:6])^2))
        plsda.pred <- as.data.frame(plsda.pred)
        sum <- 0
        for (q in 1:nrow(plsda.pred)) {
          for (w in 1:ncol(plsda.pred)) {
            #q <- as.integer(q)
            sum <- sum + (plsda.pred[q,w]-test[q,w])^2
            #print(sum)
          }
        }
        test.rmsep <- sqrt(sum/nrow(plsda.pred))
        olive.rmsep <- rbind(olive.rmsep, test.rmsep)
      } else {
        test.rmsep <- sqrt(mean((plsda.pred-test$octane)^2))
        gas.rmsep <- rbind(gas.rmsep, test.rmsep)
      }
    }
  }
}


##### Format Output #####
compare <- data.frame()
for (y in 1:3) {
  compare[y,1] <- mean(olive.rmsep[((y-1)*100+1):(y*100), 1]) #olive rmsep
  compare[y,2] <- mean(gas.rmsep[((y-1)*100+1):(y*100), 1]) #find mean adjusted r^2
}

row.names(compare) <- c('0.1', '0.3', '0.5')
colnames(compare) <- c('Olive Oil', 'Gasoline')
compare


##### Visualize Results #####
par(mfrow=c(1,2))
x=c(0.1, 0.3, 0.5)
olive <- unlist(compare[1:3,1], use.names = FALSE)
gas <- unlist(compare[1:3,2], use.names = FALSE)

plot(x, olive,
     main='Olive Oil',
     xlab='Gamma Value', ylab='RMSEP',
     col='red')
plot(x, gas,
     main='Gasoline',
     xlab='Gamma Value', ylab='RMSEP',
     col='blue')
