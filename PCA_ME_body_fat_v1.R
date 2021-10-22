# Ina Klasner
# 10/20/2021
# PCA with Measurement Error on Height and Weight Variables
# Body Fat Dataset

##### Load Libraries #####
library(pls)
library(TH.data)


##### Load Datasets #####
data('BodyFat') #get body fat dataset

#function for reading in data files
file_names <- function(gamma, rep) {
  #iris_sl_me_0.GAMMA_rep_REP.csv
  paste0('bodyfat_me_0.',
         gamma,
         '_rep_',
         rep, 
         '.csv')
}

mse <- data.frame()

##### Perform Analysis #####
#iterate through folders Height and Weight
for (x in 1:2) {
  #create string for file directory
  if(x == 1) { #go through Height ME files first
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/Body_Fat/Height'
  }
  else { # then go through Weight ME files
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/Body_Fat/Weight'
  }
  
  #iterate through the folder for each gamma and rep number
  for (i in seq(from=1, to=5, by=2)) { #iterate through gamma values 0.1, 0.3, 0.5
    for (j in 1:100) { #iterate through rep values 1, 2, ..., 100
      directory <- paste0(string, '/ME0', i, '/')
      data <- read.csv(paste0(directory, file_names(i, j))) #read in data file
      
      #create training and testing data sets
      set.seed(100)
      index <- createDataPartition(BodyFat$Age, p=0.75, list=FALSE)
      training <- data[index,] #create a training data set of 75% of the original data
      testing <- data[-index,] #create a testing set of all other data
      train <- cbind(training, BodyFat[index,1])
      test <- cbind(testing, BodyFat[-index, 1])
      colnames(train)[10] <- c('Bodyfat')
      colnames(test)[10] <- c('Bodyfat')
      
      #principal component regression on just the training set
      pcr.fit <- pcr(Bodyfat~., data=train, scale=TRUE, validation='CV')
      
      #determine the best number of principal components
      #validationplot(pcr.fit, val.type = 'MSEP') #mean squared error of prediction
      rmsep_subset <- RMSEP(pcr.fit)$val[1,,] #root mean squared error of prediction
      num_pc <- which.min(rmsep_subset) - 1 #get the best # of PCs based on cross validation
      
      #predict on testing data set
      pcr.pred <- predict(pcr.fit, test, ncomp=num_pc)
      pcr.mse <- mean((pcr.pred-test$Bodyfat)^2)
      pcr.reg <- lm(Bodyfat~., data=test)
      summary.reg <- summary.lm(pcr.reg)
      mse <- rbind(mse, c(pcr.mse, summary.reg$adj.r.squared))
    }
  }
}


##### Format Output #####
compare <- data.frame()
for (y in 1:6) {
  #compare <- rbind(compare, c(mean(mse[(y-1)*100+1:y*100, 1]), mean(mse[(y-1)*100+1:y*100, 2])))
  compare[y,1] <- mean(mse[(y-1)*100+1:y*100, 1])
  compare[y,2] <- mean(mse[(y-1)*100+1:y*100, 2])
  #print((y-1)*100+1)
  #print(mean(mse[(y-1)*100+1:y*100, 2]))
}
compare[4,1] <- mean(mse[301:400, 1])
compare[4,2] <- mean(mse[301:400, 2])
compare[5,1] <- mean(mse[401:500, 1])
compare[5,2] <- mean(mse[401:500, 2])
compare[6,1] <- mean(mse[501:600, 1])
compare[6,2] <- mean(mse[501:600, 2])
row.names(compare) <- c('Height, 0.1', 'Height, 0.3', 'Height, 0.5',
                        'Weight, 0.1', 'Weight, 0.3', 'Weight, 0.5')
colnames(compare) <- c('MSE', 'Adjusted R Squared')
compare

par(mfrow=c(1,2))
x=c(0.1, 0.3, 0.5)
h_mse <- unlist(compare[1:3,1], use.names = FALSE)
w_mse <- unlist(compare[4:6,1], use.names = FALSE)
h_rsq <- unlist(compare[1:3,2], use.names = FALSE)
w_rsq <- unlist(compare[4:6,2], use.names = FALSE)
plot(x, h_mse,
     main='MSE',
     xlab='Gamma Value', ylab='MSE',
     col='red', ylim=c(19,20.5))
points(x, w_mse, col='blue')
plot(x, h_rsq,
     main='Adjusted R Squared',
     xlab='Gamma Value', ylab='Adj R^2',
     col='red', ylim=c(0.65,0.75))
points(x, w_rsq, col='blue')
#red is ME on height
#blue is ME on weight


##### Review Original BodyFat Dataset #####
var(BodyFat$Height)
var(BodyFat$Weight)
sd(BodyFat$Height)
sd(BodyFat$Weight)
