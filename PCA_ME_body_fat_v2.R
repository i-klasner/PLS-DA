# Ina Klasner
# 10/24/2021
# PCA with Measurement Error on Height and Weight, Heterogeneous and Homogeneous 
# Body Fat Dataset

##### Load Libraries #####
library(pls)
library(caret)
library(Lock5Data)


##### Load Datasets #####
data("BodyFat") #get body fat dataset

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
num_components <- data.frame()
##### Perform Analysis #####
#iterate through folders Height and Weight
for (x in 1:4) {
  #create string for file directory
  if(x == 1) { #go through Height ME files first
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/Body_Fat/Height'
  }
  else if (x == 2) { # then go through Weight ME files
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/Body_Fat/Weight'
  }
  else if (x == 3) {
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/Body_Fat/Homo_ME_Bodyfat'
  }
  else {
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/Body_Fat/Hetero_ME_Bodyfat'
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
      ##### Component Selection #####
      # validationplot(pcr.fit, val.type = 'MSEP') #mean squared error of prediction
      # validationplot(pcr.fit, val.type='R2')
      # % validation 
      perc_comp <- cumsum(explvar(pcr.fit))
      
      # method 1: percentage rule, with cutoff at 90%
      num_comps1 <- 9 #number components by method 1
      for (comp in 1:length(perc_comp)) { #iterate through the principle components
        if (perc_comp[comp] > 90) { #stop when reach 90% mark
          num_comps1 <- comp #set number of components to use
          break
        }
      }
      
      # method 2: kaiser rule, adapt to introduced bias
      pca <- princomp(~., data=train)
      pca$loadings
      eigen <- eigenvals(pca)
      #screeplot(pca)
      num_comps2 <- 9
      for (comp in 1:(length(eigen))) {
        if(eigen[comp] < mean(eigen)*0.7) {
          num_comps2 <- comp
          break
        }
      }
      
      # method 3: slope cutoff
      num_comps3 <- 9
      for (comp in 1:(length(perc_comp)-1)) {
        diff <- abs(perc_comp[comp]-perc_comp[comp+1])
        #print(diff)
        if (diff < 4) {
          num_comps3 <- comp
          break
        }
      }
      
      # method 4: rmsep
      rmsep_subset <- RMSEP(pcr.fit)$val[1,,] #root mean squared error of prediction
      num_comps4 <- 9
      rmsep_diffs <- abs(diff(rmsep_subset))
      for (comp in 1:(length(rmsep_diffs))) {
        if (rmsep_diffs[comp] < 0.1) {
          num_comps4 <- comp-1 
          break
        }
      }
      #num_comps4 <- which.min(rmsep_subset) - 1 #get the best # of PCs based on cross validation
      #num_components <- rbind(num_components, c(num_pc))
      num_components <- rbind(num_components, c(num_comps1, num_comps2, num_comps3, num_comps4))
      num_pc <- round(mean(c(num_comps1, num_comps2, num_comps3, num_comps4)))
      
      ##### Predict on Testing Data Set #####
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
for (y in 1:12) {
  #compare <- rbind(compare, c(mean(mse[(y-1)*100+1:y*100, 1]), mean(mse[(y-1)*100+1:y*100, 2])))
  compare[y,1] <- mean(mse[((y-1)*100+1):(y*100), 1]) #find mean MSE
  compare[y,2] <- mean(mse[((y-1)*100+1):(y*100), 2]) #find mean adjusted r^2
  compare[y,3] <- mean(num_components[((y-1)*100+1):(y*100), 1]) # method 1: 90% cutoff
  compare[y,4] <- mean(num_components[((y-1)*100+1):(y*100), 2]) # method 2: Kaiser rule
  compare[y,5] <- mean(num_components[((y-1)*100+1):(y*100), 3]) # method 3: slope cutoff
  compare[y,6] <- mean(num_components[((y-1)*100+1):(y*100), 4]) # method 4: rmsep
  #print((y-1)*100+1)
  #print(mean(mse[(y-1)*100+1:y*100, 2]))
}

row.names(compare) <- c('Height, 0.1', 'Height, 0.3', 'Height, 0.5',
                        'Weight, 0.1', 'Weight, 0.3', 'Weight, 0.5',
                        'Homogeneous, 0.1', 'Homogeneous, 0.3', 'Homogeneous, 0.5',
                        'Heterogeneous, 0.1', 'Heterogeneous, 0.3', 'Heterogeneous, 0.5')
colnames(compare) <- c('MSE', 'Adjusted R Squared', 'Num Comps %', 'Num Comps Kaiser', 'Num Comps slope', 'Num Comps RMSEP')
compare

par(mfrow=c(1,2))
x=c(0.1, 0.3, 0.5)
h_mse <- unlist(compare[1:3,1], use.names = FALSE)
w_mse <- unlist(compare[4:6,1], use.names = FALSE)
hom_mse <- unlist(compare[7:9,1], use.names = FALSE)
het_mse <- unlist(compare[10:12,1], use.names = FALSE)
h_rsq <- unlist(compare[1:3,2], use.names = FALSE)
w_rsq <- unlist(compare[4:6,2], use.names = FALSE)
hom_rsq <- unlist(compare[7:9,2], use.names = FALSE)
het_rsq <- unlist(compare[10:12,2], use.names = FALSE)

plot(x, h_mse,
     main='MSE',
     xlab='Gamma Value', ylab='MSE',
     col='red', ylim=c(27,38))
points(x, w_mse, col='blue')
points(x, hom_mse, col='orange')
points(x, het_mse, col='forest green')
legend("topleft",   # Position
       inset=0.05,
       legend = c("Height", "Weight", 'Homo', 'Hetero'),
       fill = c('red', 'blue', 'orange', 'forest green'), cex=0.8 )
plot(x, h_rsq,
     main='Adjusted R Squared',
     xlab='Gamma Value', ylab='Adj R^2',
     col='red', ylim=c(0.45,0.75))
points(x, w_rsq, col='blue')
points(x, hom_rsq, col='orange')
points(x, het_rsq, col='forest green')
# red is ME on height
# blue is ME on weight with greater variance, increased MSE, relatively unchanged r^2
# orange is homogeneous case
# green is heterogeneous case


##### Review Original BodyFat Dataset #####
var(BodyFat$Height)
var(BodyFat$Weight)
sd(BodyFat$Height)
sd(BodyFat$Weight)
reg <- lm(Bodyfat~., data=BodyFat)
summary(reg) #Adjusted R-squared = 0.7332
