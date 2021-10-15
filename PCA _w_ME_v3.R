# Ina Klasner
# 10/12/2021
# PCA With Measurement Error & iris data : Sepal Length and Petal Length Measurement Error

library(caret)
library(e1071)
library(nnet)
library(datasets)

data("iris")
str(iris)

##### Read in Data Files #####
#create a function for reading in data files
file_names <- function(gamma, rep) {
  #iris_sl_me_0.GAMMA_rep_REP.csv
  name <- paste('iris_sl_me_0.',
                 gamma,
                 '_rep_',
                 rep, 
                 '.csv', sep='')
  return(name)
}


#lists for the performance measures of the different measurement error types
train_meas <- list()
test_meas <- list()

#compare with example data
sepal_01 <- c()
sepal_05 <- c()

files <- c()

##### Perform Analysis on Each File #####
for(x in 1:2) {
  #create an empty dataframe for performance measures
  perf_meas <- data.frame()
  perf_meas_test <- data.frame()
  counter <- 1
  if(x == 1) { #go through Sepal Length ME files first
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda'
  }
  else { # then go through Petal Length ME files
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/Petal_Length'
  }
  #iterate through folders, gamma and rep values
  for(i in seq(from=1, to=5, by=2)) { #iterate through gammas
    for(j in 1:100) { #iterate through files
      #i <- 1
      #j <-3
      directory <- paste0(string, '/ME0', i, '/')
      data <- read.csv(paste0(directory, file_names(i, j))) #read in data file
      files <- rbind(files, i)#paste0(directory, file_names(i, j)))
      #create training and testing sets
      set.seed(100)
      index <- createDataPartition(data$Sepal.Length, p=0.75, list=FALSE)
      training <- data[index,] #create a training data set of 75% of the original data
      testing <- data[-index,] #create a testing set of all other data
      #train_i <- iris[index,] #use species classifications from iris data
      #test_i <- iris[-index,]
      train_i <- cbind(training, iris[index,5])
      test_i <- cbind(testing, iris[-index, 5])
      
      #pca
      pca <- prcomp(training[,1:4],center=TRUE, scale.=TRUE)
      #predict <- predict(pca,training)
      predict <- predict(pca,train_i)
      train <- data.frame(predict, train_i[,5])
      #predict_1 <- predict(pca, testing)
      predict_1 <- predict(pca, test_i)
      test <- data.frame(predict_1, test_i[,5])
      
      #multinomial logistic regression
      set.seed(100)
      colnames(train) <- c('PC1', 'PC2', 'PC3', 'PC4', 'Species')
      model <- multinom(Species~PC1 + PC2, data=train)
      
      #evaluate performance, misclassification error on training
      pred <- predict(model, train)
      conf_matrix <- confusionMatrix(pred, train[,5]) #list
      stats <- data.frame(conf_matrix$byClass)
      
      #evaluate performance, misclassification error on testing
      pred_test <- predict(model, test)
      conf_matrix_test <- confusionMatrix(pred_test, test[,5]) #list
      stats_test <- data.frame(conf_matrix_test$byClass)
      
      #find misclassfication errors and calculate accuracy
      sp_acc <- c(0,0,0)
      sp_acc_test <- c(0,0,0)
      for (y in 1:3) {
        #training
        tp <- conf_matrix$table[y,y]
        fp <- sum(conf_matrix$table[y,-y])
        fn <- sum(conf_matrix$table[-y,y])
        tn <- sum(conf_matrix$table[-y,-y])
        acc <- (tp+tn)/(tp+fp+fn+tn)
        sp_acc[y] <- acc
        
        #testing
        tp_t <- conf_matrix_test$table[y,y]
        fp_t <- sum(conf_matrix_test$table[y,-y])
        fn_t <- sum(conf_matrix_test$table[-y,y])
        tn_t <- sum(conf_matrix_test$table[-y,-y])
        acc_t <- (tp_t+tn_t)/(tp_t+fp_t+fn_t+tn_t)
        sp_acc_test[y] <- acc_t
      }
      
      perf_meas <- rbind(perf_meas, c(stats$Specificity, stats$Sensitivity, stats$Precision, sp_acc))
      perf_meas_test <- rbind(perf_meas_test,
                              c(stats_test$Specificity, stats_test$Sensitivity, stats_test$Precision, sp_acc_test))
      
      counter <- counter+1
    }
  }
  
  ##### Get Statistics to Evaluate Performance ######
  
  colnames(perf_meas) <- c('set_SPE', 'ver_SPE', 'vir_SPE',
                           'set_SEN', 'ver_SEN', 'vir_SEN',
                           'set_PRE', 'ver_PRE', 'vir_PRE',
                           'set_ACC', 'ver_ACC', 'vir_ACC')
  colnames(perf_meas_test) <- c('set_SPE', 'ver_SPE', 'vir_SPE',
                                'set_SEN', 'ver_SEN', 'vir_SEN',
                                'set_PRE', 'ver_PRE', 'vir_PRE',
                                'set_ACC', 'ver_ACC', 'vir_ACC')
  
  #find average of SPE, SEN, PRE, ACC for gamma=0.1
  ###mean(perf_meas[1:100,1:3]) doesn't work
  spe1 <- (mean(perf_meas[1:100, 1]) + mean(perf_meas[1:100, 2]) + mean(perf_meas[1:100, 3]))/3
  sen1 <- (mean(perf_meas[1:100, 4]) + mean(perf_meas[1:100, 5]) + mean(perf_meas[1:100, 6]))/3
  pre1 <- (mean(perf_meas[1:100, 7]) + mean(perf_meas[1:100, 8]) + mean(perf_meas[1:100, 9]))/3
  acc1 <- (mean(perf_meas[1:100, 10]) + mean(perf_meas[1:100, 11]) + mean(perf_meas[1:100, 12]))/3
  spe1_t <- (mean(perf_meas_test[1:100, 1]) + mean(perf_meas_test[1:100, 2]) + mean(perf_meas_test[1:100, 3]))/3
  sen1_t <- (mean(perf_meas_test[1:100, 4]) + mean(perf_meas_test[1:100, 5]) + mean(perf_meas_test[1:100, 6]))/3
  pre1_t <- (mean(perf_meas_test[1:100, 7]) + mean(perf_meas_test[1:100, 8]) + mean(perf_meas_test[1:100, 9]))/3
  acc1_t <- (mean(perf_meas_test[1:100, 10]) + mean(perf_meas_test[1:100, 11]) + mean(perf_meas_test[1:100, 12]))/3
  
  #find average of SPE, SEN, PRE, ACC for gamma=0.3
  spe3 <- (mean(perf_meas[101:200, 1]) + mean(perf_meas[101:200, 2]) + mean(perf_meas[101:200, 3]))/3
  sen3 <- (mean(perf_meas[101:200, 4]) + mean(perf_meas[101:200, 5]) + mean(perf_meas[101:200, 6]))/3
  pre3 <- (mean(perf_meas[101:200, 7]) + mean(perf_meas[101:200, 8]) + mean(perf_meas[101:200, 9]))/3
  acc3 <- (mean(perf_meas[101:200, 10]) + mean(perf_meas[101:200, 11]) + mean(perf_meas[101:200, 12]))/3
  spe3_t <- (mean(perf_meas_test[101:200, 1]) + mean(perf_meas_test[101:200, 2]) + mean(perf_meas_test[101:200, 3]))/3
  sen3_t <- (mean(perf_meas_test[101:200, 4]) + mean(perf_meas_test[101:200, 5]) + mean(perf_meas_test[101:200, 6]))/3
  pre3_t <- (mean(perf_meas_test[101:200, 7]) + mean(perf_meas_test[101:200, 8]) + mean(perf_meas_test[101:200, 9]))/3
  acc3_t <- (mean(perf_meas_test[101:200, 10]) + mean(perf_meas_test[101:200, 11]) + mean(perf_meas_test[101:200, 12]))/3
  
  #find average of SPE, SEN, PRE, ACC for gamma=0.5
  spe5 <- (mean(perf_meas[201:300, 1]) + mean(perf_meas[201:300, 2]) + mean(perf_meas[201:300, 3]))/3
  sen5 <- (mean(perf_meas[201:300, 4]) + mean(perf_meas[201:300, 5]) + mean(perf_meas[201:300, 6]))/3
  pre5 <- (mean(perf_meas[201:300, 7]) + mean(perf_meas[201:300, 8]) + mean(perf_meas[201:300, 9]))/3
  acc5 <- (mean(perf_meas[201:300, 10]) + mean(perf_meas[201:300, 11]) + mean(perf_meas[201:300, 12]))/3
  spe5_t <- (mean(perf_meas_test[201:300, 1]) + mean(perf_meas_test[201:300, 2]) + mean(perf_meas_test[201:300, 3]))/3
  sen5_t <- (mean(perf_meas_test[201:300, 4]) + mean(perf_meas_test[201:300, 5]) + mean(perf_meas_test[201:300, 6]))/3
  pre5_t <- (mean(perf_meas_test[201:300, 7]) + mean(perf_meas_test[201:300, 8]) + mean(perf_meas_test[201:300, 9]))/3
  acc5_t <- (mean(perf_meas_test[201:300, 10]) + mean(perf_meas_test[201:300, 11]) + mean(perf_meas_test[201:300, 12]))/3
  
  compare <- data.frame(c('Specificity', 'Sensitivity', 'Precision', 'Accuracy'),
                        c(spe1, sen1, pre1, acc1),
                        c(spe3, sen3, pre3, acc3),
                        c(spe5, sen5, pre5, acc5))
  colnames(compare) <- c('Statistic', 'Gamma=0.1', 'Gamma=0.3', 'Gamma=0.5')
  compare_test <- data.frame(c('Specificity', 'Sensitivity', 'Precision', 'Accuracy'),
                             c(spe1_t, sen1_t, pre1_t, acc1_t),
                             c(spe3_t, sen3_t, pre3_t, acc3_t),
                             c(spe5_t, sen5_t, pre5_t, acc5_t))
  colnames(compare_test) <- c('Statistic', 'Gamma=0.1', 'Gamma=0.3', 'Gamma=0.5')
  
  train_meas <- rbind(train_meas, compare)
  test_meas <- rbind(test_meas, compare_test)
  
  if (x==1) { #if in loop for sepal length ME, record iris_sl_me_0.1_rep_3.csv and iris_sl_me_0.5_rep_3.csv results
    sepal_01 <- perf_meas_test[3,]
    sepal_05 <- perf_meas_test[203,]
  }
}


##### Print Results #####
print('Sepal Length Measurement Error Performance')
test_meas[1:4,]
print('Petal Length Measurement Error Performance')
test_meas[5:8,]

#plots
x=c(0.1, 0.3, 0.5)
y_spe =unlist(test_meas[1,-1], use.names = FALSE)
yp_spe = unlist(test_meas[5,-1], use.names = FALSE)
y_sen =unlist(test_meas[2,-1], use.names = FALSE)
yp_sen = unlist(test_meas[6,-1], use.names = FALSE)
y_pre =unlist(test_meas[3,-1], use.names = FALSE)
yp_pre = unlist(test_meas[7,-1], use.names = FALSE)
y_acc =unlist(test_meas[4,-1], use.names = FALSE)
yp_acc = unlist(test_meas[8,-1], use.names = FALSE)
par(mfrow=c(2,2))
plot(x, y_spe,
     main='Specificity Over Gamma',
     xlab='Gamma Value', ylab='Specificity',
     col='red', ylim=c(0.9, 0.95))
points(x, yp_spe, col='blue')
plot(x, y_sen,
     main='Sensitivity Over Gamma',
     xlab='Gamma Value', ylab='Sensitivity',
     col='red', ylim=c(0.8, 0.9))
points(x, yp_sen, col='blue')
plot(x, y_pre,
     main='Precision Over Gamma',
     xlab='Gamma Value', ylab='Precision',
     col='red', ylim=c(0.8, 0.9))
points(x, yp_pre, col='blue')
plot(x, y_acc,
     main='Accuracy Over Gamma',
     xlab='Gamma Value', ylab='Accuracy',
     col='red', ylim=c(0.85, 0.95))
points(x, yp_acc, col='blue')


##### Verify Accuracy By Comparing with 21108_PCA_ME_Example.R #####
#expected for gamma = 0.5:
#  specificity = [1.0000000, 0.8750000, 0.9583333]
#  sensitivity = [1.0000000, 0.9166667, 0.7500000]
#  precision = [1.0000000, 0.7857143, 0.9000000]
#actual
sepal_05

#expected for gamma = 0.1:
#  specificity = [1.0000000, 0.9166667, 0.9583333]
#  sensitivity = [1.0000000, 0.9166667, 0.8333333]
#  precision = [1.0000000, 0.8461538, 0.9090909]
#actual
sepal_01

