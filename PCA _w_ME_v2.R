# Ina Klasner
# 10/10/2021
# PCA With Measurement Error & iris data

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

#create an empty dataframe for performance measures

perf_meas <- data.frame()
perf_meas_test <- data.frame()

##### Perform Analysis on Each File #####
counter <- 1
#iterate through folders, gamma and rep values
for(i in c(1,3,5)) { #iterate through gammas
  for(j in 1:100) { #iterate through files
    directory <- paste0('/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda', '/ME0', i, '/')
    data <- read.csv(paste0(directory, file_names(i, j))) #read in data file
    
    #create training and testing sets
    set.seed(100)
    i <- createDataPartition(data$Sepal.Length, p=0.75, list=FALSE)
    training <- data[i,] #create a training data set of 75% of the original data
    testing <- data[-i,] #create a testing set of all other data
    train_i <- iris[i,] #use species classifications from iris data
    test_i <- iris[-i,]
    
    #pca
    pca <- prcomp(training[,1:4],center=TRUE, scale.=TRUE)
    #pca <- principal(training[,1:4], nfactors=2)
    predict <- predict(pca,training)
    train <- data.frame(predict, train_i[,5])
    predict_1 <- predict(pca, testing)
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
    for (i in 1:3) {
      #training
      tp <- conf_matrix$table[i,i]
      fp <- sum(conf_matrix$table[i,-i])
      fn <- sum(conf_matrix$table[-i,i])
      tn <- sum(conf_matrix$table[-i,-i])
      acc <- (tp+tn)/(tp+fp+fn+tn)
      sp_acc[i] <- acc
      
      #testing
      tp_t <- conf_matrix_test$table[i,i]
      fp_t <- sum(conf_matrix_test$table[i,-i])
      fn_t <- sum(conf_matrix_test$table[-i,i])
      tn_t <- sum(conf_matrix_test$table[-i,-i])
      acc_t <- (tp_t+tn_t)/(tp_t+fp_t+fn_t+tn_t)
      sp_acc_test[i] <- acc_t
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
  
compare
compare_test

##### Tester #####
data <- read.csv('/Users/hobbes/Desktop/michigan_tech/research/plsda/ME01/iris_sl_me_0.1_rep_1.csv')
set.seed(100)
i <- createDataPartition(data$Sepal.Length, p=0.75, list=FALSE)
training <- data[i,] #create a training data set of 75% of the original data
testing <- data[-i,] #create a testing set of all other data
train_i <- iris[i,]
test_i <- iris[-i,]
pca <- prcomp(training[,1:4],center=TRUE, scale.=TRUE)
predict <- predict(pca,training)
train <- data.frame(predict, train_i[,5])
predict_1 <- predict(pca, testing)
test <- data.frame(predict_1, test_i[,5])

set.seed(100)
model <- multinom(Species~PC1 + PC2, data=train)

pred <- predict(model, train)
conf_matrix <- confusionMatrix(pred, train[,5]) 
#conf_mat[counter] <- conf_matrix
#conf <- append(conf, conf_matrix)

stats <- data.frame(conf_matrix$byClass)
sp_acc <- c(0,0,0)
for (i in 1:3) {
  tp <- conf_matrix$table[i,i]
  fp <- sum(conf_matrix$table[i,-i])
  fn <- sum(conf_matrix$table[-i,i])
  tn <- sum(conf_matrix$table[-i,-i])
  acc <- (tp+tn)/(tp+fp+fn+tn)
  sp_acc[i] <- acc
}
