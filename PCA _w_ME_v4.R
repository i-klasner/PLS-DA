# Ina Klasner
# 10/27/2021
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

files <- c()
perf_meas <- data.frame()
perf_meas_test <- data.frame()
##### Perform Analysis on Each File #####
for(x in 1:4) {
  #create an empty dataframe for performance measures
  
  counter <- 1
  if(x == 1) { #go through Sepal Length ME files first
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/iris/Sepal_Length'
  } else if (x==2) { # then go through Petal Length ME files
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/iris/Petal_Length'
  } else if (x==3) {
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/iris/Homo_ME_IRIS'
  } else {
    string <- '/Users/hobbes/Desktop/michigan_tech/research/biostats/plsda/iris/Hetero_ME_IRIS'
  }
  #iterate through folders, gamma and rep values
  for(i in seq(from=1, to=5, by=2)) { #iterate through gammas
    for(j in 1:100) { #iterate through files
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

#find average of SPE, SEN, PRE, ACC for all gammas over all datasets
compare_t <- data.frame()
for (i in 1:12) {
  for (j in c(1,4,7,10)) {
    a <- ((i-1)*100+1)
    z <- (i*100)
    tmp_t <- (mean(perf_meas_test[a:z, j]) + mean(perf_meas_test[a:z, (j+1)]) + mean(perf_meas_test[a:z, (j+2)]))/3      
    var.name.t <- paste0('param_', j, '_t')
    assign(var.name.t, tmp_t)
  }
  compare_t <- rbind(compare_t, c(param_1_t, param_4_t, param_7_t, param_10_t))
}
colnames(compare_t) <- c('Specificity', 'Sensitivity', 'Precision', 'Accuracy')
row.names(compare_t) <- c('Sepal Gam 0.1', 'Sepal Gam 0.3', 'Sepal Gam 0.5',
                          'Petal Gam 0.1', 'Petal Gam 0.3', 'Petal Gam 0.5',
                          'Homo Gam 0.1', 'Homo Gam 0.3', 'Homo Gam 0.5',
                          'Hetero Gam 0.1', 'Hetero Gam 0.3', 'Hetero Gam 0.5')
compare_t


##### Print Results #####
print('Sepal Length Measurement Error Performance')
compare_t[1:3,]
print('Petal Length Measurement Error Performance')
compare_t[4:6,]
print('Homogeneous Measurement Error Performance')
compare_t[7:9,]
print('Heterogeneous Measurement Error Performance')
compare_t[10:12,]

##### Visualize Results #####
x=c(0.1, 0.3, 0.5)
for (x in 1:4) { #create variables for parameters on each case
  spe <- unlist(compare_t[(3*x-2):(3*x),1])   
  sen <- unlist(compare_t[(3*x-2):(3*x),2])  
  pre <- unlist(compare_t[(3*x-2):(3*x),3])  
  acc <- unlist(compare_t[(3*x-2):(3*x),4])  
  name.1 <- paste0('spe', x)
  name.2 <- paste0('sen', x)
  name.3 <- paste0('pre', x)
  name.4 <- paste0('acc', x)
  assign(name.1, spe)
  assign(name.2, sen)
  assign(name.3, pre)
  assign(name.4, acc)
}

#plotting
par(mfrow=c(2,2))
plot(x, spe1,
     main='Specificity Over Gamma',
     xlab='Gamma Value', ylab='Specificity',
     col='red', ylim=c(0.87, 0.95))
points(x, spe2, col='blue')
points(x, spe3, col='orange')
points(x, spe4, col='forest green')
legend("bottomleft",   # Position
       inset=0.05,
       legend = c("Sepal", "Petal", 'Homo', 'Hetero'),
       fill = c('red', 'blue', 'orange', 'forest green'), cex=0.5 )
plot(x, sen1,
     main='Sensitivity Over Gamma',
     xlab='Gamma Value', ylab='Sensitivity',
     col='red', ylim=c(0.77, 0.9))
points(x, sen2, col='blue')
points(x, sen3, col='orange')
points(x, sen4, col='forest green')
plot(x, pre1,
     main='Precision Over Gamma',
     xlab='Gamma Value', ylab='Precision',
     col='red', ylim=c(0.78, 0.9))
points(x, pre2, col='blue')
points(x, pre3, col='orange')
points(x, pre4, col='forest green')
plot(x, acc1,
     main='Accuracy Over Gamma',
     xlab='Gamma Value', ylab='Accuracy',
     col='red', ylim=c(0.85, 0.93))
points(x, acc2, col='blue')
points(x, acc3, col='orange')
points(x, acc4, col='forest green')

# red is ME on sepal length
# blue is ME on petal length with greater variance
# orange is homogeneous case
# green is heterogeneous case

var(iris$Sepal.Length)
var(iris$Petal.Length) #greater variance for petal length, greater effect of change in gamma
