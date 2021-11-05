# Ina Klasner
# 11/01/2021
# PLS-DA Practice with Spectral Example

##### Add Packages and Data #####
library('pls')
data('yarn')
data('oliveoil')
data('gasoline')


##### PLSR on Gasoline Data #####
gas_train <- gasoline[1:50,] #training set
gas_test <- gasoline[51:60,] # testing set

# fit plsr model with leave-one-out cross-validated predictions
gas1 <- plsr(octane ~ NIR, ncomp=10, data=gas_train, validation='LOO')
summary(gas1)

# plot root mean squared error of prediction
plot(RMSEP(gas1), legendpos='topright')
# determine 2 components is sufficient

# inspect fit
plot(gas1, ncomp=2, asp=1, line=TRUE) #prediction plot, no curvature or outliers
plot(gas1, plottype='scores', comps=1:3) #pairwise plot of scores for first 3 comps
explvar(gas1) #explained variance from each component
plot(gas1, 'loadings', comps=1:2, legendpos='topleft', labels='numbers', xlab='nm')
abline(h=0) #loading plot

# predict response variables of new observations
predict(gas1, ncomp=2, newdata=gas_test)
predplot(gas1, ncomp=2, newdata=gas_test, asp=1, line=TRUE) #plotted

# test set RMSEP
RMSEP(gas1, newdata=gas_test) #0.2445 close to 0.297 from train

# inspect the fitted model
plot(gas1, plottype='coef', ncomp=1:3, legendpos='bottomleft',
     labels='numbers', xlab='nm') #visualize regression coeffs, little diff from 2 to 3
plot(gas1, plottype='correlation', comps=1:3)

# extract regression coefficients
coef(gas1, ncomp=2)

# extract scores and % variance explained
scores(gas1)
Yscores(gas1)

# extract loadings
loadings(gas1)
Yloadings(gas1)

# extract weights
loading.weights(gas1)

# fit plsr with multiplicative scatter correction
gas2 <- plsr(octane ~ msc(NIR), ncomp=10, data=gas_train)
predict(gas2, ncomp=3, newdata=gas_test) #new spectra auto-scatter corrected

# choose number components with cross-validation
gas2.cv <- crossval(gas2, segments=10) #perform cross-validation
plot(MSEP(gas2.cv), legendpos='topright') #nearly identical
summary(gas2.cv, what='validation')


##### PLSR on Olive Oil Data #####
# fit plsr with standardized variables (divide by sd)
olive1 <- plsr(sensory~chemical, scale=TRUE, data=oliveoil)

# plot root mean squared error of prediction
plot(RMSEP(olive1), legendpos='topright')

# inspect fit
plot(olive1, ncomp=3, asp=1, line=TRUE) #prediction plot, no curvature or outliers
plot(olive1, plottype='scores', comps=1:4) #pairwise plot of scores for first 3 comps
explvar(olive1) #explained variance from each component

# inspect the fitted model
plot(olive1, plottype='correlation', comps=1:3)
