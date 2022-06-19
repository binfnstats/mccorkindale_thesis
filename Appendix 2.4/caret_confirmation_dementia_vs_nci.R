library(plyr)
library(tidyverse)
library(caret)
#library(blkbox)
library(doParallel)

registerDoParallel(30)

set.seed(123)

# file with genes and labels for dementia and NCI
dat <- readRDS("/home/DEM_NCI_genes_labels.rds")

# here $Target is a column with their clinical group - DEM or NCI
dat$Target <- dat$labels
dat$labels <- NULL


path <- "/home/yourdata"
setwd(path)

# set tune length for glmnet
tl <- 10
# set sampling for caret
ctrl <- trainControl(method = "repeatedcv",
                     repeats = 20,
                     number = 5,
                     classProbs = TRUE,
                     savePredictions = TRUE)


ranger_grid <- expand.grid(
  mtry = c(5,5000,10000,15000,20494),
  min.node.size=1,
  splitrule="gini"
)

# tuning grid for xgbtree - better to specify a grid for xgbtree as it will otherwise
# run for a very long time 
tune_grid <- expand.grid(
  nrounds = c(500,1000,1500),
  eta = c(0.1,0.3),
  max_depth = c(3,5,7,10),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)



glmnet_model <- train(x = dat[, names(dat) != "Target"],
                      y = dat$Target,
                      method = "glmnet",
                      metric = "Accuracy",
                      importance = "permutation",
                      tuneLength = tl,
                      trControl = ctrl)
save(glmnet_model, file = "glmnet_model_DEMNCI_fulldata_5fold.RData")



ranger_model <- train(x = dat[, names(dat) != "Target"],
                      y = dat$Target,
                      method = "ranger",
                      importance = "permutation",
                      tuneGrid = ranger_grid,
                      trControl = ctrl,
                      metric = "Accuracy")
save(ranger_model, file = "ranger_model_DEMNCI_fulldata_5fold.RData")

# xgboost
xgb_model <- train(x = dat[, names(dat) != "Target"],
                   y = dat$Target,
                   method = "xgbTree",
                   tuneGrid = tune_grid,
                   importance = "permutation",
                   metric = "Accuracy",
                   trControl = ctrl,
                   nthread=1)

save(xgb_model, file = "xgb_model_DEMNCI_fulldata_5fold.RData")


# load models and extract feature importance
load("ranger_model_DEMNCI_fulldata_5fold.RData")
load("glmnet_model_DEMNCI_fulldata_5fold.RData")
load("xgb_model_DEMNCI_fulldata_5fold.RData")

### dem vs nci
demnci_ranger_varimp <- varImp(ranger_model)
write.csv(demnci_ranger_varimp$importance, "demnci_ranger.csv")

demnci_xgb_varimp <- varImp(xgb_model)
write.csv(demnci_xgb_varimp$importance, "demnci_xgb.csv")

demnci_glmnet_varimp <- varImp(glmnet_model)
write.csv(demnci_glmnet_varimp$importance, "demnci_glmnet.csv")

