# Run classification analysis using 5 algorithms

library(plyr)
library(tidyverse)
library(caret)
library(Boruta)
#library(blkbox)
library(doParallel)

registerDoParallel(30)

set.seed(123)

path <- "/home/user/data"

setwd(path)


ctrl <- trainControl(method = "repeatedcv",
                     repeats = 20,
                     number = 5,
                     classProbs = TRUE,
                     savePredictions = "final",
                     allowParallel=TRUE,
                     returnData=FALSE)

tl <- 30
#tl for xgbtree shorter because it has 7 tuning parameters compared to 1-4 for others
tlx <- 10


############################################################################################


#ctrl_vs_ad
#
setwd(paste(path,"/ctrl_vs_ad", sep=""))

# labels and data files produced by the "prepare_for_m_learn" script
labels_ctrl_vs_ad = readRDS("labels_ctrl_vs_ad.rds")
data_ctrl_vs_ad = readRDS("data_ctrl_vs_ad.rds")
data_ctrl_vs_ad <- cbind(data_ctrl_vs_ad, labels_ctrl_vs_ad)

#
rf_model <- train(x = data_ctrl_vs_ad[, names(data_ctrl_vs_ad) != "labels_ctrl_vs_ad"],
                  y = data_ctrl_vs_ad$labels_ctrl_vs_ad,
                  method = "rf",
                  tuneLength = tl,
                  trControl = ctrl,
                  metric = "Kappa")
# save models
saveRDS(rf_model , "caret_rf_model_ctrl_vs_ad.rds")

rpart_model <- train(x = data_ctrl_vs_ad[, names(data_ctrl_vs_ad) != "labels_ctrl_vs_ad"],
                     y = data_ctrl_vs_ad$labels_ctrl_vs_ad,
                     method = "rpart",
                     tuneLength = tl,
                     trControl = ctrl,
                     metric = "Kappa")
# save models
saveRDS(rpart_model , "caret_rpart_model_ctrl_vs_ad.rds")

glmnet_model <- train(x = data_ctrl_vs_ad[, names(data_ctrl_vs_ad) != "labels_ctrl_vs_ad"],
                      y = data_ctrl_vs_ad$labels_ctrl_vs_ad,
                      method = "glmnet",
                      tuneLength = tl,
                      trControl = ctrl,
                      metric = "Kappa")
# save models
saveRDS(glmnet_model, "caret_glmnet_model_ctrl_vs_ad.rds")

ranger_model <- train(x = data_ctrl_vs_ad[, names(data_ctrl_vs_ad) != "labels_ctrl_vs_ad"],
                      y = data_ctrl_vs_ad$labels_ctrl_vs_ad,
                      method = "ranger",
                      tuneLength = tl,
                      trControl = ctrl,
                      metric = "Kappa")
# save models
saveRDS(ranger_model, "caret_ranger_model_ctrl_vs_ad.rds")


xgbTree_model <- train(x = data_ctrl_vs_ad[, names(data_ctrl_vs_ad) != "labels_ctrl_vs_ad"],
                       y = data_ctrl_vs_ad$labels_ctrl_vs_ad,
                       method = "xgbTree",
                       tuneLength = tlx,
                       trControl = ctrl,
                       metric = "Kappa",
                       nthread=1)
## nthread =1 necessary because xgbtree is multi-threaded by default

# save models
saveRDS(xgbTree_model, "caret_xgbtree_model_ctrl_vs_ad.rds")

############################################################################################



