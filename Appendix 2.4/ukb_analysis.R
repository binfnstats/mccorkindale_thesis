# UKB for github

library(readr)
library(data.table)
library(tidyverse)
library(ukbtools)
library(arrow)

# UKB data check

setwd("/Users/andrew/Dropbox (Sydney Uni)/Chapter4_apoe_structure/UKB")


# replace with the UKB project number eg ukb12345
# very time-consuming and requires all ukb files to be in the directory
dat_full <- ukb_df("ukbxxxxx")
# 
# Use volume scaling factor to screen for complete T1 data and to reduce dataframe size for subsequent sets
dat_t1 <- dat_full[complete.cases(dat_full$volumetric_scaling_from_t1_head_image_to_standard_space_f25000_2_0),]

# select individuals with repeat imaging (only difference is "_3_0" instead of "_2_0" at the end of the field)
# ukb_image2 <- dat_test[complete.cases(dat_test$volumetric_scaling_from_t1_head_image_to_standard_space_f25000_3_0),]
# 
# optionally save data from all fields only for individuals with neuroimaging data 
saveRDS(dat_t1,"ukbxxxxx_with_MRI.rds")

# read in apoe genotyping information and then only keep individuals with genotyping
gentype <- read.csv("full_apoe_genotype.csv")
gentype <- gentype[complete.cases(gentype$Genotype),]

# create a dummy variable for APOE4 carrier or not
# gentype$apoe4[gentype$Genotype=="e4/e4"] <- 1
# gentype$apoe4[gentype$Genotype=="e3/e4"] <- 1
# gentype$apoe4[gentype$Genotype=="e2/e4"] <- 1
# gentype$apoe4[gentype$Genotype=="e1/e4"] <- 1
# gentype$apoe4[gentype$Genotype=="e3/e3"] <- 0
# gentype$apoe4[gentype$Genotype=="e2/e3"] <- 0
# gentype$apoe4[gentype$Genotype=="e2/e2"] <- 0
# gentype$apoe4[gentype$Genotype=="e1/e2"] <- 0

#remove the "/" symbol
# gentype$apoe_txt[gentype$Genotype=="e4/e4"] <- "e4e4"
# gentype$apoe_txt[gentype$Genotype=="e3/e4"] <- "e3e4"
# gentype$apoe_txt[gentype$Genotype=="e2/e4"] <- "e2e4"
# gentype$apoe_txt[gentype$Genotype=="e1/e4"] <- "e1e4"
# gentype$apoe_txt[gentype$Genotype=="e3/e3"] <- "e3e3"
# gentype$apoe_txt[gentype$Genotype=="e2/e3"] <- "e2e3"
# gentype$apoe_txt[gentype$Genotype=="e2/e2"] <- "e2e2"
# gentype$apoe_txt[gentype$Genotype=="e1/e2"] <- "e1e2"

# now attach apoe genotyping to the two datasets

# remove the apoe genotyping subjects that do not have neuroimaging data
gentype_im1 <- gentype[gentype$eid %in% dat$eid,]

# sanity check - put the files in the same order and confirm that the eid field
# is identical in the genotype and main dataframes
# gentype_im1 <- gentype_im1[order(gentype_im1$eid),]
# ukb_image1 <- dat_t1[order(dat_t1$eid),]
# identical(dat_t1$eid,gentype_im1$eid)

# merge genotype and main dataframes together
ukb_image1_gen <- merge(gentype_im1,ukb_image1,by="eid")

# select ICD10 codes for dementia and then identify any individual that matches those fields
# for more efficiency this can be run only in the icd10 fields
icd10_list <- list_int$dementia_fields[nzchar(list_int$dementia_fields)]
df_dementia <- dat_icd %>% filter_all(any_vars(. %in% c(icd10_list)))

# 65 with a dementia diagnosis - removed

ukb_image1_gen <- ukb_image1_gen[!ukb_image1_gen$eid %in% df_dementia$eid,]

# only keep imaging variables of interest. This list was curated using all neuroimaging
# variables available in the data showcase
ukb_imonly <- ukb_image1gen[,names(ukb_image1gen) %in% vars_int$all_imaging]

# find the amount of missing data for each variable and save as a list
missing_data <- colSums(is.na(ukb_imonly))
write.csv(missing_data, "UKB_missing_imaging_data.csv")

# exclude variables with too much missing data

dat_final <- ukb_imonly[,names(ukb_imonly) %in% vars_int$vars_to_keep]

# final dataframe dimensions were 33,384 rows and 2115 columns
# 33384  2115
# removed columns as well as final columns are in supplementary tables (or available on request)

##########################################################################################
# ukb brain mri normalisation
##########################################################################################

# get list of variables to normalise
tonorm_list <- list_int$imaging_tonorm2_0
tonorm_list <- tonorm_list[nzchar(tonorm_list)]

# for loop that generates the code to normalise each variable
# for (ctr in c_list) {
#   cat(paste('dat_final$',ctr,' <- dat_final$',ctr,' * dat_final$volumetric_scaling_from_t1_head_image_to_standard_space_f25000_2_0 \n',sep=""))
# }

# run output to normalise each variable

# recode variables/convert to factors as needed
dat_fin$sex <- as.factor(dat_fin$sex_f31_0_0)
dat_fin$apoe4d <- as.factor(dat_fin$apoe4)
dat_fin$apoe_txt <- as.factor(dat_fin$apoe_txt)
dat_fin$age <- dat_fin$age_when_attended_assessment_centre_f21003_2_0
dat_fin$age_dec[dat_fin$age>40 & dat_fin$age<50] <- "40_49"
dat_fin$age_dec[dat_fin$age>49.999 & dat_fin$age<60] <- "50_59"
dat_fin$age_dec[dat_fin$age>59.999 & dat_fin$age<70] <- "60_69"
dat_fin$age_dec[dat_fin$age>69.999 & dat_fin$age<80] <- "70_79"
dat_fin$age_dec[dat_fin$age>79.999 & dat_fin$age<90] <- "80_89"
dat_fin$age_dec <- as.factor(dat_fin$age_dec)
dat_fin$education <- as.factor(dat_fin$qualifications_f6138_0_0)
dat_fin$apoe4d_txt[dat_fin$apoe4d==0] <- "Control"
dat_fin$apoe4d_txt[dat_fin$apoe4d==1] <- "APOE4"
dat_fin$apoe4d_txt <- as.factor(dat_fin$apoe4d_txt)
# get summaries of key demographic variables

# this file used as input for machine learning
saveRDS(dat_ml,"ukb_imaging_input_with_covars.rds")



#######################################################################################################
# Boruta analysis
#######################################################################################################
#

library(tidyverse)
#library(caret)
library(Boruta)
set.seed(123)

full_dat <- readRDS("ukb_input_with_freesurfer.rds")

# to do an analysis using the number of APOE4 copies as output
full_dat$e4_copies[full_dat$apoe_txt=="e2e2"] <- 0
full_dat$e4_copies[full_dat$apoe_txt=="e2e3"] <- 0
full_dat$e4_copies[full_dat$apoe_txt=="e2e4"] <- 1
full_dat$e4_copies[full_dat$apoe_txt=="e3e3"] <- 0
full_dat$e4_copies[full_dat$apoe_txt=="e3e4"] <- 1
full_dat$e4_copies[full_dat$apoe_txt=="e4e4"] <- 2
full_dat$e4_copies <- as.factor(full_dat$e4_copies)


# Exclude unwanted genotypes from analysis. In these analyses the APOE33 group
# was used as the control group 
full_dat <- full_dat[!full_dat$apoe_txt=="e2e2",]
full_dat <- full_dat[!full_dat$apoe_txt=="e2e4",]
full_dat <- full_dat[!full_dat$apoe_txt=="e2e3",]

# create a new variable to split cohort into those above and below 65 to control for age
full_dat$over_under_65[full_dat$age<65.00] <- "Under_65"
full_dat$over_under_65[full_dat$age>64.9999999] <- "Over_65"
full_dat$over_under_65 <- as.factor(full_dat$over_under_65)

# dataframes for the younger and older cohorts
dat_y <- full_dat[full_dat$over_under_65=="Under_65",]
dat_o <- full_dat[full_dat$over_under_65=="Over_65",]

#######################################################################################################
# female young
#######################################################################################################


# only younger dataset as input
dat <- dat_y

# excude males
dat <- dat[!dat$sex=="Male",]

dat$Target <- as.factor(dat$e4_copies)
# or use apoe4d variable to do a classification based on APOE4 status
# dat$Target <- as.factor(dat$apoe4d_txt)

# remove unwanted variables
dat$e4_copies <- NULL
dat$age <- NULL
dat$age_dec <- NULL
dat$apoe_txt <- NULL
dat$apoe4d <- NULL
dat$sex <- NULL
dat$apoe4d_txt <- NULL
dat$eid <- NULL
dat$volumetric_scaling_from_t1_head_image_to_standard_space_f25000_2_0 <- NULL
dat$volume_of_peripheral_cortical_grey_matter_normalised_for_head_size_f25001_2_0 <- NULL
dat$volume_of_ventricular_cerebrospinal_fluid_normalised_for_head_size_f25003_2_0 <- NULL
dat$volume_of_grey_matter_normalised_for_head_size_f25005_2_0 <- NULL
dat$volume_of_white_matter_normalised_for_head_size_f25007_2_0 <- NULL
dat$volume_of_brain_greywhite_matter_normalised_for_head_size_f25009_2_0 <- NULL

data_ukb_under65 <- dat

#covar_ukb_under65 <- readRDS("apoe_covar.rds")

labels_ukb_under65 <- dat$Target
#setwd(paste(path,"/boruta_test_ADDEM_NCI", sep=""))

rm_list <- ("Target")
data_ukb_under65 <- data_ukb_under65[,!(names(data_ukb_under65) %in% rm_list)]

#apoe
data_labels_ukb_under65 <- cbind(data_ukb_under65, labels_ukb_under65)
# altering maxRuns dramatically affects run time but lower numbers can affect stability. Make sure that the
# number of threads is set to an appropriate number for your machine otherwise Boruta will use all available cores
boruta_ukb_under65_model=Boruta(x=data_ukb_under65 , y=labels_ukb_under65, doTrace = 2 , maxRuns = 2000,num.threads=10)
saveRDS(boruta_ukb_under65_model,file="boruta_ukb_under65_model_predictE4_FEM_FS_e4_copies.rds")
final_boruta_ukb_under65_model=TentativeRoughFix(boruta_ukb_under65_model)
apoe_model_sel_feat=getSelectedAttributes(final_boruta_ukb_under65_model)
saveRDS(apoe_model_sel_feat,file="ukb_under65_model_boruta_sel_feat_predictE4_FEM_FS_e4_copies.rds")

k <-lapply(1:ncol(boruta_ukb_under65_model$ImpHistory),function(i)
  boruta_ukb_under65_model$ImpHistory[is.finite(boruta_ukb_under65_model$ImpHistory[,i]),i])
names(k) <- colnames(boruta_ukb_under65_model$ImpHistory)
boruta_ukb_under65_model_feature_rank <- sort(sapply(k,median),decreasing = TRUE)
boruta_ukb_under65_model_rank_out <- data.frame(features=names(boruta_ukb_under65_model_feature_rank), z_score=boruta_ukb_under65_model_feature_rank, row.names=NULL)

# export a csv file with all input variables ranked by z-score
write.csv(boruta_ukb_under65_model_rank_out,file="boruta_ukb_under65_model_feature_rank_predict_E4_FEM_FS_e4_copies.csv",row.names = FALSE)

#######################################################################################################
# Male old
#######################################################################################################
# 
# only older dataset as input
dat <- dat_o

# remove females
dat <- dat[!dat$sex=="Female",]

dat$Target <- as.factor(dat$e4_copies)

# remove unwanted variables
dat$e4_copies <- NULL
dat$age <- NULL
dat$age_dec <- NULL
dat$apoe_txt <- NULL
dat$apoe4d <- NULL
dat$sex <- NULL
dat$apoe4d_txt <- NULL
dat$eid <- NULL
dat$volumetric_scaling_from_t1_head_image_to_standard_space_f25000_2_0 <- NULL
dat$volume_of_peripheral_cortical_grey_matter_normalised_for_head_size_f25001_2_0 <- NULL
dat$volume_of_ventricular_cerebrospinal_fluid_normalised_for_head_size_f25003_2_0 <- NULL
dat$volume_of_grey_matter_normalised_for_head_size_f25005_2_0 <- NULL
dat$volume_of_white_matter_normalised_for_head_size_f25007_2_0 <- NULL
dat$volume_of_brain_greywhite_matter_normalised_for_head_size_f25009_2_0 <- NULL


data_ukb_over65 <- dat

#covar_ukb_over65 <- readRDS("apoe_covar.rds")

labels_ukb_over65 <- dat$Target
#setwd(paste(path,"/boruta_test_ADDEM_NCI", sep=""))

rm_list <- ("Target")
data_ukb_over65 <- data_ukb_over65[,!(names(data_ukb_over65) %in% rm_list)]

#apoe
data_labels_ukb_over65 <- cbind(data_ukb_over65, labels_ukb_over65)
boruta_ukb_over65_model=Boruta(x=data_ukb_over65 , y=labels_ukb_over65, doTrace = 2 , maxRuns = 2000,num.threads=35)
saveRDS(boruta_ukb_over65_model,file="boruta_ukb_over65_model_predictE4_MAL_FS_e4_copies.rds")
final_boruta_ukb_over65_model=TentativeRoughFix(boruta_ukb_over65_model)
apoe_model_sel_feat=getSelectedAttributes(final_boruta_ukb_over65_model)
saveRDS(apoe_model_sel_feat,file="ukb_over65_model_boruta_sel_feat_predictE4_MAL_FS_e4_copies.rds")

k <-lapply(1:ncol(boruta_ukb_over65_model$ImpHistory),function(i)
  boruta_ukb_over65_model$ImpHistory[is.finite(boruta_ukb_over65_model$ImpHistory[,i]),i])
names(k) <- colnames(boruta_ukb_over65_model$ImpHistory)
boruta_ukb_over65_model_feature_rank <- sort(sapply(k,median),decreasing = TRUE)
boruta_ukb_over65_model_rank_out <- data.frame(features=names(boruta_ukb_over65_model_feature_rank), z_score=boruta_ukb_over65_model_feature_rank, row.names=NULL)

write.csv(boruta_ukb_over65_model_rank_out,file="boruta_ukb_over65_model_feature_rank_predict_E4_MAL_FS_e4_copies.csv",row.names = FALSE)



#######################################################################################################
# caret ML classification
#######################################################################################################
#

library(plyr)
library(tidyverse)
library(caret)
library(ranger)
library(xgboost)
library(doParallel)

cl <- makePSOCKcluster(10)
registerDoParallel(cl)

set.seed(123)

tl <- 5
ctrl <- trainControl(method = "repeatedcv",
                     repeats = 3,
                     number = 5,
                     sampling = "smote",
                     classProbs = TRUE,
                     savePredictions = TRUE,
                     allowParallel= TRUE)

ranger_grid <- expand.grid(
  mtry = c(2,5,10,20,500),
  min.node.size=1,
  splitrule="gini"
)

#######################################################################################################
# Male old example
#######################################################################################################

dat <- dat_o

dat <- dat[!dat$sex=="Female",]

dat$Target <- as.factor(dat$apoe4d_txt)

# remove unwanted variables
dat$age <- NULL
dat$age_dec <- NULL
dat$apoe_txt <- NULL
dat$apoe4d <- NULL
dat$sex <- NULL
dat$apoe4d_txt <- NULL
dat$eid <- NULL
dat$volumetric_scaling_from_t1_head_image_to_standard_space_f25000_2_0 <- NULL
dat$volume_of_peripheral_cortical_grey_matter_normalised_for_head_size_f25001_2_0 <- NULL
dat$volume_of_ventricular_cerebrospinal_fluid_normalised_for_head_size_f25003_2_0 <- NULL
dat$volume_of_grey_matter_normalised_for_head_size_f25005_2_0 <- NULL
dat$volume_of_white_matter_normalised_for_head_size_f25007_2_0 <- NULL
dat$volume_of_brain_greywhite_matter_normalised_for_head_size_f25009_2_0 <- NULL
dat$over_under_65 <- NULL


glmnet_model <- train(x = dat[, names(dat) != "Target"],
                      y = dat$Target,
                      method = "glmnet",
                      metric = "Accuracy",
                      tuneLength = tl,
                      trControl = ctrl)
save(glmnet_model, file = "glmnet_model_ukb_MALE4_old.RData")

ranger_model <- train(x = dat[, names(dat) != "Target"],
                      y = dat$Target,
                      method = "ranger",
                      tuneGrid = ranger_grid,
                      importance = "permutation",
                      trControl = ctrl,
                      num.threads = 1,
                      metric = "Accuracy")
save(ranger_model, file = "ranger_model_ukb_MALE4_old.RData")

# output feature importances
# read in models if required
load("ranger_model_ukb_MALE4_old.RData")
load("glmnet_model_ukb_MALE4_old.RData")

# output csv files with feature importance
ukb_MALE4_old_ranger_varimp <- varImp(ranger_model)
write.csv(ukb_MALE4_old_ranger_varimp$importance, "ukb_MALE4_old_ranger.csv")

ukb_MALE4_old_glmnet_varimp <- varImp(glmnet_model)
write.csv(ukb_MALE4_old_glmnet_varimp$importance, "ukb_MALE4_old_glmnet.csv")


library(caret)
library(MLeval)
# create AUC plots and get AUC stats
HRLR <- list(
  glmnet = glmnet_model,
  ranger = ranger_model)

x <- evalm(HRLR, gnames = c("glmnet","ranger"))

