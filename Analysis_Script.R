#A novel combination of host protein biomarkers to distinguish bacterial from viral infections in febrile children in emergency care
#Chantal D. Tan*, Bryan van den Broek*, Rebecca S. Womersley*, Myrsini Kaforou, Nienke N. Hagedoorn, Michiel van der Flier, Heather Jackson, Henriette A. Moll, Rozemarijn Snijder, Marien I. de Jonge, Clementien L. Vermont

#Analysis Script
#Rebecca S. Womersley

data <- read.csv("Proteomics_data.csv", row.names = 1)
data$Phenotype <- as.factor(data$Phenotype)

crp_subset <- (data$CRP <= 60)

## modelling ---------------------------------------------------------
#dividing into training and testing set
set.seed(50)
part <- createDataPartition(data$Phenotype, p = 0.7, list = FALSE)
trainData <- data[part,1:8]
testData <- data[-part,1:8]


#check the phenotype distribution between those datasets
summary(trainData$Phenotype) #47 DB, 25 DV
summary(testData$Phenotype) #20 DB, 10 DV

# near zero-variance predictors? ------------------------------------------------------------
nearZeroVar(trainData) #checks uniqueness and skewness in predictors. 
nearZeroVar(testData)



# preprocessing -----------------------------------------------------------------------------
preproc <- preProcess(trainData, method = c("center", "scale"))
trainData <- predict(preproc, trainData)
testData <- predict(preproc, testData)



# 6) Building Models -----------------------------------------------------------------------------------------------------
# specifying the resampling ---------------------------------------------------------------
cv10 <- trainControl(method = "repeatedcv", 
                     number = 10, 
                     repeats = 10,
                     summaryFunction = multiClassSummary,
                     allowParallel = TRUE,
                     classProbs = TRUE,
                     savePredictions = "final")


# first - testing ability of each protein individually. -------------------------
#lets use GLM 
set.seed(50) 
ip10 <- train(Phenotype ~ IP_10, data=trainData, method='glm', metric = "ROC", trControl = cv10) #choosing our resampling method
predicted_ip10<- predict(ip10, newdata = testData, type = "prob") # get probabilities of classifications
predicted_ip10$predicted <- as.factor(ifelse(predicted_ip10$DB >= 0.5, "DB", "DV"))
predicted_ip10$correct <- as.factor(testData$Phenotype)

set.seed(50) 
trail <- train(Phenotype ~ TRAIL, data=trainData, method='glm', metric = "ROC", trControl = cv10) #choosing our resampling method
predicted_trail <- predict(trail, newdata = testData, type = "prob") # get probabilities of classifications
predicted_trail$predicted <- as.factor(ifelse(predicted_trail$DB >= 0.5, "DB", "DV"))
predicted_trail$correct <- as.factor(testData$Phenotype)

set.seed(50) 
ifng <- train(Phenotype ~ IFNG, data=trainData, method='glm', metric = "ROC", trControl = cv10) #choosing our resampling method
predicted_ifng <- predict(ifng, newdata = testData, type = "prob") # get probabilities of classifications
predicted_ifng$predicted <- as.factor(ifelse(predicted_ifng$DB >= 0.5, "DB", "DV"))
predicted_ifng$correct <- as.factor(testData$Phenotype)

set.seed(50) 
il4 <- train(Phenotype ~ IL_4, data=trainData, method='glm', metric = "ROC", trControl = cv10) #choosing our resampling method
predicted_il4 <- predict(il4, newdata = testData, type = "prob") # get probabilities of classifications
predicted_il4$predicted <- as.factor(ifelse(predicted_il4$DB >= 0.5, "DB", "DV"))
predicted_il4$correct <- as.factor(testData$Phenotype)

set.seed(50) 
pct <- train(Phenotype ~ Procalcitonin, data=trainData, method='glm', metric = "ROC", trControl = cv10) #choosing our resampling method
predicted_pct <- predict(pct, newdata = testData, type = "prob") # get probabilities of classifications
predicted_pct$predicted <- as.factor(ifelse(predicted_pct$DB >= 0.5, "DB", "DV"))
predicted_pct$correct <- as.factor(testData$Phenotype)

set.seed(50) 
lcn2 <- train(Phenotype ~ LCN2, data=trainData, method='glm', metric = "ROC", trControl = cv10) #choosing our resampling method
predicted_lcn2 <- predict(lcn2, newdata = testData, type = "prob") # get probabilities of classifications
predicted_lcn2$predicted <- as.factor(ifelse(predicted_lcn2$DB >= 0.5, "DB", "DV"))
predicted_lcn2$correct <- as.factor(testData$Phenotype)

set.seed(50) 
il6 <- train(Phenotype ~ IL_6, data=trainData, method='glm', metric = "ROC", trControl = cv10) #choosing our resampling method
predicted_il6 <- predict(il6, newdata = testData, type = "prob") # get probabilities of classifications
predicted_il6$predicted <- as.factor(ifelse(predicted_il6$DB >= 0.5, "DB", "DV"))
predicted_il6$correct <- as.factor(testData$Phenotype)



# ROC --------------------
#ROC curves of these?
#all separate
one <- roc(response = predicted_ip10$correct, predictor = predicted_ip10$DB, plot=TRUE, legacy.axes=TRUE, col="#084D49FF", lwd=4 , print.auc=TRUE, print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
one_smooth <- smooth(one, method = "density")
plot(one_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman")
plot(one,  add = TRUE,print.auc = TRUE, print.auc.x=45,print.auc.y=25, 
     legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "IP-10", line = 2.5, adj = 0, family = "Times New Roman")

two <- roc(response =predicted_trail$correct, predictor = predicted_trail$DB, 
           plot=TRUE, legacy.axes=TRUE, col="#276C69FF", lwd=4,print.auc=TRUE,
           print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
two_smooth <- smooth(two, method = "density")
plot(two_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman")
plot(two, add = TRUE,print.auc = TRUE, 
     print.auc.x=45,print.auc.y=25, legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "TRAIL", line = 2.5, adj = 0, family = "Times New Roman")

three <- roc(response =predicted_ifng$correct, predictor = predicted_ifng$DB, 
             plot=TRUE, legacy.axes=TRUE, col="#189BA0FF", lwd=4 ,print.auc=TRUE, 
             print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
three_smooth <- smooth(three, method = "density")
plot(three_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman")
plot(three, add = TRUE, lwd = 5, print.auc = TRUE, print.auc.x=45,print.auc.y=25, legacy.axes = TRUE, family = "Times New Roman")
title(main = "IFNG", line = 2.5, adj = 0, family = "Times New Roman")

four <- roc(response =predicted_il4$correct, predictor = predicted_il4$DB, 
            plot=TRUE, legacy.axes=TRUE, col="#73C1C4FF", lwd=4 , print.auc=TRUE, 
            print.auc.x=.25, print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
four_smooth <- smooth(four, method = "density")
plot(four_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman")
plot(four, add = TRUE, print.auc = TRUE, print.auc.x=45, print.auc.y=25, 
     legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "IL-4", line = 2.5, adj = 0, family = "Times New Roman")

five <- roc(response =predicted_pct$correct, predictor = predicted_pct$DB, 
            plot=TRUE, legacy.axes=TRUE, col="#BF8699FF", lwd=4 ,print.auc=TRUE, 
            print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars",percent = TRUE)
five_smooth <- smooth(five, method = "density")
plot(five_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman")
plot(five, add = TRUE, print.auc = TRUE, print.auc.x=45, print.auc.y=25, 
     legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "Procalcitonin", line = 2.5, adj = 0, family = "Times New Roman")

six <- roc(response =predicted_lcn2$correct, predictor = predicted_lcn2$DB, 
           plot=TRUE, legacy.axes=TRUE, col="#A64264FF", lwd=4 , print.auc=TRUE, 
           print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
six_smooth <- smooth(six, method = "density")
plot(six_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman")
plot(six, add = TRUE, print.auc = TRUE, print.auc.x=45, print.auc.y=25, 
     legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "LCN2", line = 2.5, adj = 0, family = "Times New Roman")

seven <- roc(response =predicted_il6$correct, predictor = predicted_il6$DB, 
             plot=TRUE, legacy.axes=TRUE, col="#830042FF", lwd=4, print.auc=TRUE, 
             print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
seven_smooth <- smooth(seven, method = "density")
plot(seven_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman")
plot(seven, add = TRUE, print.auc = TRUE, print.auc.x=45, print.auc.y=25, 
     legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "IL-6", line = 2.5, adj = 0, family = "Times New Roman")



# performing FS-PLS -----------------------------------------------------
#running FS-PLS external to this script
#results = TRAIL, LCN2, IL-6

# models for FS-PLS -----------------------------------------------------
set.seed(50)
glm_trail_lcn2 <- train(Phenotype ~ TRAIL + LCN2, data=trainData, method='glm', metric = "ROC", trControl = cv10)
predicted_trail_lcn2 <- predict(glm_trail_lcn2, newdata = testData, type = "prob") # get probabilities of classifications
predicted_trail_lcn2$predicted <- as.factor(ifelse(predicted_trail_lcn2$DB >= 0.5, "DB", "DV"))
predicted_trail_lcn2$correct <- as.factor(testData$Phenotype)
confusionMatrix(predicted_trail_lcn2$predicted, predicted_trail_lcn2$correct)

set.seed(50)
glm_fspls <- train(Phenotype ~ TRAIL + LCN2 + IL_6, data=trainData, method='glm', metric = "ROC", trControl = cv10) #choosing our resampling method
predicted_fspls <- predict(glm_fspls, newdata = testData, type = "prob") # get probabilities of classifications
predicted_fspls$predicted <- as.factor(ifelse(predicted_fspls$DB >= 0.5, "DB", "DV"))
predicted_fspls$correct <- as.factor(testData$Phenotype)
confusionMatrix(predicted_fspls$predicted, predicted_fspls$correct)


# ROC --------------------------------------------------------------------------------------------------------------------------------------------
#ROC curves
roc(response = predicted_trail$correct, predictor = predicted_trail$DB, 
    plot=TRUE, legacy.axes=TRUE, col="#276C69FF", lwd=4 , print.auc=TRUE, 
    print.auc.x=50,print.auc.y=45, ci=TRUE, ci.type="bars", percent = TRUE, family = "Times New Roman")
plot.roc(x = predicted_trail_lcn2$correct, predictor = predicted_trail_lcn2$DB, 
         col="#e99787", lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=50,
         print.auc.y=40, ci=TRUE, percent = TRUE, family = "Times New Roman")
plot.roc(x = predicted_fspls$correct, predictor = predicted_fspls$DB, 
         col="#F34A4A", lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=50, 
         print.auc.y=35, ci=TRUE, percent = TRUE, family = "Times New Roman")

legend("bottomright", legend=c("TRAIL","TRAIL + LCN2", "TRAIL + LCN2 + IL-6"), col=c("#276C69FF", "#e99787","#F34A4A"), lwd=4, text.font = "serif")
title(main = "FS-PLS Signature Performance in Testing Set", adj = 0, line = 2.5)


#roc.test
roc.one <- roc(response = predicted_trail$correct, predictor = predicted_trail$DB, 
               plot=TRUE, legacy.axes=TRUE, col="#FDBE83", lwd=4 , print.auc=TRUE, 
               print.auc.x=50,print.auc.y=45, ci=TRUE, ci.type="bars", percent = TRUE)
roc.two <- roc(response = predicted_trail_lcn2$correct, predictor = predicted_trail_lcn2$DB, 
               col="#C8A3B5", lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=50,
               print.auc.y=40, ci=TRUE, percent = TRUE)
roc.three <- roc(response = predicted_fspls$correct, predictor = predicted_fspls$DB, 
                 col="#2F4E68", lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=50, 
                 print.auc.y=35, ci=TRUE, percent = TRUE)

roc.test(roc.one, roc.two, method = "delong", paired = TRUE)
roc.test(roc.two, roc.three, method = "delong", paired = TRUE)
roc.test(roc.one, roc.three, method = "delong", paired = TRUE)

## put together a master list of model prediciotn for roc plot
predict_master_2 <- data.frame(correct = testData$Phenotype, TRAIL = predicted_trail$DB, TRAIL_LCN2 = predicted_trail_lcn2$DB, FSPLS = predicted_fspls$DB)
roc.list_2 <- roc(correct ~ TRAIL + TRAIL_LCN2 + FSPLS, data = predict_master_2, 
                  ci = TRUE, legacy.axes = TRUE)
ci.list_2 <- lapply(roc.list_2, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list_2 <- lapply(ci.list_2, function(ciobj) 
  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3]))

#plot roc --------
#defining characteristics of the plot first 
cols_2 <- c("#276C69FF", "#e99787","#F34A4A") #define colours
models_2 <- c("TRAIL", "TRAIL + LCN2", "TRAIL + LCN2 + IL-6") #text for the models

#making labels
cis_2 <- NULL
aucs_2 <- NULL
for (i in 1:3) { 
  aucs_2[i] <- paste(as.character(round(roc.list_2[[i]]$ci, 3))[2])
  cis_2[i] <- paste0("(", paste(as.character(round(roc.list_2[[i]]$ci, 3))[c(1,3)], collapse = " - "), ")")
}
auc_ci_2 <- data.frame(models = models_2, auc = aucs_2, ci = cis_2, label = paste(models_2, ";  ", "AUC: ", aucs_2," ", cis_2, sep = "")) #making a df with labels together in it

#generate rocs
p_2 <- ggroc(roc.list_2, size = 1.5) + theme_minimal(base_size = 15) + 
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=1, color = "grey") + 
  coord_equal() + 
  scale_colour_manual(labels = auc_ci_2$label,values = cols_2) +
  labs(y= "Sensitivity", x= "Specificity", color = " ") + 
  theme(legend.position = c(0.8, 0.3), text = element_text(family = "Times New Roman"))


#add CI
for(i in 1:3) {
  p_2 <- p_2 + geom_ribbon(
    data = dat.ci.list_2[[i]],
    aes(x = x, ymin = lower, ymax = upper),
    fill = cols_2[i],
    alpha = 0.15,
    inherit.aes = F) + theme(text = element_text(family = "Times New Roman"))
} 

p_2



# Using a subset of CRP<60 DBs and testing performance of FS-PLS signature ----------------------------------------------------------------------
data_withcrp <- na.omit(data[crp_subset, ])
crp_train <- predict(preproc, data_withcrp)
crp_test <- predict(preproc, data_withcrp)


# using FS-PLS -----------------------------------------------------

set.seed(50) 
predicted_trail_crp <- predict(trail, newdata = crp_test, type = "prob") # get probabilities of classifications
predicted_trail_crp$predicted <- as.factor(ifelse(predicted_trail_crp$DB >= 0.5, "DB", "DV"))
predicted_trail_crp$correct <- as.factor(crp_test$Phenotype)

set.seed(50)
predicted_trail_lcn2_crp <- predict(glm_trail_lcn2, newdata = crp_test, type = "prob") # get probabilities of classifications
predicted_trail_lcn2_crp$predicted <- as.factor(ifelse(predicted_trail_lcn2_crp$DB >= 0.5, "DB", "DV"))
predicted_trail_lcn2_crp$correct <- as.factor(crp_test$Phenotype)

set.seed(50)
predicted_fspls_crp <- predict(glm_fspls, newdata = crp_test, type = "prob") # get probabilities of classifications
predicted_fspls_crp$predicted <- as.factor(ifelse(predicted_fspls_crp$DB >= 0.5, "DB", "DV"))
predicted_fspls_crp$correct <- as.factor(crp_test$Phenotype)
confusionMatrix(predicted_fspls_crp$predicted, predicted_fspls_crp$correct)

# ROC --------------------------------------------------------------------------------------------------------------------------------------------
#ROC curves of these?
roc(response =predicted_trail_crp$correct, predictor = predicted_trail_crp$DB, plot=TRUE, legacy.axes=TRUE, col="#276C69FF", lwd=4 , print.auc= TRUE, print.auc.x=50,print.auc.y=45, ci=TRUE, ci.type="bars", percent = TRUE)
plot.roc(x = predicted_trail_lcn2_crp$correct, predictor = predicted_trail_lcn2_crp$DB, col="#e99787", lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=50, print.auc.y=40, ci=TRUE, percent = TRUE)
plot.roc(x = predicted_fspls_crp$correct, predictor = predicted_fspls_crp$DB, col="#F34A4A", lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=50, print.auc.y=35, ci=TRUE, percent = TRUE)

legend("bottomright", legend=c("TRAIL","TRAIL + LCN2",  "TRAIL + LCN2 + IL-6"), col=c("#276C69FF", "#e99787","#F34A4A"), lwd=4)
title("FS-PLS Signature Performance in Patients with CRP < 60 mg/L", line = 2.5, adj = 0)

#plot roc --------
## put together a master list of model prediciotn for roc plot
predict_master_3 <- data.frame(correct = crp_test$Phenotype, TRAIL = predicted_trail_crp$DB, TRAIL_LCN2 = predicted_trail_lcn2_crp$DB, FSPLS = predicted_fspls_crp$DB)
roc.list_3 <- roc(correct ~ TRAIL + TRAIL_LCN2 + FSPLS, data = predict_master_3, ci = TRUE)
ci.list_3 <- lapply(roc.list_3, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list_3 <- lapply(ci.list_3, function(ciobj) 
  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3]))
#defining characteristics of the plot first 
cols_3 <- c("#276C69FF", "#e99787","#F34A4A") #define colours
models_3 <- c("TRAIL", "TRAIL + LCN2", "TRAIL + LCN2 + IL-6") #text for the models

#making labels
cis_3 <- NULL
aucs_3 <- NULL
for (i in 1:3) { 
  aucs_3[i] <- paste(as.character(round(roc.list_3[[i]]$ci, 3))[2])
  cis_3[i] <- paste0("(", paste(as.character(round(roc.list_3[[i]]$ci, 3))[c(1,3)], collapse = " - "), ")")
}
auc_ci_3 <- data.frame(models = models_3, auc = aucs_3, ci = cis_3, label = paste(models_3, ";  ", "AUC: ", aucs_3," ", cis_3, sep = "")) #making a df with labels together in it

#generate rocs
p_3 <- ggroc(roc.list_3, size = 1.5) + theme_minimal(base_size = 15) + 
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=1, color = "grey") + 
  coord_equal() + 
  scale_colour_manual(labels = auc_ci_3$label,values = cols_3) +
  labs(y= "Sensitivity", x= "Specificity", colour = " ", title = "ROC for performance of FS-PLS protein signature in CRP < 60 subgroup") + theme(legend.position = "none") 


#add CI
for(i in 1:3) {
  p_3 <- p_3 + geom_ribbon(
    data = dat.ci.list_3[[i]],
    aes(x = x, ymin = lower, ymax = upper),
    fill = cols_3[i],
    alpha = 0.15,
    inherit.aes = F) 
} 

p_3


# Testing the MeMed BV(R) signature in CRP <60 group -------------------------------------
# MeMed signature = CRP + TRAIL + IP-10
memed_train <- predict(preproc, data_withcrp)
memed_test <- predict(preproc, data_withcrp)

set.seed(50)
glm_memed <- train(Phenotype ~ TRAIL + IP_10 + CRP, data=memed_train, method='glm', metric = "ROC", trControl = cv10) #choosing our resampling method
predicted_memed <- predict(glm_memed, newdata = memed_test, type = "prob") # get probabilities of classifications
predicted_memed$predicted <- as.factor(ifelse(predicted_memed$DB >= 0.5, "DB", "DV"))
predicted_memed$correct <- as.factor(memed_test$Phenotype)
confusionMatrix(predicted_memed$predicted, predicted_memed$correct)

#roc
roc(response =predicted_memed$correct, predictor = predicted_memed$DB, plot=TRUE, legacy.axes=TRUE, col="#276C69FF", lwd=4 , print.auc=TRUE, print.auc.x=50,print.auc.y=40, ci=TRUE, ci.type="bars", percent = TRUE)
plot.roc(x = predicted_fspls_crp$correct, predictor = predicted_fspls_crp$DB, col="#e99787", lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=50, print.auc.y=35, ci=TRUE, percent = TRUE)
title("Performance of FS-PLS and MeMed Three-protein Signatures in patients with CRP < 60 mg/L", adj = 0, line = 2.5)
legend("bottomright", legend=c("MeMed Signature", "FS-PLS Signature"), col=c( "#276C69FF","#e99787"), lwd=4)

#plot roc --------
## put together a master list of model prediciotn for roc plot
predict_master_4 <- data.frame(correct = crp_test$Phenotype, MeMed = predicted_memed$DB, FSPLS = predicted_fspls_crp$DB)
roc.list_4 <- roc(correct ~ MeMed + FSPLS, data = predict_master_4, ci = TRUE)
ci.list_4 <- lapply(roc.list_4, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list_4 <- lapply(ci.list_4, function(ciobj) 
  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3]))


#defining characteristics of the plot first 
cols_4 <- c("#276C69FF", "#e99787") #define colours
models_4 <- c("MeMed Signature", "FS-PLS Signature") #text for the models

#making labels
cis_4 <- NULL
aucs_4 <- NULL
for (i in 1:2) { 
  aucs_4[i] <- paste(as.character(round(roc.list_4[[i]]$ci, 3))[2])
  cis_4[i] <- paste0("(", paste(as.character(round(roc.list_4[[i]]$ci, 3))[c(1,3)], collapse = " - "), ")")
}
auc_ci_4 <- data.frame(models = models_4, auc = aucs_4, ci = cis_4, label = paste(models_4, ";  ", "AUC: ", aucs_4," ", cis_4, sep = "")) #making a df with labels together in it

#generate rocs
p_4 <- ggroc(roc.list_4, size = 1.5) + theme_minimal(base_size = 15) + 
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=1, color = "grey") + 
  coord_equal() + 
  scale_colour_manual(labels = auc_ci_4$label,values = cols_4) +
  labs(y= "Sensitivity", x= "Specificity", color = " ") + theme(legend.position = "none")


#add CI
for(i in 1:2) {
  p_4 <- p_4 + geom_ribbon(
    data = dat.ci.list_4[[i]],
    aes(x = x, ymin = lower, ymax = upper),
    fill = cols_4[i],
    alpha = 0.15,
    inherit.aes = F) 
} 

p_4



#test significane between fs-pls signature and memed using roc.test
one <- roc(response = predicted_memed$correct, predictor = predicted_memed$DB, plot=TRUE, legacy.axes=TRUE, col="#276C69FF", lwd=4 , print.auc=TRUE, print.auc.x=50,print.auc.y=40, ci=TRUE, ci.type="bars", percent = TRUE)
two <- roc(response = predicted_fspls_crp$correct, predictor = predicted_fspls_crp$DB,plot=TRUE, legacy.axes=TRUE, col="#e99787", lwd=4,  print.auc=TRUE, print.auc.x=50,print.auc.y=40, ci=TRUE, ci.type="bars", percent = TRUE)

#paired because they are the same patients in both tests. 
roc.test(one, two, paired = TRUE, method = "delong")



       
       

