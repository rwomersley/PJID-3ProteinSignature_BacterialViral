# final script - proteomics
library(readxl)
library(dplyr)
library(conflicted)
conflict_prefer("filter", "dplyr")
library(caret)
library(pROC)
library(FactoMineR)
library(factoextra)


rawdata <- as.data.frame(read_xlsx("Proteomics_with_Meta_final_June-22.xlsx"))

row.names(rawdata) <- rawdata$PatientID
rawdata$Phenotype <- as.factor(rawdata$Phenotype_PERFORM)
levels(rawdata$Phenotype) <- c("DB", "DV")
crp <- (rawdata$CRP_database <= 60)
colnames(rawdata)[1:7] <- c("IL_6", "IP_10", "TRAIL", "Procalcitonin", "LCN2", "IFNG", "IL_4")
data <- as.data.frame(rawdata[,c(1:7,31)])


final_data <- data[,c(2,3,6,7,4,5,1,8)]

## modelling ---------------------------------------------------------
#dividing into training and testing set
set.seed(50)
part <- createDataPartition(data$Phenotype, p = 0.7, list = FALSE)
trainData <- final_data[part,]
testData <- final_data[-part,]


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


# are variables correlated?-----------------------------------------------------------------------------------------------------
#function for testing significance for correlation
library(psych)
corr <- corr.test(trainData[,1:7], adjust = "bonferroni") #bonferroni correction

colnames(corr$r) <- c('IP-10', "TRAIL","IFNG",'IL-4', "Procalcitonin", "LCN2",'IL-6')
colnames(corr$p) <- c('IP-10', "TRAIL","IFNG",'IL-4', "Procalcitonin", "LCN2",'IL-6')

rownames(corr$r) <- c('IP-10', "TRAIL","IFNG",'IL-4', "Procalcitonin", "LCN2",'IL-6')
rownames(corr$p) <- c('IP-10', "TRAIL","IFNG",'IL-4', "Procalcitonin", "LCN2",'IL-6')

corr
print(corr, short = FALSE)

library("corrplot")
corrplot(corr$r, method = "square", type = "upper" )
col <- colorRampPalette(c("#4477AA", "#77AADD","#FFFFFF","#EE9988","#BB4444"))
corrplot(corr$r, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = corr$p, sig.level = 1, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE, 
)


# and a quick explore of what we've got -----------------------------------------------------------
#boxplot
#training
featurePlot(x = trainData[, -8], 
            y = trainData$Phenotype, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))
#testing set
featurePlot(x = testData[, -8], 
            y = testData$Phenotype, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))
#density plot
#training set
featurePlot(x = (trainData[, -8]), 
            y = trainData$Phenotype, 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))
#testing set
featurePlot(x = testData[, -8], 
            y = testData$Phenotype, 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

# 6) FITTING MODELS -----------------------------------------------------------------------------------------------------
# hypothesis ---- #yes - sequential build. #no - direct build. #not sure? - stepwise regression (forwards selection or barckwards elimination) or best subsets regression. 
#not sure - lets do a stepwise regression. 

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
plot(one_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman", cex.lab = 1.5, cex.axis = 1.5)
plot(one,  add = TRUE,print.auc = TRUE, print.auc.x=60,print.auc.y=25, print.auc.cex = 1.3, 
     legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "IP-10", line = 2.5, adj = 0, family = "Times New Roman", cex.main = 2)

two <- roc(response =predicted_trail$correct, predictor = predicted_trail$DB, 
           plot=TRUE, legacy.axes=TRUE, col="#276C69FF", lwd=4,print.auc=TRUE,
           print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
two_smooth <- smooth(two, method = "density")
plot(two_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman", cex.lab = 1.5, cex.axis = 1.5)
plot(two, add = TRUE,print.auc = TRUE, print.auc.cex = 1.5, 
     print.auc.x=70,print.auc.y=25, legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "TRAIL", line = 2.5, adj = 0, family = "Times New Roman", cex.main = 2)

three <- roc(response =predicted_ifng$correct, predictor = predicted_ifng$DB, 
             plot=TRUE, legacy.axes=TRUE, col="#189BA0FF", lwd=4 ,print.auc=TRUE, 
             print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
three_smooth <- smooth(three, method = "density")
plot(three_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman", cex.lab = 1.5, cex.axis = 1.5)
plot(three, add = TRUE, lwd = 5, print.auc = TRUE, print.auc.x=70,print.auc.y=25, legacy.axes = TRUE, print.auc.cex = 1.5, family = "Times New Roman")
title(main = "IFNG", line = 2.5, adj = 0, family = "Times New Roman", cex.main = 2)

four <- roc(response =predicted_il4$correct, predictor = predicted_il4$DB, 
            plot=TRUE, legacy.axes=TRUE, col="#73C1C4FF", lwd=4 , print.auc=TRUE, 
            print.auc.x=.25, print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
four_smooth <- smooth(four, method = "density")
plot(four_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman", cex.lab = 1.5, cex.axis = 1.5)
plot(four, add = TRUE, print.auc = TRUE,  print.auc.x=70,print.auc.y=25, print.auc.cex = 1.5, 
     legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "IL-4", line = 2.5, adj = 0, family = "Times New Roman", cex.main = 2)

five <- roc(response =predicted_pct$correct, predictor = predicted_pct$DB, 
            plot=TRUE, legacy.axes=TRUE, col="#BF8699FF", lwd=4 ,print.auc=TRUE, 
            print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars",percent = TRUE)
five_smooth <- smooth(five, method = "density")
plot(five_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman", cex.lab = 1.5, cex.axis = 1.5)
plot(five, add = TRUE, print.auc = TRUE,  print.auc.x=70,print.auc.y=25, print.auc.cex = 1.5, 
     legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "Procalcitonin", line = 2.5, adj = 0, family = "Times New Roman", cex.main = 2)

six <- roc(response =predicted_lcn2$correct, predictor = predicted_lcn2$DB, 
           plot=TRUE, legacy.axes=TRUE, col="#A64264FF", lwd=4 , print.auc=TRUE, 
           print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
six_smooth <- smooth(six, method = "density")
plot(six_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman", cex.lab = 1.5, cex.axis = 1.5)
plot(six, add = TRUE, print.auc = TRUE, print.auc.x=70,print.auc.y=25, print.auc.cex = 1.5, 
     legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "LCN2", line = 2.5, adj = 0, family = "Times New Roman", cex.main = 2)

seven <- roc(response =predicted_il6$correct, predictor = predicted_il6$DB, 
             plot=TRUE, legacy.axes=TRUE, col="#830042FF", lwd=4, print.auc=TRUE, 
             print.auc.x=.25,print.auc.y=.4, ci=TRUE, ci.type="bars", percent = TRUE)
seven_smooth <- smooth(seven, method = "density")
plot(seven_smooth, col = "grey",  legacy.axes = TRUE, percent = TRUE, lwd = 5, family = "Times New Roman", cex.lab = 1.5, cex.axis = 1.5)
plot(seven, add = TRUE, print.auc = TRUE, print.auc.x=70,print.auc.y=25, print.auc.cex = 1.5, 
     legacy.axes = TRUE, lwd = 5, family = "Times New Roman")
title(main = "IL-6", line = 2.5, adj = 0, family = "Times New Roman", cex.main = 2)




# building models--------------------------------------------------------------------------
# performing FS-PLS -----------------------------------------------------
#fspls_example (1).R 

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
    print.auc.cex = 1.5, cex.lab = 1.5, cex.axis = 1.5,
    print.auc.x=35, print.auc.y=30, ci=TRUE, ci.type="bars", percent = TRUE, family = "Times New Roman")
plot.roc(x = predicted_trail_lcn2$correct, predictor = predicted_trail_lcn2$DB, 
         col="#e99787", lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=35, print.auc.cex = 1.5,
         print.auc.y=25, ci=TRUE, percent = TRUE, family = "Times New Roman")
plot.roc(x = predicted_fspls$correct, predictor = predicted_fspls$DB, 
         col="#F34A4A", lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=35, print.auc.cex = 1.5,
         print.auc.y=20, ci=TRUE, percent = TRUE, family = "Times New Roman")

legend("bottomright", legend=c("TRAIL","TRAIL + LCN2", "TRAIL + LCN2 + IL-6"), col=c("#276C69FF", "#e99787","#F34A4A"), lwd=4, cex = 1.5, text(vfont = "serif")) #font = "serif")
title(main = "FS-PLS Signature Performance in Testing Set", adj = 0, line = 2.5, cex.main = 1.5)


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

## put together a master list of model prediction for roc plot
predict_master_2 <- data.frame(correct = testData$Phenotype, TRAIL = predicted_trail$DB, TRAIL_LCN2 = predicted_trail_lcn2$DB, FSPLS = predicted_fspls$DB)
roc.list_2 <- roc(correct ~ TRAIL + TRAIL_LCN2 + FSPLS, data = predict_master_2, percent = TRUE,
                  ci = TRUE, legacy.axes = TRUE)
ci.list_2 <- lapply(roc.list_2, ci.se, specificities = seq(0, 1, l = 25))
for (i in 1:3) {
  rownames(ci.list_2[[i]]) <- sub("%", "", rownames(ci.list_2[[i]]))
  rownames(ci.list_2[[i]]) <- as.numeric(rownames(ci.list_2[[i]]))*100
}

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
  scale_colour_manual(labels = auc_ci_2$label, values = cols_2) +
  labs(y= "Sensitivity", x= "Specificity", color = " ") + 
  theme(legend.position = c(0.7, 0.2), text = element_text(family = "Times New Roman", size = 20))


#add CI
for(i in 1:3) {
  p_2 <- p_2 + geom_ribbon(
    data = dat.ci.list_2[[i]],
    aes(x = x, ymin = lower, ymax = upper),
    fill = cols_2[i],
    alpha = 0.15, 
    inherit.aes = F) + theme(text = element_text(family = "Times New Roman"), )
} 

p_2



# Using a subset of CRP<60 DBs and testing performance of FS-PLS signature ----------------------------------------------------------------------
data_withcrp <- na.omit(rawdata[crp, c(1:7,10, 31)])
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
roc(response =predicted_trail_crp$correct, predictor = predicted_trail_crp$DB, plot=TRUE, legacy.axes=TRUE, col="#276C69FF", lwd=4 , 
    print.auc= TRUE, print.auc.x=35,print.auc.y=30, ci=TRUE, ci.type="bars", percent = TRUE,
    print.auc.cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, family = "Times New Roman")
plot.roc(x = predicted_trail_lcn2_crp$correct, predictor = predicted_trail_lcn2_crp$DB, col="#e99787", family = "Times New Roman",
         print.auc.cex = 1.5,lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=35, print.auc.y=25, ci=TRUE, percent = TRUE)
plot.roc(x = predicted_fspls_crp$correct, predictor = predicted_fspls_crp$DB, col="#F34A4A", family = "Times New Roman", 
         print.auc.cex = 1.5, lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=35, print.auc.y=20, ci=TRUE, percent = TRUE)

legend("bottomright", legend=c("TRAIL","TRAIL + LCN2",  "TRAIL + LCN2 + IL-6"), col=c("#276C69FF", "#e99787","#F34A4A"), lwd=4, cex = 1.5)# family = windowsFont("Times New Roman"))
title("FS-PLS Signature Performance in Patients with CRP < 60 mg/L", line = 2.5, adj = 0, cex.main = 1.5)

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


# Testing the MeMed signature in CRP <60 group -------------------------------------
# MeMed signature = CRP + TRAIL + IP-10
memed_train <- predict(preproc, data_withcrp)
memed_test <- predict(preproc, data_withcrp)

set.seed(50)
glm_memed <- train(Phenotype ~ TRAIL + IP_10 + CRP_database, data=memed_train, method='glm', metric = "ROC", trControl = cv10) #choosing our resampling method
predicted_memed <- predict(glm_memed, newdata = memed_test, type = "prob") # get probabilities of classifications
predicted_memed$predicted <- as.factor(ifelse(predicted_memed$DB >= 0.5, "DB", "DV"))
predicted_memed$correct <- as.factor(memed_test$Phenotype)
confusionMatrix(predicted_memed$predicted, predicted_memed$correct)

#roc
roc(response =predicted_memed$correct, predictor = predicted_memed$DB, plot=TRUE, legacy.axes=TRUE, col="#276C69FF", lwd=4 , print.auc=TRUE, print.auc.x=35,print.auc.y=10, ci=TRUE, ci.type="bars", percent = TRUE, print.auc.cex = 1.5, cex.axis = 1.5, cex.lab = 1.5, family = "Times New Roman")
plot.roc(x = predicted_fspls_crp$correct, predictor = predicted_fspls_crp$DB, col="#e99787", lwd=4, add=TRUE, print.auc=TRUE, print.auc.x=35, print.auc.y=5, ci=TRUE, percent = TRUE, print.auc.cex = 1.5, family = "Times New Roman")
title("Performance of FS-PLS and MeMed Three-protein Signatures in patients with CRP < 60 mg/L", adj = 0, line = 2.5, cex.main = 1.2, family = "Times New Roman")
legend(x = 27, y = 10, legend=c("MeMed Signature", "FS-PLS Signature"), col=c( "#276C69FF","#e99787"), lwd=4, cex = 1.2)

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

# supplementary figures -----------------------------------------
# PCA ------------------------------------------------
# -  - - - - 
#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
# -  - - - - 

normalised <- data
normalised[,1:7] <- lapply(normalised[,1:7], min_max_norm)
colnames(normalised)[8] <- "Phenotype"
levels(normalised$Phenotype) <- c("Definite Bacterial", "Definite Viral")

p <- p <- PCA(normalised, quali.sup = 8)

fviz_pca_ind(p, geom.ind = "point", # show points only (nbut not "text")
             palette = c("#276C69FF",  "#BF8699FF"),
             legend.title = "Phenotype", habillage = 8, mean.point = FALSE, addEllipses = TRUE) + 
  theme(text = element_text(family = "Times New Roman", size = 20)) + labs(x = paste("PC1 (", floor(p$eig[1,2]), "%)", sep = ""), y = paste("PC2 (", floor(p$eig[2,2]), "%)", sep = ""))



#pairs plot
library(GGally)
ggpairs(normalised[,1:8], mapping = ggplot2::aes(color = Phenotype)) + 
  scale_colour_manual(values = c("#ad2d8d", "#58a16f")) + 
  scale_fill_manual(values = c("#ad2d8d", "#58a16f"))


# boxplots ------------------------------------------------------
boxplot <- final_data
boxplot[,1:7] <- lapply(boxplot[,1:7], min_max_norm)
boxplot$fspls <- boxplot$LCN2 + boxplot$IL_6 - boxplot$TRAIL

ggplot(boxplot) +
  aes(x = Phenotype, y = fspls, fill = Phenotype) +
  geom_boxplot() + geom_jitter(width = 0.2) +
  scale_fill_manual(
    values = c(DB = "#AD2D8D",
               DV = "#58A16F")
    ) +
  theme_minimal() + labs(title = "DRS using FS-PLS selected proteins", y = "DRS")


ggplot(boxplot) +
  aes(x = Phenotype, y = fspls, fill = Phenotype) +
  geom_boxplot() + 
  scale_fill_manual(
    values = c(DB = "#AD2D8D",
               DV = "#58A16F")
  ) +
  theme_minimal() + labs(title = "DRS using FS-PLS selected proteins", y = "DRS")



#boxplot for supplementary 
library(tidyr)
library(ggpubr)
gathered <- gather(data = normalised, 'IL_6', 'IP_10', 'TRAIL', 'Procalcitonin', 'LCN2', 'IFNG', 'IL_4', key = "Protein", value = "Normalised Expression")
gathered$Protein <- factor(gathered$Protein, levels = c("IP_10", "TRAIL", "IFNG", "IL_4", "Procalcitonin", "LCN2", "IL_6"))
levels(gathered$Protein) <- c("IP-10", "TRAIL", "IFNG", "IL-4", "Procalcitonin", "LCN2", "IL-6")

#with jitter?
plot <- ggplot(gathered, aes(Phenotype, `Normalised Expression`, fill = Phenotype)) + 
  geom_boxplot(outlier.size = 2, outlier.shape = 24) + 
  geom_jitter(width = 0.1, size = 0.4, alpha = 0.5) + 
  facet_wrap(vars(Protein)) + theme_classic() + 
  scale_fill_manual(values = c("#276C69FF",  "#BF8699FF")) +
  theme(text = element_text(size = 15, family = "Times New Roman"))

#add stats to compare between groups. 
plot + stat_compare_means(method = "wilcox.test", 
                          aes(group = Phenotype,  label = ..p.signif..), 
                          label.x = 1.5, label.y = 0.9, size = 7)

#individual
#with jitter?
plot <- gathered %>%
  filter(Protein %in% "IL-6") %>%
  ggplot( aes(Phenotype, `Normalised Expression`, fill = Phenotype)) + 
  geom_boxplot(outlier.size = 5, outlier.shape = 24) + 
  geom_jitter(width = 0.1, size = 3, alpha = 0.5) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + title(main = "IL-6") +
  scale_fill_manual(values = c("#276C69FF",  "#BF8699FF")) +
  theme(text = element_text(size = 25, family = "Times New Roman"), legend.position = "none") +
  labs(title = "IL-6", element_text(size = 7))

#add stats to compare between groups. 
plot + stat_compare_means(method = "wilcox.test", 
                          aes(group = Phenotype,  label = ..p.signif..), 
                          label.x = 1.5, label.y = 0.9, size = 15)

       
       


       
       

