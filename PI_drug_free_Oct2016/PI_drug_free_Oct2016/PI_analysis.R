# PI full task analysis

library(foreign)
library(psych)
library(reshape2)
library(ggplot2)
library(car)
library(Rmisc)
library(gmodels)
library("wesanderson", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")

sessList <- list(27, 28, 29, 30, 31)#list session nubmers here
ratList <- list(25,28,35,36,46,49,50,53,57,63,64)#list rat number here
x <- length(sessList)
y <- length(ratList)
rowsNeeded <- (x*y)

#to indicate lesion condition of animals
sapGroup = list(28, 36,50,63,64)
conGroup = list(25,35,46,49,53,57)

summary <- matrix(data = NA, nrow = rowsNeeded, ncol=5)
colnames(summary)<- c('ratID', 'session', 'AB_acc', 'AC_acc', 'DE_acc')
criteriaOverview <- matrix(data = NA, nrow = rowsNeeded, ncol=2)

i = 1
for (sessNum in sessList){
  for (rat in ratList) {
    filename <- paste("U3_",rat,"_",sessNum,".csv", sep="")
    criteriaFilename <- paste("U3_CriteriaMet",rat,"_",sessNum,".csv", sep="")
    print(filename)
    currentData <- read.csv(filename, header=TRUE, sep="\t")
    currentData$resp2 <- ifelse(currentData$Response == "correct", 1,0)
    currentAB <- subset(currentData, trialType == "AB")
    ABacc <- mean(currentAB$resp2)
    currentDE <- subset(currentData, trialType == "DE")
    DEacc <- mean(currentDE$resp2)
    currentAC <- subset(currentData, trialType == "CA")
    ACacc <- mean(currentAC$resp2)
    #sliding window criteria?
    summary[i,0:5] <- c(rat, sessNum, ABacc, ACacc, DEacc)
    #criteriaCurrentFile
    i=i+1

  }
}


#prepare wide and long form data frames
finalSummaryW <- data.frame(summary)
finalSummaryW$id <- finalSummaryW$rat
finalSummaryW$condition <- ifelse(is.element(finalSummaryW$ratID, sapGroup) ,'sap','control')


z = rowsNeeded*3 
finalSummaryL <- reshape(finalSummaryW, varying = c('AB_acc', 'AC_acc', 'DE_acc'), v.names="accuracy",
timevar='trialType',times=c('AB', 'CA', 'DE'), new.row.names=1:z, direction = 'long')

#data analysis
summaryTable <- aggregate(accuracy~trialType+session+condition, data=finalSummaryL, FUN=(mean))
tempTable <- aggregate(accuracy~trialType+session, data=finalSummaryL, FUN=(sd))
summaryTable<- as.data.frame(summaryTable)
summaryTable$sd <- tempTable$accuracy*100
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
summaryTable$se <- tempTable$accuracy*100
summaryTable$accuracy <- summaryTable$accuracy * 100

finalSummaryW$AB_CA_diff <- finalSummaryW$AB_acc - finalSummaryW$AC_acc
finalSummaryW$DE_CA_diff <- finalSummaryW$DE_acc - finalSummaryW$AC_acc
finalSummaryW$AB_DE_diff <- finalSummaryW$AB_acc - finalSummaryW$DE_acc


#plotting all together
finalSummaryL$session <- as.factor(finalSummaryL$session)
library(ggplot2)
library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
ggplot(summaryTable, aes(x=session, y=accuracy, fill=trialType)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  #geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se), width=.2,
                #position=position_dodge(.9)) +
  scale_fill_brewer(palette = "Dark2")

#plotting seperated by condition
library(ggplot2)
library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
ggplot(summaryTable, aes(x=session, y=accuracy, fill=trialType)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=wes_palette("Royal1"))+
  facet_grid(~condition)+
  ggtitle("PI accuracy with SEM error bars")

#inferential stats
#mixed ANOVA 3X5X2 trial type by accuracy by lesion type


#repeated measures ANOVA 3x5 for trial type X accuracy within each group
finalSummaryL$ratID <- as.factor(finalSummaryL$ratID)
finalSummaryL_AC_DE_only <- subset(finalSummaryL, trialType != "AB")
# 
# #mixed anova
# anovaBetween <- aov(accuracy~condition*trialType*session + Error(ratID/(trialType*session))+condition, data=finalSummaryL_AC_DE_only)
# summary(anovaBetween)
# 
# #within ANOVa for sap rats only
# sapL <- subset(finalSummaryL_AC_DE_only, condition="sap")
# anovaWithinSap <- aov(accuracy~trialType*session +Error(ratID/(trialType*session)), data=sapL)
# summary(anovaWithinSap)
# 
# #within ANOVA for con rats only
# conL <- subset(finalSummaryL_AC_DE_only, condition="con")
# anovaWithinCon <- aov(accuracy~trialType*session +Error(ratID/(trialType*session)), data=conL)
# summary(anovaWithinCon)
# 
# 
# #anova or trialType by condition
# model1 <- lm(cbind(AB_acc, AC_acc, DE_acc) ~ condition, data = finalSummaryW)
# TrialType <- factor(rep(1:3))
# idata <- data.frame(TrialType)
# test <- Anova(model1, idata=idata, idesign = ~TrialType, type=3)
# summary(test, multivariate = FALSE)
# finalSummaryW$session <- as.factor(finalSummaryW$session)
# 
# model2 <- lm(cbind(AB_acc, AC_acc, DE_acc) ~ session, data = finalSummaryW)
# TrialType <- factor(rep(1:3))
# idata <- data.frame(TrialType)
# test <- Anova(model2, idata=idata, idesign = ~TrialType, type=3)
# summary(test, multivariate = FALSE)
# 
# #to do two within subjects factors need to add columns to wide format so that there's trialtypeXsession
# #then need to include all those levels in the model and make a factor for each
# 
# model <- lm(cbind(recall_10, recall_11, recall_20, recall_21, recall_30, recall_31) ~
#               1, data = rcwide)
# # generate factor design matrix
# stimulus <- factor(rep(1:3, each = 2))
# drunk <- factor(rep(0:1, times = 3))
# idata <- data.frame(stimulus, drunk)
# # run ws ANOVA using the Anova command in package car
# test <- Anova(model, idata = idata, idesign = ~stimulus * drunk, type = 3)