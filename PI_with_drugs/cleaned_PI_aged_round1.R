# PI full task analysis

library(foreign)
library(psych)
library(reshape2)
library(ggplot2)
library(car)
library(Rmisc)
library(gmodels)
library(wesanderson)

# for AB acquisition
sessList <- list(23,24,25)#list session nubmers here
ratList <- list(25,28,35,36,46,49,50,53,57,63,64)#list rat number here
x <- length(sessList)
y <- length(ratList)
rowsNeeded <- (x*y)

#to indicate lesion condition of animals
sapGroup = list(28, 36,50,63,64)
conGroup = list(25,35,46,49,53,57)

ABsummary <- matrix(data = NA, nrow = rowsNeeded, ncol=3)
colnames(ABsummary)<- c('ratID', 'session', 'AB_acc')

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
    ABsummary[i,0:5] <- c(rat, sessNum, ABacc, ACacc, DEacc)
    i=i+1
    
  }
}

ABsummaryTable <- aggregate(accuracy~trialType+session+condition, data=ABfinalSummaryL, FUN=(mean))
tempTable <- aggregate(accuracy~trialType+session, data=ABfinalSummaryL, FUN=(sd))
ABsummaryTable<- as.data.frame(ABsummaryTable)
ABsummaryTable$sd <- tempTable$accuracy*100
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
ABsummaryTable$se <- tempTable$accuracy*100
ABsummaryTable$accuracy <- ABsummaryTable$accuracy * 100

#for full task
sessList <- list(27, 28, 29, 30)#list session nubmers here
ratList <- list(25,28,35,36,46,49,50,53,57,63,64)#list rat number here
x <- length(sessList)
y <- length(ratList)
rowsNeeded <- (x*y)

#to indicate lesion condition of animals
sapGroup = list(28, 36,50,63,64)
conGroup = list(25,35,46,49,53,57)

FULLsummary <- matrix(data = NA, nrow = rowsNeeded, ncol=5)
colnames(FULLsummary)<- c('ratID', 'session', 'AB_acc', 'AC_acc', 'DE_acc')

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
    FULLsummary[i,0:5] <- c(rat, sessNum, ABacc, ACacc, DEacc)
    i=i+1
    
  }
}

#prepare wide and long form data frames for full task
finalSummaryW <- data.frame(FULLsummary)
finalSummaryW$condition <- ifelse(is.element(finalSummaryW$ratID, sapGroup) ,'sap','control')
z = rowsNeeded*3 
FULLfinalSummaryL <- reshape(finalSummaryW, varying = c('AB_acc', 'AC_acc', 'DE_acc'), v.names="accuracy",
                         timevar='trialType',times=c('AB', 'CA', 'DE'), new.row.names=1:z, direction = 'long')
FULLfinalSummaryL$ratID <-as.factor(FULLfinalSummaryL$ratID)
FULLfinalSummaryL$session <-as.factor(FULLfinalSummaryL$session)
finalSummaryW$ratID <-as.factor(finalSummaryW$ratID)
finalSummaryW$session <-as.factor(finalSummaryW$session)
#making summary table for full task
FULLsummaryTable <- aggregate(accuracy~trialType+session+condition, data=FULLfinalSummaryL, FUN=(mean))
tempTable <- aggregate(accuracy~trialType+session, data=FULLfinalSummaryL, FUN=(sd))
FULLsummaryTable<- as.data.frame(FULLsummaryTable)
FULLsummaryTable$sd <- tempTable$accuracy*100
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
FULLsummaryTable$se <- tempTable$accuracy*100
FULLsummaryTable$accuracy <- FULLsummaryTable$accuracy * 100
