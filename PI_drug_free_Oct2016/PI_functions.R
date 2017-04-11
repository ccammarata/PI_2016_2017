# PI full task analysis

library(foreign)
library(psych)
library(reshape2)
library(ggplot2)
library(car)
library(Rmisc)
library(gmodels)
library(wesanderson)

# functions for full task
PI_fulltask_df <- function(sessList,ratList,sapGroup,conGroup){
  library(foreign)
  library(psych)
  library(reshape2)
  x <- length(sessList)
  y <- length(ratList)
  rowsNeeded <- (x*y)
  summary <- matrix(data = NA, nrow = rowsNeeded, ncol=5)
  colnames(summary)<- c('ratID', 'session', 'AB_acc', 'DE_acc','CA_acc')
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
      summary[i,0:5] <- c(rat, sessNum, ABacc, ACacc, DEacc)
      i=i+1
    }
  }
  finalSummaryW <- data.frame(summary)
  finalSummaryW$id <- finalSummaryW$rat
  finalSummaryW$condition <- ifelse(is.element(finalSummaryW$ratID, sapGroup) ,'sap','control')
  z = rowsNeeded*3 
  finalSummaryL <- reshape(finalSummaryW, varying = c('AB_acc', 'AC_acc', 'DE_acc'), v.names="accuracy",
                           timevar='trialType',times=c('AB', 'CA', 'DE'), new.row.names=1:z, direction = 'long')
  return(finalSummaryL)
}

PI_fulltask_df_summary <- function(finalsummaryL){
  summaryTable <- aggregate(accuracy~trialType+session+condition, data=finalSummaryL, FUN=(mean))
  tempTable <- aggregate(accuracy~trialType+session, data=finalSummaryL, FUN=(sd))
  summaryTable<- as.data.frame(summaryTable)
  summaryTable$sd <- tempTable$accuracy*100
  tempTable$accuracy <- tempTable$accuracy/sqrt(y)
  summaryTable$se <- tempTable$accuracy*100
  summaryTable$accuracy <- summaryTable$accuracy * 100
  return(summaryTable)
}

#### functions for AB only

PI_AB_df <- function(sessList,ratList,sapGroup,conGroup){
  library(foreign)
  library(psych)
  library(reshape2)
  x <- length(sessList)
  y <- length(ratList)
  rowsNeeded <- (x*y)
  summary <- matrix(data = NA, nrow = rowsNeeded, ncol=3)
  colnames(summary)<- c('ratID', 'session', 'AB_acc')
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
      summary[i,0:5] <- c(rat, sessNum, ABacc)
      i=i+1
    }
  }
  finalSummaryW <- data.frame(summary)
  finalSummaryW$id <- finalSummaryW$rat
  finalSummaryW$condition <- ifelse(is.element(finalSummaryW$ratID, sapGroup) ,'sap','control')
  return(finalSummaryW)
}

PI_AB_df_summary <- function(finalsummaryW){
  summaryTable <- aggregate(accuracy~trialType+session+condition, data=finalSummaryW, FUN=(mean))
  tempTable <- aggregate(accuracy~trialType+session, data=finalSummaryW, FUN=(sd))
  summaryTable<- as.data.frame(summaryTable)
  summaryTable$sd <- tempTable$accuracy*100
  tempTable$accuracy <- tempTable$accuracy/sqrt(y)
  summaryTable$se <- tempTable$accuracy*100
  summaryTable$accuracy <- summaryTable$accuracy * 100
  return(summaryTable)
}