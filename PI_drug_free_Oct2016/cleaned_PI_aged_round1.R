options(contrasts=c("contr.sum", "contr.poly"))
setwd("~/Documents/GitHub/PI_2016_2017/PI_drug_free_Oct2016")
library(foreign)
library(psych)
library(reshape2)
library(ggplot2)
library(car)
library(ez)
library(Rmisc)
library(gmodels)
library(wesanderson)

# for AB acquisition
sessList <- list(23,24,25)#list session nubmers here
ratList <- list(25,35,46,49,53,57)#list rat number here
x <- length(sessList)
y <- length(ratList)
rowsNeeded <- (x*y)

#to indicate lesion condition of animals
sapGroup = list(28,36,50,63,64)
conGroup = list(25,35,46,49,53,57)

ABsummary <- matrix(data = NA, nrow = rowsNeeded, ncol=3)
colnames(ABsummary)<- c('ratID', 'session', 'accuracy')

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
    ABsummary[i,0:3] <- c(rat, sessNum, ABacc)
    i=i+1
    
  }
}
ABsummary <- data.frame(ABsummary)
ABsummary$ratID <- as.factor(ABsummary$ratID)
ABsummary$session <- as.factor(ABsummary$session)
ABsummary$condition <- ifelse(is.element(ABsummary$ratID, sapGroup) ,'sap','control')

ABsummaryTable <- aggregate(accuracy~session+condition, data=ABsummary, FUN=(mean))
tempTable <- aggregate(accuracy~session, data=ABsummary, FUN=(sd))
ABsummaryTable<- as.data.frame(ABsummaryTable)
ABsummaryTable$sd <- tempTable$accuracy*100
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
ABsummaryTable$se <- tempTable$accuracy*100
ABsummaryTable$accuracy <- ABsummaryTable$accuracy * 100


#for full task
sessList <- list(27, 28, 29, 30)#list session nubmers here
ratList <- list(25,35,46,49,53,57)#list rat number here
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
FULLfinalSummaryW <- data.frame(FULLsummary)
FULLfinalSummaryW$condition <- ifelse(is.element(FULLfinalSummaryW$ratID, sapGroup) ,'sap','control')
z = rowsNeeded*3 
FULLfinalSummaryL <- reshape(FULLfinalSummaryW, varying = c('AB_acc', 'AC_acc', 'DE_acc'), v.names="accuracy",
                         timevar='trialType',times=c('AB', 'CA', 'DE'), new.row.names=1:z, direction = 'long')
FULLfinalSummaryL$ratID <-as.factor(FULLfinalSummaryL$ratID)
FULLfinalSummaryL$session <-as.factor(FULLfinalSummaryL$session)
FULLfinalSummaryL$trialType <- as.factor(FULLfinalSummaryL$trialType)
FULLfinalSummaryW$ratID <-as.factor(FULLfinalSummaryW$ratID)
FULLfinalSummaryW$session <-as.factor(FULLfinalSummaryW$session)

#making summary table for full task
FULLsummaryTable <- aggregate(accuracy~trialType+session+condition, data=FULLfinalSummaryL, FUN=(mean))
tempTable <- aggregate(accuracy~trialType+session, data=FULLfinalSummaryL, FUN=(sd))
FULLsummaryTable<- as.data.frame(FULLsummaryTable)
FULLsummaryTable$sd <- tempTable$accuracy*100
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
FULLsummaryTable$se <- tempTable$accuracy*100
FULLsummaryTable$accuracy <- FULLsummaryTable$accuracy * 100


#part 1 - AB acquisition
#check distribution
describe(ABsummary)
hist(ABsummary$accuracy)
#distribution is actually a bit strange, though all above chance - also pretty small number of rats
plot(accuracy~session, data=ABsummary)
mean(ABsummary$accuracy)
# no interesting results, though by the end of the first session the group average was almost 80% accuracy

#part 2 - full task
describe(FULLfinalSummaryW)
describe(FULLfinalSummaryL)
hist(FULLfinalSummaryL$accuracy)
#accuracy for the full task is normally distributed 
plot(accuracy~session, data=FULLfinalSummaryL)
#looks like session didn't ahve much effect
plot(accuracy~trialType, data=FULLfinalSummaryL)
#weak but present effect for trial type

#inferential tests of full task dataaov1 <- aov_ez("ratID","accuracy",FULLfinalSummaryL,within = c("trialType", "session"))
#not sure if above will let me check assumptions
testAOV <- ezANOVA(data=FULLfinalSummaryL, dv=accuracy, wid = ratID,within=c(trialType,session),detailed=TRUE,return_aov = TRUE)