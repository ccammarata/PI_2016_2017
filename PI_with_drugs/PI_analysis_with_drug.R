# PI full task analysis

library(foreign)
library(psych)
library(reshape2)
library(ggplot2)
#library(car) might not need this if I'm using afex
library(Rmisc)
library(gmodels)
library(wesanderson)
library(afex)
library(car)
library(phia)

#list rat number here
ratList <- list(28, 36,50,63,64,25,35,46,49,53,57)
#to indicate lesion condition of animals
sapGroup = list(28, 36,50,63,64)
conGroup = list(25,35,46,49,53,57)



#make a data frame for the AB acquisition days
sessList <- list(35,36,37,43,44,45)#list session nubmers here
x <- length(sessList)
y <- length(ratList)
rowsNeeded <- (x*y)

#making an empty matrix
AB_acquisition <- matrix(data = NA, nrow = rowsNeeded, ncol=4)
colnames(AB_acquisition)<- c('ratID', 'session', 'task',' AB_acc')

i = 1
for (sessNum in sessList){
  for (rat in ratList) {
    filename <- paste("U3_",rat,"_",sessNum,".csv", sep="")
    print(filename)
    currentData <- read.csv(filename, header=TRUE, sep="\t")
    currentData$resp2 <- ifelse(currentData$Response == "correct", 1,0)
    currentAB <- subset(currentData, trialType == "AB")
    ABacc <- mean(currentAB$resp2)
    #sliding window criteria?
    AB_acquisition[i,0:4] <- c(rat, sessNum, 'AB_aq', ABacc)
    i=i+1
  }
}
#making the matrix into a data frame and adding lesion condition
AB_acquisitionW <- data.frame(AB_acquisition) #making it into a data frame
AB_acquisitionW$ratID <- factor(AB_acquisitionW$ratID) #making rat number a factor
AB_acquisitionW$condition <- ifelse(is.element(AB_acquisitionW$ratID, sapGroup) ,'sap','control') #adding lesion condition

###!!!### make the accuracy values character and then numeric

#make a data frame for the full task days
sessList <- list(38,39,40,41,46,47,48,49)#list session nubmers here
##!!## note that I am ignoring the 5th day of round 1
x <- length(sessList)
y <- length(ratList)
rowsNeeded <- (x*y)

#to indicate drug condition by day
drug1First <- list(25,28,46,50,53,64)#here I'll enter the ID numbers of the animlas that had drug 1 first; all others had drug 2 first
firstWeek <- list(38,39,40,41)

#making an empty matrix
summary <- matrix(data = NA, nrow = rowsNeeded, ncol=6)
colnames(summary)<- c('ratID', 'session', 'task', 'AB_acc', 'AC_acc', 'DE_acc')
criteriaOverview <- matrix(data = NA, nrow = rowsNeeded, ncol=2)

i = 1
for (sessNum in sessList){
  for (rat in ratList) {
    filename <- paste("U3_",rat,"_",sessNum,".csv", sep="")
    print(filename)
    currentData <- read.csv(filename, header=TRUE, sep="\t")
    currentData$resp2 <- ifelse(currentData$Response == "correct", 1,0)
    currentAB <- subset(currentData, trialType == "AB")
    ABacc <- mean(currentAB$resp2)
    currentDE <- subset(currentData, trialType == "DE")
    DEacc <- mean(currentDE$resp2)
    currentAC <- subset(currentData, trialType == "CA")
    ACacc <- mean(currentAC$resp2)
    summary[i,0:6] <- c(rat, sessNum, 'fullTask',ABacc, ACacc, DEacc)
    i=i+1
    
  }
}

#prepare wide form data frames
fullTaskW <- data.frame(summary) #making it into a data frame

fullTaskW$AB_acc <- as.character(fullTaskW$AB_acc)
fullTaskW$AB_acc <- as.numeric(fullTaskW$AB_acc)
fullTaskW$AC_acc <- as.character(fullTaskW$AC_acc)
fullTaskW$AC_acc <- as.numeric(fullTaskW$AC_acc)
fullTaskW$DE_acc <- as.character(fullTaskW$DE_acc)
fullTaskW$DE_acc <- as.numeric(fullTaskW$DE_acc)
fullTaskW$condition <- ifelse(is.element(fullTaskW$ratID, sapGroup) ,'sap','control') #adding lesion condition
fullTaskW$round <- recode(fullTaskW$session, "c(38,39,40,41)='round1';else='round2'")
fullTaskW$drug <- ifelse(is.element(fullTaskW$ratID,drug1First),
                         ifelse(fullTaskW$round=='round1','drug1','drug2'),
                         ifelse(fullTaskW$round=='round1','drug2','drug1'))
fullTaskW$drug <- factor(fullTaskW$drug)
fullTaskW$day <- recode(fullTaskW$session, "c(38,46)='full1';c(39,47)='full2';c(40,48)='full3';c(41,49)='full4'") 
fullTaskW$condition <- factor(fullTaskW$condition)

#making a long form data set
z = rowsNeeded*3 
fullTaskL <- reshape(fullTaskW, varying = c('AB_acc', 'AC_acc', 'DE_acc'), v.names="accuracy",
timevar='trialType',times=c('AB', 'CA', 'DE'), new.row.names=1:z, direction = 'long')

#calculating difference scores in wide data frame
fullTaskW$PI_effect <- fullTaskW$DE_acc - fullTaskW$AC_acc
fullTaskW$AB_CA_diff <- fullTaskW$AB_acc - fullTaskW$AC_acc
#fullTaskW$AB_DE_diff <- fullTaskW$AB_acc - fullTaskW$DE_acc




#confirm data frames look right
str(AB_acquisitionW)
str(fullTaskW)
str(fullTaskL)
xtabs(~ drug+round, data=fullTaskW)
aggregate(accuracy~trialType+drug+condition, data=fullTaskL, FUN=(mean))


#summary tables for making certain plots
summaryTable <- aggregate(accuracy~trialType+day+condition, data=fullTaskL, FUN=(mean))
tempTable <- aggregate(accuracy~trialType+day, data=fullTaskL, FUN=(sd))
summaryTable<- as.data.frame(summaryTable)
summaryTable$sd <- tempTable$accuracy*100
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
summaryTable$se <- tempTable$accuracy*100
summaryTable$accuracy <- summaryTable$accuracy * 100


#checking how accuracies are distributed
hist(fullTaskL$accuracy)

#exploratory plots
p1= ggplot(data=fullTaskL, aes(day, accuracy))
p1+geom_boxplot(aes(col=trialType))

p2= ggplot(data=fullTaskL, aes(drug, accuracy))
p2+geom_boxplot(aes(col=trialType))

p3= ggplot(data=fullTaskL, aes(condition, accuracy))
p3+geom_boxplot(aes(col=trialType))

p4= ggplot(data=fullTaskL, aes(round, accuracy))
p4+geom_boxplot(aes(col=trialType))

p4.1= ggplot(data=fullTaskL, aes(round, accuracy))
p4.1+geom_boxplot(aes(col=trialType))+facet_grid(drug ~ .)

p5= ggplot(data=fullTaskL, aes(day, accuracy))
p5+geom_boxplot(aes(col=trialType))+facet_grid(drug ~ .)

p6= ggplot(data=fullTaskL, aes(day, accuracy))
p6+geom_boxplot(aes(col=trialType))+facet_grid(condition ~ drug)


#inferential stats

anova1 <- aov_ez("ratID","accuracy",fullTaskL, within = c("day","drug","trialType"),
       between = c("condition"))
summary(anova1)

AC_DE_only <- subset(fullTaskL, trialType != "AB")

a2 <- aov_ez("ratID","accuracy",AC_DE_only, within = c("day","drug","trialType"),
             between = c("condition"))
summary(a2)
#even if I'm only looking at AC vs. DE, there is a significant effect of trialType

#shows a significant main effect for trial type, and a significant 4 way interaction
#no 3 way interaction (not significant if you don't count the day)
#should do simple main effects - I first am intersted in whether there is a
# 3 way interaction in any of the days



#to explore the significant main effect of trialType, pairwise comparisons to see what's 
#driving the difference - using adjusted p values
lsmeans(anova1, "trialType", contr = "pairwise", adjust = "holm")
#shows that CA is significantly lower than either of the other two - in particular it's much
#lower then AB, but there is no significant difference between AB and DE
lsmeans(anova1, "day", contr = "pairwise", adjust = "holm")

#probing the interaction in a few ways
#coarse simple main effect of trialType at each level of condition
testInteractions(anova1$lm, across="trialType",fixed=c("condition"),idata=anova1$Anova$idata)
#suggests that the effect of trialType was only significant for the controls, marginally
#significant for the sap rats

#finer simple main effect by condition and drug
testInteractions(anova1$lm, across="trialType",fixed=c("condition","drug"),idata=anova1$Anova$idata)
#suggests that trialType only had a significan effect during drug 1, but both lesion condition
#did have it then

