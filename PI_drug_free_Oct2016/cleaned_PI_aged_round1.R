options(contrasts=c("contr.sum", "contr.poly"))
library(foreign)
library(psych)
library(reshape2)
library(ggplot2)
library(car)
library(ez)
library(Rmisc)
library(gmodels)
library(wesanderson)
library(afex)

# for AB acquisition
sessList <- list(23,24,25)#list session nubmers here
ratList <- list(25,35,46,49,53,57,28,36,50,63,64)#list rat number here
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
#ABsummary$condition <- ifelse(is.element(ABsummary$ratID, sapGroup) ,'sap','control')

ABsummaryTable <- aggregate(accuracy~session, data=ABsummary, FUN=(mean))
tempTable <- aggregate(accuracy~session, data=ABsummary, FUN=(sd))
ABsummaryTable<- as.data.frame(ABsummaryTable)
ABsummaryTable$sd <- tempTable$accuracy*100
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
ABsummaryTable$se <- tempTable$accuracy*100
ABsummaryTable$accuracy <- ABsummaryTable$accuracy * 100


#for full task
sessList <- list(27, 28, 29, 30,31)#list session nubmers here
ratList <- list(25,35,46,49,53,57,28,36,50,63,64)#list rat number here
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
#FULLfinalSummaryW$condition <- ifelse(is.element(FULLfinalSummaryW$ratID, sapGroup) ,'sap','control')
z = rowsNeeded*3 
FULLfinalSummaryL <- reshape(FULLfinalSummaryW, varying = c('AB_acc', 'AC_acc', 'DE_acc'), v.names="accuracy",
                         timevar='trialType',times=c('AB', 'CA', 'DE'), new.row.names=1:z, direction = 'long')
FULLfinalSummaryL$ratID <-as.factor(FULLfinalSummaryL$ratID)
FULLfinalSummaryL$session <-as.factor(FULLfinalSummaryL$session)
FULLfinalSummaryL$trialType <- as.factor(FULLfinalSummaryL$trialType)
FULLfinalSummaryW$ratID <-as.factor(FULLfinalSummaryW$ratID)
FULLfinalSummaryW$session <-as.factor(FULLfinalSummaryW$session)

#making summary table for full task
FULLsummaryTable <- aggregate(accuracy~trialType, data=FULLfinalSummaryL, FUN=(mean))
tempTable <- aggregate(accuracy~trialType, data=FULLfinalSummaryL, FUN=(sd))
FULLsummaryTable<- as.data.frame(FULLsummaryTable)
FULLsummaryTable$sd <- tempTable$accuracy
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
FULLsummaryTable$se <- tempTable$accuracy
FULLsummaryTable$accuracy <- FULLsummaryTable$accuracy


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
p1= ggplot(data=FULLfinalSummaryL, aes(trialType,accuracy))
p1+geom_boxplot(aes(fill=trialType)) +
  coord_cartesian(ylim=c(0,.8))+
  ggtitle("Accuracy by Trial type") +
  labs(x="Trial Type",y="Accuracy")+
  scale_fill_manual(values=wes_palette("Cavalcanti"))+
  guides(fill=guide_legend(title="trial type"))


#inferential tests of full task
qqnorm(FULLfinalSummaryL$accuracy)
qqline(FULLfinalSummaryL$accuracy) #to check for normality of residuals
aov1 <- ezANOVA(data=FULLfinalSummaryL, dv=accuracy, wid = ratID,within=c(trialType,session),detailed=TRUE,return_aov = TRUE)
aov1$`Mauchly's Test for Sphericity` #to quckly check sphericity, which is fine in this case
#I used ezANOVA to get the sphericity test but I can't figure out how to do post hoc analyses with that output
aov2 <- aov_ez("ratID","accuracy",data=FULLfinalSummaryL,within = c("trialType","session"))
summary(aov2)
lsmeans(aov2, "trialType", contr = "pairwise", adjust = "holm")

ggplot(FULLsummaryTable, aes(x=trialType, y=accuracy, fill=trialType)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  coord_cartesian(ylim=c(0,.8))+
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=wes_palette("Royal1"))+
  ggtitle("Accuracy by trial trype")+
  guides(fill=guide_legend(title="trial type")) +
  labs(x="Trial Type",y="Accuracy")

FULLsummaryTable <- aggregate(accuracy~trialType+session, data=FULLfinalSummaryL, FUN=(mean))
tempTable <- aggregate(accuracy~trialType+session, data=FULLfinalSummaryL, FUN=(sd))
FULLsummaryTable<- as.data.frame(FULLsummaryTable)
FULLsummaryTable$sd <- tempTable$accuracy
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
FULLsummaryTable$se <- tempTable$accuracy
FULLsummaryTable$accuracy <- FULLsummaryTable$accuracy

subSummaryTable <-subset(FULLsummaryTable, trialType != 'AB')

ggplot(subSummaryTable, aes(x=session,y=accuracy,colour = trialType))+
  geom_line(aes(group=trialType))+geom_point()+
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se),width=.05)

ggplot(FULLsummaryTable, aes(x=session,y=accuracy,colour = trialType))+
  geom_line(aes(group=trialType))+
  geom_point()+geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se),width=.1,position=position_dodge(.05))+
  coord_cartesian(ylim=c(.4,.8))+
  guides(colour=guide_legend(title="trial type"))+
  scale_color_manual(values=wes_palette("Royal1"))
  



