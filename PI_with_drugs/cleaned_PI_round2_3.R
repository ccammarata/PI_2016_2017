options(contrasts=c("contr.sum", "contr.poly"))
setwd("~/GitHub/PI_2016_2017/PI_with_drugs")
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
colnames(AB_acquisition)<- c('ratID', 'session', 'task','accuracy')

i = 1
for (sessNum in sessList){
  for (rat in ratList) {
    filename <- paste("U3_",rat,"_",sessNum,".csv", sep="")
    print(filename)
    currentData <- read.csv(filename, header=TRUE, sep="\t")
    currentData$resp2 <- ifelse(currentData$Response == "correct", 1,0)
    currentAB <- subset(currentData, trialType == "AB")
    ABacc <- mean(currentAB$resp2)
    AB_acquisition[i,0:4] <- c(rat, sessNum, 'AB_aq', ABacc)
    i=i+1
  }
}
#making the matrix into a data frame and adding lesion condition
AB_acquisitionW <- data.frame(AB_acquisition) #making it into a data frame
AB_acquisitionW$ratID <- factor(AB_acquisitionW$ratID) #making rat number a factor
#AB_acquisitionW$condition <- ifelse(is.element(AB_acquisitionW$ratID, sapGroup) ,'sap','control') #adding lesion condition
AB_acquisitionW$accuracy <- as.character(AB_acquisitionW$accuracy)
AB_acquisitionW$accuracy <- as.numeric(AB_acquisitionW$accuracy)

ABsummaryTable <- aggregate(accuracy~session, data=AB_acquisitionW, FUN=(mean))
tempTable <- aggregate(accuracy~session, data=AB_acquisition, FUN=(sd))
ABsummaryTable<- as.data.frame(ABsummaryTable)
ABsummaryTable$sd <- tempTable$accuracy*100
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
ABsummaryTable$se <- tempTable$accuracy*100
ABsummaryTable$accuracy <- ABsummaryTable$accuracy * 100

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
#fullTaskW$condition <- ifelse(is.element(fullTaskW$ratID, sapGroup) ,'sap','control') #adding lesion condition
fullTaskW$round <- recode(fullTaskW$session, "c(38,39,40,41)='round1';else='round2'")
fullTaskW$drug <- ifelse(is.element(fullTaskW$ratID,drug1First),
                         ifelse(fullTaskW$round=='round1','drug1','drug2'),
                         ifelse(fullTaskW$round=='round1','drug2','drug1'))
fullTaskW$drug <- factor(fullTaskW$drug)
fullTaskW$day <- recode(fullTaskW$session, "c(38,46)='full1';c(39,47)='full2';c(40,48)='full3';c(41,49)='full4'") 
#fullTaskW$condition <- factor(fullTaskW$condition)

#making a long form data set
z = rowsNeeded*3 
fullTaskL <- reshape(fullTaskW, varying = c('AB_acc', 'AC_acc', 'DE_acc'), v.names="accuracy",
                     timevar='trialType',times=c('AB', 'CA', 'DE'), new.row.names=1:z, direction = 'long')

summaryTable <- aggregate(accuracy~trialType+drug, data=fullTaskL, FUN=(mean))
tempTable <- aggregate(accuracy~trialType, data=fullTaskL, FUN=(sd))
summaryTable<- as.data.frame(summaryTable)
summaryTable$sd <- tempTable$accuracy
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
summaryTable$se <- tempTable$accuracy
summaryTable$accuracy <- summaryTable$accuracy

# analysis if AB acquisition
hist(AB_acquisitionW$accuracy)
describe(AB_acquisitionW)
plot(accuracy~session, data=AB_acquisitionW)
#doesn't indicate systamatic change over the days
#in both rounds they learned the pair withing the fist dat, although they 
#didn't do as well on the first day of the first pair

#anlaysis of full task
hist(fullTaskL$accuracy)
describe(fullTaskL)
describe(fullTaskW)

qqnorm(fullTaskL$accuracy)
qqline(fullTaskL$accuracy)
#normality of residuals seems fine
anova1 <- aov_ez("ratID","accuracy",fullTaskL, within = c("day","round","trialType"))
p1= ggplot(data=fullTaskL, aes(trialType, accuracy))
p1+geom_boxplot(aes(col=trialType))+facet_grid(round ~ .)

#doing an ANOVA to check sphericity/HOV
aov2 <- ezANOVA(data=fullTaskL, dv=accuracy, wid = ratID,within=c(trialType,day,drug), detailed=TRUE,return_aov = TRUE)
aov2$`Mauchly's Test for Sphericity` #sphericity not significant
anova2 <- aov_ez("ratID","accuracy",fullTaskL, within = c("day","drug","trialType"))
summary(anova2)#only significant effect is trial type; highly significant
lsmeans(anova2, "trialType", contr = "pairwise", adjust = "holm")#significant differences between AB vs. Ca and DE vs. Ca
p2= ggplot(data=fullTaskL, aes(trialType, accuracy))
p2+geom_boxplot(aes(fill=trialType))+facet_grid(drug ~ .)+
  coord_cartesian(ylim=c(0,1))+
  ggtitle("Accuracy by Trial type and drug") +
  labs(x="Trial Type",y="Accuracy")+
  scale_fill_manual(values=wes_palette("Cavalcanti"))+
  guides(fill=guide_legend(title="trial type"))

ggplot(summaryTable, aes(x=trialType, y=accuracy, fill=trialType)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  facet_wrap(drug ~ .)+
  coord_cartesian(ylim=c(0,1))+
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=wes_palette("Royal1"))+
  ggtitle("Accuracy by trial trype")+
  guides(fill=guide_legend(title="trial type")) +
  labs(x="Trial Type",y="Accuracy")
