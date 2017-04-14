options(contrasts=c("contr.sum", "contr.poly"))
library(foreign)
library(psych)
library(reshape2)
library(ggplot2)
library(wesanderson)
library(afex)


# reading in the data from individual rat/day CSVs 
sessList <- list(27, 28, 29, 30,31)
ratList <- list(25,35,46,49,53,57,28,36,50,63,64)
x <- length(sessList)
y <- length(ratList)
rowsNeeded <- (x*y)
sapGroup = list(28, 36,50,63,64)
conGroup = list(25,35,46,49,53,57)

fullTask <- matrix(data = NA, nrow = rowsNeeded, ncol=5)
colnames(fullTask)<- c('ratID', 'session', 'AB_acc', 'AC_acc', 'DE_acc')

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
    fullTask[i,0:5] <- c(rat, sessNum, ABacc, ACacc, DEacc)
    i=i+1
    
  }
}

#prepare wide and long form data frames for full task
fullTaskWidedf <- data.frame(fullTask)
z = rowsNeeded*3 
fullTaskLongdf <- reshape(fullTaskWidedf, varying = c('AB_acc', 'AC_acc', 'DE_acc'), v.names="accuracy",
                             timevar='trialType',times=c('AB', 'CA', 'DE'), new.row.names=1:z, direction = 'long')
fullTaskLongdf$ratID <-as.factor(fullTaskLongdf$ratID)
fullTaskLongdf$session <-as.factor(fullTaskLongdf$session)
fullTaskLongdf$trialType <- as.factor(fullTaskLongdf$trialType)
fullTaskWidedf$ratID <-as.factor(fullTaskWidedf$ratID)
fullTaskWidedf$session <-as.factor(fullTaskWidedf$session)


#boxplot of accuracy by trial type, uses long data frame
p1= ggplot(data=fullTaskLongdf, aes(trialType,accuracy))
p1+geom_boxplot(aes(fill=trialType)) +
  coord_cartesian(ylim=c(0,.8))+
  ggtitle("Accuracy by Trial type") +
  labs(x="Trial Type",y="Accuracy")+
  scale_fill_manual(values=wes_palette("Cavalcanti"))+
  guides(fill=guide_legend(title="trial type"))

# making a summary table for the line graphs
FULLsummaryTable <- aggregate(accuracy~trialType+session, data=fullTaskLongdf, FUN=(mean))
tempTable <- aggregate(accuracy~trialType+session, data=fullTaskLongdf, FUN=(sd))
FULLsummaryTable<- as.data.frame(FULLsummaryTable)
FULLsummaryTable$sd <- tempTable$accuracy
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
FULLsummaryTable$se <- tempTable$accuracy
FULLsummaryTable$accuracy <- FULLsummaryTable$accuracy

ggplot(FULLsummaryTable, aes(x=session,y=accuracy,colour = trialType))+
  geom_line(aes(group=trialType))+
  geom_point()+geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se),width=.1,position=position_dodge(.05))+
  coord_cartesian(ylim=c(.4,.8))+
  guides(colour=guide_legend(title="trial type"))


# making a summary table without session
### for some reason if I do a barchart from a summary table that includes sesison, but I don't have session as one of the variables
#### in my barchart (i.e. there's a variable in my summary table that I collapse over), then my plots come out weird because
#### it trieds to sort of make the graph twice and smoosh them on top of one another... so I've just taken to making a version that
#### has only the variable I need

NoSESSsummaryTable <- aggregate(accuracy~trialType, data=fullTaskLongdf, FUN=(mean))
tempTable <- aggregate(accuracy~trialType, data=fullTaskLongdf, FUN=(sd))
NoSESSsummaryTable<- as.data.frame(NoSESSsummaryTable)
NoSESSsummaryTable$sd <- tempTable$accuracy
tempTable$accuracy <- tempTable$accuracy/sqrt(y)
NoSESSsummaryTable$se <- tempTable$accuracy




ggplot(NoSESSsummaryTable, aes(x=trialType, y=accuracy, fill=trialType)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  coord_cartesian(ylim=c(0,.8))+
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values=wes_palette("Royal1"))+
  ggtitle("Accuracy by trial trype")+
  guides(fill=guide_legend(title="trial type")) +
  labs(x="Trial Type",y="Accuracy")


