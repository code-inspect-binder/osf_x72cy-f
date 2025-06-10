####################################################################
## Preprocessing script for data from the Qualtrics survey format ##
####################################################################

rm(list = ls())

###########################################################
## 0. Set Working Directory And Load Relevant R Packages ##
###########################################################
# setwd("...")
library(qualtRics) # For loading data directly from the qualtrics exported format


#################################################################################
## 1. Load the existing data for the selected studies (from Camerer and Klein) ##
#################################################################################
replication.results     <- na.omit(read.csv2(file = "Data/OverviewDescriptions.csv",stringsAsFactors = FALSE))[1:27,] # get the results for the 27 studies 
replication.outcomes    <- as.numeric(replication.results$RS.outcome)                               # 1 = replicated, 0 = did not replicate
replication.effectsizes <- as.numeric(replication.results$RS.effect.size..correlation.coefficient.) # effect sizes as correlation coefficient (r)


########################################################
## 2. Load the data from Qualtrics for the lay people ##
########################################################
raw_data <- readSurvey("Data/QualtricsTest.csv")

# Specify some exclusion criteria:
prediction.data <- raw_data[raw_data$DistributionChannel=="anonymous",]


# data structure: 
# if condition = 1:  description only [block labeL:des.]
# if condition = 2:  description + BF [block label: bf.]
# understanding      = studynumber_label_0 
# replication belief = studynumber_questionlabel_belief
# confidence         = studynumber_questionlabel_conf_1

nstudies   <- 27                    # we have 26 studies + bogus study 
nsubjects  <- nrow(prediction.data) # the number of participants

# extract relevant data & calculate Brier score for each participant
lay.data     <- NULL
for (i in 1:nsubjects) {
  study.order  <- as.numeric(unlist(strsplit(prediction.data$preorder[i], split="|", fixed = TRUE)))[2:28] 
  label        <- ifelse(prediction.data$Conditie[i]==1,"des.","bf.")
  subject.data <- NULL
  for (j in 1:nstudies){
    understanding          <- eval(parse(text = paste0("prediction.data$`",study.order[j],"_",label,"0`"  )))[i]           # 1 = did not understand, otherwise: NA 
    replication.belief     <- eval(parse(text = paste0("prediction.data$`",study.order[j],"_",label,"belief`"  )))[i] - 1  # now 0 = will not be replicated, 1 = will be replicated 
    confidence.rating      <- eval(parse(text = paste0("prediction.data$`",study.order[j],"_",label,"conf_1`")))[i]        # on a scale from 0 - 100
    confidence.rating      <- ifelse(replication.belief == 0, confidence.rating*-1, confidence.rating)                     # make confidence in replication failure negative
    confidence.rating      <- confidence.rating / 200 + .5                                                                 # convert to 0-1 scale
    study                  <- j
    condition              <- ifelse(label == "des.","DescriptionOnly","DescriptionPlusStatistics")
    replication.outcome    <- replication.outcomes[j]
    replication.effectsize <- replication.effectsizes[j]
    ind.subject.data       <- cbind(study,condition,understanding,replication.belief,confidence.rating,
                                    replication.outcome,replication.effectsize)
    subject.data           <- rbind(subject.data,ind.subject.data)
  }
  subject      <- i
  ind.lay.data <- cbind(subject,subject.data)
  lay.data     <- rbind(lay.data,ind.lay.data)
}
rm(study.order,label,condition,subject.data,understanding,replication.belief,confidence.rating,
   i,j,study,replication.outcome,replication.effectsize,ind.subject.data,subject,ind.lay.data)


##############################################
## 3. Data exclusion based on set criteria  ##
##############################################
# criteria:  
# 1. if participants failed the attention check (i.e., did not press 'NO' and 75% (range 70-80 is allowed))
# 2. if a study description is not understood, exlcude this study for this participant
# 3. if a study is not understood by > 50% of the participants, exclude this study
# 4. if a participant does not understand > 50% of the studies, exclude this participant

bogus.study     <- 27
clean.data      <- as.data.frame(lay.data, stringsAsFactors = FALSE)
correct.range   <- clean.data[clean.data$study==bogus.study,]$confidence.rating
correct.range   <- rep(correct.range >= .1 & correct.range <= .15, each = nstudies) # correct range NO and 70-80% --> .1-.15 on the confidence scale
clean.data.1    <- clean.data[correct.range,]                                       # apply 1.
clean.data.2    <- clean.data.1[is.na(clean.data.1$understanding),]                 # apply 2.
remove.studies  <- which(table(clean.data.2$study)<nsubjects/2)                     # indicate which studies are understood by less than half of the people
clean.data.3    <- clean.data.2[! clean.data.2$study %in% remove.studies,]          # apply 3.
remove.subjects <- which(table(clean.data.3$subject)<nstudies/2)                    # indicate which subjects understood less than half of the studies
clean.data.4    <- clean.data.3[! clean.data.3$subject %in% remove.subjects,]       # apply 4. 

full.data               <- clean.data.4
full.data$understanding <- NULL         # delete empty column 

# set classes of columns in dataframe, because they got messed up 
factor.cols             <- c("subject","study","condition")
numeric.cols            <- c("replication.belief","confidence.rating","replication.outcome","replication.effectsize")
full.data[factor.cols]  <- lapply(full.data[factor.cols], as.factor)
full.data[numeric.cols] <- lapply(full.data[numeric.cols], as.numeric)

rm(clean.data,clean.data.1,clean.data.2,clean.data.3, factor.cols, numeric.cols)

# recount remaining number of studies and subjects 
nstudies  <- length(unique(full.data$study))
nsubjects <- length(unique(full.data$subject))


###################################################
## 4. Calculate the Brier score per participant  ##
###################################################
for (i in 1:nsubjects){
  subject.data   <- subset(full.data, subject == i)
  full.data$brier.score[full.data$subject == i] <- BrierScore(subject.data$replication.outcome, subject.data$confidence.rating)
}
rm(i,subject.data)

out.full                        <- full.data
out.brierscores                 <- aggregate(brier.score ~ condition + subject, data = full.data, FUN = mean)
out.brierscores$log.brier.score <- log(out.brierscores$brier.score)


##########################
## 4. Export Data Files ##
##########################
write.csv(out.brierscores, file = "TestData_BrierScores.csv", row.names = FALSE) # Export data for t-test in JASP (RQ 1)
write.csv(out.full,        file = "TestData_Full.csv", row.names = FALSE)        # Export data for R analysis (quality check, RQ 2, RQ 3)






           