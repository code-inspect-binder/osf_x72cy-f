########################################################################
##   Preprocessing script for data from the Qualtrics survey format   ##
##                                                                    ##
##    0. Set Working Directory And Load Relevant R Packages           ##
########################################################################

  # rm(list = ls())
  # setwd('...') # set working directory to 'Analysis_Code' folder
  library(DescTools) # to calculate the Brier score

#################################################################################
## 1. Load the existing data for the selected studies (from Camerer and Klein) ##
#################################################################################

  replication.results     <- na.omit(read.csv(file = '../Data/Processed_Data/DescriptionAndStatisticsPerStudy.csv',stringsAsFactors = FALSE))[1:27,] # get the results for the 27 studies 
  replication.outcomes    <- as.numeric(replication.results$RS.outcome)                               # 1 = replicated, 0 = did not replicate
  replication.effectsizes <- as.numeric(replication.results$RS.effect.size..correlation.coefficient.) # effect sizes as correlation coefficient (r)

########################################################
## 2. Load the data from Qualtrics for the lay people ##
########################################################

all.data <- NULL
for (s in 1:3){
  if (s == 1){
    raw_data <- read.csv('../Data/Raw_Data/Data_UvA_Lab_AN.csv', header = TRUE, stringsAsFactors = FALSE) #for the UvA lab
    subject.offset <- 0 
  }
  if (s == 2){
    raw_data <- read.csv('../Data/Raw_Data/Data_MTurk_AN.csv', header = TRUE, stringsAsFactors = FALSE) #for MTurk
    subject.offset <- nsubjects 
  }
  if (s == 3){
    raw_data <- read.csv('../Data/Raw_Data/Data_Social_Media_AN.csv', header = TRUE, stringsAsFactors = FALSE) # for social media 
    subject.offset <- last.subject
  }
  # Specify some exclusion criteria:
  prediction.data <- raw_data[raw_data$DistributionChannel=='anonymous',]           # Exclude previews and test generations from Qualtrics 
  prediction.data <- prediction.data[!is.na(prediction.data$knowledge),]            # Include: finished 
  prediction.data <- prediction.data[prediction.data$knowledge==5,]                 # Include: not familiar with studies 

  # data structure: 
  # if condition = 1:  description only [block labeL: des.]
  # if condition = 2:  description + BF [block label: bf.]
  # understanding      = studynumber_conditionlabel_0 
  # replication belief = studynumber_conditionlabel_belief
  # confidence         = studynumber_conditionlabel_conf_1
  # for lab website and social media version there are separate Dutch blocks 1=English, 2=Dutch
  
  nstudies   <- 28                    # we have 27 studies + bogus study 
  nsubjects  <- nrow(prediction.data) # the number of participants
  
  # extract relevant data for each participant
  lay.data     <- NULL
  for (i in 1:nsubjects) {
    if(s == 1){
      studied.psy <- "TRUE_UvA"
    } else {
      studied.psy <- ifelse(prediction.data$psych[i] == 1, "TRUE", "FALSE")
    }
    study.order  <- as.numeric(unlist(strsplit(prediction.data$preorder[i], split='|', fixed = TRUE)))[2:29] 
    label        <- ifelse(prediction.data$Conditie[i]==1,'des.','bf.')                                             # between-subjects condition (description only or BF)
    language     <- ifelse(prediction.data$language[i]==2,'NL','')                                                  # display language, if 2=Dutch, 'NL' is added to the question
    duration     <- prediction.data$Duration..in.seconds.[i]
    first.study  <- strsplit(prediction.data$preorder[i], "[|]")[[1]][2]                                            # gets number of first study (2nd object, because sequence starts with '|')
    if (first.study == "28") {first.study <- strsplit(prediction.data$preorder[i], "[|]")[[1]][3]}
    last.study   <- tail(strsplit(prediction.data$preorder[i], "[|]")[[1]], n = 1)
    if (last.study == "28") {last.study <- tail(strsplit(prediction.data$preorder[i], "[|]")[[1]], n = 2)[1]}
    subject.data <- NULL
    for (j in 1:nstudies){
      understanding          <- eval(parse(text = paste0('prediction.data$X',j,'_',label,'0',language)))[i]         # 1 = did not understand, otherwise: NA 
      replication.belief     <- eval(parse(text = paste0('prediction.data$X',j,'_',label,'belief',language)))[i]    # now 0 = will not be replicated, 1 = will be replicated 
      replication.belief     <- ifelse(replication.belief==2,0,replication.belief)                                  # change 2 to 0 --> belief in failure
      confidence.rating      <- eval(parse(text = paste0('prediction.data$X',j,'_',label,'conf',language,'_1')))[i] # on a scale from 0 - 100
      confidence.rating      <- ifelse(replication.belief == 0, confidence.rating*-1, confidence.rating)            # make confidence in replication failure negative
      confidence.rating      <- confidence.rating / 200 + .5                                                        # convert to 0-1 scale
      study                  <- j
      condition              <- ifelse(label == 'des.','DescriptionOnly','DescriptionPlusEvidence')
      replication.outcome    <- replication.outcomes[j]
      replication.effectsize <- replication.effectsizes[j]
      ind.subject.data       <- cbind(study,condition,understanding,studied.psy,replication.belief,confidence.rating,
                                      replication.outcome,replication.effectsize)
      subject.data           <- rbind(subject.data,ind.subject.data)
    }
    subject      <- subject.offset + i 
    ind.lay.data <- cbind(subject,subject.data, first.study, last.study, duration)
    lay.data     <- rbind(lay.data,ind.lay.data)
  }
  last.subject <- subject
  rm(study.order,label,condition,subject.data,understanding,replication.belief,confidence.rating,
     i,j,study,replication.outcome,replication.effectsize,ind.subject.data,subject,ind.lay.data, language, first.study, last.study)
  
  all.data  <- rbind(all.data,lay.data)
}
rm(last.subject,s, lay.data, prediction.data,raw_data, subject.offset)

##############################################
## 3. Data exclusion based on set criteria  ##
##############################################

# criteria:  
# 1. if participants failed the attention check (i.e., did not press 'NO' and 75% (range 70-80 is allowed))
# 2. if a study description is not understood, exlcude this study for this participant
# 3. if a study is not understood by > 50% of the participants, exclude this study
# 4. if a participant does not understand > 50% of the studies, exclude this participant
  
  bogus.study              <- 28
  clean.data               <- as.data.frame(all.data, stringsAsFactors = FALSE)
  correct.range            <- clean.data[clean.data$study==bogus.study,]$confidence.rating
  correct.range            <- rep(correct.range >= .1 & correct.range <= .15, each = nstudies)    # correct range NO and 70-80% --> .1-.15 on the confidence scale
  clean.data.1             <- clean.data[correct.range,]                                          # apply 1.
  incomprehensible.studies <- sum(as.numeric(clean.data.1$understanding), na.rm = TRUE)           # count how many studies were not understood in total across all participants 
  subs.incomp.studies      <- length(unique(clean.data.1$subject[clean.data.1$understanding==1])) # count number of subjects who indicated not to understand at least one description
  clean.data.2             <- clean.data.1[is.na(clean.data.1$understanding),]                    # apply 2.
  remove.studies           <- which(table(clean.data.2$study)<nsubjects/2)                        # indicate which studies are understood by less than half of the people
  clean.data.3             <- clean.data.2[! clean.data.2$study %in% remove.studies,]             # apply 3.
  remove.subjects          <- which(table(clean.data.3$subject)<nstudies/2)                       # indicate which subjects understood less than half of the studies
  clean.data.4             <- clean.data.3[! clean.data.3$subject %in% remove.subjects,]          # apply 4. 
  
  full.data               <- clean.data.4
  full.data$understanding <- NULL         # delete empty column 
  
  # set classes of columns in dataframe, because they got messed up 
  factor.cols             <- c('condition', 'studied.psy')
  numeric.cols            <- c('subject','study','replication.belief','confidence.rating','replication.outcome','replication.effectsize','first.study','last.study','duration')
  full.data[factor.cols]  <- lapply(full.data[factor.cols], as.factor)
  full.data[numeric.cols] <- lapply(full.data[numeric.cols], as.numeric)
  
  rm(clean.data,clean.data.1,clean.data.2,clean.data.3, clean.data.4,factor.cols, numeric.cols, bogus.study,correct.range)
  
  # recount remaining number of studies and subjects 
  nstudies  <- length(unique(full.data$study))
  nsubjects <- length(unique(as.numeric(full.data$subject)))
  
  # remove rows for the bogus study
  full.data <- full.data[!full.data$study==28,]
  full.data$guessed.correctly <- ifelse(full.data$replication.belief == full.data$replication.outcome, 1, 0)
  
  # supplemental analysis: remove psychology students in MTurk and Social Media subset
  full.data.no.psy <- full.data[full.data$studied.psy != "TRUE", ]

###################################################
## 4. Calculate the Brier score per participant  ##
###################################################

# add consecutive number to the data
  subs <- sort(unique(full.data$subject))
  for (i in seq_along(subs)){
    full.data$sub.number[full.data$subject==subs[i]] <- i
  }
  
  for (i in 1:nsubjects){
    subject.data   <- subset(full.data, sub.number == i)
    full.data$brier.score[  full.data$sub.number == i] <- BrierScore(subject.data$replication.outcome, subject.data$confidence.rating)
    full.data$total.correct[full.data$sub.number == i] <- sum(subject.data$guessed.correctly)
  }
  rm(i,subject.data)
  
  out.full                        <- full.data
  out.brierscores                 <- aggregate(brier.score ~ condition + subject, data = full.data, FUN = mean)
  total.correct                   <- aggregate(total.correct ~ subject, data = full.data, FUN = mean) 
  out.brierscores$total.correct   <- total.correct$total.correct 
  out.brierscores$log.brier.score <- log(out.brierscores$brier.score)
  rm(total.correct)
  
# add consecutive number to the data for full.data.no.psy data set
  subs <- sort(unique(full.data.no.psy$subject))
  for (i in seq_along(subs)){
    full.data.no.psy$sub.number[full.data.no.psy$subject==subs[i]] <- i
  }
  
  for (i in 1:nsubjects){
    subject.data   <- subset(full.data.no.psy, sub.number == i)
    full.data.no.psy$brier.score[  full.data.no.psy$sub.number == i] <- BrierScore(subject.data$replication.outcome, subject.data$confidence.rating)
    full.data.no.psy$total.correct[full.data.no.psy$sub.number == i] <- sum(subject.data$guessed.correctly)
  }
  rm(i,subject.data)
  
  out.full.no.psy                        <- full.data.no.psy
  out.brierscores.no.psy                 <- aggregate(brier.score ~ condition + subject, data = full.data.no.psy, FUN = mean)
  total.correct                          <- aggregate(total.correct ~ subject, data = full.data.no.psy, FUN = mean) 
  out.brierscores.no.psy$total.correct   <- total.correct$total.correct 
  out.brierscores.no.psy$log.brier.score <- log(out.brierscores.no.psy$brier.score)
  rm(total.correct)

##############################################################
## 5.  Data for first vs. last study comparison & duration  ##
##############################################################
  
# for every subject, get predictions for the first study and for the last study 
# and add the overall completion time per subject 
  order.data <- matrix(nrow = nsubjects, ncol = 7)
  colnames(order.data) <- c('subject','condition','correct.first','correct.last','confidence.first','confidence.last','duration(sec)')
  for (i in 1:nsubjects){
    subject.data     <- subset(full.data, sub.number == i)
    correct.first    <- subject.data$guessed.correctly[subject.data$study == subject.data$first.study]
    if (length(correct.first) == 0) {correct.first <- NA}
    correct.last     <- subject.data$guessed.correctly[subject.data$study == subject.data$last.study ]
    if (length(correct.last) == 0) {correct.last <- NA}
    confidence.first <- subject.data$confidence.rating[subject.data$study == subject.data$first.study]
    if (length(confidence.first) == 0) {confidence.first <- NA}
    confidence.last  <- subject.data$confidence.rating[subject.data$study == subject.data$last.study ]
    if (length(confidence.last) == 0) {confidence.last <- NA}
    order.data[i,]   <- cbind(subject.data$subject[1], subject.data$condition[1], 
                              correct.first, correct.last, confidence.first, confidence.last,subject.data$duration[1])
  }
  order.data <- as.data.frame(order.data)
  order.data$condition <- ifelse(order.data$condition == 1, 'DescriptionOnly','DescriptionPlusEvidence')
  order.data$confidence.first <- abs((order.data$confidence.first -.5) * 2) # transform scale of the confidence ratings from 0,1 to -1,1 and makes absolute to assess strength
  order.data$confidence.last  <- abs((order.data$confidence.last  -.5) * 2)
  out.order.duration <- order.data
  
#############################################################################
## 6. Did MTurk workers and participants in Social Media Study Psychology? ##
#############################################################################
  
  subj.index <- unique(full.data$subject)
  studied.psy.data <- data.frame(subject     = unique(full.data$subject), 
                                 studied.psy = NA)
  for(i in 1:nsubjects){
    studied.psy.data$studied.psy[i] <- full.data[full.data$subject == subj.index[i], ]$studied.psy[1]
  }
  table(studied.psy.data$studied.psy) 
  # 22.5% (i.e., 24 out of 107) of the subjects studied psychology, but did not finish their phD
  
  rm(subj.index, studied.psy.data)

##########################
## 7. Export Data Files ##
##########################

  # write.csv(out.brierscores,        file = '../Data/Processed_Data/FinalData_BrierScores.csv',      row.names = FALSE) # Export data for t-test in JASP (RQ 1)
  # write.csv(out.full,               file = '../Data/Processed_Data/FinalData_Full.csv',             row.names = FALSE) # Export data for R analysis (quality check, RQ 2, RQ 3)

# exporting data for analyses in supplements
  # write.csv(out.order.duration,     file = '../Data/Processed_Data/FinalData_OrderAndDuration.csv', row.names = FALSE) # Export data for extra analysis of first vs. last prediction
  # write.csv(out.brierscores.no.psy, file = '../Data/Processed_Data/FinalData_BrierScoresNoPsy.csv', row.names = FALSE) # Export data subset analysis: t-test in JASP (RQ 1)
  # write.csv(out.full.no.psy,        file = '../Data/Processed_Data/FinalData_FullNoPsy.csv',        row.names = FALSE) # Export data for subset analysis: R analysis (quality check, RQ 2, RQ 3)
  
######################################
## 8. Remove Unnecessary Data files ##
######################################
  
  # delete unnecessary files
  rm(list=ls()[-which(ls() %in% ls(pattern = "out"))])
  
