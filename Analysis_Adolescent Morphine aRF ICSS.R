############### ADOLESCENT MORPHINE ADMINISTRATION ANALYSIS ####################
# This program takes in data processed by the script "RF Theta Calculation"    #
# and run statistics to analyze findings. Data from each experiment are loaded,#
# analyzed, and visualized.                                                    #

# PREPARE R -----
## Set up paths
computerusergithubpath <- c('C:/Users/rmdon/OneDrive/Desktop/GitHub_Repositories/') # Path for Rachel's laptop to github repositories
computeruserboxpath <- c('C:/Users/rmdon/Box/') # Path for Rachel's laptop to Box drive

## Load required libraries
source(paste(computerusergithubpath,"Analysis-PBAdolescentMorphineReward/Functions/Functions_LoadPackages.R",sep=''))

## Load all custom functions
source.all(paste(computerusergithubpath,'JRoitmanICSS/Functions/Abridged RF ICSS/',sep=''))
source.all(paste(computerusergithubpath,'Analysis-PBAdolescentMorphineReward/Functions/',sep=''))

## Set working directory to the raw data folder
setwd(paste(computeruserboxpath,'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Raw Data',sep=''))

## Set up paths to export figures and analysis to
### Overall paths
path_analysis <- c(paste(computeruserboxpath, 'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Analysis/',sep=''))
path_figures <- c(paste(computeruserboxpath, 'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Figures/Figure Panels/',sep=''))

## Specific paths
path_analysis_RRFtest <- paste(path_analysis,"Adolescent Reduced Rate Frequency Morphine Test/",sep="")
path_analysis_RRFtest_means <- paste(path_analysis,"Adolescent Reduced Rate Frequency Morphine Test/Means Tables/",sep="")
path_figures_RRFtest <- paste(path_analysis_RRFtest,"Figures/",sep="")
path_subjectfigures_RRFtest <- paste(path_figures_RRFtest,"Subject Figures/",sep="") # Subject shaping figures path


## Read in subject key
subjectkey <- read_csv("SubjectKey_Adolescent.csv")
subjectexc <- c()


# PREPARE DATA -----
## Load raw data -----
RFadolescentdata_raw <- read_csv("RawData_RRFICSSMorphineTest_Adolescent.csv") # Read in raw data
RFadolescentdata_raw <- RFadolescentdata_raw %>% filter (Phase != "RETEST", !(SubjectID==576 & UniqueSessionID %in% c("BL_1_AM","BL_1_PM"))) # Remove Retest sessions and first two baselines for 576 (wrong amplitude set)

## Prepare variables -----
### Set column names containing lever presses by frequency
freqcolnames <- c('Freq_141', 'Freq_100', 'Freq_89', 'Freq_79', 'Freq_71', 'Freq_63', 'Freq_50') # Set column names of lever presses for each frequency
lococolnames <- c('Loc_141', 'Loc_100', 'Loc_89', 'Loc_79', 'Loc_71', 'Loc_63', 'Loc_50') # Set column names of bream break counts for each frequency

### Add Sex and Group to raw
RFadolescentdata_raw$Sex <- subjectkey$Sex[match(RFadolescentdata_raw$SubjectID, subjectkey$SubjectID)]
RFadolescentdata_raw$Group <- subjectkey$ExperimentalGroup[match(RFadolescentdata_raw$SubjectID, subjectkey$SubjectID)]
RFadolescentdata_raw$ShapingInclude <- subjectkey$ShapingInclude[match(RFadolescentdata_raw$SubjectID, subjectkey$SubjectID)]

### Factor variables
RFadolescentdata_raw$Sex <- factor(RFadolescentdata_raw$Sex, levels = c("M","F"))
RFadolescentdata_raw$Group <- factor(RFadolescentdata_raw$Group, levels = c("AdolS","AdolM"))
RFadolescentdata_raw$SessionTime <- factor(RFadolescentdata_raw$SessionTime, levels = c("AM","PM"))

RFadolescentdata_raw$Phase_raw <- RFadolescentdata_raw$Phase
RFadolescentdata_raw$Phase <- factor(if_else(RFadolescentdata_raw$Phase=="POST" & RFadolescentdata_raw$Day <=3,"POST_EARLY", 
                                                      if_else(RFadolescentdata_raw$Phase=="POST" & RFadolescentdata_raw$Day >=4,"POST_LATE", 
                                                              RFadolescentdata_raw$Phase)), levels = c("BL","INJ","POST_EARLY","POST_LATE"))
### Prepare Session variables
RFadolescentdata_raw$Day_raw <- factor(RFadolescentdata_raw$Day) # Factor day by phase
RFadolescentdata_raw$Day <- ifelse(RFadolescentdata_raw$Phase == "BL", 0, RFadolescentdata_raw$Day) # Create collapsed Day variable where Baseline=0
RFadolescentdata_raw$Day <- factor(RFadolescentdata_raw$Day, levels = unique(RFadolescentdata_raw$Day)) # Factor collapsed Day variable

allsessions <- c(unique(RFadolescentdata_raw$UniqueSessionID))
sessionIDs <- c(0,0,0,0,1:34)
RFadolescentdata_raw$Session <- sessionIDs[match(RFadolescentdata_raw$UniqueSessionID, allsessions)]
RFadolescentdata_raw$Session <- factor(RFadolescentdata_raw$Session)


### Add total locomotion, max lever presses, and total lever presses per pass
RFadolescentdata_raw$sumLoco <- rowSums(as.matrix(RFadolescentdata_raw[lococolnames])) # Add total locomotion counts
RFadolescentdata_raw$maxLP <- rowMaxs(as.matrix(RFadolescentdata_raw[freqcolnames])) # Max lever pressing per pass
RFadolescentdata_raw$sumLP <- rowSums(as.matrix(RFadolescentdata_raw[freqcolnames])) # Max lever pressing per pass

### Add lever presses as a percent of max pass lever presses
freqpercmaxcolnames <- c('PercMax_Freq_141', 'PercMax_Freq_100', 'PercMax_Freq_89', 'PercMax_Freq_79', 'PercMax_Freq_71', 'PercMax_Freq_63', 'PercMax_Freq_50') # Set column names of lever presses for each frequency

RFadolescentdata_raw <- RFadolescentdata_raw %>% mutate(PercMax_Freq_141 = findPercMax(Freq_141,maxLP), PercMax_Freq_100 = findPercMax(Freq_100,maxLP), PercMax_Freq_89 = findPercMax(Freq_89,maxLP),
                                                                PercMax_Freq_79 = findPercMax(Freq_79,maxLP), PercMax_Freq_71 = findPercMax(Freq_71,maxLP), PercMax_Freq_63 = findPercMax(Freq_63,maxLP), 
                                                                PercMax_Freq_50 = findPercMax(Freq_50,maxLP)) # Find percent of max lever presses for each frequency
### Check Variables
unique(RFadolescentdata_raw$SubjectID)
unique(RFadolescentdata_raw$Sex)
unique(RFadolescentdata_raw[grepl('F',RFadolescentdata_raw$Sex),"SubjectID"])
unique(RFadolescentdata_raw[grepl('M',RFadolescentdata_raw$Sex),"SubjectID"])

unique(RFadolescentdata_raw$Group)
unique(RFadolescentdata_raw[grepl('AdolM',RFadolescentdata_raw$Group),"SubjectID"])
unique(RFadolescentdata_raw[grepl('AdolS',RFadolescentdata_raw$Group),"SubjectID"])

unique(RFadolescentdata_raw$ShapingInclude)
unique(RFadolescentdata_raw$Phase) 
unique(RFadolescentdata_raw$Day)
unique(RFadolescentdata_raw$Session)
unique(RFadolescentdata_raw$SessionTime)
#view(RFadolescentdata_raw)

## Plot individual subject passes
# RRFICSS_ploteachRF(RFadolescentdata_raw,"all",freqcolnames,path_figures_subject,"png")

### Clean data -----
RFadolescentdata_cleaned <- RRFICSS_cleandata(RFadolescentdata_raw,freqcolnames) # Clean data


### Fit RF Model - Linear Model -----
RFadolescentdata_LP_raw <- RRFICSS_LM_mediancenter(RFadolescentdata_cleaned,freqcolnames,"raw") # Fit linear model to lever press values

RFadolescentsubject_LP_lm <- RFadolescentdata_LP_raw %>% # Prep data for statistics - find subject average for each day
  filter(Int > 0, Pass != 1, sumLP >=10, ShapingInclude=="INCLUDE", !SubjectID %in% subjectexc) %>%
  group_by(SubjectID, Sex, Group, ShapingInclude, UniqueSessionID, Phase, Phase_raw, Day, Day_raw, Session, SessionTime, Box, Amplitude) %>%
  dplyr::summarize(npass = n_distinct(Pass),
                   mean_Int = mean(Int, na.rm=TRUE),
                   mean_Slope = mean(Coef, na.rm=TRUE))

### Fit RF Model - logistic model -----
#RFadolescentdata_LGM_raw <- RRFICSS_LGM_fitmodel(RFadolescentdata_cleaned,freqcolnames) # Fit logistic model to lever press values

#RFadolescentsubject_LGM <- RFadolescentdata_LGM_raw %>% # Prep data for statistics - find subject average for each day
#   filter(!SubjectID %in% subjectexc, ShapingInclude=="INCLUDE", Pass > 1, sumLP >=10, Slope > 0, RFExclude == 0) %>%
#   group_by(SubjectID, Sex, Group, ShapingInclude, UniqueSessionID, Phase, Phase_raw, Day, Day_collapsed, SessionTime, Box, Amplitude) %>%
#   dplyr::summarize(npass = n_distinct(Pass),
#                    mean_Theta = mean(Theta, na.rm=TRUE),
#                    mean_M50 = mean(M50, na.rm=TRUE),
#                    mean_Alpha = mean(Alpha, na.rm=TRUE))
# 
# RFadolescentsubject_LGM <- RRFICSS_LGM_findChange(RFadolescentsubject_LGM,"BL") # Add percent change from baseline


### Prepare Motor Behavior Variables -----
RFadolescentsubject_loco <- RFadolescentdata_raw %>% # Prep data for statistics - find subject average for each day
  filter(Pass > 1, ShapingInclude=="INCLUDE", !SubjectID %in% subjectexc) %>%
  group_by(SubjectID, Sex, Group, ShapingInclude, UniqueSessionID, Phase, Phase_raw, Day, Day_raw, Session, SessionTime, Box, Amplitude) %>%
  dplyr::summarize(npass = n_distinct(Pass),
                   session_sumLP = sum(sumLP, na.rm=TRUE),
                   pass_sumLP = mean(sumLP, na.rm=TRUE),
                   pass_maxLP = mean(maxLP, na.rm=TRUE),
                   session_Loco = sum(sumLoco, na.rm=TRUE),
                   pass_Loco = mean(sumLoco, na.rm=TRUE),
                   finalblock_Loco = mean(Loc_50, na.rm=TRUE))

RFadolescentsubject_loco <- RRFICSS_loco_findChange(RFadolescentsubject_loco,"BL") # Add percent change from baseline

### Add Percent Change from Baseline -----
#RFadolescentsubject_lm <- RRFICSS_lm_findChange(RFadolescentsubject_LP_lm,"BL")
#RFadolescentsubject_loco <- RRFICSS_loco_findChange(RFadolescentsubject_loco,"BL")
RFadolescentsubject_LP_lm <- RRFICSS_lm_findChange(RFadolescentsubject_LP_lm,"BL")
RFadolescentsubject_loco <- RRFICSS_loco_findChange(RFadolescentsubject_loco,"BL")

## Export processed data -----
write.csv(RFadolescentsubject_LP_lm, paste(path_analysis,"aRFadolescent_LMdata.csv",sep=""))
write.csv(RFadolescentsubject_loco, paste(path_analysis,"aRFadolescent_Locodata.csv",sep=""))


# PREPARE MEANS -----
## Session Means -----
RFadolmeans_LP_lm <- RFadolescentsubject_LP_lm %>%
  group_by(Group, Phase, UniqueSessionID, Day, Session, SessionTime) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   Int_mean = mean(mean_Int, na.rm=TRUE), Int_sd = sd(mean_Int, na.rm=TRUE), Int_se = Int_sd / sqrt(n),
                   IntC_mean = mean(mean_Int_C, na.rm=TRUE), IntC_sd = sd(mean_Int_C, na.rm=TRUE), IntC_se = IntC_sd / sqrt(n),
                   
                   Slope_mean = mean(mean_Slope, na.rm=TRUE), Slope_sd = sd(mean_Slope, na.rm=TRUE), Slope_se = Slope_sd / sqrt(n),
                   SlopeC_mean = mean(mean_Slope_C, na.rm=TRUE), SlopeC_sd = sd(mean_Slope_C, na.rm=TRUE), SlopeC_se = SlopeC_sd / sqrt(n))
RFadolmeans_LP_lm$sessionlabel = factor(RFadolmeans_LP_lm$UniqueSessionID, labels = c(1:36))
RFadolmeans_LP_lm


RFadolmeans_loco <- RFadolescentsubject_loco %>%
  group_by(Group, Phase, UniqueSessionID, Day, Session, SessionTime) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   
                   pass_maxLP_mean = mean(pass_maxLP, na.rm=TRUE), pass_maxLP_sd = sd(pass_maxLP,na.rm=TRUE), pass_maxLP_se = pass_maxLP_sd / sqrt(n),
                   pass_maxLPC_mean = mean(pass_maxLP_C, na.rm=TRUE), pass_maxLPC_sd = sd(pass_maxLP_C,na.rm=TRUE), pass_maxLPC_se = pass_maxLPC_sd / sqrt(n),
                   
                   session_sumLP_mean = mean(session_sumLP, na.rm=TRUE), session_sumLP_sd = sd(session_sumLP,na.rm=TRUE), session_sumLP_se = session_sumLP_sd / sqrt(n),
                   session_sumLPC_mean = mean(session_sumLP_C, na.rm=TRUE), session_sumLPC_sd = sd(session_sumLP_C,na.rm=TRUE), session_sumLPC_se = session_sumLPC_sd / sqrt(n),
                   
                   pass_sumLP_mean = mean(pass_sumLP, na.rm=TRUE), pass_sumLP_sd = sd(pass_sumLP,na.rm=TRUE), pass_sumLP_se = pass_sumLP_sd / sqrt(n),
                   pass_sumLPC_mean = mean(pass_sumLP_C, na.rm=TRUE), pass_sumLPC_sd = sd(pass_sumLP_C,na.rm=TRUE), pass_sumLPC_se = pass_sumLPC_sd / sqrt(n),
                   
                   finalblock_Loco_mean = mean(finalblock_Loco, na.rm=TRUE), finalblock_Loco_sd = sd(finalblock_Loco,na.rm=TRUE), finalblock_Loco_se = finalblock_Loco_sd / sqrt(n),
                   finalblock_LocoC_mean = mean(finalblock_Loco_C, na.rm=TRUE), finalblock_LocoC_sd = sd(finalblock_Loco_C,na.rm=TRUE), finalblock_LocoC_se = finalblock_LocoC_sd / sqrt(n))
RFadolmeans_loco$sessionlabel = factor(RFadolmeans_loco$UniqueSessionID, labels = c(1:36))
RFadolmeans_loco


## Overall Phase Means -----
RFadolmeans_Phase_LP_lm <- RFadolescentsubject_LP_lm %>%
  group_by(Group,Phase,SessionTime) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   Int_mean = mean(mean_Int, na.rm=TRUE), Int_sd = sd(mean_Int, na.rm=TRUE), Int_se = Int_sd / sqrt(n),
                   IntC_mean = mean(mean_Int_C, na.rm=TRUE), IntC_sd = sd(mean_Int_C, na.rm=TRUE), IntC_se = IntC_sd / sqrt(n),
                   IntPC_mean = mean(mean_Int_PC, na.rm=TRUE), IntPC_sd = sd(mean_Int_PC, na.rm=TRUE), IntPC_se = IntPC_sd / sqrt(n),
                   
                   
                   Slope_mean = mean(mean_Slope, na.rm=TRUE), Slope_sd = sd(mean_Slope, na.rm=TRUE), Slope_se = Slope_sd / sqrt(n),
                   SlopeC_mean = mean(mean_Slope_C, na.rm=TRUE), SlopeC_sd = sd(mean_Slope_C, na.rm=TRUE), SlopeC_se = SlopeC_sd / sqrt(n),
                   SlopePC_mean = mean(mean_Slope_PC, na.rm=TRUE), SlopePC_sd = sd(mean_Slope_PC, na.rm=TRUE), SlopePC_se = SlopePC_sd / sqrt(n))



RFadolmeans_Phase_loco <- RFadolescentsubject_loco %>%
  group_by(Group,Phase,SessionTime) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   
                   pass_maxLP_mean = mean(pass_maxLP, na.rm=TRUE), pass_maxLP_sd = sd(pass_maxLP,na.rm=TRUE), pass_maxLP_se = pass_maxLP_sd / sqrt(n),
                   pass_maxLPC_mean = mean(pass_maxLP_C, na.rm=TRUE), pass_maxLPC_sd = sd(pass_maxLP_C,na.rm=TRUE), pass_maxLPC_se = pass_maxLPC_sd / sqrt(n),
                   pass_maxLPPC_mean = mean(pass_maxLP_PC, na.rm=TRUE), pass_maxLPPC_sd = sd(pass_maxLP_PC,na.rm=TRUE), pass_maxLPPC_se = pass_maxLPPC_sd / sqrt(n),
                   
                   session_sumLP_mean = mean(session_sumLP, na.rm=TRUE), session_sumLP_sd = sd(session_sumLP,na.rm=TRUE), session_sumLP_se = session_sumLP_sd / sqrt(n),
                   session_sumLPC_mean = mean(session_sumLP_C, na.rm=TRUE), session_sumLPC_sd = sd(session_sumLP_C,na.rm=TRUE), session_sumLPC_se = session_sumLPC_sd / sqrt(n),
                   session_sumLPPC_mean = mean(session_sumLP_PC, na.rm=TRUE), session_sumLPPC_sd = sd(session_sumLP_PC,na.rm=TRUE), session_sumLPPC_se = session_sumLPPC_sd / sqrt(n),
                   
                   pass_sumLP_mean = mean(pass_sumLP, na.rm=TRUE), pass_sumLP_sd = sd(pass_sumLP,na.rm=TRUE), pass_sumLP_se = pass_sumLP_sd / sqrt(n),
                   pass_sumLPC_mean = mean(pass_sumLP_C, na.rm=TRUE), pass_sumLPC_sd = sd(pass_sumLP_C,na.rm=TRUE), pass_sumLPC_se = pass_sumLPC_sd / sqrt(n),
                   pass_sumLPPC_mean = mean(pass_sumLP_PC, na.rm=TRUE), pass_sumLPPC_sd = sd(pass_sumLP_PC,na.rm=TRUE), pass_sumLPPC_se = pass_sumLPPC_sd / sqrt(n),
                   
                   finalblock_Loco_mean = mean(finalblock_Loco, na.rm=TRUE), finalblock_Loco_sd = sd(finalblock_Loco,na.rm=TRUE), finalblock_Loco_se = finalblock_Loco_sd / sqrt(n),
                   finalblock_LocoC_mean = mean(finalblock_Loco_C, na.rm=TRUE), finalblock_LocoC_sd = sd(finalblock_Loco_C,na.rm=TRUE), finalblock_LocoC_se = finalblock_LocoC_sd / sqrt(n),
                   finalblock_LocoPC_mean = mean(finalblock_Loco_PC, na.rm=TRUE), finalblock_LocoPC_sd = sd(finalblock_Loco_PC,na.rm=TRUE), finalblock_LocoPC_se = finalblock_LocoPC_sd / sqrt(n))

view(RFadolmeans_Phase_LP_lm)
view(RFadolmeans_Phase_loco)


## Subject Phase Means -----
RFadolmeans_SubjectPhase_LP_lm <- RFadolescentsubject_LP_lm %>%
  group_by(SubjectID,Group,Sex,Phase,SessionTime) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   Int_mean = mean(mean_Int, na.rm=TRUE), Int_sd = sd(mean_Int, na.rm=TRUE), Int_se = Int_sd / sqrt(n),
                   IntC_mean = mean(mean_Int_C, na.rm=TRUE), IntC_sd = sd(mean_Int_C, na.rm=TRUE), IntC_se = IntC_sd / sqrt(n),
                   
                   Slope_mean = mean(mean_Slope, na.rm=TRUE), Slope_sd = sd(mean_Slope, na.rm=TRUE), Slope_se = Slope_sd / sqrt(n),
                   SlopeC_mean = mean(mean_Slope_C, na.rm=TRUE), SlopeC_sd = sd(mean_Slope_C, na.rm=TRUE), SlopeC_se = SlopeC_sd / sqrt(n))

RFadolmeans_SubjectPhase_loco <- RFadolescentsubject_loco %>%
  group_by(SubjectID,Group,Sex,Phase,SessionTime) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   
                   pass_maxLP_mean = mean(pass_maxLP, na.rm=TRUE), pass_maxLP_sd = sd(pass_maxLP,na.rm=TRUE), pass_maxLP_se = pass_maxLP_sd / sqrt(n),
                   pass_maxLPC_mean = mean(pass_maxLP_C, na.rm=TRUE), pass_maxLPC_sd = sd(pass_maxLP_C,na.rm=TRUE), pass_maxLPC_se = pass_maxLPC_sd / sqrt(n),
                   
                   session_sumLP_mean = mean(session_sumLP, na.rm=TRUE), session_sumLP_sd = sd(session_sumLP,na.rm=TRUE), session_sumLP_se = session_sumLP_sd / sqrt(n),
                   session_sumLPC_mean = mean(session_sumLP_C, na.rm=TRUE), session_sumLPC_sd = sd(session_sumLP_C,na.rm=TRUE), session_sumLPC_se = session_sumLPC_sd / sqrt(n),
                   
#                   pass_sumLP_mean = mean(pass_sumLP, na.rm=TRUE), pass_sumLP_sd = sd(pass_sumLP,na.rm=TRUE), pass_sumLP_se = pass_sumLP_sd / sqrt(n),
#                   pass_sumLPC_mean = mean(pass_sumLP_C, na.rm=TRUE), pass_sumLPC_sd = sd(pass_sumLP_C,na.rm=TRUE), pass_sumLPC_se = pass_sumLPC_sd / sqrt(n),
                   
                   finalblock_Loco_mean = mean(finalblock_Loco, na.rm=TRUE), finalblock_Loco_sd = sd(finalblock_Loco,na.rm=TRUE), finalblock_Loco_se = finalblock_Loco_sd / sqrt(n),
                   finalblock_LocoC_mean = mean(finalblock_Loco_C, na.rm=TRUE), finalblock_LocoC_sd = sd(finalblock_Loco_C,na.rm=TRUE), finalblock_LocoC_se = finalblock_LocoC_sd / sqrt(n))


## Export Means Tables -----
write.csv(RFadolmeans_Phase_LP_lm, paste(path_analysis_RRFtest_means,"RRFadolescent_Means_Phase_RRFvariables.csv",sep=""))
write.csv(RFadolmeans_Phase_loco, paste(path_analysis_RRFtest_means,"RRFadolescent_Means_Phase_Motorvariables.csv",sep=""))

write.csv(RFadolmeans_LP_lm, paste(path_analysis_RRFtest_means,"RRFadolescent_Means_PhaseandSession_RRFvariables.csv",sep=""))
write.csv(RFadolmeans_loco, paste(path_analysis_RRFtest_means,"RRFadolescent_Means_PhaseandSession_Motorvariables.csv",sep=""))


# STATISTICS -----
## Descriptive Statistics -----
### Subject ns -----
RFadol_Int_ns <- RFadolescentsubject_LP_lm %>%
  group_by(Sex, Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID))
RFadol_Int_ns

RFadol_excluded_ns <- RFadolescentdata_raw %>% filter(ShapingInclude!="INCLUDE") %>% group_by(Sex) %>% dplyr::summarize(n = n_distinct(SubjectID))
RFadol_excluded_ns # Excluded subject ns

amps <- RFadolescentsubject_LP_lm %>% group_by(SubjectID, Amplitude) %>% dplyr::summarize() # Included subject amplitudes

min(amps$Amplitude) # Overall group minimum amplitude
max(amps$Amplitude) # Overall group maximum amplitude

### Overall Phase Means -----
RFadolmeans_Phase_lm_overall <- RFadolescentsubject_LP_lm %>%
  group_by(Group,Phase) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   Int_mean = mean(mean_Int, na.rm=TRUE), Int_sd = sd(mean_Int, na.rm=TRUE), Int_se = Int_sd / sqrt(n),
                   IntC_mean = mean(mean_Int_C, na.rm=TRUE), IntC_sd = sd(mean_Int_C, na.rm=TRUE), IntC_se = IntC_sd / sqrt(n),
                   
                   Slope_mean = mean(mean_Slope, na.rm=TRUE), Slope_sd = sd(mean_Slope, na.rm=TRUE), Slope_se = Slope_sd / sqrt(n),
                   SlopeC_mean = mean(mean_Slope_C, na.rm=TRUE), SlopeC_sd = sd(mean_Slope_C, na.rm=TRUE), SlopeC_se = SlopeC_sd / sqrt(n))

RFadolmeans_Phase_loco_overall <- RFadolescentsubject_loco %>%
  group_by(Group,Phase) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   
                   pass_maxLP_mean = mean(pass_maxLP, na.rm=TRUE), pass_maxLP_sd = sd(pass_maxLP,na.rm=TRUE), pass_maxLP_se = pass_maxLP_sd / sqrt(n),
                   pass_maxLPC_mean = mean(pass_maxLP_C, na.rm=TRUE), pass_maxLPC_sd = sd(pass_maxLP_C,na.rm=TRUE), pass_maxLPC_se = pass_maxLPC_sd / sqrt(n),
                   
                   session_sumLP_mean = mean(session_sumLP, na.rm=TRUE), session_sumLP_sd = sd(session_sumLP,na.rm=TRUE), session_sumLP_se = session_sumLP_sd / sqrt(n),
                   session_sumLPC_mean = mean(session_sumLP_C, na.rm=TRUE), session_sumLPC_sd = sd(session_sumLP_C,na.rm=TRUE), session_sumLPC_se = session_sumLPC_sd / sqrt(n),
                   
                   pass_sumLP_mean = mean(pass_sumLP, na.rm=TRUE), pass_sumLP_sd = sd(pass_sumLP,na.rm=TRUE), pass_sumLP_se = pass_sumLP_sd / sqrt(n),
                   pass_sumLPC_mean = mean(pass_sumLP_C, na.rm=TRUE), pass_sumLPC_sd = sd(pass_sumLP_C,na.rm=TRUE), pass_sumLPC_se = pass_sumLPC_sd / sqrt(n),
                   
                   finalblock_Loco_mean = mean(finalblock_Loco, na.rm=TRUE), finalblock_Loco_sd = sd(finalblock_Loco,na.rm=TRUE), finalblock_Loco_se = finalblock_Loco_sd / sqrt(n),
                   finalblock_LocoC_mean = mean(finalblock_Loco_C, na.rm=TRUE), finalblock_LocoC_sd = sd(finalblock_Loco_C,na.rm=TRUE), finalblock_LocoC_se = finalblock_LocoC_sd / sqrt(n))

view(RFadolmeans_Phase_loco_overall)

RFadolmeans_Phase_loco_overall_sex <- RFadolescentsubject_loco %>%
  group_by(Group,Phase,Sex) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   
                   pass_maxLP_mean = mean(pass_maxLP, na.rm=TRUE), pass_maxLP_sd = sd(pass_maxLP,na.rm=TRUE), pass_maxLP_se = pass_maxLP_sd / sqrt(n),
                   pass_maxLPC_mean = mean(pass_maxLP_C, na.rm=TRUE), pass_maxLPC_sd = sd(pass_maxLP_C,na.rm=TRUE), pass_maxLPC_se = pass_maxLPC_sd / sqrt(n),
                   
                   session_sumLP_mean = mean(session_sumLP, na.rm=TRUE), session_sumLP_sd = sd(session_sumLP,na.rm=TRUE), session_sumLP_se = session_sumLP_sd / sqrt(n),
                   session_sumLPC_mean = mean(session_sumLP_C, na.rm=TRUE), session_sumLPC_sd = sd(session_sumLP_C,na.rm=TRUE), session_sumLPC_se = session_sumLPC_sd / sqrt(n),
                   
                   pass_sumLP_mean = mean(pass_sumLP, na.rm=TRUE), pass_sumLP_sd = sd(pass_sumLP,na.rm=TRUE), pass_sumLP_se = pass_sumLP_sd / sqrt(n),
                   pass_sumLPC_mean = mean(pass_sumLP_C, na.rm=TRUE), pass_sumLPC_sd = sd(pass_sumLP_C,na.rm=TRUE), pass_sumLPC_se = pass_sumLPC_sd / sqrt(n),
                   
                   finalblock_Loco_mean = mean(finalblock_Loco, na.rm=TRUE), finalblock_Loco_sd = sd(finalblock_Loco,na.rm=TRUE), finalblock_Loco_se = finalblock_Loco_sd / sqrt(n),
                   finalblock_LocoC_mean = mean(finalblock_Loco_C, na.rm=TRUE), finalblock_LocoC_sd = sd(finalblock_Loco_C,na.rm=TRUE), finalblock_LocoC_se = finalblock_LocoC_sd / sqrt(n))

view(RFadolmeans_Phase_loco_overall_sex)

## Intercept -----
### Phase Analysis by Group -----
#### Adolescent Saline
##### Mixed effects model by group, phase, session time, and sex with random intercept by subject
RFadol_Int_S_model <- glmmTMB(mean_Int ~ Phase*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Group=="AdolS"), family=stats::gaussian)
summary(RFadol_Int_S_model)

#### Compare full model to the null model
RFadol_Int_S_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_S_nullmodel, RFadol_Int_S_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadol_Int_S_emmeans <-emmeans(RFadol_Int_S_model, pairwise ~ Phase|SessionTime, adjust = 'none')
summary(RFadol_Int_S_emmeans)


#### Adolescent Morphine
##### Mixed effects model by group, phase, session time, and sex with random intercept by subject
RFadol_Int_M_model <- glmmTMB(mean_Int ~ Phase*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Group=="AdolM"), family=stats::gaussian)
summary(RFadol_Int_M_model)

#### Compare full model to the null model
RFadol_Int_M_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_M_nullmodel, RFadol_Int_M_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadol_Int_M_emmeans <-emmeans(RFadol_Int_M_model, pairwise ~ Phase|SessionTime, adjust = 'none')
summary(RFadol_Int_M_emmeans)


### Session Analysis by Group -----
#### Adolescent Saline
##### Baseline
##### Mixed effects model by day with random intercept by subject
RFadol_Int_S_BL_model <- glmmTMB(mean_Int ~ Day_raw*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="BL"&Group=="AdolS"),
                                     family = stats::gaussian())
summary(RFadol_Int_S_BL_model)

##### Compare full model to the null model
RFadol_Int_S_BL_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="BL"&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_S_BL_model, RFadol_Int_S_BL_nullmodel, test = "LRT") # Full model is a better fit than the null mode

###### Overall, there is no significant effect of Session in the baseline phase for the saline group.


#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_Int_S_INJ_model <- glmmTMB(mean_Int ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","INJ")&Group=="AdolS"),
                               family = stats::gaussian())
summary(RFadol_Int_S_INJ_model)

##### Compare full model to the null model
RFadol_Int_S_INJ_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","INJ")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_S_INJ_nullmodel, RFadol_Int_S_INJ_model, test = "LRT") # Full model is a better fit than the null mode

###### Overall, intercept decreases in the AM sessions in the injection phase for the saline group.


#### Post-Injection - Early Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_Int_S_POSTE_model <- glmmTMB(mean_Int ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","POST_EARLY")&Group=="AdolS"),
                                family = stats::gaussian())
summary(RFadol_Int_S_POSTE_model)

##### Compare full model to the null model
RFadol_Int_S_POSTE_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","POST_EARLY")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_S_POSTE_nullmodel, RFadol_Int_S_POSTE_model, test = "LRT") # Full model is not better fit than the null mode

###### Overall, no difference across sessions in post withdrawal early phase


#### Post-Injection - Late Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_Int_S_POSTL_model <- glmmTMB(mean_Int ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","POST_LATE")&Group=="AdolS"),
                                  family = stats::gaussian())
summary(RFadol_Int_S_POSTL_model)

##### Compare full model to the null model
RFadol_Int_S_POSTL_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","POST_LATE")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_S_POSTL_nullmodel, RFadol_Int_S_POSTL_model, test = "LRT") # Full model is a better fit than the null model
###### Overall, while some sessions are lower than baseline the null model is not significantly different than the saturated model.


#### Adolescent Morphine
##### Baseline
##### Mixed effects model by day with random intercept by subject
RFadol_Int_M_BL_model <- glmmTMB(mean_Int ~ Day_raw*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="BL"&Group=="AdolM"),
                                 family = stats::gaussian())
summary(RFadol_Int_M_BL_model)

##### Compare full model to the null model
RFadol_Int_M_BL_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="BL"&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_M_BL_model, RFadol_Int_M_BL_nullmodel, test = "LRT") # Full model is a better fit than the null mode

###### Overall, there is no significant effect of Session in the baseline phase for the morphine group.


#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_Int_M_INJ_model <- glmmTMB(mean_Int ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","INJ")&Group=="AdolM"),
                                  family = stats::gaussian())
summary(RFadol_Int_M_INJ_model)

##### Compare full model to the null model
RFadol_Int_M_INJ_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","INJ")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_M_INJ_nullmodel, RFadol_Int_M_INJ_model, test = "LRT") # Full model is a better fit than the null mode

###### Overall, session 13 is significantly higher than baseline for the morphine group.


#### Post-Injection - Early Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_Int_M_POSTE_model <- glmmTMB(mean_Int ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","POST_EARLY")&Group=="AdolM"),
                                    family = stats::gaussian())
summary(RFadol_Int_M_POSTE_model)

##### Compare full model to the null model
RFadol_Int_M_POSTE_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","POST_EARLY")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_M_POSTE_nullmodel, RFadol_Int_M_POSTE_model, test = "LRT") # Full model is not a better fit than the null mode


#### Post-Injection - Late Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_Int_M_POSTL_model <- glmmTMB(mean_Int ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","POST_LATE")&Group=="AdolM"),
                                    family = stats::gaussian())
summary(RFadol_Int_M_POSTL_model)

##### Compare full model to the null model
RFadol_Int_M_POSTL_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase %in% c("BL","POST_LATE")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_M_POSTL_nullmodel, RFadol_Int_M_POSTL_model, test = "LRT") # Full model is a better fit than the null mode

###### Overall, morning sessions are significantly lower than baseline for the morphine group.


### Between Group Analysis by Phase -----
#### Baseline
RFadol_Int_BL_model <- glmmTMB(mean_Int ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="BL"), family=stats::gaussian)
summary(RFadol_Int_BL_model)

##### Compare full model to the null model
RFadol_Int_BL_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_BL_model, RFadol_Int_BL_nullmodel, test = "LRT") # Full model is a better fit than the null mode

#### Injection
RFadol_Int_INJ_model <- glmmTMB(mean_Int ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="INJ"), family=stats::gaussian)
summary(RFadol_Int_INJ_model)

##### Compare full model to the null model
RFadol_Int_INJ_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="INJ"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_INJ_nullmodel, RFadol_Int_INJ_model, test = "LRT") # Full model is a better fit than the null mode

##### Post hoc tests with emmeans and pairs
RFadol_Int_INJ_emmeans <-emmeans(RFadol_Int_INJ_model, pairwise ~ Group|SessionTime, adjust = 'none')
summary(RFadol_Int_INJ_emmeans)


#### Post-Injection Early Phase
RFadol_Int_POSTE_model <- glmmTMB(mean_Int ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="POST_EARLY"), family=stats::gaussian)
summary(RFadol_Int_POSTE_model)

##### Compare full model to the null model
RFadol_Int_POSTE_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="POST_EARLY"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_POSTE_nullmodel, RFadol_Int_POSTE_model, test = "LRT") # Full model is a better fit than the null mode

#### Post-Injection Late Phase
RFadol_Int_POSTL_model <- glmmTMB(mean_Int ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="POST_LATE"), family=stats::gaussian)
summary(RFadol_Int_POSTL_model)

##### Compare full model to the null model
RFadol_Int_POSTL_nullmodel <- glmmTMB(mean_Int ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase=="POST_LATE"), family = stats::gaussian()) # Make the null model
anova(RFadol_Int_POSTL_nullmodel, RFadol_Int_POSTL_model, test = "LRT") # Full model is a better fit than the null mode

### Percent Change Group Analysis -----
#### Mixed effects model by phase and sex with random intercept by subject
RFadolescent_Int_PC_model <- glmmTMB(mean_Int_PC ~ Phase*Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase!="BL"), family=stats::gaussian)
summary(RFadolescent_Int_PC_model)

#### Compare full model to the null model
RFadolescent_Int_PC_nullmodel <- glmmTMB(mean_Int_PC ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_LP_lm, Phase!="BL"), family = stats::gaussian()) # Make the null model
anova(RFadolescent_Int_PC_nullmodel, RFadolescent_Int_PC_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadolescent_Int_PC_within_emmeans <-emmeans(RFadolescent_Int_PC_model, pairwise ~ Phase|Group|SessionTime, adjust = 'none')
summary(RFadolescent_Int_PC_within_emmeans)

RFadolescent_Int_PC_between_emmeans <-emmeans(RFadolescent_Int_PC_model, pairwise ~ Group|Phase|SessionTime, adjust = 'none')
summary(RFadolescent_Int_PC_between_emmeans)


### Export Analysis Tables -----
#### Phase Tables
RFadol_Int_M_model_paramtable <- mixedmodel_paramtable(RFadol_Int_M_model)
RFadol_Int_M_model_fittable <- mixedmodel_fittable(RFadol_Int_M_model, RFadol_Int_M_nullmodel)
RFadol_Int_M_model_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTime(RFadol_Int_M_emmeans,RFadol_Int_M_model)

write.csv(RFadol_Int_M_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Phase_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_M_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Phase_GroupM_model_fittable.csv",sep=""))
write.csv(RFadol_Int_M_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Phase_GroupM_emmeans_contraststable.csv",sep=""))

RFadol_Int_S_model_paramtable <- mixedmodel_paramtable(RFadol_Int_S_model)
RFadol_Int_S_model_fittable <- mixedmodel_fittable(RFadol_Int_S_model, RFadol_Int_S_nullmodel)
RFadol_Int_S_model_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTime(RFadol_Int_S_emmeans,RFadol_Int_S_model)

write.csv(RFadol_Int_S_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Phase_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_S_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Phase_GroupS_model_fittable.csv",sep=""))
write.csv(RFadol_Int_S_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Phase_GroupS_emmeans_contraststable.csv",sep=""))

#### Phase percent change
RFadolescent_Int_PC_model_paramtable <- mixedmodel_paramtable(RFadolescent_Int_PC_model)
RFadolescent_Int_PC_model_fittable <- mixedmodel_fittable(RFadolescent_Int_PC_model, RFadolescent_Int_PC_nullmodel)
RFadolescent_Int_PC_between_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTimePhase(RFadolescent_Int_PC_between_emmeans,RFadolescent_Int_PC_model)

write.csv(RFadolescent_Int_PC_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_PhasePC_model_paramtable.csv",sep=""))
write.csv(RFadolescent_Int_PC_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_PhasePC_model_fittable.csv",sep=""))
write.csv(RFadolescent_Int_PC_between_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_Int_PhasePC_emmeans_contraststable.csv",sep=""))


#### Session Tables
##### Baseline
RFadol_Int_M_BL_model_paramtable <- mixedmodel_paramtable(RFadol_Int_M_BL_model)
RFadol_Int_M_BL_model_fittable <- mixedmodel_fittable(RFadol_Int_M_BL_model, RFadol_Int_M_BL_nullmodel)

write.csv(RFadol_Int_M_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_BL_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_M_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_BL_GroupM_model_fittable.csv",sep=""))

RFadol_Int_S_BL_model_paramtable <- mixedmodel_paramtable(RFadol_Int_S_BL_model)
RFadol_Int_S_BL_model_fittable <- mixedmodel_fittable(RFadol_Int_S_BL_model, RFadol_Int_S_BL_nullmodel)

write.csv(RFadol_Int_S_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_BL_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_S_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_BL_GroupS_model_fittable.csv",sep=""))

##### Injection
RFadol_Int_M_INJ_model_paramtable <- mixedmodel_paramtable(RFadol_Int_M_INJ_model)
RFadol_Int_M_INJ_model_fittable <- mixedmodel_fittable(RFadol_Int_M_INJ_model, RFadol_Int_M_INJ_nullmodel)

write.csv(RFadol_Int_M_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_INJ_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_M_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_INJ_GroupM_model_fittable.csv",sep=""))

RFadol_Int_S_INJ_model_paramtable <- mixedmodel_paramtable(RFadol_Int_S_INJ_model)
RFadol_Int_S_INJ_model_fittable <- mixedmodel_fittable(RFadol_Int_S_INJ_model, RFadol_Int_S_INJ_nullmodel)

write.csv(RFadol_Int_S_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_INJ_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_S_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_INJ_GroupS_model_fittable.csv",sep=""))

##### Post Injection - Early
RFadol_Int_M_POSTE_model_paramtable <- mixedmodel_paramtable(RFadol_Int_M_POSTE_model)
RFadol_Int_M_POSTE_model_fittable <- mixedmodel_fittable(RFadol_Int_M_POSTE_model, RFadol_Int_M_POSTE_nullmodel)

write.csv(RFadol_Int_M_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_POSTE_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_M_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_POSTE_GroupM_model_fittable.csv",sep=""))

RFadol_Int_S_POSTE_model_paramtable <- mixedmodel_paramtable(RFadol_Int_S_POSTE_model)
RFadol_Int_S_POSTE_model_fittable <- mixedmodel_fittable(RFadol_Int_S_POSTE_model, RFadol_Int_S_POSTE_nullmodel)

write.csv(RFadol_Int_S_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_POSTE_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_S_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_POSTE_GroupS_model_fittable.csv",sep=""))

##### Post Injection - Late
RFadol_Int_M_POSTL_model_paramtable <- mixedmodel_paramtable(RFadol_Int_M_POSTL_model)
RFadol_Int_M_POSTL_model_fittable <- mixedmodel_fittable(RFadol_Int_M_POSTL_model, RFadol_Int_M_POSTL_nullmodel)

write.csv(RFadol_Int_M_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_POSTL_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_M_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_POSTL_GroupM_model_fittable.csv",sep=""))

RFadol_Int_S_POSTL_model_paramtable <- mixedmodel_paramtable(RFadol_Int_S_POSTL_model)
RFadol_Int_S_POSTL_model_fittable <- mixedmodel_fittable(RFadol_Int_S_POSTL_model, RFadol_Int_S_POSTL_nullmodel)

write.csv(RFadol_Int_S_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_POSTL_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_S_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_Session_POSTL_GroupS_model_fittable.csv",sep=""))

#### Between Group Tables
##### Baseline
RFadol_Int_BL_model_paramtable <- mixedmodel_paramtable(RFadol_Int_BL_model)
RFadol_Int_BL_model_fittable <- mixedmodel_fittable(RFadol_Int_BL_model, RFadol_Int_BL_nullmodel)

write.csv(RFadol_Int_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_BetweenGroups_BL_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_BetweenGroups_BL_model_fittable.csv",sep=""))

##### Injection
RFadol_Int_INJ_model_paramtable <- mixedmodel_paramtable(RFadol_Int_INJ_model)
RFadol_Int_INJ_model_fittable <- mixedmodel_fittable(RFadol_Int_INJ_model, RFadol_Int_INJ_nullmodel)
RFadol_Int_INJ_model_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTime(RFadol_Int_INJ_emmeans,RFadol_Int_INJ_model)

write.csv(RFadol_Int_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_BetweenGroups_INJ_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_BetweenGroups_INJ_model_fittable.csv",sep=""))
write.csv(RFadol_Int_INJ_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_Int_BetweenGroups_INJ_emmeans_contraststable.csv",sep=""))

##### Post Injection - Early
RFadol_Int_POSTE_model_paramtable <- mixedmodel_paramtable(RFadol_Int_POSTE_model)
RFadol_Int_POSTE_model_fittable <- mixedmodel_fittable(RFadol_Int_POSTE_model, RFadol_Int_POSTE_nullmodel)

write.csv(RFadol_Int_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_BetweenGroups_POSTE_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_BetweenGroups_POSTE_model_fittable.csv",sep=""))

##### Post Injection - Late
RFadol_Int_POSTL_model_paramtable <- mixedmodel_paramtable(RFadol_Int_POSTL_model)
RFadol_Int_POSTL_model_fittable <- mixedmodel_fittable(RFadol_Int_POSTL_model, RFadol_Int_POSTL_nullmodel)

write.csv(RFadol_Int_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_Int_BetweenGroups_POSTL_model_paramtable.csv",sep=""))
write.csv(RFadol_Int_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_Int_BetweenGroups_POSTL_model_fittable.csv",sep=""))




## Block Locomotion -----
### Phase Analysis by Group -----
#### Adolescent Saline
##### Mixed effects model by group, phase, session time, and sex with random intercept by subject
RFadol_blockloco_S_model <- glmmTMB(finalblock_Loco ~ Phase*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolS"), family=stats::gaussian)
summary(RFadol_blockloco_S_model)

#### Compare full model to the null model
RFadol_blockloco_S_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_S_nullmodel, RFadol_blockloco_S_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadol_blockloco_S_emmeans <-emmeans(RFadol_blockloco_S_model, pairwise ~ Phase|SessionTime, adjust = 'none')
summary(RFadol_blockloco_S_emmeans)


#### Adolescent Morphine
##### Mixed effects model by group, phase, session time, and sex with random intercept by subject
RFadol_blockloco_M_model <- glmmTMB(finalblock_Loco ~ Phase*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolM"), family=stats::gaussian)
summary(RFadol_blockloco_M_model)

#### Compare full model to the null model
RFadol_blockloco_M_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_M_nullmodel, RFadol_blockloco_M_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadol_blockloco_M_emmeans <-emmeans(RFadol_blockloco_M_model, pairwise ~ Phase|SessionTime, adjust = 'none')
summary(RFadol_blockloco_M_emmeans)



### Session Analysis by Group -----
#### Adolescent Saline
##### Baseline
##### Mixed effects model by day with random intercept by subject
RFadol_blockloco_S_BL_model <- glmmTMB(finalblock_Loco ~ Day_raw*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolS"),
                                   family = stats::gaussian())
summary(RFadol_blockloco_S_BL_model)

##### Compare full model to the null model
RFadol_blockloco_S_BL_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_S_BL_model, RFadol_blockloco_S_BL_nullmodel, test = "LRT") # Full model is a better fit than the null mode
###### Overall, there is no significant effect of Session in the baseline phase for the saline group.


#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_blockloco_S_INJ_model <- glmmTMB(finalblock_Loco ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolS"),
                                    family = stats::gaussian())
summary(RFadol_blockloco_S_INJ_model)

##### Compare full model to the null model
RFadol_blockloco_S_INJ_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_S_INJ_nullmodel, RFadol_blockloco_S_INJ_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, blockloco doesn't change across the injection phase for the saline group.


#### Post-Injection - Early Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_blockloco_S_POSTE_model <- glmmTMB(finalblock_Loco ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolS"),
                                      family = stats::gaussian())
summary(RFadol_blockloco_S_POSTE_model)

##### Compare full model to the null model
RFadol_blockloco_S_POSTE_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_S_POSTE_nullmodel, RFadol_blockloco_S_POSTE_model, test = "LRT") # Full model is a better fit than the null mode

###### Overall, blockloco is doesn't change for the saline group.


#### Post-Injection - Late Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_blockloco_S_POSTL_model <- glmmTMB(finalblock_Loco ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolS"),
                                      family = stats::gaussian())
summary(RFadol_blockloco_S_POSTL_model)

##### Compare full model to the null model
RFadol_blockloco_S_POSTL_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_S_POSTL_nullmodel, RFadol_blockloco_S_POSTL_model, test = "LRT") # Full model is a better fit than the null mode

###### Overall, blockloco is decreases at the end of the post late phase for the saline group.


#### Adolescent Morphine
##### Baseline
##### Mixed effects model by day with random intercept by subject
RFadol_blockloco_M_BL_model <- glmmTMB(finalblock_Loco ~ Day_raw*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolM"),
                                       family = stats::gaussian())
summary(RFadol_blockloco_M_BL_model)

##### Compare full model to the null model
RFadol_blockloco_M_BL_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_M_BL_model, RFadol_blockloco_M_BL_nullmodel, test = "LRT") # Full model is a better fit than the null mode
###### Overall, there is no significant effect of Session but there is a significant effect of session time in the baseline phase for the morphine group.


#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_blockloco_M_INJ_model <- glmmTMB(finalblock_Loco ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolM"),
                                        family = stats::gaussian())
summary(RFadol_blockloco_M_INJ_model)

##### Compare full model to the null model
RFadol_blockloco_M_INJ_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_M_INJ_nullmodel, RFadol_blockloco_M_INJ_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, there are all AM sessions have an increase in locomotion


#### Post-Injection - Early Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_blockloco_M_POSTE_model <- glmmTMB(finalblock_Loco ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolM"),
                                          family = stats::gaussian())
summary(RFadol_blockloco_M_POSTE_model)

##### Compare full model to the null model
RFadol_blockloco_M_POSTE_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_M_POSTE_nullmodel, RFadol_blockloco_M_POSTE_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, sessions  16 and 18 are significantly lower than baseline for the morphine group.


#### Post-Injection - Late Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_blockloco_M_POSTL_model <- glmmTMB(finalblock_Loco ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolM"),
                                          family = stats::gaussian())
summary(RFadol_blockloco_M_POSTL_model)

##### Compare full model to the null model
RFadol_blockloco_M_POSTL_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_M_POSTL_nullmodel, RFadol_blockloco_M_POSTL_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, all late phase sessions are significantly lower than baseline for the morphine group.



### Between Group Analysis by Phase -----
#### Baseline
RFadol_blockloco_BL_model <- glmmTMB(finalblock_Loco ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"), family=stats::gaussian)
summary(RFadol_blockloco_BL_model)

##### Compare full model to the null model
RFadol_blockloco_BL_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_BL_model, RFadol_blockloco_BL_nullmodel, test = "LRT") # Full model is not a better fit than the null model


#### Injection
RFadol_blockloco_INJ_model <- glmmTMB(finalblock_Loco ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="INJ"), family=stats::gaussian)
summary(RFadol_blockloco_INJ_model)

##### Compare full model to the null model
RFadol_blockloco_INJ_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="INJ"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_INJ_nullmodel, RFadol_blockloco_INJ_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadol_blockloco_INJ_emmeans <-emmeans(RFadol_blockloco_INJ_model, pairwise ~ Group|SessionTime, adjust = 'none')
summary(RFadol_blockloco_INJ_emmeans)


#### Post-Injection Early Phase
RFadol_blockloco_POSTE_model <- glmmTMB(finalblock_Loco ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_EARLY"), family=stats::gaussian)
summary(RFadol_blockloco_POSTE_model)

##### Compare full model to the null model
RFadol_blockloco_POSTE_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_EARLY"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_POSTE_nullmodel, RFadol_blockloco_POSTE_model, test = "LRT") # Full model is a better fit than the null model

#### Post-Injection Late Phase
RFadol_blockloco_POSTL_model <- glmmTMB(finalblock_Loco ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_LATE"), family=stats::gaussian)
summary(RFadol_blockloco_POSTL_model)

##### Compare full model to the null model
RFadol_blockloco_POSTL_nullmodel <- glmmTMB(finalblock_Loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_LATE"), family = stats::gaussian()) # Make the null model
anova(RFadol_blockloco_POSTL_nullmodel, RFadol_blockloco_POSTL_model, test = "LRT") # Full model is a better fit than the null model


### Percent Change Group Analysis -----
#### Mixed effects model by phase and sex with random finalblock_Locoercept by subject
RFadolescent_finalblock_Loco_PC_model <- glmmTMB(finalblock_Loco_PC ~ Phase*Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase!="BL"), family=stats::gaussian)
summary(RFadolescent_finalblock_Loco_PC_model)

#### Compare full model to the null model
RFadolescent_finalblock_Loco_PC_nullmodel <- glmmTMB(finalblock_Loco_PC ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase!="BL"), family = stats::gaussian()) # Make the null model
anova(RFadolescent_finalblock_Loco_PC_nullmodel, RFadolescent_finalblock_Loco_PC_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadolescent_finalblock_Loco_PC_within_emmeans <-emmeans(RFadolescent_finalblock_Loco_PC_model, pairwise ~ Phase|Group|SessionTime, adjust = 'none')
summary(RFadolescent_finalblock_Loco_PC_within_emmeans)

RFadolescent_finalblock_Loco_PC_between_emmeans <-emmeans(RFadolescent_finalblock_Loco_PC_model, pairwise ~ Group|Phase|SessionTime, adjust = 'none')
summary(RFadolescent_finalblock_Loco_PC_between_emmeans)



### Export Analysis Tables -----
#### Phase Tables
RFadol_blockloco_M_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_M_model)
RFadol_blockloco_M_model_fittable <- mixedmodel_fittable(RFadol_blockloco_M_model, RFadol_blockloco_M_nullmodel)
RFadol_blockloco_M_model_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTime(RFadol_blockloco_M_emmeans,RFadol_blockloco_M_model)

write.csv(RFadol_blockloco_M_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Phase_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_M_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Phase_GroupM_model_fittable.csv",sep=""))
write.csv(RFadol_blockloco_M_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Phase_GroupM_emmeans_contraststable.csv",sep=""))

RFadol_blockloco_S_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_S_model)
RFadol_blockloco_S_model_fittable <- mixedmodel_fittable(RFadol_blockloco_S_model, RFadol_blockloco_S_nullmodel)
RFadol_blockloco_S_model_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTime(RFadol_blockloco_S_emmeans,RFadol_blockloco_S_model)

write.csv(RFadol_blockloco_S_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Phase_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_S_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Phase_GroupS_model_fittable.csv",sep=""))
write.csv(RFadol_blockloco_S_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Phase_GroupS_emmeans_contraststable.csv",sep=""))

#### Phase percent change
RFadolescent_finalblock_Loco_PC_model_paramtable <- mixedmodel_paramtable(RFadolescent_finalblock_Loco_PC_model)
RFadolescent_finalblock_Loco_PC_model_fittable <- mixedmodel_fittable(RFadolescent_finalblock_Loco_PC_model, RFadolescent_finalblock_Loco_PC_nullmodel)
RFadolescent_finalblock_Loco_PC_between_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTimePhase(RFadolescent_finalblock_Loco_PC_between_emmeans,RFadolescent_finalblock_Loco_PC_model)

write.csv(RFadolescent_finalblock_Loco_PC_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_PhasePC_model_paramtable.csv",sep=""))
write.csv(RFadolescent_finalblock_Loco_PC_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_PhasePC_model_fittable.csv",sep=""))
write.csv(RFadolescent_finalblock_Loco_PC_between_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_PhasePC_emmeans_contraststable.csv",sep=""))


#### Session Tables
##### Baseline
RFadol_blockloco_M_BL_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_M_BL_model)
RFadol_blockloco_M_BL_model_fittable <- mixedmodel_fittable(RFadol_blockloco_M_BL_model, RFadol_blockloco_M_BL_nullmodel)

write.csv(RFadol_blockloco_M_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_BL_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_M_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_BL_GroupM_model_fittable.csv",sep=""))

RFadol_blockloco_S_BL_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_S_BL_model)
RFadol_blockloco_S_BL_model_fittable <- mixedmodel_fittable(RFadol_blockloco_S_BL_model, RFadol_blockloco_S_BL_nullmodel)

write.csv(RFadol_blockloco_S_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_BL_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_S_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_BL_GroupS_model_fittable.csv",sep=""))

##### Injection
RFadol_blockloco_M_INJ_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_M_INJ_model)
RFadol_blockloco_M_INJ_model_fittable <- mixedmodel_fittable(RFadol_blockloco_M_INJ_model, RFadol_blockloco_M_INJ_nullmodel)

write.csv(RFadol_blockloco_M_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_INJ_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_M_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_INJ_GroupM_model_fittable.csv",sep=""))

RFadol_blockloco_S_INJ_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_S_INJ_model)
RFadol_blockloco_S_INJ_model_fittable <- mixedmodel_fittable(RFadol_blockloco_S_INJ_model, RFadol_blockloco_S_INJ_nullmodel)

write.csv(RFadol_blockloco_S_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_INJ_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_S_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_INJ_GroupS_model_fittable.csv",sep=""))

##### Post Injection - Early
RFadol_blockloco_M_POSTE_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_M_POSTE_model)
RFadol_blockloco_M_POSTE_model_fittable <- mixedmodel_fittable(RFadol_blockloco_M_POSTE_model, RFadol_blockloco_M_POSTE_nullmodel)

write.csv(RFadol_blockloco_M_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_POSTE_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_M_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_POSTE_GroupM_model_fittable.csv",sep=""))

RFadol_blockloco_S_POSTE_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_S_POSTE_model)
RFadol_blockloco_S_POSTE_model_fittable <- mixedmodel_fittable(RFadol_blockloco_S_POSTE_model, RFadol_blockloco_S_POSTE_nullmodel)

write.csv(RFadol_blockloco_S_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_POSTE_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_S_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_POSTE_GroupS_model_fittable.csv",sep=""))

##### Post Injection - Late
RFadol_blockloco_M_POSTL_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_M_POSTL_model)
RFadol_blockloco_M_POSTL_model_fittable <- mixedmodel_fittable(RFadol_blockloco_M_POSTL_model, RFadol_blockloco_M_POSTL_nullmodel)

write.csv(RFadol_blockloco_M_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_POSTL_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_M_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_POSTL_GroupM_model_fittable.csv",sep=""))

RFadol_blockloco_S_POSTL_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_S_POSTL_model)
RFadol_blockloco_S_POSTL_model_fittable <- mixedmodel_fittable(RFadol_blockloco_S_POSTL_model, RFadol_blockloco_S_POSTL_nullmodel)

write.csv(RFadol_blockloco_S_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_POSTL_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_S_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_Session_POSTL_GroupS_model_fittable.csv",sep=""))

#### Between Group Tables
##### Baseline
RFadol_blockloco_BL_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_BL_model)
RFadol_blockloco_BL_model_fittable <- mixedmodel_fittable(RFadol_blockloco_BL_model, RFadol_blockloco_BL_nullmodel)

write.csv(RFadol_blockloco_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_BetweenGroups_BL_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_BetweenGroups_BL_model_fittable.csv",sep=""))

##### Injection
RFadol_blockloco_INJ_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_INJ_model)
RFadol_blockloco_INJ_model_fittable <- mixedmodel_fittable(RFadol_blockloco_INJ_model, RFadol_blockloco_INJ_nullmodel)
RFadol_blockloco_INJ_model_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTime(RFadol_blockloco_INJ_emmeans,RFadol_blockloco_INJ_model)

write.csv(RFadol_blockloco_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_BetweenGroups_INJ_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_BetweenGroups_INJ_model_fittable.csv",sep=""))
write.csv(RFadol_blockloco_INJ_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_BetweenGroups_INJ_emmeans_contraststable.csv",sep=""))

##### Post Injection - Early
RFadol_blockloco_POSTE_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_POSTE_model)
RFadol_blockloco_POSTE_model_fittable <- mixedmodel_fittable(RFadol_blockloco_POSTE_model, RFadol_blockloco_POSTE_nullmodel)

write.csv(RFadol_blockloco_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_BetweenGroups_POSTE_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_BetweenGroups_POSTE_model_fittable.csv",sep=""))

##### Post Injection - Late
RFadol_blockloco_POSTL_model_paramtable <- mixedmodel_paramtable(RFadol_blockloco_POSTL_model)
RFadol_blockloco_POSTL_model_fittable <- mixedmodel_fittable(RFadol_blockloco_POSTL_model, RFadol_blockloco_POSTL_nullmodel)

write.csv(RFadol_blockloco_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_BetweenGroups_POSTL_model_paramtable.csv",sep=""))
write.csv(RFadol_blockloco_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_blockloco_BetweenGroups_POSTL_model_fittable.csv",sep=""))



# FIGURES -----
## Example Subject Plots -----
# Prep data
examplesubjectdata <- RFadolescentdata_raw %>% dplyr::filter(SubjectID == 536) %>% 
  pivot_longer(cols = c('Freq_141','Freq_100','Freq_89','Freq_79','Freq_71','Freq_63','Freq_50'), names_to = 'Frequency', values_to='LeverPresses')

examplesubjectdata$Pass <- factor(examplesubjectdata$Pass)

# Prep plotting variables
xorder <- c('Freq_50', 'Freq_63', 'Freq_71', 'Freq_79', 'Freq_89', 'Freq_100', 'Freq_141')
xlabels <- c('50','63','71','79','89','100','141')


AICSS_ymax <- 55
AICSS_ymin <- 0

# Plot Reduced RF Example - Baseline
MBLplot <- ggplot(subset(examplesubjectdata, UniqueSessionID == "BL_1_AM"), aes(x = Frequency, y = LeverPresses, group = Pass)) +
  geom_line(linewidth = mylinewidth, aes(color = Pass), show.legend = FALSE) + scale_x_discrete(limits = xorder, labels = xlabels) + 
  geom_point(aes(x = Frequency, y = LeverPresses, color = Pass), size = mypointsize, show.legend = FALSE) +
  scale_color_manual(values=c("#8d8d8d", "#0000FF", "#950EF2", "#FF0080", "#0FC127")) +  # Old values: values=c("#DC8B06", "#0E6BF2", "#950EF2", "#0FC127", "#D4212F")
  scale_y_continuous(limits = c(AICSS_ymin, AICSS_ymax), breaks=seq(0,50,10), expand=c(0,0)) +
  mytheme +
  ylab("Lever Presses") +
  ggtitle("Baseline (PN43)")

MBLplot

# Plot Reduced RF Example - Morphine Administration
MAplot <- ggplot(subset(examplesubjectdata, UniqueSessionID == "INJ_6_AM"), aes(x = Frequency, y = LeverPresses, group = Pass)) +
  geom_line(linewidth = mylinewidth, aes(color = Pass), show.legend = FALSE) + scale_x_discrete(limits = xorder, labels = xlabels) + 
  geom_point(aes(x = Frequency, y = LeverPresses, color = Pass), size = mypointsize, show.legend = FALSE) +
  scale_color_manual(values=c("#8d8d8d", "#0000FF", "#950EF2", "#FF0080", "#0FC127")) +  # Old values: values=c("#DC8B06", "#0E6BF2", "#950EF2", "#0FC127", "#D4212F")
  scale_y_continuous(limits = c(AICSS_ymin, AICSS_ymax), breaks=seq(0,50,10), expand=c(0,0)) +
  mytheme +
  ylab("Lever Presses") +
  ggtitle("Morphine (PN50)")

MAplot

# Plot Reduced RF Example - Morphine Withdrawal
MWplot <- ggplot(subset(examplesubjectdata, UniqueSessionID == "POST_2_AM"), aes(x = Frequency, y = LeverPresses, group = Pass)) +
  geom_line(linewidth = mylinewidth, aes(color = Pass), show.legend = FALSE) + scale_x_discrete(limits = xorder, labels = xlabels) + 
  geom_point(aes(x = Frequency, y = LeverPresses, color = Pass), size = mypointsize, show.legend = FALSE) +
  scale_color_manual(values=c("#8d8d8d", "#0000FF", "#950EF2", "#FF0080", "#0FC127")) +  # Old values: values=c("#DC8B06", "#0E6BF2", "#950EF2", "#0FC127", "#D4212F")
  scale_y_continuous(limits = c(AICSS_ymin, AICSS_ymax), breaks=seq(0,50,10), expand=c(0,0)) +
  mytheme +
  ylab("Lever Presses") +
  ggtitle("Withdrawal (PN53)")

MWplot

# Make empty tile for key
keyplot <- ggplot() + mytheme

# Create full panel
exampleRRFplots <- ggarrange(MBLplot,MAplot, MWplot,keyplot, ncol=4, nrow=1, widths=c(1,1,1,.8), labels = c("A", "B", "C", ""))
exampleRRFplots

## Figure Panels -----
### Prepare plot colors and variables -----
# Morphine admin color "AdolM" = "#FEDB16"
groupcolors <- c("AdolS" = AdolScolor, "AdolM" = AdolMcolor)
Mphasecolors <- c("BL" = ICSSBLcolor, "INJ"=AdolMcolor,"POST_EARLY"=AdolEPcolor, "POST_LATE"=AdolLPcolor)
Sphasecolors <- c("BL" = ICSSBLcolor, "INJ"=AdolScolor,"POST_EARLY"=ICSSBLcolor,"POST_LATE"=ICSSBLcolor)

sexshapes <- c("M" = 0, "F" = 1)

# Percent change colors
PCphasecolors <- c("AdolSINJ"=AdolScolor,"AdolMINJ"=AdolMcolor,"AdolSPOST_EARLY"=ICSSBLcolor,"AdolMPOST_EARLY"=AdolEPcolor,"AdolSPOST_LATE"=ICSSBLcolor,"AdolMPOST_LATE"=AdolLPcolor)
RFadolmeans_Phase_LP_lm$PCfill <- paste(RFadolmeans_Phase_LP_lm$Group,RFadolmeans_Phase_LP_lm$Phase,sep="")
RFadolmeans_Phase_loco$PCfill <- paste(RFadolmeans_Phase_loco$Group,RFadolmeans_Phase_loco$Phase,sep="")


# Set up x axis ticks and labels - label day (every other session)
sessionticks <- c(5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35)
sessionlabels <- c('45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60')
#sessionlabels <- c('1','2','3','4','5','6','7','1','2','3','4','5','6','7','8','9')

injsession <- c(1,3,5,7,9,11,13)

### Intercept -----
ymax_Int <- 20
ymin_Int <- -10

ymax_Int_bar <- 30
ymin_Int_bar <- 0

injmarker_Int <- c(rep(ymin_Int,7))
injticks_Int <- data.frame(injsession,injmarker_Int)

# Intercept Change Line
IntC_LP_line <- ggplot(subset(RFadolmeans_LP_lm,Phase!="BL"), aes(x=sessionlabel, y=IntC_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadolmeans_LP_lm,Phase!="BL"), aes(ymin=IntC_mean-IntC_se, ymax=IntC_mean+IntC_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name ="Postnatal Day", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0, 0), limits = c(ymin_Int, ymax_Int), breaks = seq(-10, ymax_Int, 10)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("Presses at Median Frequency", sep = "")) +
  geom_point(data=injticks_Int, aes(x=injsession, y=injmarker_Int, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE) + # Add injection ticks
  annotate("text", x = 13, y = 17, label = "*", size = mystarsize) + # Add asterisks for morphine sessions 
  annotate("text", x = 14, y = -4.25, label = "*", size = mystarsize)+
  annotate("text", x = 21, y = -7.5, label = "*", size = mystarsize)+
  annotate("text", x = 23, y = -6, label = "*", size = mystarsize)+
  annotate("text", x = 25, y = -8, label = "*", size = mystarsize)+
  annotate("text", x = 26, y = -5, label = "*", size = mystarsize)+
  annotate("text", x = 27, y = -7.5, label = "*", size = mystarsize)+
  annotate("text", x = 29, y = -8.75, label = "*", size = mystarsize)+
  annotate("text", x = 31, y = -7, label = "*", size = mystarsize) +
  annotate("text", x = 7, y = -4.75, label = "#", size = 2)+ # Add hashtags for saline sessions
  annotate("text", x = 9, y = -4.25, label = "#", size = 2)+
  annotate("text", x = 11, y = -4.75, label = "#", size = 2)+
  annotate("text", x = 13, y = -3.75, label = "#", size = 2)

IntC_LP_line

# Saline Group Intercept by Phase Bar
Int_LP_S_bar <- ggplot(subset(RFadolmeans_Phase_LP_lm,Group=="AdolS" & SessionTime == "AM"), aes(x=Phase, y=Int_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_LP_lm, Group=="AdolS" & SessionTime == "AM"), mapping = aes(x = Phase, y = Int_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_LP_lm,Group=="AdolS" & SessionTime == "AM"), aes(ymin=Int_mean-Int_se, ymax=Int_mean+Int_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0, 0), limits = c(ymin_Int_bar, ymax_Int_bar), breaks = seq(ymin_Int_bar, ymax_Int_bar, 10)) + 
  scale_color_manual("legend", values = Sphasecolors) +
  scale_fill_manual("legend", values = Sphasecolors) +
  ylab("Presses at Median Frequency") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Control", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(14) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2) +
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(17.5) ,label.size = 2, label = "##",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(21) ,label.size = 2, label = "##",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

Int_LP_S_bar

# Morphine Group Intercept by Phase Bar
Int_LP_M_bar <- ggplot(subset(RFadolmeans_Phase_LP_lm,Group=="AdolM" & SessionTime == "AM"), aes(x=Phase, y=Int_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_LP_lm, Group=="AdolM" & SessionTime == "AM"), mapping = aes(x = Phase, y = Int_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_LP_lm,Group=="AdolM" & SessionTime == "AM"), aes(ymin=Int_mean-Int_se, ymax=Int_mean+Int_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0, 0), limits = c(ymin_Int_bar, ymax_Int_bar), breaks = seq(ymin_Int_bar, ymax_Int_bar, 10)) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  ylab("Presses at Median Frequency") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Morphine", sep = "")) +
  geom_hline(yintercept = 0, col = 'black', linewidth = 1)+
  geom_bracket(xmin=c(2), xmax = c(3), y.position=c(17) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2), xmax = c(4), y.position=c(20.5) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(24.5) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(28) ,label.size = mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

Int_LP_M_bar

# Create full panel
Int_LP <- ggarrange(Int_LP_S_bar, Int_LP_M_bar, IntC_LP_line, ncol=3, nrow=1, widths=c(1.6,1.6,4), labels = c("D", "E", "F"))
Int_LP


### Final Pass Locomotion -----
ymax_blockLocoC <- 32
ymin_blockLocoC <- -18

ymax_blockLoco_bar <- 82
ymin_blockLoco_bar <- 0

injmarker_blockLoco <- c(rep(ymin_blockLocoC,7))
injticks_blockLoco <- data.frame(injsession,injmarker_blockLoco)

# Final Pass Locomotion Change Line
blockLoco_LP_line <- ggplot(subset(RFadolmeans_loco,Phase!="BL"), aes(x=sessionlabel, y=finalblock_LocoC_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadolmeans_loco,Phase!="BL"), aes(ymin=finalblock_LocoC_mean-finalblock_LocoC_se, ymax=finalblock_LocoC_mean+finalblock_LocoC_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name ="Postnatal Day", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blockLocoC, ymax_blockLocoC), breaks = seq(-15, ymax_blockLocoC, 15)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("Last Block Locomotion", sep = "")) +
  geom_point(data=injticks_blockLoco, aes(x=injsession, y=injmarker_blockLoco, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE)+ # Add injection ticks
  annotate("text", x = 3, y = 23, label = "*", size = mystarsize) + # Add asterisks for morphine sessions 
  annotate("text", x = 4, y = 4, label = "*", size = mystarsize)+
  annotate("text", x = 5, y = 21, label = "*", size = mystarsize)+
  annotate("text", x = 7, y = 24, label = "*", size = mystarsize)+
  annotate("text", x = 8, y = 4, label = "*", size = mystarsize)+
  annotate("text", x = 9, y = 20, label = "*", size = mystarsize)+
  annotate("text", x = 10, y = 4, label = "*", size = mystarsize)+
  annotate("text", x = 11, y = 28, label = "*", size = mystarsize)+
  annotate("text", x = 12, y = 4, label = "*", size = mystarsize)+
  annotate("text", x = 13, y = 29, label = "*", size = mystarsize)+
  annotate("text", x = 14, y = 4, label = "*", size = mystarsize)

blockLoco_LP_line

# Saline Group Final Pass Locomotion by Phase Bar
blockLoco_LP_S_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "AM"), aes(x=Phase, y=finalblock_Loco_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolS" & SessionTime == "AM"), mapping = aes(x = Phase, y = finalblock_Loco_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "AM"), aes(ymin=finalblock_Loco_mean-finalblock_Loco_se, ymax=finalblock_Loco_mean+finalblock_Loco_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blockLoco_bar, ymax_blockLoco_bar), breaks = seq(ymin_blockLoco_bar, ymax_blockLoco_bar,25)) + 
  scale_color_manual("legend", values = Sphasecolors) +
  scale_fill_manual("legend", values = Sphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Block Locomotion") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Control", sep = ""))

blockLoco_LP_S_bar

# Morphine Group Final Pass Locomotion by Phase Bar
blockLoco_LP_M_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "AM"), aes(x=Phase, y=finalblock_Loco_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolM" & SessionTime == "AM"), mapping = aes(x = Phase, y = finalblock_Loco_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "AM"), aes(ymin=finalblock_Loco_mean-finalblock_Loco_se, ymax=finalblock_Loco_mean+finalblock_Loco_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blockLoco_bar, ymax_blockLoco_bar), breaks = seq(ymin_blockLoco_bar, ymax_blockLoco_bar, 25)) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Block Locomotion") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Morphine", sep = "")) +   
  geom_bracket(xmin=c(1), xmax = c(1.9), y.position=c(65) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.1), xmax = c(3), y.position=c(65) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.1), xmax = c(4), y.position=c(74) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

blockLoco_LP_M_bar

# Create full panel
blockLoco <- ggarrange(blockLoco_LP_S_bar,blockLoco_LP_M_bar, blockLoco_LP_line, ncol=3, nrow=1, widths=c(1.6,1.6,4), labels = c("G", "H", "I"))
blockLoco


### Intercept Percent Change -----
ymax_IntPC <- 120
ymin_IntPC <- -60

Int_PC_bar <- ggplot(subset(RFadolmeans_Phase_LP_lm,Phase!="BL" & SessionTime=="AM"), aes(x=Phase, y=IntPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadolmeans_Phase_LP_lm,Phase!="BL" & SessionTime=="AM"), aes(ymin=IntPC_mean-IntPC_se, ymax=IntPC_mean+IntPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  #geom_point(subset(PCmeans_lm_subject, Phase!="BL"), mapping = aes(x = Phase, y = IntPC_mean), position = position_dodge(width = 0.95), size = mysubjectpointsize, stroke = mypointstroke, shape = 1, color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_IntPC, ymax_IntPC), breaks = seq(-60, ymax_IntPC, 60)) + 
  scale_x_discrete(expand = c(0.25, 0), labels=c('INJ','POST EARLY','POST LATE')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change") + xlab('Phase') +
  mytheme +
  theme(axis.title.x = element_blank())+
  ggtitle(paste("Presses at Median Frequency", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)+
  geom_bracket(xmin=c(.8), xmax = c(1.2), y.position=c(109) ,label.size = 2, label = "&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

Int_PC_bar


### Final Block Locomotion Percent Change -----
ymax_blocklocoPC <- 120
ymin_blocklocoPC <- -60

blockloco_PC_bar <- ggplot(subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="AM"), aes(x=Phase, y=finalblock_LocoPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="AM"), aes(ymin=finalblock_LocoPC_mean-finalblock_LocoPC_se, ymax=finalblock_LocoPC_mean+finalblock_LocoPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  #geom_point(subset(PCmeans_lm_subject, Phase!="BL"), mapping = aes(x = Phase, y = blocklocoPC_mean), position = position_dodge(width = 0.95), size = mysubjectpointsize, stroke = mypointstroke, shape = 1, color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blocklocoPC, ymax_blocklocoPC), breaks = seq(-60, ymax_blocklocoPC, 60)) + 
  scale_x_discrete(expand = c(0.25, 0), labels=c('INJ','POST EARLY','POST LATE')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change") + xlab('Phase') +
  mytheme +
  theme(axis.title.x = element_blank())+
  ggtitle(paste("Last Block Locomotion", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)+
  geom_bracket(xmin=c(.8), xmax = c(1.2), y.position=c(113) ,label.size = 2, label = "&&&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(2.8), xmax = c(3.2), y.position=c(50) ,label.size = 2, label = "&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

blockloco_PC_bar

# Create full panel
IntLoco_PC <- ggarrange(Int_PC_bar, blockloco_PC_bar, ncol=2, nrow=1, widths=c(4,4), labels = c("J", "K"))
IntLoco_PC


## Figure 4 -----
# Create full panel
figure4 <- ggarrange(exampleRRFplots,Int_LP,blockLoco,IntLoco_PC, ncol=1, nrow=4)
figure4
ggsave(plot=figure4, path=path_figures, device="pdf", filename="F4_AllPanels.pdf", width=RFpanelwidth, height=RFpanelheight*4)




# SUPPLEMENTARY -----
## STATISTICS -----

## Max Lever Presses -----
### Phase Analysis by Group -----
#### Adolescent Saline
##### Mixed effects model by group, phase, session time, and sex with random intercept by subject
RFadol_maxLP_S_model <- glmmTMB(pass_maxLP ~ Phase*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolS"), family=stats::gaussian)
summary(RFadol_maxLP_S_model)

#### Compare full model to the null model
RFadol_maxLP_S_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_S_nullmodel, RFadol_maxLP_S_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadol_maxLP_S_emmeans <-emmeans(RFadol_maxLP_S_model, pairwise ~ Phase|SessionTime, adjust = 'none')
summary(RFadol_maxLP_S_emmeans)

#### Adolescent Morphine
##### Mixed effects model by group, phase, session time, and sex with random intercept by subject
RFadol_maxLP_M_model <- glmmTMB(pass_maxLP ~ Phase*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolM"), family=stats::gaussian)
summary(RFadol_maxLP_M_model)

#### Compare full model to the null model
RFadol_maxLP_M_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_M_nullmodel, RFadol_maxLP_M_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadol_maxLP_M_emmeans <-emmeans(RFadol_maxLP_M_model, pairwise ~ Phase|SessionTime, adjust = 'none')
summary(RFadol_maxLP_M_emmeans)


### Session Analysis by Group -----
#### Adolescent Saline
##### Baseline
##### Mixed effects model by day with random intercept by subject
RFadol_maxLP_S_BL_model <- glmmTMB(pass_maxLP ~ Day_raw*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolS"),
                                   family = stats::gaussian())
summary(RFadol_maxLP_S_BL_model)

##### Compare full model to the null model
RFadol_maxLP_S_BL_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_S_BL_model, RFadol_maxLP_S_BL_nullmodel, test = "LRT") # Full model is a better fit than the null mode

###### Overall, there is no significant effect of Session in the baseline phase for the saline group.


#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_maxLP_S_INJ_model <- glmmTMB(pass_maxLP ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolS"),
                                    family = stats::gaussian())
summary(RFadol_maxLP_S_INJ_model)

##### Compare full model to the null model
RFadol_maxLP_S_INJ_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_S_INJ_nullmodel, RFadol_maxLP_S_INJ_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, maxLP in the morning sessions decreases significantly across the injection phase for the saline group.


#### Post-Injection - Early Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_maxLP_S_POSTE_model <- glmmTMB(pass_maxLP ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolS"),
                                      family = stats::gaussian())
summary(RFadol_maxLP_S_POSTE_model)

##### Compare full model to the null model
RFadol_maxLP_S_POSTE_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_S_POSTE_nullmodel, RFadol_maxLP_S_POSTE_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, maxLP is lower than baseline for the saline group.


#### Post-Injection - Late Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_maxLP_S_POSTL_model <- glmmTMB(pass_maxLP ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolS"),
                                      family = stats::gaussian())
summary(RFadol_maxLP_S_POSTL_model)

##### Compare full model to the null model
RFadol_maxLP_S_POSTL_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_S_POSTL_nullmodel, RFadol_maxLP_S_POSTL_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, maxLP is lower than baseline in post late phase for the saline group.


#### Adolescent Morphine
##### Baseline
##### Mixed effects model by day with random intercept by subject
RFadol_maxLP_M_BL_model <- glmmTMB(pass_maxLP ~ Day_raw*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolM"),
                                   family = stats::gaussian())
summary(RFadol_maxLP_M_BL_model)

##### Compare full model to the null model
RFadol_maxLP_M_BL_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_M_BL_model, RFadol_maxLP_M_BL_nullmodel, test = "LRT") # Full model is a better fit than the null mode
###### Overall, there is no significant effect of Session in the baseline phase for the morphine group.


#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_maxLP_M_INJ_model <- glmmTMB(pass_maxLP ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolM"),
                                    family = stats::gaussian())
summary(RFadol_maxLP_M_INJ_model)

##### Compare full model to the null model
RFadol_maxLP_M_INJ_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_M_INJ_nullmodel, RFadol_maxLP_M_INJ_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, there is no significant effect of session for the morphine group.


#### Post-Injection - Early Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_maxLP_M_POSTE_model <- glmmTMB(pass_maxLP ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolM"),
                                      family = stats::gaussian())
summary(RFadol_maxLP_M_POSTE_model)

##### Compare full model to the null model
RFadol_maxLP_M_POSTE_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_M_POSTE_nullmodel, RFadol_maxLP_M_POSTE_model, test = "LRT") # Full model is close to a better fit than the null mode


#### Post-Injection - Late Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_maxLP_M_POSTL_model <- glmmTMB(pass_maxLP ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolM"),
                                      family = stats::gaussian())
summary(RFadol_maxLP_M_POSTL_model)

##### Compare full model to the null model
RFadol_maxLP_M_POSTL_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_M_POSTL_nullmodel, RFadol_maxLP_M_POSTL_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, morning late phase sessions have higher maxLP for morphine group


### Between Group Analysis by Phase -----
#### Baseline
RFadol_maxLP_BL_model <- glmmTMB(pass_maxLP ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"), family=stats::gaussian)
summary(RFadol_maxLP_BL_model)

##### Compare full model to the null model
RFadol_maxLP_BL_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_BL_model, RFadol_maxLP_BL_nullmodel, test = "LRT") # Full model is a better fit than the null mode

#### Injection
RFadol_maxLP_INJ_model <- glmmTMB(pass_maxLP ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="INJ"), family=stats::gaussian)
summary(RFadol_maxLP_INJ_model)

##### Compare full model to the null model
RFadol_maxLP_INJ_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="INJ"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_INJ_nullmodel, RFadol_maxLP_INJ_model, test = "LRT") # Full model is a better fit than the null mode


#### Post-Injection Early Phase
RFadol_maxLP_POSTE_model <- glmmTMB(pass_maxLP ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_EARLY"), family=stats::gaussian)
summary(RFadol_maxLP_POSTE_model)

##### Compare full model to the null model
RFadol_maxLP_POSTE_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_EARLY"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_POSTE_nullmodel, RFadol_maxLP_POSTE_model, test = "LRT") # Full model is a better fit than the null mode

#### Post hoc tests with emmeans and pairs
RFadol_maxLP_POSTE_emmeans <-emmeans(RFadol_maxLP_POSTE_model, pairwise ~ Group|SessionTime, adjust = 'none')
summary(RFadol_maxLP_POSTE_emmeans)


#### Post-Injection Late Phase
RFadol_maxLP_POSTL_model <- glmmTMB(pass_maxLP ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_LATE"), family=stats::gaussian)
summary(RFadol_maxLP_POSTL_model)

##### Compare full model to the null model
RFadol_maxLP_POSTL_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_LATE"), family = stats::gaussian()) # Make the null model
anova(RFadol_maxLP_POSTL_nullmodel, RFadol_maxLP_POSTL_model, test = "LRT") # Full model is a better fit than the null mode\


### Percent Change Group Analysis -----
#### Mixed effects model by phase and sex with random pass_maxLPercept by subject
RFadolescent_pass_maxLP_PC_model <- glmmTMB(pass_maxLP_PC ~ Phase*Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase!="BL"), family=stats::gaussian)
summary(RFadolescent_pass_maxLP_PC_model)

#### Compare full model to the null model
RFadolescent_pass_maxLP_PC_nullmodel <- glmmTMB(pass_maxLP_PC ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase!="BL"), family = stats::gaussian()) # Make the null model
anova(RFadolescent_pass_maxLP_PC_nullmodel, RFadolescent_pass_maxLP_PC_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadolescent_pass_maxLP_PC_within_emmeans <-emmeans(RFadolescent_pass_maxLP_PC_model, pairwise ~ Phase|Group|SessionTime, adjust = 'none')
summary(RFadolescent_pass_maxLP_PC_within_emmeans)

RFadolescent_pass_maxLP_PC_between_emmeans <-emmeans(RFadolescent_pass_maxLP_PC_model, pairwise ~ Group|Phase|SessionTime, adjust = 'none')
summary(RFadolescent_pass_maxLP_PC_between_emmeans)



### Export Analysis Tables -----
#### Phase Tables
RFadol_maxLP_M_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_M_model)
RFadol_maxLP_M_model_fittable <- mixedmodel_fittable(RFadol_maxLP_M_model, RFadol_maxLP_M_nullmodel)
RFadol_maxLP_M_model_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTime(RFadol_maxLP_M_emmeans,RFadol_maxLP_M_model)

write.csv(RFadol_maxLP_M_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Phase_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_M_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Phase_GroupM_model_fittable.csv",sep=""))
write.csv(RFadol_maxLP_M_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Phase_GroupM_emmeans_contraststable.csv",sep=""))

RFadol_maxLP_S_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_S_model)
RFadol_maxLP_S_model_fittable <- mixedmodel_fittable(RFadol_maxLP_S_model, RFadol_maxLP_S_nullmodel)
RFadol_maxLP_S_model_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTime(RFadol_maxLP_S_emmeans,RFadol_maxLP_S_model)

write.csv(RFadol_maxLP_S_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Phase_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_S_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Phase_GroupS_model_fittable.csv",sep=""))
write.csv(RFadol_maxLP_S_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Phase_GroupS_emmeans_contraststable.csv",sep=""))

#### Phase percent change
RFadolescent_pass_maxLP_PC_model_paramtable <- mixedmodel_paramtable(RFadolescent_pass_maxLP_PC_model)
RFadolescent_pass_maxLP_PC_model_fittable <- mixedmodel_fittable(RFadolescent_pass_maxLP_PC_model, RFadolescent_pass_maxLP_PC_nullmodel)
RFadolescent_pass_maxLP_PC_between_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTimePhase(RFadolescent_pass_maxLP_PC_between_emmeans,RFadolescent_pass_maxLP_PC_model)

write.csv(RFadolescent_pass_maxLP_PC_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_PhasePC_model_paramtable.csv",sep=""))
write.csv(RFadolescent_pass_maxLP_PC_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_PhasePC_model_fittable.csv",sep=""))
write.csv(RFadolescent_pass_maxLP_PC_between_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_PhasePC_emmeans_contraststable.csv",sep=""))

#### Session Tables
##### Baseline
RFadol_maxLP_M_BL_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_M_BL_model)
RFadol_maxLP_M_BL_model_fittable <- mixedmodel_fittable(RFadol_maxLP_M_BL_model, RFadol_maxLP_M_BL_nullmodel)

write.csv(RFadol_maxLP_M_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_BL_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_M_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_BL_GroupM_model_fittable.csv",sep=""))

RFadol_maxLP_S_BL_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_S_BL_model)
RFadol_maxLP_S_BL_model_fittable <- mixedmodel_fittable(RFadol_maxLP_S_BL_model, RFadol_maxLP_S_BL_nullmodel)

write.csv(RFadol_maxLP_S_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_BL_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_S_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_BL_GroupS_model_fittable.csv",sep=""))

##### Injection
RFadol_maxLP_M_INJ_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_M_INJ_model)
RFadol_maxLP_M_INJ_model_fittable <- mixedmodel_fittable(RFadol_maxLP_M_INJ_model, RFadol_maxLP_M_INJ_nullmodel)

write.csv(RFadol_maxLP_M_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_INJ_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_M_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_INJ_GroupM_model_fittable.csv",sep=""))

RFadol_maxLP_S_INJ_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_S_INJ_model)
RFadol_maxLP_S_INJ_model_fittable <- mixedmodel_fittable(RFadol_maxLP_S_INJ_model, RFadol_maxLP_S_INJ_nullmodel)

write.csv(RFadol_maxLP_S_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_INJ_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_S_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_INJ_GroupS_model_fittable.csv",sep=""))

##### Post Injection - Early
RFadol_maxLP_M_POSTE_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_M_POSTE_model)
RFadol_maxLP_M_POSTE_model_fittable <- mixedmodel_fittable(RFadol_maxLP_M_POSTE_model, RFadol_maxLP_M_POSTE_nullmodel)

write.csv(RFadol_maxLP_M_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_POSTE_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_M_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_POSTE_GroupM_model_fittable.csv",sep=""))

RFadol_maxLP_S_POSTE_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_S_POSTE_model)
RFadol_maxLP_S_POSTE_model_fittable <- mixedmodel_fittable(RFadol_maxLP_S_POSTE_model, RFadol_maxLP_S_POSTE_nullmodel)

write.csv(RFadol_maxLP_S_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_POSTE_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_S_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_POSTE_GroupS_model_fittable.csv",sep=""))

##### Post Injection - Late
RFadol_maxLP_M_POSTL_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_M_POSTL_model)
RFadol_maxLP_M_POSTL_model_fittable <- mixedmodel_fittable(RFadol_maxLP_M_POSTL_model, RFadol_maxLP_M_POSTL_nullmodel)

write.csv(RFadol_maxLP_M_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_POSTL_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_M_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_POSTL_GroupM_model_fittable.csv",sep=""))

RFadol_maxLP_S_POSTL_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_S_POSTL_model)
RFadol_maxLP_S_POSTL_model_fittable <- mixedmodel_fittable(RFadol_maxLP_S_POSTL_model, RFadol_maxLP_S_POSTL_nullmodel)

write.csv(RFadol_maxLP_S_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_POSTL_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_S_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_Session_POSTL_GroupS_model_fittable.csv",sep=""))

#### Between Group Tables
##### Baseline
RFadol_maxLP_BL_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_BL_model)
RFadol_maxLP_BL_model_fittable <- mixedmodel_fittable(RFadol_maxLP_BL_model, RFadol_maxLP_BL_nullmodel)

write.csv(RFadol_maxLP_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_BetweenGroups_BL_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_BetweenGroups_BL_model_fittable.csv",sep=""))

##### Injection
RFadol_maxLP_INJ_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_INJ_model)
RFadol_maxLP_INJ_model_fittable <- mixedmodel_fittable(RFadol_maxLP_INJ_model, RFadol_maxLP_INJ_nullmodel)

write.csv(RFadol_maxLP_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_BetweenGroups_INJ_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_BetweenGroups_INJ_model_fittable.csv",sep=""))

##### Post Injection - Early
RFadol_maxLP_POSTE_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_POSTE_model)
RFadol_maxLP_POSTE_model_fittable <- mixedmodel_fittable(RFadol_maxLP_POSTE_model, RFadol_maxLP_POSTE_nullmodel)
RFadol_maxLP_POSTE_model_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTime(RFadol_maxLP_POSTE_emmeans,RFadol_maxLP_POSTE_model)

write.csv(RFadol_maxLP_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_BetweenGroups_POSTE_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_BetweenGroups_POSTE_model_fittable.csv",sep=""))
write.csv(RFadol_maxLP_POSTE_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_BetweenGroups_POSTE_emmeans_contraststable.csv",sep=""))

##### Post Injection - Late
RFadol_maxLP_POSTL_model_paramtable <- mixedmodel_paramtable(RFadol_maxLP_POSTL_model)
RFadol_maxLP_POSTL_model_fittable <- mixedmodel_fittable(RFadol_maxLP_POSTL_model, RFadol_maxLP_POSTL_nullmodel)

write.csv(RFadol_maxLP_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_BetweenGroups_POSTL_model_paramtable.csv",sep=""))
write.csv(RFadol_maxLP_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_maxLP_BetweenGroups_POSTL_model_fittable.csv",sep=""))



## Total Session Lever Presses -----
### Scale the data -----
#### NOTE: session_sumLP must be scaled to prevent convergence issues in the later within phase day models. 
#### All models are fit to the scaled and centered sumLP data, and tables will be back-transformed for interpretation at the end.
RFadolescentsubject_loco$session_sumLP_scaled <- scale(RFadolescentsubject_loco$session_sumLP)

# Save scaling parameters for back-transformation
mean_session_sumLP <- mean(RFadolescentsubject_loco$session_sumLP)
sd_session_sumLP <- sd(RFadolescentsubject_loco$session_sumLP)

### Phase Analysis by Group -----
#### Adolescent Saline
##### Mixed effects model by group, phase, session time, and sex with random intercept by subject
RFadol_sumLP_S_model <- glmmTMB(session_sumLP_scaled ~ Phase*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolS"), family=stats::gaussian)
summary(RFadol_sumLP_S_model)

# Convert random-effect variance and SD back to raw scale
sumLP_S_var_scaled <- 0.1040
sumLP_S_sd_scaled  <- sqrt(sumLP_S_var_scaled)

sumLP_S_var_raw <- sumLP_S_var_scaled * (sd_session_sumLP^2)
sumLP_S_sd_raw  <- sumLP_S_sd_scaled  * sd_session_sumLP

sumLP_S_var_raw
sumLP_S_sd_raw


#### Compare full model to the null model
RFadol_sumLP_S_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_S_nullmodel, RFadol_sumLP_S_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadol_sumLP_S_emmeans <-emmeans(RFadol_sumLP_S_model, pairwise ~ Phase|SessionTime, adjust = 'none')
summary(RFadol_sumLP_S_emmeans)


#### Adolescent Morphine
##### Mixed effects model by group, phase, session time, and sex with random intercept by subject
RFadol_sumLP_M_model <- glmmTMB(session_sumLP_scaled ~ Phase*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolM"), family=stats::gaussian)
summary(RFadol_sumLP_M_model)

# Convert random-effect variance and SD back to raw scale
sumLP_M_var_scaled <- 0.4852
sumLP_M_sd_scaled  <- sqrt(sumLP_M_var_scaled)

sumLP_M_var_raw <- sumLP_M_var_scaled * (sd_session_sumLP^2)
sumLP_M_sd_raw  <- sumLP_M_sd_scaled  * sd_session_sumLP

sumLP_M_var_raw
sumLP_M_sd_raw


#### Compare full model to the null model
RFadol_sumLP_M_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_M_nullmodel, RFadol_sumLP_M_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadol_sumLP_M_emmeans <-emmeans(RFadol_sumLP_M_model, pairwise ~ Phase|SessionTime, adjust = 'none')
summary(RFadol_sumLP_M_emmeans)

### Session Analysis by Group -----
#### Adolescent Morphine
##### Baseline
##### Mixed effects model by day with random intercept by subject
RFadol_sumLP_M_BL_model <- glmmTMB(session_sumLP_scaled ~ Day_raw*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolM"),
                                   family = stats::gaussian())
summary(RFadol_sumLP_M_BL_model)

# Convert random-effect variance and SD back to raw scale
sumLP_M_BL_var_scaled <- 1.1973
sumLP_M_BL_sd_scaled  <- sqrt(sumLP_M_BL_var_scaled)

sumLP_M_BL_var_raw <- sumLP_M_BL_var_scaled * (sd_session_sumLP^2)
sumLP_M_BL_sd_raw  <- sumLP_M_BL_sd_scaled  * sd_session_sumLP

sumLP_M_BL_var_raw
sumLP_M_BL_sd_raw

##### Compare full model to the null model
RFadol_sumLP_M_BL_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_M_BL_model, RFadol_sumLP_M_BL_nullmodel, test = "LRT") # Full model is a better fit than the null mode

###### Overall, there is no significant effect of Session in the baseline phase for the morphine group.


#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_sumLP_M_INJ_model <- glmmTMB(session_sumLP_scaled ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolM"),
                                    family = stats::gaussian())
summary(RFadol_sumLP_M_INJ_model)

# Convert random-effect variance and SD back to raw scale
sumLP_M_INJ_var_scaled <- 0.950
sumLP_M_INJ_sd_scaled  <- sqrt(sumLP_M_INJ_var_scaled)

sumLP_M_INJ_var_raw <- sumLP_M_INJ_var_scaled * (sd_session_sumLP^2)
sumLP_M_INJ_sd_raw  <- sumLP_M_INJ_sd_scaled  * sd_session_sumLP

sumLP_M_INJ_var_raw
sumLP_M_INJ_sd_raw

##### Compare full model to the null model
RFadol_sumLP_M_INJ_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_M_INJ_nullmodel, RFadol_sumLP_M_INJ_model, test = "LRT") # Full model is a better fit than the null mode

###### Overall, there are some effects by session of morphine admin with increases in AM and decreases in PM sessions


#### Post-Injection - Early Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_sumLP_M_POSTE_model <- glmmTMB(session_sumLP_scaled ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolM"),
                                      family = stats::gaussian())
summary(RFadol_sumLP_M_POSTE_model)

# Convert random-effect variance and SD back to raw scale
sumLP_M_POSTE_var_scaled <- 0.4308
sumLP_M_POSTE_sd_scaled  <- sqrt(sumLP_M_POSTE_var_scaled)

sumLP_M_POSTE_var_raw <- sumLP_M_POSTE_var_scaled * (sd_session_sumLP^2)
sumLP_M_POSTE_sd_raw  <- sumLP_M_POSTE_sd_scaled  * sd_session_sumLP

sumLP_M_POSTE_var_raw
sumLP_M_POSTE_sd_raw

##### Compare full model to the null model
RFadol_sumLP_M_POSTE_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_M_POSTE_nullmodel, RFadol_sumLP_M_POSTE_model, test = "LRT") # Full model is a better fit than the null mode

###### Overall, sessions 17-20 are significantly lower than baseline for the morphine group.


#### Post-Injection - Late Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_sumLP_M_POSTL_model <- glmmTMB(session_sumLP_scaled ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolM"),
                                      family = stats::gaussian())
summary(RFadol_sumLP_M_POSTL_model)


# Convert random-effect variance and SD back to raw scale
sumLP_M_POSTL_var_scaled <- 0.3407
sumLP_M_POSTL_sd_scaled  <- sqrt(sumLP_M_POSTL_var_scaled)

sumLP_M_POSTL_var_raw <- sumLP_M_POSTL_var_scaled * (sd_session_sumLP^2)
sumLP_M_POSTL_sd_raw  <- sumLP_M_POSTL_sd_scaled  * sd_session_sumLP

sumLP_M_POSTL_var_raw
sumLP_M_POSTL_sd_raw

##### Compare full model to the null model
RFadol_sumLP_M_POSTL_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolM"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_M_POSTL_nullmodel, RFadol_sumLP_M_POSTL_model, test = "LRT") # Full model is a better fit than the null mode

###### Overall, all late phase sessions are significantly lower than baseline for the morphine group.



#### Adolescent Saline
##### Baseline
##### Mixed effects model by day with random intercept by subject
RFadol_sumLP_S_BL_model <- glmmTMB(session_sumLP_scaled ~ Day_raw*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolS"),
                                   family = stats::gaussian())
summary(RFadol_sumLP_S_BL_model)

# Convert random-effect variance and SD back to raw scale
sumLP_S_BL_var_scaled <- 0.1513
sumLP_S_BL_sd_scaled  <- sqrt(sumLP_S_BL_var_scaled)

sumLP_S_BL_var_raw <- sumLP_S_BL_var_scaled * (sd_session_sumLP^2)
sumLP_S_BL_sd_raw  <- sumLP_S_BL_sd_scaled  * sd_session_sumLP

sumLP_S_BL_var_raw
sumLP_S_BL_sd_raw

##### Compare full model to the null model
RFadol_sumLP_S_BL_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_S_BL_model, RFadol_sumLP_S_BL_nullmodel, test = "LRT") # Full model is a better fit than the null mode
###### Overall, there is no significant effect of Session in the baseline phase for the saline group.


#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_sumLP_S_INJ_model <- glmmTMB(session_sumLP_scaled ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolS"),
                                    family = stats::gaussian())
summary(RFadol_sumLP_S_INJ_model)

# Convert random-effect variance and SD back to raw scale
sumLP_S_INJ_var_scaled <- 0.08909
sumLP_S_INJ_sd_scaled  <- sqrt(sumLP_S_INJ_var_scaled)

sumLP_S_INJ_var_raw <- sumLP_S_INJ_var_scaled * (sd_session_sumLP^2)
sumLP_S_INJ_sd_raw  <- sumLP_S_INJ_sd_scaled  * sd_session_sumLP

sumLP_S_INJ_var_raw
sumLP_S_INJ_sd_raw

##### Compare full model to the null model
RFadol_sumLP_S_INJ_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","INJ")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_S_INJ_nullmodel, RFadol_sumLP_S_INJ_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, sumLP decreases significantly across the injection phase for the saline group.


#### Post-Injection - Early Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_sumLP_S_POSTE_model <- glmmTMB(session_sumLP_scaled ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolS"),
                                      family = stats::gaussian())
summary(RFadol_sumLP_S_POSTE_model)

# Convert random-effect variance and SD back to raw scale
sumLP_S_POSTE_var_scaled <- 0.07623
sumLP_S_POSTE_sd_scaled  <- sqrt(sumLP_S_POSTE_var_scaled)

sumLP_S_POSTE_var_raw <- sumLP_S_POSTE_var_scaled * (sd_session_sumLP^2)
sumLP_S_POSTE_sd_raw  <- sumLP_S_POSTE_sd_scaled  * sd_session_sumLP

sumLP_S_POSTE_var_raw
sumLP_S_POSTE_sd_raw

##### Compare full model to the null model
RFadol_sumLP_S_POSTE_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_EARLY")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_S_POSTE_nullmodel, RFadol_sumLP_S_POSTE_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, sumLP is lower than baseline for the saline group.


#### Post-Injection - Late Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadol_sumLP_S_POSTL_model <- glmmTMB(session_sumLP_scaled ~ Day*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolS"),
                                      family = stats::gaussian())
summary(RFadol_sumLP_S_POSTL_model)

# Convert random-effect variance and SD back to raw scale
sumLP_S_POSTL_var_scaled <- 0.1236
sumLP_S_POSTL_sd_scaled  <- sqrt(sumLP_S_POSTL_var_scaled)

sumLP_S_POSTL_var_raw <- sumLP_S_POSTL_var_scaled * (sd_session_sumLP^2)
sumLP_S_POSTL_sd_raw  <- sumLP_S_POSTL_sd_scaled  * sd_session_sumLP

sumLP_S_POSTL_var_raw
sumLP_S_POSTL_sd_raw

##### Compare full model to the null model
RFadol_sumLP_S_POSTL_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase %in% c("BL","POST_LATE")&Group=="AdolS"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_S_POSTL_nullmodel, RFadol_sumLP_S_POSTL_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, sumLP is lower than baseline in post late phase for the saline group.


### Between Group Analysis by Phase -----
#### Baseline
RFadol_sumLP_BL_model <- glmmTMB(session_sumLP_scaled ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"), family=stats::gaussian)
summary(RFadol_sumLP_BL_model)

##### Compare full model to the null model
RFadol_sumLP_BL_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_BL_model, RFadol_sumLP_BL_nullmodel, test = "LRT") # Full model is a better fit than the null model


#### Injection
RFadol_sumLP_INJ_model <- glmmTMB(session_sumLP_scaled ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="INJ"), family=stats::gaussian)
summary(RFadol_sumLP_INJ_model)

##### Compare full model to the null model
RFadol_sumLP_INJ_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="INJ"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_INJ_nullmodel, RFadol_sumLP_INJ_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadol_sumLP_INJ_emmeans <-emmeans(RFadol_sumLP_INJ_model, pairwise ~ Group|SessionTime, adjust = 'none')
summary(RFadol_sumLP_INJ_emmeans)

#### Post-Injection Early Phase
RFadol_sumLP_POSTE_model <- glmmTMB(session_sumLP_scaled ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_EARLY"), family=stats::gaussian)
summary(RFadol_sumLP_POSTE_model)

##### Compare full model to the null model
RFadol_sumLP_POSTE_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_EARLY"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_POSTE_nullmodel, RFadol_sumLP_POSTE_model, test = "LRT") # Full model is not a better fit than the null model


#### Post-Injection Late Phase
RFadol_sumLP_POSTL_model <- glmmTMB(session_sumLP_scaled ~ Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_LATE"), family=stats::gaussian)
summary(RFadol_sumLP_POSTL_model)

##### Compare full model to the null model
RFadol_sumLP_POSTL_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase=="POST_LATE"), family = stats::gaussian()) # Make the null model
anova(RFadol_sumLP_POSTL_nullmodel, RFadol_sumLP_POSTL_model, test = "LRT") # Full model is not a better fit than the null model

### Percent Change Group Analysis -----
#### Mixed effects model by phase and sex with random session_sumLPercept by subject
RFadolescent_session_sumLP_PC_model <- glmmTMB(session_sumLP_PC ~ Phase*Group*SessionTime + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase!="BL"), family=stats::gaussian)
summary(RFadolescent_session_sumLP_PC_model)

#### Compare full model to the null model
RFadolescent_session_sumLP_PC_nullmodel <- glmmTMB(session_sumLP_PC ~ 1 + Sex + (1|SubjectID), data = subset(RFadolescentsubject_loco, Phase!="BL"), family = stats::gaussian()) # Make the null model
anova(RFadolescent_session_sumLP_PC_nullmodel, RFadolescent_session_sumLP_PC_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadolescent_session_sumLP_PC_within_emmeans <-emmeans(RFadolescent_session_sumLP_PC_model, pairwise ~ Phase|Group|SessionTime, adjust = 'none')
summary(RFadolescent_session_sumLP_PC_within_emmeans)

RFadolescent_session_sumLP_PC_between_emmeans <-emmeans(RFadolescent_session_sumLP_PC_model, pairwise ~ Group|Phase|SessionTime, adjust = 'none')
summary(RFadolescent_session_sumLP_PC_between_emmeans)



### Export Analysis Tables -----
#### Phase Tables
RFadol_sumLP_M_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_M_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_M_model_fittable <- mixedmodel_fittable(RFadol_sumLP_M_model, RFadol_sumLP_M_nullmodel)
RFadol_sumLP_M_model_emmeans_contraststable <- mixedmodel_rescalecontrasts(mixedmodel_emmeanstable_SessionTime(RFadol_sumLP_M_emmeans,RFadol_sumLP_M_model),sd_session_sumLP)

write.csv(RFadol_sumLP_M_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Phase_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_M_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Phase_GroupM_model_fittable.csv",sep=""))
write.csv(RFadol_sumLP_M_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Phase_GroupM_emmeans_contraststable.csv",sep=""))

RFadol_sumLP_S_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_S_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_S_model_fittable <- mixedmodel_fittable(RFadol_sumLP_S_model, RFadol_sumLP_S_nullmodel)
RFadol_sumLP_S_model_emmeans_contraststable <- mixedmodel_rescalecontrasts(mixedmodel_emmeanstable_SessionTime(RFadol_sumLP_S_emmeans,RFadol_sumLP_S_model),sd_session_sumLP)

write.csv(RFadol_sumLP_S_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Phase_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_S_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Phase_GroupS_model_fittable.csv",sep=""))
write.csv(RFadol_sumLP_S_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Phase_GroupS_emmeans_contraststable.csv",sep=""))

#### Phase percent change
RFadolescent_session_sumLP_PC_model_paramtable <- mixedmodel_paramtable(RFadolescent_session_sumLP_PC_model)
RFadolescent_session_sumLP_PC_model_fittable <- mixedmodel_fittable(RFadolescent_session_sumLP_PC_model, RFadolescent_session_sumLP_PC_nullmodel)
RFadolescent_session_sumLP_PC_between_emmeans_contraststable <- mixedmodel_emmeanstable_SessionTimePhase(RFadolescent_session_sumLP_PC_between_emmeans,RFadolescent_session_sumLP_PC_model)

write.csv(RFadolescent_session_sumLP_PC_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_PhasePC_model_paramtable.csv",sep=""))
write.csv(RFadolescent_session_sumLP_PC_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_PhasePC_model_fittable.csv",sep=""))
write.csv(RFadolescent_session_sumLP_PC_between_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_PhasePC_emmeans_contraststable.csv",sep=""))

#### Session Tables
##### Baseline
RFadol_sumLP_M_BL_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_M_BL_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_M_BL_model_fittable <- mixedmodel_fittable(RFadol_sumLP_M_BL_model, RFadol_sumLP_M_BL_nullmodel)

write.csv(RFadol_sumLP_M_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_BL_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_M_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_BL_GroupM_model_fittable.csv",sep=""))

RFadol_sumLP_S_BL_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_S_BL_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_S_BL_model_fittable <- mixedmodel_fittable(RFadol_sumLP_S_BL_model, RFadol_sumLP_S_BL_nullmodel)

write.csv(RFadol_sumLP_S_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_BL_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_S_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_BL_GroupS_model_fittable.csv",sep=""))

##### Injection
RFadol_sumLP_M_INJ_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_M_INJ_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_M_INJ_model_fittable <- mixedmodel_fittable(RFadol_sumLP_M_INJ_model, RFadol_sumLP_M_INJ_nullmodel)

write.csv(RFadol_sumLP_M_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_INJ_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_M_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_INJ_GroupM_model_fittable.csv",sep=""))

RFadol_sumLP_S_INJ_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_S_INJ_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_S_INJ_model_fittable <- mixedmodel_fittable(RFadol_sumLP_S_INJ_model, RFadol_sumLP_S_INJ_nullmodel)

write.csv(RFadol_sumLP_S_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_INJ_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_S_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_INJ_GroupS_model_fittable.csv",sep=""))

##### Post Injection - Early
RFadol_sumLP_M_POSTE_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_M_POSTE_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_M_POSTE_model_fittable <- mixedmodel_fittable(RFadol_sumLP_M_POSTE_model, RFadol_sumLP_M_POSTE_nullmodel)

write.csv(RFadol_sumLP_M_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_POSTE_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_M_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_POSTE_GroupM_model_fittable.csv",sep=""))

RFadol_sumLP_S_POSTE_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_S_POSTE_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_S_POSTE_model_fittable <- mixedmodel_fittable(RFadol_sumLP_S_POSTE_model, RFadol_sumLP_S_POSTE_nullmodel)

write.csv(RFadol_sumLP_S_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_POSTE_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_S_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_POSTE_GroupS_model_fittable.csv",sep=""))

##### Post Injection - Late
RFadol_sumLP_M_POSTL_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_M_POSTL_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_M_POSTL_model_fittable <- mixedmodel_fittable(RFadol_sumLP_M_POSTL_model, RFadol_sumLP_M_POSTL_nullmodel)

write.csv(RFadol_sumLP_M_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_POSTL_GroupM_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_M_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_POSTL_GroupM_model_fittable.csv",sep=""))

RFadol_sumLP_S_POSTL_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_S_POSTL_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_S_POSTL_model_fittable <- mixedmodel_fittable(RFadol_sumLP_S_POSTL_model, RFadol_sumLP_S_POSTL_nullmodel)

write.csv(RFadol_sumLP_S_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_POSTL_GroupS_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_S_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_Session_POSTL_GroupS_model_fittable.csv",sep=""))

#### Between Group Tables
##### Baseline
RFadol_sumLP_BL_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_BL_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_BL_model_fittable <- mixedmodel_fittable(RFadol_sumLP_BL_model, RFadol_sumLP_BL_nullmodel)

write.csv(RFadol_sumLP_BL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_BetweenGroups_BL_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_BL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_BetweenGroups_BL_model_fittable.csv",sep=""))

##### Injection
RFadol_sumLP_INJ_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_INJ_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_INJ_model_fittable <- mixedmodel_fittable(RFadol_sumLP_INJ_model, RFadol_sumLP_INJ_nullmodel)
RFadol_sumLP_INJ_model_emmeans_contraststable <- mixedmodel_rescalecontrasts(mixedmodel_emmeanstable_SessionTime(RFadol_sumLP_INJ_emmeans,RFadol_sumLP_INJ_model),sd_session_sumLP)

write.csv(RFadol_sumLP_INJ_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_BetweenGroups_INJ_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_INJ_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_BetweenGroups_INJ_model_fittable.csv",sep=""))
write.csv(RFadol_sumLP_INJ_model_emmeans_contraststable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_BetweenGroups_INJ_emmeans_contraststable.csv",sep=""))

##### Post Injection - Early
RFadol_sumLP_POSTE_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_POSTE_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_POSTE_model_fittable <- mixedmodel_fittable(RFadol_sumLP_POSTE_model, RFadol_sumLP_POSTE_nullmodel)

write.csv(RFadol_sumLP_POSTE_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_BetweenGroups_POSTE_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_POSTE_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_BetweenGroups_POSTE_model_fittable.csv",sep=""))

##### Post Injection - Late
RFadol_sumLP_POSTL_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadol_sumLP_POSTL_model),sd_session_sumLP,mean_session_sumLP)
RFadol_sumLP_POSTL_model_fittable <- mixedmodel_fittable(RFadol_sumLP_POSTL_model, RFadol_sumLP_POSTL_nullmodel)

write.csv(RFadol_sumLP_POSTL_model_paramtable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_BetweenGroups_POSTL_model_paramtable.csv",sep=""))
write.csv(RFadol_sumLP_POSTL_model_fittable, paste(path_analysis_RRFtest,"RRFadolescent_sumLP_BetweenGroups_POSTL_model_fittable.csv",sep=""))



## FIGURES -----
### Supplementary Figure 3 Panels -----
#### Max Lever Presses -----
ymax_maxLPC <- 10
ymin_maxLPC <- -18

ymax_maxLP_bar <- 60
ymin_maxLP_bar <- 0

injmarker_maxLP <- c(rep(ymin_maxLPC,7))
injticks_maxLP <- data.frame(injsession,injmarker_maxLP)

# Max LP Change Line
maxLPC_LP_line <- ggplot(subset(RFadolmeans_loco,Phase!="BL"), aes(x=sessionlabel, y=pass_maxLPC_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadolmeans_loco,Phase!="BL"), aes(ymin=pass_maxLPC_mean-pass_maxLPC_se, ymax=pass_maxLPC_mean+pass_maxLPC_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name ="Postnatal Day", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLPC, ymax_maxLPC), breaks = seq(-16, ymax_maxLPC, 8)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("Max Lever Presses", sep = "")) +
  geom_point(data=injticks_maxLP, aes(x=injsession, y=injmarker_maxLP, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE)+ # Add injection ticks
  annotate("text", x = 17, y = 1, label = "*", size = mystarsize) + # Add asterisks for morphine sessions 
  annotate("text", x = 19, y = 1, label = "*", size = mystarsize)+
  annotate("text", x = 21, y = -2.25, label = "*", size = mystarsize)+
  annotate("text", x = 25, y = -2.25, label = "*", size = mystarsize)+
  annotate("text", x = 27, y = 1, label = "*", size = mystarsize)+
  annotate("text", x = 29, y = 1, label = "*", size = mystarsize)+
  annotate("text", x = 31, y = 1, label = "*", size = mystarsize)+
  annotate("text", x = 6, y = -9.5, label = "#", size = 2) +
  annotate("text", x = 7, y = -10.75, label = "#", size = 2)+ # Add hashtags for saline sessions
  annotate("text", x = 9, y = -10.5, label = "#", size = 2)+
  annotate("text", x = 11, y = -11.75, label = "#", size = 2)+
  annotate("text", x = 13, y = -10.75, label = "#", size = 2)+
  annotate("text", x = 17, y = -12.75, label = "#", size = 2)+
  annotate("text", x = 21, y = -14, label = "#", size = 2)+
  annotate("text", x = 25, y = -13.5, label = "#", size = 2)+
  annotate("text", x = 27, y = -12.5, label = "#", size = 2)+
  annotate("text", x = 29, y = -13.5, label = "#", size = 2)+
  annotate("text", x = 31, y = -14.25, label = "#", size = 2)

maxLPC_LP_line

# Saline Group Max LP by Phase Bar
maxLP_LP_S_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "AM"), aes(x=Phase, y=pass_maxLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolS" & SessionTime == "AM"), mapping = aes(x = Phase, y = pass_maxLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "AM"), aes(ymin=pass_maxLP_mean-pass_maxLP_se, ymax=pass_maxLP_mean+pass_maxLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLP_bar, ymax_maxLP_bar), breaks = seq(ymin_maxLP_bar, ymax_maxLP_bar,20)) + 
  scale_color_manual("legend", values = Sphasecolors) +
  scale_fill_manual("legend", values = Sphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Max Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Control", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(1.9), y.position=c(38) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2) +
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(47) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(54) ,label.size = 2, label = "###",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(2.1), xmax = c(4), y.position=c(38) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

maxLP_LP_S_bar

# Morphine Group Max LP by Phase Bar
maxLP_LP_M_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "AM"), aes(x=Phase, y=pass_maxLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolM" & SessionTime == "AM"), mapping = aes(x = Phase, y = pass_maxLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "AM"), aes(ymin=pass_maxLP_mean-pass_maxLP_se, ymax=pass_maxLP_mean+pass_maxLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLP_bar, ymax_maxLP_bar), breaks = seq(ymin_maxLP_bar, ymax_maxLP_bar, 20)) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Max Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Morphine", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(50) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(56) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2), xmax = c(4), y.position=c(42) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

maxLP_LP_M_bar

# Create full panel
maxLP <- ggarrange(maxLP_LP_S_bar,maxLP_LP_M_bar, maxLPC_LP_line, ncol=3, nrow=1, widths=c(1.6,1.6,4), labels = c("A", "B", "C"))
maxLP

#### Total Session Lever Presses -----
ymax_sessionsumLPC <- 350
ymin_sessionsumLPC <- -250

ymax_sessionLP_bar <- 850
ymin_sessionLP_bar <- 0

injmarker_sessionsumLP <- c(rep(ymin_sessionsumLPC,7))
injticks_sessionsumLP <- data.frame(injsession,injmarker_sessionsumLP)

# Session Lever Presses Change Line
sessionsumLPC_LP_line <- ggplot(subset(RFadolmeans_loco,Phase!="BL"), aes(x=sessionlabel, y=session_sumLPC_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadolmeans_loco,Phase!="BL"), aes(ymin=session_sumLPC_mean-session_sumLPC_se, ymax=session_sumLPC_mean+session_sumLPC_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name ="Postnatal Day", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sessionsumLPC, ymax_sessionsumLPC), breaks = seq(-150, ymax_sessionsumLPC, 150)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("Total Session Lever Presses", sep = "")) +
  geom_point(data=injticks_sessionsumLP, aes(x=injsession, y=injmarker_sessionsumLP, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE)+ # Add injection ticks
  annotate("text", x = 13, y = 300, label = "*", size = mystarsize) + # Add asterisks for morphine sessions 
  annotate("text", x = 14, y = -190, label = "*", size = mystarsize)+
  annotate("text", x = 17, y = -220, label = "*", size = mystarsize)+
  annotate("text", x = 19, y = -200, label = "*", size = mystarsize)+
  annotate("text", x = 21, y = -220, label = "*", size = mystarsize)+
  annotate("text", x = 23, y = -185, label = "*", size = mystarsize)+
  annotate("text", x = 25, y = -230, label = "*", size = mystarsize)+
  annotate("text", x = 27, y = -220, label = "*", size = mystarsize)+
  annotate("text", x = 29, y = -235, label = "*", size = mystarsize)+
  annotate("text", x = 31, y = -195, label = "*", size = mystarsize)+
  annotate("text", x = 7, y = -130, label = "#", size = 2)+ # Add hashtags for saline sessions
  annotate("text", x = 9, y = -120, label = "#", size = 2)+
  annotate("text", x = 11, y = -125, label = "#", size = 2)+
  annotate("text", x = 13, y = -120, label = "#", size = 2)+
  annotate("text", x = 15, y = 75, label = "#", size = 2)+
  annotate("text", x = 17, y = -20, label = "#", size = 2)+
  annotate("text", x = 19, y = 25, label = "#", size = 2)+
  annotate("text", x = 21, y = -25, label = "#", size = 2)+
  annotate("text", x = 23, y = 25, label = "#", size = 2)+
  annotate("text", x = 25, y = 25, label = "#", size = 2)+
  annotate("text", x = 27, y = 25, label = "#", size = 2)+
  annotate("text", x = 29, y = 25, label = "#", size = 2)+
  annotate("text", x = 31, y = 25, label = "#", size = 2)

sessionsumLPC_LP_line

# Saline Group Session Lever Pressing by Phase Bar
sessionLP_LP_S_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "AM"), aes(x=Phase, y=session_sumLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolS" & SessionTime == "AM"), mapping = aes(x = Phase, y = session_sumLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "AM"), aes(ymin=session_sumLP_mean-session_sumLP_se, ymax=session_sumLP_mean+session_sumLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sessionLP_bar, ymax_sessionLP_bar), breaks = seq(ymin_sessionLP_bar, ymax_sessionLP_bar,250)) + 
  scale_color_manual("legend", values = Sphasecolors) +
  scale_fill_manual("legend", values = Sphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Session Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Control", sep = ""))+
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(440) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(530) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(620) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

sessionLP_LP_S_bar

# Morphine Group Session Lever Pressing by Phase Bar
sessionLP_LP_M_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "AM"), aes(x=Phase, y=session_sumLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolM" & SessionTime == "AM"), mapping = aes(x = Phase, y = session_sumLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "AM"), aes(ymin=session_sumLP_mean-session_sumLP_se, ymax=session_sumLP_mean+session_sumLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sessionLP_bar, ymax_sessionLP_bar), breaks = seq(ymin_sessionLP_bar, ymax_sessionLP_bar, 250)) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Session Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Morphine", sep = "")) +   
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(700) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(800) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2), xmax = c(3), y.position=c(500) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2), xmax = c(4), y.position=c(600) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

sessionLP_LP_M_bar

# Create full panel
sessionLP <- ggarrange(sessionLP_LP_S_bar, sessionLP_LP_M_bar, sessionsumLPC_LP_line, ncol=3, nrow=1, widths=c(1.6,1.6,4), labels = c("D", "E", "F"))
sessionLP


#### Max Lever Presses Percent Change -----
ymax_maxLPPC <- 80
ymin_maxLPPC <- -80

maxLP_PC_bar <- ggplot(subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="AM"), aes(x=Phase, y=pass_maxLPPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="AM"), aes(ymin=pass_maxLPPC_mean-pass_maxLPPC_se, ymax=pass_maxLPPC_mean+pass_maxLPPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLPPC, ymax_maxLPPC), breaks = seq(-80, ymax_maxLPPC, 40)) + 
  scale_x_discrete(expand = c(0.25, 0), name ="Phase", labels=c('INJ','POST EARLY','POST LATE')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change from Baseline") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Max Lever Presses", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)

maxLP_PC_bar

#### Total Session Lever Presses Percent Change -----
ymax_sumLPPC <- 80
ymin_sumLPPC <- -80

sumLP_PC_bar <- ggplot(subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="AM"), aes(x=Phase, y=session_sumLPPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="AM"), aes(ymin=session_sumLPPC_mean-session_sumLPPC_se, ymax=session_sumLPPC_mean+session_sumLPPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  #geom_point(subset(PCmeans_lm_subject, Phase!="BL"), mapping = aes(x = Phase, y = sumLPPC_mean), position = position_dodge(width = 0.95), size = mysubjectpointsize, stroke = mypointstroke, shape = 1, color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sumLPPC, ymax_sumLPPC), breaks = seq(-80, ymax_sumLPPC, 40)) + 
  scale_x_discrete(expand = c(0.25, 0), name ="Phase", labels=c('INJ','POST EARLY','POST LATE')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change from Baseline") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Total Session Lever Presses", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)#+

sumLP_PC_bar

# Create full panel
maxsumLP_PC <- ggarrange(maxLP_PC_bar, sumLP_PC_bar, ncol=2, nrow=1, widths=c(2,2), labels = c("G", "H"))
maxsumLP_PC


### Supplementary Figure 3 -----
suppfigure3 <- ggarrange(maxLP,sessionLP,maxsumLP_PC, ncol=1, nrow=3)
suppfigure3
ggsave(plot=suppfigure3, path=path_figures, device="pdf", filename="SupplementaryF3_AllPanels.pdf", width=RFpanelwidth, height=RFpanelheight*3)


### Supplementary Figure 4 Panels -----


#### PM Intercept -----
ymax_Int_bar <- 30
ymin_Int_bar <- 0

# Saline Group Intercept by Phase Bar
Int_PM_S_bar <- ggplot(subset(RFadolmeans_Phase_LP_lm,Group=="AdolS" & SessionTime == "PM"), aes(x=Phase, y=Int_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_LP_lm, Group=="AdolS" & SessionTime == "PM"), mapping = aes(x = Phase, y = Int_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_LP_lm,Group=="AdolS" & SessionTime == "PM"), aes(ymin=Int_mean-Int_se, ymax=Int_mean+Int_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0, 0), limits = c(ymin_Int_bar, ymax_Int_bar), breaks = seq(ymin_Int_bar, ymax_Int_bar, 10)) + 
  scale_color_manual("legend", values = Sphasecolors) +
  scale_fill_manual("legend", values = Sphasecolors) +
  ylab("Presses at Median Frequency") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Control", sep = "")) + 
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(14) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2) +
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(17.5) ,label.size = 2, label = "##",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(21) ,label.size = 2, label = "##",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

Int_PM_S_bar

# Morphine Group Intercept by Phase Bar
Int_PM_M_bar <- ggplot(subset(RFadolmeans_Phase_LP_lm,Group=="AdolM" & SessionTime == "PM"), aes(x=Phase, y=Int_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_LP_lm, Group=="AdolM" & SessionTime == "PM"), mapping = aes(x = Phase, y = Int_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_LP_lm,Group=="AdolM" & SessionTime == "PM"), aes(ymin=Int_mean-Int_se, ymax=Int_mean+Int_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0, 0), limits = c(ymin_Int_bar, ymax_Int_bar), breaks = seq(ymin_Int_bar, ymax_Int_bar, 10)) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  ylab("Presses at Median Frequency") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Morphine", sep = "")) +
  geom_hline(yintercept = 0, col = 'black', linewidth = 1)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(28) ,label.size = mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

Int_PM_M_bar

# PM Intercept Percent Change
ymax_IntPC <- 120
ymin_IntPC <- -60

Int_PC_PM_bar <- ggplot(subset(RFadolmeans_Phase_LP_lm,Phase!="BL" & SessionTime=="PM"), aes(x=Phase, y=IntPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadolmeans_Phase_LP_lm,Phase!="BL" & SessionTime=="PM"), aes(ymin=IntPC_mean-IntPC_se, ymax=IntPC_mean+IntPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_IntPC, ymax_IntPC), breaks = seq(-60, ymax_IntPC, 60)) + 
  scale_x_discrete(expand = c(0.25, 0), labels=c('INJ','POST EARLY','POST LATE')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change from BL") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Presses at Median Frequency", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)

Int_PC_PM_bar


# Create full panel
Int_PM <- ggarrange(Int_PM_S_bar, Int_PM_M_bar, Int_PC_PM_bar, ncol=3, nrow=1, widths=c(1.6,1.6,4), labels = c("A", "B", "C"))
Int_PM


#### PM Final Pass Locomotion -----
ymax_blockLoco_bar <- 82
ymin_blockLoco_bar <- 0

# Saline Group Final Pass Locomotion by Phase Bar
blockLoco_PM_S_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "PM"), aes(x=Phase, y=finalblock_Loco_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolS" & SessionTime == "PM"), mapping = aes(x = Phase, y = finalblock_Loco_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "PM"), aes(ymin=finalblock_Loco_mean-finalblock_Loco_se, ymax=finalblock_Loco_mean+finalblock_Loco_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blockLoco_bar, ymax_blockLoco_bar), breaks = seq(ymin_blockLoco_bar, ymax_blockLoco_bar,25)) + 
  scale_color_manual("legend", values = Sphasecolors) +
  scale_fill_manual("legend", values = Sphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Block Locomotion") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Control", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(55) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(64) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(73) ,label.size = 2, label = "##",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

blockLoco_PM_S_bar

# Morphine Group Final Pass Locomotion by Phase Bar
blockLoco_PM_M_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "PM"), aes(x=Phase, y=finalblock_Loco_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolM" & SessionTime == "PM"), mapping = aes(x = Phase, y = finalblock_Loco_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "PM"), aes(ymin=finalblock_Loco_mean-finalblock_Loco_se, ymax=finalblock_Loco_mean+finalblock_Loco_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blockLoco_bar, ymax_blockLoco_bar), breaks = seq(ymin_blockLoco_bar, ymax_blockLoco_bar, 25)) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Block Locomotion") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Morphine", sep = ""))
blockLoco_PM_M_bar

# PM Final Block Locomotion Percent Change
ymax_blocklocoPC <- 120
ymin_blocklocoPC <- -60

blockloco_PC_PM_bar <- ggplot(subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="PM"), aes(x=Phase, y=finalblock_LocoPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="PM"), aes(ymin=finalblock_LocoPC_mean-finalblock_LocoPC_se, ymax=finalblock_LocoPC_mean+finalblock_LocoPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blocklocoPC, ymax_blocklocoPC), breaks = seq(-60, ymax_blocklocoPC, 60)) + 
  scale_x_discrete(expand = c(0.25, 0), labels=c('INJ','POST EARLY','POST LATE')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change from BL") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Last Block Locomotion", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)

blockloco_PC_PM_bar


# Create full panel
blockLoco_PM <- ggarrange(blockLoco_PM_S_bar,blockLoco_PM_M_bar, blockloco_PC_PM_bar, ncol=3, nrow=1, widths=c(1.6,1.6,4), labels = c("D", "E", "F"))
blockLoco_PM


#### PM Max Lever Presses -----
ymax_maxLP_bar <- 60
ymin_maxLP_bar <- 0

# Saline Group Max LP by Phase Bar
maxLP_PM_S_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "PM"), aes(x=Phase, y=pass_maxLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolS" & SessionTime == "PM"), mapping = aes(x = Phase, y = pass_maxLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "PM"), aes(ymin=pass_maxLP_mean-pass_maxLP_se, ymax=pass_maxLP_mean+pass_maxLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLP_bar, ymax_maxLP_bar), breaks = seq(ymin_maxLP_bar, ymax_maxLP_bar,20)) + 
  scale_color_manual("legend", values = Sphasecolors) +
  scale_fill_manual("legend", values = Sphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Max Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Control", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(1.9), y.position=c(37) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2) +
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(44) ,label.size = 2, label = "###",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(51) ,label.size = 2, label = "###",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(2.1), xmax = c(3), y.position=c(30) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(2.1), xmax = c(4), y.position=c(37) ,label.size = 2, label = "##",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

maxLP_PM_S_bar

# Morphine Group Max LP by Phase Bar
maxLP_PM_M_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "PM"), aes(x=Phase, y=pass_maxLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolM" & SessionTime == "PM"), mapping = aes(x = Phase, y = pass_maxLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "PM"), aes(ymin=pass_maxLP_mean-pass_maxLP_se, ymax=pass_maxLP_mean+pass_maxLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLP_bar, ymax_maxLP_bar), breaks = seq(ymin_maxLP_bar, ymax_maxLP_bar, 20)) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Max Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Morphine", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(50) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(56) ,label.size = mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2), xmax = c(4), y.position=c(42) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

maxLP_PM_M_bar

# Max Lever Presses Percent Change
ymax_maxLPPC <- 80
ymin_maxLPPC <- -80

maxLP_PC_PM_bar <- ggplot(subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="PM"), aes(x=Phase, y=pass_maxLPPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="PM"), aes(ymin=pass_maxLPPC_mean-pass_maxLPPC_se, ymax=pass_maxLPPC_mean+pass_maxLPPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLPPC, ymax_maxLPPC), breaks = seq(-80, ymax_maxLPPC, 40)) + 
  scale_x_discrete(expand = c(0.25, 0), name ="Phase", labels=c('INJ','POST EARLY','POST LATE')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change from BL") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Max Lever Presses", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)

maxLP_PC_PM_bar

# Create full panel
maxLP_PM <- ggarrange(maxLP_PM_S_bar,maxLP_PM_M_bar, maxLP_PC_PM_bar, ncol=3, nrow=1, widths=c(1.6,1.6,4), labels = c("G", "H", "I"))
maxLP_PM

#### PM Total Session Lever Presses -----
ymax_sessionLP_bar <- 850
ymin_sessionLP_bar <- 0

# Saline Group Session Lever Pressing by Phase Bar
sessionLP_PM_S_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "PM"), aes(x=Phase, y=session_sumLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolS" & SessionTime == "PM"), mapping = aes(x = Phase, y = session_sumLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolS" & SessionTime == "PM"), aes(ymin=session_sumLP_mean-session_sumLP_se, ymax=session_sumLP_mean+session_sumLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sessionLP_bar, ymax_sessionLP_bar), breaks = seq(ymin_sessionLP_bar, ymax_sessionLP_bar,250)) + 
  scale_color_manual("legend", values = Sphasecolors) +
  scale_fill_manual("legend", values = Sphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Session Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Control", sep = ""))+
  geom_bracket(xmin=c(1), xmax = c(1.9), y.position=c(400) ,label.size = 2, label = "##",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(490) ,label.size = 2, label = "###",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(650) ,label.size = 2, label = "###",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(2.1), xmax = c(3), y.position=c(400) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(2), xmax = c(4), y.position=c(580) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

sessionLP_PM_S_bar

# Morphine Group Session Lever Pressing by Phase Bar
sessionLP_PM_M_bar <- ggplot(subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "PM"), aes(x=Phase, y=session_sumLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadolmeans_SubjectPhase_loco, Group=="AdolM" & SessionTime == "PM"), mapping = aes(x = Phase, y = session_sumLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Group=="AdolM" & SessionTime == "PM"), aes(ymin=session_sumLP_mean-session_sumLP_se, ymax=session_sumLP_mean+session_sumLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sessionLP_bar, ymax_sessionLP_bar), breaks = seq(ymin_sessionLP_bar, ymax_sessionLP_bar, 250)) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Session Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Morphine", sep = "")) +   
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(620) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(710) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(800) ,label.size = mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

sessionLP_PM_M_bar

# Total Session Lever Presses Percent Change
ymax_sumLPPC <- 80
ymin_sumLPPC <- -80

sumLP_PC_PM_bar <- ggplot(subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="PM"), aes(x=Phase, y=session_sumLPPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadolmeans_Phase_loco,Phase!="BL" & SessionTime=="PM"), aes(ymin=session_sumLPPC_mean-session_sumLPPC_se, ymax=session_sumLPPC_mean+session_sumLPPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sumLPPC, ymax_sumLPPC), breaks = seq(-80, ymax_sumLPPC, 40)) + 
  scale_x_discrete(expand = c(0.25, 0), name ="Phase", labels=c('INJ','POST EARLY','POST LATE')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change from BL") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Total Session Lever Presses", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)#+

sumLP_PC_PM_bar

# Create full panel
sessionLP_PM <- ggarrange(sessionLP_PM_S_bar, sessionLP_PM_M_bar, sumLP_PC_PM_bar, ncol=3, nrow=1, widths=c(1.6,1.6,4), labels = c("J", "K", "L"))
sessionLP_PM


### Supplementary Figure 4 -----
# Create full panel
suppfigure4 <- ggarrange(Int_PM,blockLoco_PM,maxLP_PM,sessionLP_PM, ncol=1, nrow=4)
suppfigure4
ggsave(plot=suppfigure4, path=path_figures, device="pdf", filename="SupplementaryF4_AllPanels.pdf", width=RFpanelwidth, height=RFpanelheight*4)

