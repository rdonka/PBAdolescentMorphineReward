############# ADULT MORPHINE RE-EXPOSURE ADMINISTRATION ANALYSIS ################
# This program takes in data processed by the script "RF Theta Calculation"    #
# and run statistics to analyze findings. Data from each experiment are loaded,#
# analyzed, and visualized.                                                    #


# Prepare R for Analysis -----
## Set up paths
computerusergithubpath <- c('C:/Users/rmdon/OneDrive/Desktop/GitHub_Repositories/') # Path for Rachel's laptop to github repositories
computeruserboxpath <- c('C:/Users/rmdon/Box/') # Path for Rachel's laptop to Box drive

# Load required libraries
source(paste(computerusergithubpath,"Analysis-PBAdolescentMorphineReward/Functions/Functions_LoadPackages.R",sep=''))

## Load all custom functions
source.all(paste(computerusergithubpath,'JRoitmanICSS/Functions/Standard RF ICSS/',sep=''))
source.all(paste(computerusergithubpath,'Analysis-PBAdolescentMorphineReward/Functions/',sep=''))

# Set working directory to the raw data folder
setwd(paste(computeruserboxpath,'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Raw Data',sep=''))

##### Set up paths to export figures and analysis to
### Overall paths
path_analysis <- c(paste(computeruserboxpath, 'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Analysis/',sep=''))
path_figures <- c(paste(computeruserboxpath, 'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Figures/Figure Panels/',sep=''))

# Specific paths
path_analysis_RFtest <- paste(path_analysis,"Adult Re-Exposure Rate Frequency Morphine Test/",sep="")
path_analysis_RFtest_means <- paste(path_analysis,"Adult Re-Exposure Rate Frequency Morphine Test/Means Tables/",sep="")
path_figures_RFtest <- paste(path_analysis_RFtest,"Figures/",sep="")
path_subjectfigures_RFtest <- paste(path_analysis_RFtest,"Subject Figures/",sep="") # Subject shaping figures path


##### Read in subject key
subjectkey <- read_csv("SubjectKey_Adolescent.csv")
subjectexc <- c()

# PREPARE DATA -----
## Read in data and prepare variables -----
### Prep col names
freqcolnames <- c('141', '126', '112', '100', '89', '79', '71', '63', '56', '50', '45', '40','35','32','28') # Set column names of lever presses for each frequency

### Read in data and fit model
RFadultexpdata_raw <- RFICSS_fitmodel("RawData_RFICSSMorphineTest_AdultExperienced.csv")

### Add Sex and Group to raw
RFadultexpdata_raw$Sex <- subjectkey$Sex[match(RFadultexpdata_raw$SubjectID, subjectkey$SubjectID)]
RFadultexpdata_raw$Group <- subjectkey$ExperimentalGroup[match(RFadultexpdata_raw$SubjectID, subjectkey$SubjectID)]
RFadultexpdata_raw$ShapingInclude <- subjectkey$ShapingInclude[match(RFadultexpdata_raw$SubjectID, subjectkey$SubjectID)]

### Prepare categorical variables and factor
RFadultexpdata_raw$Sex <- factor(RFadultexpdata_raw$Sex, levels = c("M","F"))
RFadultexpdata_raw$Group <- factor(RFadultexpdata_raw$Group, levels=c('AdolS','AdolM'))
RFadultexpdata_raw$Phase_raw <- factor(RFadultexpdata_raw$Phase, levels = c("BL","INJ","POST"))
RFadultexpdata_raw$Phase <- factor(if_else(RFadultexpdata_raw$Phase=="POST" & RFadultexpdata_raw$Day <=3,"POST_EARLY", 
                                                      if_else(RFadultexpdata_raw$Phase=="POST" & RFadultexpdata_raw$Day >=4,"POST_LATE", 
                                                              RFadultexpdata_raw$Phase)), levels = c("BL","INJ","POST_EARLY","POST_LATE"))

RFadultexpdata_raw$Session <- factor(RFadultexpdata_raw$Day) # Factor day by phase
RFadultexpdata_raw$Session_collapsed <- ifelse(RFadultexpdata_raw$Phase == "BL", 0, RFadultexpdata_raw$Day) # Create collapsed Day variable where Baseline=0
RFadultexpdata_raw$Session_collapsed <- factor(RFadultexpdata_raw$Session_collapsed, levels = unique(RFadultexpdata_raw$Session_collapsed)) # Factor collapsed Day variable

### Filter data
RFadultexpdata_raw <- filter(RFadultexpdata_raw,!is.na(Phase), !Condition %in% c("MW10","MW11","MW12"), ShapingInclude=="INCLUDE")

### Add max lever presses, and total lever presses per pass
RFadultexpdata_raw$maxLP <- rowMaxs(as.matrix(RFadultexpdata_raw[freqcolnames])) # Max lever pressing per pass
RFadultexpdata_raw$sumLP <- rowSums(as.matrix(RFadultexpdata_raw[freqcolnames])) # Max lever pressing per pass

### Check Variables
unique(RFadultexpdata_raw$SubjectID)
unique(RFadultexpdata_raw$Sex)
unique(RFadultexpdata_raw$Group)
unique(RFadultexpdata_raw$ShapingInclude)

unique(RFadultexpdata_raw$Day)
unique(RFadultexpdata_raw$Session)
unique(RFadultexpdata_raw$Session_collapsed)

unique(RFadultexpdata_raw$Phase) 
unique(RFadultexpdata_raw$Phase_raw)
#view(RFadultexpdata_raw)

## Prepare subject averages by day -----
### Find subject averages by day
RFadultexpdata_lm <- RFadultexpdata_raw %>% filter(!is.na(Phase), FinalInclusion == "INCLUDE", Pass > 1, !SubjectID %in% subjectexc, ShapingInclude=="INCLUDE") %>%
  group_by(SubjectID, Sex, Group, ShapingInclude, Condition, UniqueSessionID, Phase, Phase_raw, Day, Session, Session_collapsed, Box, Amplitude) %>%
  dplyr::summarize(npass = n_distinct(Pass),
                   mean_theta = mean(theta, na.rm=TRUE), # find the means for theta, M50, alpha, sumLP, and maxLP
                   mean_M50 = mean(M50, na.rm=TRUE),
                   mean_alpha = mean(alpha, na.rm=TRUE),
                   mean_alpha95 = mean(alpha95, na.rm=TRUE))

RFadultexpdata_loco <- RFadultexpdata_raw %>% filter(!is.na(Phase), Pass > 1, !SubjectID %in% subjectexc, ShapingInclude=="INCLUDE") %>%
  group_by(SubjectID, Sex, Group, ShapingInclude, Condition, UniqueSessionID, Phase, Phase_raw, Day, Session, Session_collapsed, Box, Amplitude) %>%
  dplyr::summarize(npass = n_distinct(Pass),
                   pass_sumLP = mean(TotalLeverPresses, na.rm=TRUE),
                   session_sumLP = sum(TotalLeverPresses, na.rm=TRUE),
                   pass_maxLP = mean(maxLP, na.rm=TRUE),
                   finalblock_loco = mean(BeamBreaks, na.rm=TRUE))


## Add percent change from baseline -----
RFadultexpsubject_lm <- RFICSS_LM_findChange(RFadultexpdata_lm,"BL")
RFadultexpsubject_loco <- RFICSS_loco_findChange(RFadultexpdata_loco,"BL")

## Export processed data -----
write.csv(RFadultexpsubject_lm, paste(path_analysis,"RFadult_reexposure_LMdata.csv",sep=""))
write.csv(RFadultexpsubject_loco, paste(path_analysis,"RFadult_reexposure_Locodata.csv",sep=""))


# PREPARE MEANS -----
## Individual Session Means -----
sessionlabels <- (c(1:15))

RFadultexpmeans_lm <- RFadultexpsubject_lm %>%
  group_by(UniqueSessionID,Phase,Phase_raw,Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   theta_mean = mean(mean_theta, na.rm=TRUE), theta_sd = sd(mean_theta, na.rm=TRUE), theta_se = theta_sd / sqrt(n),
                   thetaC_mean = mean(mean_theta_C, na.rm=TRUE), thetaC_sd = sd(mean_theta_C, na.rm=TRUE), thetaC_se = thetaC_sd / sqrt(n),
                   
                   M50_mean = mean(mean_M50, na.rm=TRUE), M50_sd = sd(mean_M50, na.rm=TRUE), M50_se = M50_sd / sqrt(n),
                   M50C_mean = mean(mean_M50_C, na.rm=TRUE), M50C_sd = sd(mean_M50_C, na.rm=TRUE), M50C_se = M50C_sd / sqrt(n),
                   
                   alpha_mean = mean(mean_alpha, na.rm=TRUE), alpha_sd = sd(mean_alpha, na.rm=TRUE), alpha_se = alpha_sd / sqrt(n),
                   alphaC_mean = mean(mean_alpha_C, na.rm=TRUE), alphaC_sd = sd(mean_alpha_C, na.rm=TRUE), alphaC_se = alphaC_sd / sqrt(n))
RFadultexpmeans_lm$sessionlabel = factor(RFadultexpmeans_lm$UniqueSessionID, labels = c(1:20))
RFadultexpmeans_lm


RFadultexpmeans_loco <- RFadultexpsubject_loco %>%
  group_by(UniqueSessionID,Phase,Phase_raw,Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   
                   pass_maxLP_mean = mean(pass_maxLP, na.rm=TRUE), pass_maxLP_sd = sd(pass_maxLP,na.rm=TRUE), pass_maxLP_se = pass_maxLP_sd / sqrt(n),
                   pass_maxLPC_mean = mean(pass_maxLP_C, na.rm=TRUE), pass_maxLPC_sd = sd(pass_maxLP_C,na.rm=TRUE), pass_maxLPC_se = pass_maxLPC_sd / sqrt(n),
                   
                   session_sumLP_mean = mean(session_sumLP, na.rm=TRUE), session_sumLP_sd = sd(session_sumLP,na.rm=TRUE), session_sumLP_se = session_sumLP_sd / sqrt(n),
                   session_sumLPC_mean = mean(session_sumLP_C, na.rm=TRUE), session_sumLPC_sd = sd(session_sumLP_C,na.rm=TRUE), session_sumLPC_se = session_sumLPC_sd / sqrt(n),
                   
                   pass_sumLP_mean = mean(pass_sumLP, na.rm=TRUE), pass_sumLP_sd = sd(pass_sumLP,na.rm=TRUE), pass_sumLP_se = pass_sumLP_sd / sqrt(n),
                   pass_sumLPC_mean = mean(pass_sumLP_C, na.rm=TRUE), pass_sumLPC_sd = sd(pass_sumLP_C,na.rm=TRUE), pass_sumLPC_se = pass_sumLPC_sd / sqrt(n),
                   
                   finalblock_loco_mean = mean(finalblock_loco, na.rm=TRUE), finalblock_loco_sd = sd(finalblock_loco,na.rm=TRUE), finalblock_loco_se = finalblock_loco_sd / sqrt(n),
                   finalblock_locoC_mean = mean(finalblock_loco_C, na.rm=TRUE), finalblock_locoC_sd = sd(finalblock_loco_C,na.rm=TRUE), finalblock_locoC_se = finalblock_locoC_sd / sqrt(n))
RFadultexpmeans_loco$sessionlabel = factor(RFadultexpmeans_loco$UniqueSessionID, labels = c(1:20))
RFadultexpmeans_loco


## Overall Phase Means -----
RFadultexpmeans_Phase_lm <- RFadultexpsubject_lm %>%
  group_by(Phase,Phase_raw,Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   theta_mean = mean(mean_theta, na.rm=TRUE), theta_sd = sd(mean_theta, na.rm=TRUE), theta_se = theta_sd / sqrt(n),
                   thetaC_mean = mean(mean_theta_C, na.rm=TRUE), thetaC_sd = sd(mean_theta_C, na.rm=TRUE), thetaC_se = thetaC_sd / sqrt(n),
                   thetaPC_mean = mean(mean_theta_PC, na.rm=TRUE), thetaPC_sd = sd(mean_theta_PC, na.rm=TRUE), thetaPC_se = thetaPC_sd / sqrt(n),
                   
                   M50_mean = mean(mean_M50, na.rm=TRUE), M50_sd = sd(mean_M50, na.rm=TRUE), M50_se = M50_sd / sqrt(n),
                   M50C_mean = mean(mean_M50_C, na.rm=TRUE), M50C_sd = sd(mean_M50_C, na.rm=TRUE), M50C_se = M50C_sd / sqrt(n),
                   M50PC_mean = mean(mean_M50_PC, na.rm=TRUE), M50PC_sd = sd(mean_M50_PC, na.rm=TRUE), M50PC_se = M50PC_sd / sqrt(n),
                   
                   
                   alpha_mean = mean(mean_alpha, na.rm=TRUE), alpha_sd = sd(mean_alpha, na.rm=TRUE), alpha_se = alpha_sd / sqrt(n),
                   alphaC_mean = mean(mean_alpha_C, na.rm=TRUE), alphaC_sd = sd(mean_alpha_C, na.rm=TRUE), alphaC_se = alphaC_sd / sqrt(n),
                   alphaPC_mean = mean(mean_alpha_PC, na.rm=TRUE), alphaPC_sd = sd(mean_alpha_PC, na.rm=TRUE), alphaPC_se = alphaPC_sd / sqrt(n))


RFadultexpmeans_Phase_loco <- RFadultexpsubject_loco %>%
  group_by(Phase,Phase_raw,Group) %>%
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
                   
                   finalblock_loco_mean = mean(finalblock_loco, na.rm=TRUE), finalblock_loco_sd = sd(finalblock_loco,na.rm=TRUE), finalblock_loco_se = finalblock_loco_sd / sqrt(n),
                   finalblock_locoC_mean = mean(finalblock_loco_C, na.rm=TRUE), finalblock_locoC_sd = sd(finalblock_loco_C,na.rm=TRUE), finalblock_locoC_se = finalblock_locoC_sd / sqrt(n),
                   finalblock_locoPC_mean = mean(finalblock_loco_PC, na.rm=TRUE), finalblock_locoPC_sd = sd(finalblock_loco_PC,na.rm=TRUE), finalblock_locoPC_se = finalblock_locoPC_sd / sqrt(n))


## Subject Phase Means -----
RFadultexpmeans_SubjectPhase_lm <- RFadultexpsubject_lm %>%
  group_by(SubjectID,Sex,Phase,Phase_raw,Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   
                   theta_mean = mean(mean_theta, na.rm=TRUE), theta_sd = sd(mean_theta, na.rm=TRUE), theta_se = theta_sd / sqrt(n),
                   thetaC_mean = mean(mean_theta_C, na.rm=TRUE), thetaC_sd = sd(mean_theta_C, na.rm=TRUE), thetaC_se = thetaC_sd / sqrt(n),
                   thetaPC_mean = mean(mean_theta_PC, na.rm=TRUE), thetaPC_sd = sd(mean_theta_PC, na.rm=TRUE), thetaPC_se = thetaPC_sd / sqrt(n),
                   
                   M50_mean = mean(mean_M50, na.rm=TRUE), M50_sd = sd(mean_M50, na.rm=TRUE), M50_se = M50_sd / sqrt(n),
                   M50C_mean = mean(mean_M50_C, na.rm=TRUE), M50C_sd = sd(mean_M50_C, na.rm=TRUE), M50C_se = M50C_sd / sqrt(n),
                   M50PC_mean = mean(mean_M50_PC, na.rm=TRUE), M50PC_sd = sd(mean_M50_PC, na.rm=TRUE), M50PC_se = M50PC_sd / sqrt(n),
                   
                   alpha_mean = mean(mean_alpha, na.rm=TRUE), alpha_sd = sd(mean_alpha, na.rm=TRUE), alpha_se = alpha_sd / sqrt(n),
                   alphaC_mean = mean(mean_alpha_C, na.rm=TRUE), alphaC_sd = sd(mean_alpha_C, na.rm=TRUE), alphaC_se = alphaC_sd / sqrt(n),
                   alphaPC_mean = mean(mean_alpha_PC, na.rm=TRUE), alphaPC_sd = sd(mean_alpha_PC, na.rm=TRUE), alphaPC_se = alphaPC_sd / sqrt(n))


RFadultexpmeans_SubjectPhase_loco <- RFadultexpsubject_loco %>%
  group_by(SubjectID,Sex,Phase,Phase_raw,Group) %>%
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
                   
                   finalblock_loco_mean = mean(finalblock_loco, na.rm=TRUE), finalblock_loco_sd = sd(finalblock_loco,na.rm=TRUE), finalblock_loco_se = finalblock_loco_sd / sqrt(n),
                   finalblock_locoC_mean = mean(finalblock_loco_C, na.rm=TRUE), finalblock_locoC_sd = sd(finalblock_loco_C,na.rm=TRUE), finalblock_locoC_se = finalblock_locoC_sd / sqrt(n),
                   finalblock_locoPC_mean = mean(finalblock_loco_PC, na.rm=TRUE), finalblock_locoPC_sd = sd(finalblock_loco_PC,na.rm=TRUE), finalblock_locoPC_se = finalblock_locoPC_sd / sqrt(n))

### Export Means Tables -----
write.csv(RFadultexpmeans_Phase_lm, paste(path_analysis_RFtest_means,"RFadultexp_Means_Phase_RFvariables.csv",sep=""))
write.csv(RFadultexpmeans_Phase_loco, paste(path_analysis_RFtest_means,"RFadultexp_Means_Phase_Motorvariables.csv",sep=""))

write.csv(RFadultexpmeans_lm, paste(path_analysis_RFtest_means,"RFadultexp_Means_PhaseandSession_RFvariables.csv",sep=""))
write.csv(RFadultexpmeans_loco, paste(path_analysis_RFtest_means,"RFadultexp_Means_PhaseandSession_Motorvariables.csv",sep=""))


# STATISTICS -----

## Subject ns -----
RFadultexp__ns <- RFadultexpsubject_lm %>%
  group_by(Sex,Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID))
RFadultexp__ns


## M50 -----
### Exploratory Analysis: Sex -----
#### Mixed effects model with sex interaction and random intercept by subject
RFadultexp_M50_sexInt_model <- glmmTMB(mean_M50 ~ Phase*Group*Sex + (1|SubjectID), data = RFadultexpsubject_lm, family=stats::gaussian)
summary(RFadultexp_M50_sexInt_model)
Anova(RFadultexp_M50_sexInt_model, type=2)

##### Overall, there is an interaction of sex and phase but it's not group specific or consistent across phases.
##### As the interaction is lacking and there was not an a priori hypothesis of sex differences, sex is included as a
##### fixed effect in final models.

### Phase Analysis -----
#### Mixed effects model by phase and sex with random intercept by subject
RFadultexp_M50_model <- glmmTMB(mean_M50 ~ Phase*Group + Sex + (1|SubjectID), data = RFadultexpsubject_lm, family=stats::gaussian)
summary(RFadultexp_M50_model)

#### Compare full model to the null model
RFadultexp_M50_nullmodel <- glmmTMB(mean_M50 ~ 1 + Sex + (1|SubjectID), data = RFadultexpsubject_lm, family = stats::gaussian()) # Make the null model
anova(RFadultexp_M50_nullmodel, RFadultexp_M50_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultexp_M50_within_emmeans <-emmeans(RFadultexp_M50_model, pairwise ~ Phase|Group, adjust = 'none')
summary(RFadultexp_M50_within_emmeans)

RFadultexp_M50_between_emmeans <-emmeans(RFadultexp_M50_model, pairwise ~ Group|Phase, adjust = 'none')
summary(RFadultexp_M50_between_emmeans)

### Session Analysis -----
#### Baseline
RFadultexp_M50_BL_model <- glmmTMB(mean_M50 ~ Session*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_lm, Phase=="BL"),
                                     family = stats::gaussian())
summary(RFadultexp_M50_BL_model)

##### Compare full model to the null model
RFadultexp_M50_BL_nullmodel <- glmmTMB(mean_M50 ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_lm, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultexp_M50_BL_nullmodel, RFadultexp_M50_BL_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, there is no significant effect of day in either group in the baseline phase.

#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultexp_M50_INJ_model <- glmmTMB(mean_M50 ~ Session_collapsed*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_lm, Phase %in% c("BL","INJ")),
                                      family = stats::gaussian())
summary(RFadultexp_M50_INJ_model)

##### Compare full model to the null model
RFadultexp_M50_INJ_nullmodel <- glmmTMB(mean_M50 ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_lm, Phase %in% c("BL","INJ")), family = stats::gaussian()) # Make the null model
anova(RFadultexp_M50_INJ_nullmodel, RFadultexp_M50_INJ_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultexp_M50_INJ_within_emmeans <- emmeans(RFadultexp_M50_INJ_model, pairwise ~ Session_collapsed|Group, adjust = 'none')
summary(RFadultexp_M50_INJ_within_emmeans)

RFadultexp_M50_INJ_between_emmeans <- emmeans(RFadultexp_M50_INJ_model, pairwise ~ Group|Session_collapsed, adjust = 'none')
summary(RFadultexp_M50_INJ_between_emmeans)


#### Post-Injection
##### Mixed effects model by day with random intercept by subject
RFadultexp_M50_POST_model <- glmmTMB(mean_M50 ~ Session_collapsed*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_lm, Phase %in% c("BL","POST_EARLY","POST_LATE")),
                                        family = stats::gaussian())
summary(RFadultexp_M50_POST_model)

##### Compare full model to the null model
RFadultexp_M50_POST_nullmodel <- glmmTMB(mean_M50 ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_lm, Phase %in% c("BL","POST_EARLY","POST_LATE")), family = stats::gaussian()) # Make the null model
anova(RFadultexp_M50_POST_nullmodel, RFadultexp_M50_POST_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultexp_M50_POST_within_emmeans <- emmeans(RFadultexp_M50_POST_model, pairwise ~ Session_collapsed|Group, adjust = 'none')
summary(RFadultexp_M50_POST_within_emmeans)

RFadultexp_M50_POST_between_emmeans <- emmeans(RFadultexp_M50_POST_model, pairwise ~ Group|Session_collapsed, adjust = 'none')
summary(RFadultexp_M50_POST_between_emmeans)


### Percent Change Group Analysis -----
#### Mixed effects model by phase and sex with random intercept by subject
RFadultexp_M50_PC_model <- glmmTMB(mean_M50_PC ~ Phase*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_lm, Phase!="BL"), family=stats::gaussian)
summary(RFadultexp_M50_PC_model)

#### Compare full model to the null model
RFadultexp_M50_PC_nullmodel <- glmmTMB(mean_M50_PC ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_lm, Phase!="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultexp_M50_PC_nullmodel, RFadultexp_M50_PC_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultexp_M50_PC_within_emmeans <-emmeans(RFadultexp_M50_PC_model, pairwise ~ Phase|Group, adjust = 'none')
summary(RFadultexp_M50_PC_within_emmeans)

RFadultexp_M50_PC_between_emmeans <-emmeans(RFadultexp_M50_PC_model, pairwise ~ Group|Phase, adjust = 'none')
summary(RFadultexp_M50_PC_between_emmeans)


### Export Analysis Tables -----
#### Phase Tables
RFadultexp_M50_model_paramtable <- mixedmodel_paramtable(RFadultexp_M50_model)
RFadultexp_M50_model_fittable <- mixedmodel_fittable(RFadultexp_M50_model, RFadultexp_M50_nullmodel)
RFadultexp_M50_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_M50_within_emmeans,RFadultexp_M50_model)
RFadultexp_M50_between_emmeans_contraststable <- mixedmodel_emmeanstable_Phase(RFadultexp_M50_between_emmeans,RFadultexp_M50_model)

write.csv(RFadultexp_M50_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_M50_model_paramtable.csv",sep=""))
write.csv(RFadultexp_M50_model_fittable, paste(path_analysis_RFtest,"RFadultexp_M50_model_fittable.csv",sep=""))
write.csv(RFadultexp_M50_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_M50_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_M50_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_M50_emmeans_between_contraststable.csv",sep=""))

#### Session Tables
##### Baseline
RFadultexp_M50_BL_model_paramtable <- mixedmodel_paramtable(RFadultexp_M50_BL_model)
RFadultexp_M50_BL_model_fittable <- mixedmodel_fittable(RFadultexp_M50_BL_model, RFadultexp_M50_BL_nullmodel)

write.csv(RFadultexp_M50_BL_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_M50_Session_BL_model_paramtable.csv",sep=""))
write.csv(RFadultexp_M50_BL_model_fittable, paste(path_analysis_RFtest,"RFadultexp_M50_Session_BL_model_fittable.csv",sep=""))

##### Injection
RFadultexp_M50_INJ_model_paramtable <- mixedmodel_paramtable(RFadultexp_M50_INJ_model)
RFadultexp_M50_INJ_model_fittable <- mixedmodel_fittable(RFadultexp_M50_INJ_model, RFadultexp_M50_INJ_nullmodel)
RFadultexp_M50_INJ_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_M50_INJ_within_emmeans,RFadultexp_M50_INJ_model)
RFadultexp_M50_INJ_between_emmeans_contraststable <- mixedmodel_emmeanstable_Session(RFadultexp_M50_INJ_between_emmeans,RFadultexp_M50_INJ_model)

write.csv(RFadultexp_M50_INJ_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_M50_Session_INJ_model_paramtable.csv",sep=""))
write.csv(RFadultexp_M50_INJ_model_fittable, paste(path_analysis_RFtest,"RFadultexp_M50_Session_INJ_model_fittable.csv",sep=""))
write.csv(RFadultexp_M50_INJ_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_M50_Session_INJ_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_M50_INJ_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_M50_Session_INJ_emmeans_between_contraststable.csv",sep=""))


##### Post-Injection Early
RFadultexp_M50_POST_model_paramtable <- mixedmodel_paramtable(RFadultexp_M50_POST_model)
RFadultexp_M50_POST_model_fittable <- mixedmodel_fittable(RFadultexp_M50_POST_model, RFadultexp_M50_POST_nullmodel)
RFadultexp_M50_POST_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_M50_POST_within_emmeans,RFadultexp_M50_POST_model)
RFadultexp_M50_POST_between_emmeans_contraststable <- mixedmodel_emmeanstable_Session(RFadultexp_M50_POST_between_emmeans,RFadultexp_M50_POST_model)

write.csv(RFadultexp_M50_POST_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_M50_Session_POST_model_paramtable.csv",sep=""))
write.csv(RFadultexp_M50_POST_model_fittable, paste(path_analysis_RFtest,"RFadultexp_M50_Session_POST_model_fittable.csv",sep=""))
write.csv(RFadultexp_M50_POST_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_M50_Session_POST_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_M50_POST_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_M50_Session_POST_emmeans_between_contraststable.csv",sep=""))


#### Percent Change Group Tables
RFadultexp_M50_PC_model_paramtable <- mixedmodel_paramtable(RFadultexp_M50_PC_model)
RFadultexp_M50_PC_model_fittable <- mixedmodel_fittable(RFadultexp_M50_PC_model, RFadultexp_M50_PC_nullmodel)
RFadultexp_M50_PC_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_M50_PC_within_emmeans,RFadultexp_M50_PC_model)
RFadultexp_M50_PC_between_emmeans_contraststable <- mixedmodel_emmeanstable_Phase(RFadultexp_M50_PC_between_emmeans,RFadultexp_M50_PC_model)

write.csv(RFadultexp_M50_PC_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_M50_PC_model_paramtable.csv",sep=""))
write.csv(RFadultexp_M50_PC_model_fittable, paste(path_analysis_RFtest,"RFadultexp_M50_PC_model_fittable.csv",sep=""))
write.csv(RFadultexp_M50_PC_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_M50_PC_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_M50_PC_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_M50_PC_emmeans_between_contraststable.csv",sep=""))

## Block Locomotion -----
### Phase Analysis -----
#### Mixed effects model by phase and sex with random intercept by subject
RFadultexp_blockloco_model <- glmmTMB(finalblock_loco ~ Phase*Group + Sex + (1|SubjectID), data = RFadultexpsubject_loco, family=stats::gaussian)
summary(RFadultexp_blockloco_model)

#### Compare full model to the null model
RFadultexp_blockloco_nullmodel <- glmmTMB(finalblock_loco ~ 1 + (1|SubjectID), data = RFadultexpsubject_loco, family = stats::gaussian()) # Make the null model
anova(RFadultexp_blockloco_nullmodel, RFadultexp_blockloco_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultexp_blockloco_within_emmeans <-emmeans(RFadultexp_blockloco_model, pairwise ~ Phase|Group, adjust = 'none')
summary(RFadultexp_blockloco_within_emmeans)

RFadultexp_blockloco_between_emmeans <-emmeans(RFadultexp_blockloco_model, pairwise ~ Group|Phase, adjust = 'none')
summary(RFadultexp_blockloco_between_emmeans)


### Session Analysis -----
#### Baseline
RFadultexp_blockloco_BL_model <- glmmTMB(finalblock_loco ~ Session*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase=="BL"),
                                     family = stats::gaussian())
summary(RFadultexp_blockloco_BL_model)

##### Compare full model to the null model
RFadultexp_blockloco_BL_nullmodel <- glmmTMB(finalblock_loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultexp_blockloco_BL_nullmodel, RFadultexp_blockloco_BL_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, there is no significant effect of day in either group in the baseline phase.

#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultexp_blockloco_INJ_model <- glmmTMB(finalblock_loco ~ Session_collapsed*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","INJ")),
                                      family = stats::gaussian())
summary(RFadultexp_blockloco_INJ_model)

##### Compare full model to the null model
RFadultexp_blockloco_INJ_nullmodel <- glmmTMB(finalblock_loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","INJ")), family = stats::gaussian()) # Make the null model
anova(RFadultexp_blockloco_INJ_nullmodel, RFadultexp_blockloco_INJ_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultexp_blockloco_INJ_within_emmeans <- emmeans(RFadultexp_blockloco_INJ_model, pairwise ~ Session_collapsed|Group, adjust = 'none')
summary(RFadultexp_blockloco_INJ_within_emmeans)

RFadultexp_blockloco_INJ_between_emmeans <- emmeans(RFadultexp_blockloco_INJ_model, pairwise ~ Group|Session_collapsed, adjust = 'none')
summary(RFadultexp_blockloco_INJ_between_emmeans)


#### Post-Injection
##### Mixed effects model by day with random intercept by subject
RFadultexp_blockloco_POST_model <- glmmTMB(finalblock_loco ~ Session_collapsed*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","POST_EARLY","POST_LATE")),
                                       family = stats::gaussian())
summary(RFadultexp_blockloco_POST_model)

##### Compare full model to the null model
RFadultexp_blockloco_POST_nullmodel <- glmmTMB(finalblock_loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","POST_EARLY","POST_LATE")), family = stats::gaussian()) # Make the null model
anova(RFadultexp_blockloco_POST_nullmodel, RFadultexp_blockloco_POST_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultexp_blockloco_POST_within_emmeans <- emmeans(RFadultexp_blockloco_POST_model, pairwise ~ Session_collapsed|Group, adjust = 'none')
summary(RFadultexp_blockloco_POST_within_emmeans)

RFadultexp_blockloco_POST_between_emmeans <- emmeans(RFadultexp_blockloco_POST_model, pairwise ~ Group|Session_collapsed, adjust = 'none')
summary(RFadultexp_blockloco_POST_between_emmeans)

### Percent Change Group Analysis -----
#### Mixed effects model by phase and sex with random intercept by subject
RFadultexp_blockloco_PC_model <- glmmTMB(finalblock_loco_PC ~ Phase*Group + Sex, data = subset(RFadultexpsubject_loco, Phase!="BL"), family=stats::gaussian)
summary(RFadultexp_blockloco_PC_model)

#### Compare full model to the null model
RFadultexp_blockloco_PC_nullmodel <- glmmTMB(finalblock_loco_PC ~ 1 + Sex, data = subset(RFadultexpsubject_loco, Phase!="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultexp_blockloco_PC_nullmodel, RFadultexp_blockloco_PC_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultexp_blockloco_PC_within_emmeans <-emmeans(RFadultexp_blockloco_PC_model, pairwise ~ Phase|Group, adjust = 'none')
summary(RFadultexp_blockloco_PC_within_emmeans)

RFadultexp_blockloco_PC_between_emmeans <-emmeans(RFadultexp_blockloco_PC_model, pairwise ~ Group|Phase, adjust = 'none')
summary(RFadultexp_blockloco_PC_between_emmeans)


### Export Analysis Tables -----
#### Phase Tables
RFadultexp_blockloco_model_paramtable <- mixedmodel_paramtable(RFadultexp_blockloco_model)
RFadultexp_blockloco_model_fittable <- mixedmodel_fittable(RFadultexp_blockloco_model, RFadultexp_blockloco_nullmodel)
RFadultexp_blockloco_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_blockloco_within_emmeans,RFadultexp_blockloco_model)
RFadultexp_blockloco_between_emmeans_contraststable <- mixedmodel_emmeanstable_Phase(RFadultexp_blockloco_between_emmeans,RFadultexp_blockloco_model)

write.csv(RFadultexp_blockloco_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_blockloco_model_paramtable.csv",sep=""))
write.csv(RFadultexp_blockloco_model_fittable, paste(path_analysis_RFtest,"RFadultexp_blockloco_model_fittable.csv",sep=""))
write.csv(RFadultexp_blockloco_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_blockloco_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_blockloco_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_blockloco_emmeans_between_contraststable.csv",sep=""))

#### Session Tables
##### Baseline
RFadultexp_blockloco_BL_model_paramtable <- mixedmodel_paramtable(RFadultexp_blockloco_BL_model)
RFadultexp_blockloco_BL_model_fittable <- mixedmodel_fittable(RFadultexp_blockloco_BL_model, RFadultexp_blockloco_BL_nullmodel)

write.csv(RFadultexp_blockloco_BL_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_blockloco_Session_BL_model_paramtable.csv",sep=""))
write.csv(RFadultexp_blockloco_BL_model_fittable, paste(path_analysis_RFtest,"RFadultexp_blockloco_Session_BL_model_fittable.csv",sep=""))

##### Injection
RFadultexp_blockloco_INJ_model_paramtable <- mixedmodel_paramtable(RFadultexp_blockloco_INJ_model)
RFadultexp_blockloco_INJ_model_fittable <- mixedmodel_fittable(RFadultexp_blockloco_INJ_model, RFadultexp_blockloco_INJ_nullmodel)
RFadultexp_blockloco_INJ_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_blockloco_INJ_within_emmeans,RFadultexp_blockloco_INJ_model)
RFadultexp_blockloco_INJ_between_emmeans_contraststable <- mixedmodel_emmeanstable_Session(RFadultexp_blockloco_INJ_between_emmeans,RFadultexp_blockloco_INJ_model)

write.csv(RFadultexp_blockloco_INJ_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_blockloco_Session_INJ_model_paramtable.csv",sep=""))
write.csv(RFadultexp_blockloco_INJ_model_fittable, paste(path_analysis_RFtest,"RFadultexp_blockloco_Session_INJ_model_fittable.csv",sep=""))
write.csv(RFadultexp_blockloco_INJ_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_blockloco_Session_INJ_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_blockloco_INJ_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_blockloco_Session_INJ_emmeans_between_contraststable.csv",sep=""))


##### Post-Injection Early
RFadultexp_blockloco_POST_model_paramtable <- mixedmodel_paramtable(RFadultexp_blockloco_POST_model)
RFadultexp_blockloco_POST_model_fittable <- mixedmodel_fittable(RFadultexp_blockloco_POST_model, RFadultexp_blockloco_POST_nullmodel)
RFadultexp_blockloco_POST_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_blockloco_POST_within_emmeans,RFadultexp_blockloco_POST_model)
RFadultexp_blockloco_POST_between_emmeans_contraststable <- mixedmodel_emmeanstable_Session(RFadultexp_blockloco_POST_between_emmeans,RFadultexp_blockloco_POST_model)

write.csv(RFadultexp_blockloco_POST_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_blockloco_Session_POST_model_paramtable.csv",sep=""))
write.csv(RFadultexp_blockloco_POST_model_fittable, paste(path_analysis_RFtest,"RFadultexp_blockloco_Session_POST_model_fittable.csv",sep=""))
write.csv(RFadultexp_blockloco_POST_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_blockloco_Session_POST_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_blockloco_POST_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_blockloco_Session_POST_emmeans_between_contraststable.csv",sep=""))


#### Percent Change Group Tables
RFadultexp_blockloco_PC_model_paramtable <- mixedmodel_paramtable(RFadultexp_blockloco_PC_model)
#RFadultexp_blockloco_PC_model_fittable <- mixedmodel_fittable(RFadultexp_blockloco_PC_model, RFadultexp_blockloco_PC_nullmodel)
RFadultexp_blockloco_PC_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_blockloco_PC_within_emmeans,RFadultexp_blockloco_PC_model)
RFadultexp_blockloco_PC_between_emmeans_contraststable <- mixedmodel_emmeanstable_Phase(RFadultexp_blockloco_PC_between_emmeans,RFadultexp_blockloco_PC_model)

write.csv(RFadultexp_blockloco_PC_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_blockloco_PC_model_paramtable.csv",sep=""))
#write.csv(RFadultexp_blockloco_PC_model_fittable, paste(path_analysis_RFtest,"RFadultexp_blockloco_PC_model_fittable.csv",sep=""))
write.csv(RFadultexp_blockloco_PC_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_blockloco_PC_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_blockloco_PC_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_blockloco_PC_emmeans_between_contraststable.csv",sep=""))



# FIGURES -----
## Example Subject Plots -----
### Prep data
examplesubjectdata <- RFadultexpdata_raw %>% dplyr::filter(SubjectID == 538) %>% 
  pivot_longer(cols = all_of(freqcolnames), names_to = 'Frequency', values_to='LeverPresses')

examplesubjectdata$Pass <- factor(examplesubjectdata$Pass)

# Prep plotting variables
xorder <- c('28','32','35','40','45','50','56','63','71','79','89','100','112','126','141')
#xlabels <- c('28','32','35','40','45','50','56','63','71','79','89','100','112','126','141')
xlabels <- c('28','','35','','45','','56','','71','','89','','112','','141')

RFexample_ymax <- 60
RFexample_ymin <- 0

# Plot Reduced RF Example - Baseline
MBLplot <- ggplot(subset(examplesubjectdata, UniqueSessionID == "BL_3"), aes(x = Frequency, y = LeverPresses, group = Pass)) +
  geom_line(linewidth = mylinewidth, aes(color = Pass), show.legend = FALSE) + scale_x_discrete(limits = xorder, labels = xlabels) + 
  geom_point(aes(x = Frequency, y = LeverPresses, color = Pass), size = mypointsize, show.legend = FALSE) +
  scale_color_manual(values=c("#8d8d8d", "#0000FF", "#950EF2", "#FF0080", "#0FC127")) +  # Old values: values=c("#DC8B06", "#0E6BF2", "#950EF2", "#0FC127", "#D4212F")
  scale_y_continuous(limits = c(RFexample_ymin, RFexample_ymax), breaks=seq(RFexample_ymin,RFexample_ymax,20), expand=c(0,0)) +
  mytheme +
  ylab("Lever Presses") +
  ggtitle("Baseline")

MBLplot

# Plot Reduced RF Example - Morphine Administration
MAplot <- ggplot(subset(examplesubjectdata, UniqueSessionID == "INJ_6"), aes(x = Frequency, y = LeverPresses, group = Pass)) +
  geom_line(linewidth = mylinewidth, aes(color = Pass), show.legend = FALSE) + scale_x_discrete(limits = xorder, labels = xlabels) + 
  geom_point(aes(x = Frequency, y = LeverPresses, color = Pass), size = mypointsize, show.legend = FALSE) +
  scale_color_manual(values=c("#8d8d8d", "#0000FF", "#950EF2", "#FF0080", "#0FC127")) +  # Old values: values=c("#DC8B06", "#0E6BF2", "#950EF2", "#0FC127", "#D4212F")
  scale_y_continuous(limits = c(RFexample_ymin, RFexample_ymax), breaks=seq(RFexample_ymin,RFexample_ymax,20), expand=c(0,0)) +
  mytheme +
  ylab("Lever Presses") +
  ggtitle("Morphine (Session 7)")

MAplot

# Plot Reduced RF Example - Morphine Withdrawal
MWplot <- ggplot(subset(examplesubjectdata, UniqueSessionID == "POST_2"), aes(x = Frequency, y = LeverPresses, group = Pass)) +
  geom_line(linewidth = mylinewidth, aes(color = Pass), show.legend = FALSE) + scale_x_discrete(limits = xorder, labels = xlabels) + 
  geom_point(aes(x = Frequency, y = LeverPresses, color = Pass), size = mypointsize, show.legend = FALSE) +
  scale_color_manual(values=c("#8d8d8d", "#0000FF", "#950EF2", "#FF0080", "#0FC127")) +  # Old values: values=c("#DC8B06", "#0E6BF2", "#950EF2", "#0FC127", "#D4212F")
  scale_y_continuous(limits = c(RFexample_ymin, RFexample_ymax), breaks=seq(RFexample_ymin,RFexample_ymax,20), expand=c(0,0)) +
  mytheme +
  ylab("Lever Presses") +
  ggtitle("Withdrawal (Session 9)")

MWplot

# Make empty tile for key
keyplot <- ggplot() + mytheme

# Create full panel
exampleRFplots <- ggarrange(MBLplot,MAplot, MWplot,keyplot, ncol=4, nrow=1, widths=c(1,1,1,.8), labels = c("A", "B", "C", ""))
exampleRFplots



## Figure 5 Panels -----
### Prepare plot colors and variables -----
groupcolors <- c("AdolS" = AdolSAdultcolor, "AdolM" = AdolMAdultcolor)
AdolSMphasecolors <- c("BL" = ICSSBLcolor, "INJ"=AdolSAdultcolor,"POST_EARLY"=AdolSAdultEPcolor,"POST_LATE"=AdolSAdultLPcolor)
AdolMMphasecolors <- c("BL" = ICSSBLcolor, "INJ"=AdolMAdultcolor,"POST_EARLY"=AdolMAdultEPcolor,"POST_LATE"=AdolMAdultLPcolor)

sexshapes <- c("M" = 0, "F" = 1)
# Set up x axis ticks and labels - label day (every other session)
sessionticks <- c(5:20)
#sessionlabels <- c('5','6','7','8','9','10','11','12','13','14','15','16','17','18','19')
sessionlabels <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16')

injsession <- c(1,2,3,4,5,6,7)

### Prepare plot colors and variables -----
PCphasecolors <- c("AdolSINJ"=AdolSAdultcolor,"AdolMINJ"=AdolMAdultcolor,"AdolSPOST_EARLY"=AdolSAdultEPcolor,"AdolMPOST_EARLY"=AdolMAdultEPcolor,"AdolSPOST_LATE"=AdolSAdultLPcolor,"AdolMPOST_LATE"=AdolMAdultLPcolor)

# Add fill as a variable to means
RFadultexpmeans_Phase_lm$PCfill <- paste(RFadultexpmeans_Phase_lm$Group,RFadultexpmeans_Phase_lm$Phase,sep="")
RFadultexpmeans_SubjectPhase_lm$PCfill <- paste(RFadultexpmeans_SubjectPhase_lm$Group,RFadultexpmeans_SubjectPhase_lm$Phase,sep="")

RFadultexpmeans_Phase_loco$PCfill <- paste(RFadultexpmeans_Phase_loco$Group,RFadultexpmeans_Phase_loco$Phase,sep="")
RFadultexpmeans_SubjectPhase_loco$PCfill <- paste(RFadultexpmeans_SubjectPhase_loco$Group,RFadultexpmeans_SubjectPhase_loco$Phase,sep="")



### M50 -----
ymax_M50C <- 9
ymin_M50C <- -38

ymax_M50_bar <- 143
ymin_M50_bar <- 0

injmarker_M50 <- c(rep(ymin_M50C,7))
injticks_M50 <- data.frame(injsession,injmarker_M50)

# M50 Change Line
M50C_line <- ggplot(subset(RFadultexpmeans_lm,Phase!="BL"), aes(x=sessionlabel, y=M50C_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadultexpmeans_lm,Phase!="BL"), aes(ymin=M50C_mean-M50C_se, ymax=M50C_mean+M50C_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name ="Session", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_M50C, ymax_M50C), breaks = seq(-36, ymax_M50C, 12)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("M50", sep = "")) +
  geom_point(data=injticks_M50, aes(x=injsession, y=injmarker_M50, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE)+ # Add injection ticks
  annotate("text", x = 1, y = -20.5, label = "*", size = mystarsize) + # Add asterisks for morphine group session significance
  annotate("text", x = 2, y = -25, label = "**", size = mystarsize)+
  annotate("text", x = 3, y = -19.5, label = "*", size = mystarsize)+
  annotate("text", x = 4, y = -25, label = "***", size = mystarsize)+
  annotate("text", x = 5, y = -29, label = "***", size = mystarsize)+
  annotate("text", x = 6, y = -23, label = "***", size = mystarsize)+
  annotate("text", x = 7, y = -27, label = "***", size = mystarsize)+
  annotate("text", x = 8, y = -16.5, label = "#", size = 2)+ # Add hashtags for control group session significance
  annotate("text", x = 9, y = -19, label = "##", size = 2)+
  annotate("text", x = 2, y = -27, label = "###", size = 2)+
  annotate("text", x = 5, y = -31, label = "&&", size = 2)+ # Add & sign for significance between control and morphine groups
  annotate("text", x = 8, y = -19.5, label = "&", size = 2)+
  annotate("text", x = 9, y = -22, label = "&&", size = 2)

M50C_line

# Saline Group Intercept by Phase Bar
M50_AdolS_bar <- ggplot(subset(RFadultexpmeans_Phase_lm,Group=="AdolS"), aes(x=Phase, y=M50_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadultexpmeans_SubjectPhase_lm, Group=="AdolS"), mapping = aes(x = Phase, y = M50_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadultexpmeans_Phase_lm,Group=="AdolS"), aes(ymin=M50_mean-M50_se, ymax=M50_mean+M50_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_M50_bar, ymax_M50_bar), breaks = seq(ymin_M50_bar, ymax_M50_bar, 40)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = AdolSMphasecolors) +
  scale_fill_manual("legend", values = AdolSMphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("M50") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Adol Control", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(1.9), y.position=c(119) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2) +
  geom_bracket(xmin=c(1), xmax = c(2.9), y.position=c(133) ,label.size = 2, label = "##",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2) +
  geom_bracket(xmin=c(2.1), xmax = c(4), y.position=c(119) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2) +
  geom_bracket(xmin=c(3.1), xmax = c(4), y.position=c(133) ,label.size = 2, label = "#",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

M50_AdolS_bar

# Morphine Group Intercept by Phase Bar
M50_AdolM_bar <- ggplot(subset(RFadultexpmeans_Phase_lm,Group=="AdolM"), aes(x=Phase, y=M50_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadultexpmeans_SubjectPhase_lm, Group=="AdolM"), mapping = aes(x = Phase, y = M50_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadultexpmeans_Phase_lm,Group=="AdolM"), aes(ymin=M50_mean-M50_se, ymax=M50_mean+M50_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_M50_bar, ymax_M50_bar), breaks = seq(ymin_M50_bar, ymax_M50_bar, 40)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = AdolMMphasecolors) +
  scale_fill_manual("legend", values = AdolMMphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("M50") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Adol Morphine", sep = ""))+
  geom_bracket(xmin=c(1), xmax = c(1.9), y.position=c(125) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.1), xmax = c(2.9), y.position=c(125) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3.1), xmax = c(4), y.position=c(125) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

M50_AdolM_bar

# Create full panel
M50 <- ggarrange(M50_AdolS_bar,M50_AdolM_bar, M50C_line, ncol=3, nrow=1, widths=c(1.8,1.8,4), labels = c("D", "E", "F"))
M50
ggsave(plot=M50, path=path_figures, device="pdf", filename="F5_ReExposureAdultRF_MorphineTest_M50.pdf", width=RFpanelwidth, height=RFpanelheight)



### Final Block Locomotion -----
ymax_blocklocoC <- 26
ymin_blocklocoC <- -26

ymax_blockloco_bar <- 85
ymin_blockloco_bar <- 0

injmarker_blockloco <- c(rep(ymin_blocklocoC,7))
injticks_blockloco <- data.frame(injsession,injmarker_blockloco)

# blockloco Change Line
blockloco_line <- ggplot(subset(RFadultexpmeans_loco,Phase!="BL"), aes(x=sessionlabel, y=finalblock_locoC_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadultexpmeans_loco,Phase!="BL"), aes(ymin=finalblock_locoC_mean-finalblock_locoC_se, ymax=finalblock_locoC_mean+finalblock_locoC_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name ="Session", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blocklocoC, ymax_blocklocoC), breaks = seq(-24, ymax_blocklocoC, 12)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("Last Block Locomotion", sep = "")) +
  geom_point(data=injticks_blockloco, aes(x=injsession, y=injmarker_blockloco, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE)+ # Add injection ticks
  annotate("text", x = 5, y = 18, label = "*", size = mystarsize) + # Add asterisks for morphine group session significance
  annotate("text", x = 6, y = 23, label = "*", size = mystarsize)+
  annotate("text", x = 7, y = 18.5, label = "*", size = mystarsize) + # Add asterisks for morphine group session significance
  annotate("text", x = 8, y = -15.5, label = "**", size = mystarsize)+
  annotate("text", x = 9, y = -12, label = "*", size = mystarsize)+
  annotate("text", x = 11, y = -15, label = "*", size = mystarsize)+
  annotate("text", x = 1, y = -19, label = "#", size = 2)+ # Add hashtags for control group session significance
  annotate("text", x = 11, y = -17, label = "#", size = 2)+
  annotate("text", x = 13, y = -12.5, label = "#", size = 2)

blockloco_line

# Adolescent Saline Group loco by Phase Bar
blockloco_AdolS_bar <- ggplot(subset(RFadultexpmeans_Phase_loco,Group=="AdolS"), aes(x=Phase, y=finalblock_loco_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadultexpmeans_SubjectPhase_loco, Group=="AdolS"), mapping = aes(x = Phase, y = finalblock_loco_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadultexpmeans_Phase_loco,Group=="AdolS"), aes(ymin=finalblock_loco_mean-finalblock_loco_se, ymax=finalblock_loco_mean+finalblock_loco_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blockloco_bar, ymax_blockloco_bar), breaks = seq(ymin_blockloco_bar, ymax_blockloco_bar, 25)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = AdolSMphasecolors) +
  scale_fill_manual("legend", values = AdolSMphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Block Locomotion") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Adol Control", sep = ""))

blockloco_AdolS_bar

# Adolescent Morphine Group loco by Phase Bar
blockloco_AdolM_bar <- ggplot(subset(RFadultexpmeans_Phase_loco,Group=="AdolM"), aes(x=Phase, y=finalblock_loco_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadultexpmeans_SubjectPhase_loco, Group=="AdolM"), mapping = aes(x = Phase, y = finalblock_loco_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadultexpmeans_Phase_loco,Group=="AdolM"), aes(ymin=finalblock_loco_mean-finalblock_loco_se, ymax=finalblock_loco_mean+finalblock_loco_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blockloco_bar, ymax_blockloco_bar), breaks = seq(ymin_blockloco_bar, ymax_blockloco_bar, 25)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = AdolMMphasecolors) +
  scale_fill_manual("legend", values = AdolMMphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Block Locomotion") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Adol Morphine", sep = ""))+
  geom_bracket(xmin=c(1), xmax = c(1.9), y.position=c(50) ,label.size = mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.1), xmax = c(3), y.position=c(50) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.1), xmax = c(4), y.position=c(62) ,label.size = mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

blockloco_AdolM_bar

# Create full panel
blockloco <- ggarrange(blockloco_AdolS_bar,blockloco_AdolM_bar,blockloco_line, ncol=3, nrow=1, widths=c(1.8,1.8,4), labels = c("G","H","I"))
blockloco

### M50 Percent Change -----
ymax_M50PC <- 10
ymin_M50PC <- -22

M50_PC_bar <- ggplot(subset(RFadultexpmeans_Phase_lm,Phase!="BL"), aes(x=Phase, y=M50PC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadultexpmeans_Phase_lm,Phase!="BL"), aes(ymin=M50PC_mean-M50PC_se, ymax=M50PC_mean+M50PC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  #geom_point(subset(PCmeans_lm_subject, Phase!="BL"), mapping = aes(x = Phase, y = M50PC_mean), position = position_dodge(width = 0.95), size = mysubjectpointsize, stroke = mypointstroke, shape = 1, color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_M50PC, ymax_M50PC), breaks = seq(-20, ymax_M50PC, 10)) + 
  scale_x_discrete(expand = c(0.25, 0), name ="Phase", labels=c('INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change from BL") + xlab('Phase') +
  mytheme +
  ggtitle(paste("M50", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)+
  geom_bracket(xmin=c(0.83), xmax = c(1.17), y.position=c(-19) ,label.size = 2, label = "&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(1.83), xmax = c(2.17), y.position=c(-12) ,label.size = 2, label = "&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

M50_PC_bar


### Final Block Locomotion Percent Change -----
ymax_blocklocoPC <- 120
ymin_blocklocoPC <- -60

blockloco_PC_bar <- ggplot(subset(RFadultexpmeans_Phase_loco,Phase!="BL"), aes(x=Phase, y=finalblock_locoPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadultexpmeans_Phase_loco,Phase!="BL"), aes(ymin=finalblock_locoPC_mean-finalblock_locoPC_se, ymax=finalblock_locoPC_mean+finalblock_locoPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  #geom_point(subset(PCmeans_lm_subject, Phase!="BL"), mapping = aes(x = Phase, y = blocklocoPC_mean), position = position_dodge(width = 0.95), size = mysubjectpointsize, stroke = mypointstroke, shape = 1, color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blocklocoPC, ymax_blocklocoPC), breaks = seq(-60, ymax_blocklocoPC, 60)) + 
  scale_x_discrete(expand = c(0.25, 0), name ="Phase", labels=c('INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change from BL") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Last Block Locomotion", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)+
  geom_bracket(xmin=c(0.83), xmax = c(1.17), y.position=c(112) ,label.size = 2, label = "&&&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

blockloco_PC_bar

# Create full panel
M50locoPC <- ggarrange(M50_PC_bar,blockloco_PC_bar, ncol=2, nrow=1, widths=c(1,1), labels = c("J","K"))
M50locoPC



## Create Figure 5 -----
# Create full panel
figure5 <- ggarrange(exampleRFplots,M50,blockloco,M50locoPC, ncol=1, nrow=4)
figure5
ggsave(plot=figure5, path=path_figures, device="pdf", filename="F5_AllPanels.pdf", width=RFpanelwidth, height=RFpanelheight*4)




# CITE PACKAGES -----
library(grateful)
cite_packages(path = paste(computerusergithubpath,"Analysis-PBAdolescentMorphineReward/Functions/Functions_LoadPackages.R", sep=""),
              out.dir = path_analysis, out.format = c("docx")
)





# SUPPLEMENTAL -----

## STATISTICS -----
## Max Lever Presses -----
### Phase Analysis -----
#### Mixed effects model by phase and sex with random intercept by subject
RFadultexp_maxLP_model <- glmmTMB(pass_maxLP ~ Phase*Group + Sex + (1|SubjectID), data = RFadultexpsubject_loco, family=stats::gaussian)
summary(RFadultexp_maxLP_model)

#### Compare full model to the null model
RFadultexp_maxLP_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = RFadultexpsubject_loco, family = stats::gaussian()) # Make the null model
anova(RFadultexp_maxLP_nullmodel, RFadultexp_maxLP_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultexp_maxLP_within_emmeans <-emmeans(RFadultexp_maxLP_model, pairwise ~ Phase|Group, adjust = 'none')
summary(RFadultexp_maxLP_within_emmeans)

RFadultexp_maxLP_between_emmeans <-emmeans(RFadultexp_maxLP_model, pairwise ~ Group|Phase, adjust = 'none')
summary(RFadultexp_maxLP_between_emmeans)



### Session Analysis -----
#### Baseline
RFadultexp_maxLP_BL_model <- glmmTMB(pass_maxLP ~ Session*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase=="BL"),
                                     family = stats::gaussian())
summary(RFadultexp_maxLP_BL_model)

##### Compare full model to the null model
RFadultexp_maxLP_BL_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultexp_maxLP_BL_nullmodel, RFadultexp_maxLP_BL_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, there is no significant effect of day in either group in the baseline phase.

#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultexp_maxLP_INJ_model <- glmmTMB(pass_maxLP ~ Session_collapsed*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","INJ")),
                                      family = stats::gaussian())
summary(RFadultexp_maxLP_INJ_model)

##### Compare full model to the null model
RFadultexp_maxLP_INJ_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","INJ")), family = stats::gaussian()) # Make the null model
anova(RFadultexp_maxLP_INJ_nullmodel, RFadultexp_maxLP_INJ_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultexp_maxLP_INJ_within_emmeans <- emmeans(RFadultexp_maxLP_INJ_model, pairwise ~ Session_collapsed|Group, adjust = 'none')
summary(RFadultexp_maxLP_INJ_within_emmeans)

RFadultexp_maxLP_INJ_between_emmeans <- emmeans(RFadultexp_maxLP_INJ_model, pairwise ~ Group|Session_collapsed, adjust = 'none')
summary(RFadultexp_maxLP_INJ_between_emmeans)


#### Post-Injection
##### Mixed effects model by day with random intercept by subject
RFadultexp_maxLP_POST_model <- glmmTMB(pass_maxLP ~ Session_collapsed*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","POST_EARLY","POST_LATE")),
                                       family = stats::gaussian())
summary(RFadultexp_maxLP_POST_model)

##### Compare full model to the null model
RFadultexp_maxLP_POST_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","POST_EARLY","POST_LATE")), family = stats::gaussian()) # Make the null model
anova(RFadultexp_maxLP_POST_nullmodel, RFadultexp_maxLP_POST_model, test = "LRT") # Full model is a not better fit than the null mode

### Percent Change Group Analysis -----
#### Mixed effects model by phase and sex with random intercept by subject
RFadultexp_maxLP_PC_model <- glmmTMB(pass_maxLP_PC ~ Phase*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase!="BL"), family=stats::gaussian)
summary(RFadultexp_maxLP_PC_model)

#### Compare full model to the null model
RFadultexp_maxLP_PC_nullmodel <- glmmTMB(pass_maxLP_PC ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase!="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultexp_maxLP_PC_nullmodel, RFadultexp_maxLP_PC_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultexp_maxLP_PC_within_emmeans <-emmeans(RFadultexp_maxLP_PC_model, pairwise ~ Phase|Group, adjust = 'none')
summary(RFadultexp_maxLP_PC_within_emmeans)

RFadultexp_maxLP_PC_between_emmeans <-emmeans(RFadultexp_maxLP_PC_model, pairwise ~ Group|Phase, adjust = 'none')
summary(RFadultexp_maxLP_PC_between_emmeans)


### Export Analysis Tables -----
#### Phase Tables
RFadultexp_maxLP_model_paramtable <- mixedmodel_paramtable(RFadultexp_maxLP_model)
RFadultexp_maxLP_model_fittable <- mixedmodel_fittable(RFadultexp_maxLP_model, RFadultexp_maxLP_nullmodel)
RFadultexp_maxLP_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_maxLP_within_emmeans,RFadultexp_maxLP_model)
RFadultexp_maxLP_between_emmeans_contraststable <- mixedmodel_emmeanstable_Phase(RFadultexp_maxLP_between_emmeans,RFadultexp_maxLP_model)

write.csv(RFadultexp_maxLP_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_maxLP_model_paramtable.csv",sep=""))
write.csv(RFadultexp_maxLP_model_fittable, paste(path_analysis_RFtest,"RFadultexp_maxLP_model_fittable.csv",sep=""))
write.csv(RFadultexp_maxLP_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_maxLP_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_maxLP_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_maxLP_emmeans_between_contraststable.csv",sep=""))

#### Session Tables
##### Baseline
RFadultexp_maxLP_BL_model_paramtable <- mixedmodel_paramtable(RFadultexp_maxLP_BL_model)
RFadultexp_maxLP_BL_model_fittable <- mixedmodel_fittable(RFadultexp_maxLP_BL_model, RFadultexp_maxLP_BL_nullmodel)

write.csv(RFadultexp_maxLP_BL_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_maxLP_Session_BL_model_paramtable.csv",sep=""))
write.csv(RFadultexp_maxLP_BL_model_fittable, paste(path_analysis_RFtest,"RFadultexp_maxLP_Session_BL_model_fittable.csv",sep=""))

##### Injection
RFadultexp_maxLP_INJ_model_paramtable <- mixedmodel_paramtable(RFadultexp_maxLP_INJ_model)
RFadultexp_maxLP_INJ_model_fittable <- mixedmodel_fittable(RFadultexp_maxLP_INJ_model, RFadultexp_maxLP_INJ_nullmodel)
RFadultexp_maxLP_INJ_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_maxLP_INJ_within_emmeans,RFadultexp_maxLP_INJ_model)
RFadultexp_maxLP_INJ_between_emmeans_contraststable <- mixedmodel_emmeanstable_Session(RFadultexp_maxLP_INJ_between_emmeans,RFadultexp_maxLP_INJ_model)

write.csv(RFadultexp_maxLP_INJ_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_maxLP_Session_INJ_model_paramtable.csv",sep=""))
write.csv(RFadultexp_maxLP_INJ_model_fittable, paste(path_analysis_RFtest,"RFadultexp_maxLP_Session_INJ_model_fittable.csv",sep=""))
write.csv(RFadultexp_maxLP_INJ_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_maxLP_Session_INJ_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_maxLP_INJ_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_maxLP_Session_INJ_emmeans_between_contraststable.csv",sep=""))


##### Post-Injection
RFadultexp_maxLP_POST_model_paramtable <- mixedmodel_paramtable(RFadultexp_maxLP_POST_model)
RFadultexp_maxLP_POST_model_fittable <- mixedmodel_fittable(RFadultexp_maxLP_POST_model, RFadultexp_maxLP_POST_nullmodel)

write.csv(RFadultexp_maxLP_POST_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_maxLP_Session_POST_model_paramtable.csv",sep=""))
write.csv(RFadultexp_maxLP_POST_model_fittable, paste(path_analysis_RFtest,"RFadultexp_maxLP_Session_POST_model_fittable.csv",sep=""))


#### Percent Change Group Tables
RFadultexp_maxLP_PC_model_paramtable <- mixedmodel_paramtable(RFadultexp_maxLP_PC_model)
RFadultexp_maxLP_PC_model_fittable <- mixedmodel_fittable(RFadultexp_maxLP_PC_model, RFadultexp_maxLP_PC_nullmodel)
RFadultexp_maxLP_PC_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_maxLP_PC_within_emmeans,RFadultexp_maxLP_PC_model)
RFadultexp_maxLP_PC_between_emmeans_contraststable <- mixedmodel_emmeanstable_Phase(RFadultexp_maxLP_PC_between_emmeans,RFadultexp_maxLP_PC_model)

write.csv(RFadultexp_maxLP_PC_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_maxLP_PC_model_paramtable.csv",sep=""))
write.csv(RFadultexp_maxLP_PC_model_fittable, paste(path_analysis_RFtest,"RFadultexp_maxLP_PC_model_fittable.csv",sep=""))
write.csv(RFadultexp_maxLP_PC_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_maxLP_PC_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_maxLP_PC_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_maxLP_PC_emmeans_between_contraststable.csv",sep=""))


## Total Session Lever Presses -----
### Scale the data -----
#### NOTE: session_sumLP must be scaled to prevent convergence issues in the later within phase day models. 
#### All models are fit to the scaled and centered sumLP data, and tables will be back-transformed for interpretation at the end.
RFadultexpsubject_loco$session_sumLP_scaled <- scale(RFadultexpsubject_loco$session_sumLP)

# Save scaling parameters for back-transformation
mean_session_sumLP <- mean(RFadultexpdata_loco$session_sumLP)
sd_session_sumLP <- sd(RFadultexpdata_loco$session_sumLP)


### Phase Analysis -----
#### Mixed effects model by phase and sex with random intercept by subject
RFadultexp_sumLP_model <- glmmTMB(session_sumLP_scaled ~ Phase*Group + Sex + (1|SubjectID), data = RFadultexpsubject_loco, family=stats::gaussian)
summary(RFadultexp_sumLP_model)

# Convert random-effect variance and SD back to raw scale
sumLP_var_scaled <- 0.5513
sumLP_sd_scaled  <- sqrt(sumLP_var_scaled)

sumLP_var_raw <- sumLP_var_scaled * (sd_session_sumLP^2)
sumLP_sd_raw  <- sumLP_sd_scaled  * sd_session_sumLP

sumLP_var_raw
sumLP_sd_raw

#### Compare full model to the null model
RFadultexp_sumLP_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = RFadultexpsubject_loco, family = stats::gaussian()) # Make the null model
anova(RFadultexp_sumLP_nullmodel, RFadultexp_sumLP_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultexp_sumLP_within_emmeans <-emmeans(RFadultexp_sumLP_model, pairwise ~ Phase|Group, adjust = 'none')
summary(RFadultexp_sumLP_within_emmeans)

RFadultexp_sumLP_between_emmeans <-emmeans(RFadultexp_sumLP_model, pairwise ~ Group|Phase, adjust = 'none')
summary(RFadultexp_sumLP_between_emmeans)

### Session Analysis -----
#### Baseline
RFadultexp_sumLP_BL_model <- glmmTMB(session_sumLP_scaled ~ Session*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase=="BL"),
                                     family = stats::gaussian())
summary(RFadultexp_sumLP_BL_model)


# Convert random-effect variance and SD back to raw scale
sumLP_BL_var_scaled <- 0.53718
sumLP_BL_sd_scaled  <- sqrt(sumLP_BL_var_scaled)

sumLP_BL_var_raw <- sumLP_BL_var_scaled * (sd_session_sumLP^2)
sumLP_BL_sd_raw  <- sumLP_BL_sd_scaled  * sd_session_sumLP

sumLP_BL_var_raw
sumLP_BL_sd_raw


##### Compare full model to the null model
RFadultexp_sumLP_BL_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultexp_sumLP_BL_nullmodel, RFadultexp_sumLP_BL_model, test = "LRT") # Full model is a better fit than the null mode
###### Overall, there is no significant effect of day in either group in the baseline phase.

#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultexp_sumLP_INJ_model <- glmmTMB(session_sumLP_scaled ~ Session_collapsed*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","INJ")),
                                      family = stats::gaussian())
summary(RFadultexp_sumLP_INJ_model)

# Convert random-effect variance and SD back to raw scale
sumLP_INJ_var_scaled <- 0.8084
sumLP_INJ_sd_scaled  <- sqrt(sumLP_INJ_var_scaled)

sumLP_INJ_var_raw <- sumLP_INJ_var_scaled * (sd_session_sumLP^2)
sumLP_INJ_sd_raw  <- sumLP_INJ_sd_scaled  * sd_session_sumLP

sumLP_INJ_var_raw
sumLP_INJ_sd_raw

##### Compare full model to the null model
RFadultexp_sumLP_INJ_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","INJ")), family = stats::gaussian()) # Make the null model
anova(RFadultexp_sumLP_INJ_nullmodel, RFadultexp_sumLP_INJ_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultexp_sumLP_INJ_within_emmeans <- emmeans(RFadultexp_sumLP_INJ_model, pairwise ~ Session_collapsed|Group, adjust = 'none')
summary(RFadultexp_sumLP_INJ_within_emmeans)

RFadultexp_sumLP_INJ_between_emmeans <- emmeans(RFadultexp_sumLP_INJ_model, pairwise ~ Group|Session_collapsed, adjust = 'none')
summary(RFadultexp_sumLP_INJ_between_emmeans)


#### Post-Injection
##### Mixed effects model by day with random intercept by subject
RFadultexp_sumLP_POST_model <- glmmTMB(session_sumLP_scaled ~ Session_collapsed*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","POST_EARLY","POST_LATE")),
                                       family = stats::gaussian())
summary(RFadultexp_sumLP_POST_model)

# Convert random-effect variance and SD back to raw scale
sumLP_POST_var_scaled <- 0.3700
sumLP_POST_sd_scaled  <- sqrt(sumLP_POST_var_scaled)

sumLP_POST_var_raw <- sumLP_POST_var_scaled * (sd_session_sumLP^2)
sumLP_POST_sd_raw  <- sumLP_POST_sd_scaled  * sd_session_sumLP

sumLP_POST_var_raw
sumLP_POST_sd_raw

##### Compare full model to the null model
RFadultexp_sumLP_POST_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase %in% c("BL","POST_EARLY","POST_LATE")), family = stats::gaussian()) # Make the null model
anova(RFadultexp_sumLP_POST_nullmodel, RFadultexp_sumLP_POST_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultexp_sumLP_POST_within_emmeans <- emmeans(RFadultexp_sumLP_POST_model, pairwise ~ Session_collapsed|Group, adjust = 'none')
summary(RFadultexp_sumLP_POST_within_emmeans)

RFadultexp_sumLP_POST_between_emmeans <- emmeans(RFadultexp_sumLP_POST_model, pairwise ~ Group|Session_collapsed, adjust = 'none')
summary(RFadultexp_sumLP_POST_between_emmeans)



### Percent Change Group Analysis -----
#### Mixed effects model by phase and sex with random intercept by subject
RFadultexp_sumLP_PC_model <- glmmTMB(session_sumLP_PC ~ Phase*Group + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase!="BL"), family=stats::gaussian)
summary(RFadultexp_sumLP_PC_model)

#### Compare full model to the null model
RFadultexp_sumLP_PC_nullmodel <- glmmTMB(session_sumLP_PC ~ 1 + Sex + (1|SubjectID), data = subset(RFadultexpsubject_loco, Phase!="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultexp_sumLP_PC_nullmodel, RFadultexp_sumLP_PC_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultexp_sumLP_PC_within_emmeans <-emmeans(RFadultexp_sumLP_PC_model, pairwise ~ Phase|Group, adjust = 'none')
summary(RFadultexp_sumLP_PC_within_emmeans)

RFadultexp_sumLP_PC_between_emmeans <-emmeans(RFadultexp_sumLP_PC_model, pairwise ~ Group|Phase, adjust = 'none')
summary(RFadultexp_sumLP_PC_between_emmeans)


### Export Analysis Tables -----
#### Phase Tables
RFadultexp_sumLP_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadultexp_sumLP_model),sd_session_sumLP,mean_session_sumLP)
RFadultexp_sumLP_model_fittable <- mixedmodel_fittable(RFadultexp_sumLP_model, RFadultexp_sumLP_nullmodel)
RFadultexp_sumLP_within_emmeans_contraststable <- mixedmodel_rescalecontrasts(mixedmodel_emmeanstable_Group(RFadultexp_sumLP_within_emmeans,RFadultexp_sumLP_model),sd_session_sumLP)
RFadultexp_sumLP_between_emmeans_contraststable <- mixedmodel_rescalecontrasts(mixedmodel_emmeanstable_Phase(RFadultexp_sumLP_between_emmeans,RFadultexp_sumLP_model),sd_session_sumLP)

write.csv(RFadultexp_sumLP_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_sumLP_model_paramtable.csv",sep=""))
write.csv(RFadultexp_sumLP_model_fittable, paste(path_analysis_RFtest,"RFadultexp_sumLP_model_fittable.csv",sep=""))
write.csv(RFadultexp_sumLP_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_sumLP_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_sumLP_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_sumLP_emmeans_between_contraststable.csv",sep=""))

#### Session Tables
##### Baseline
RFadultexp_sumLP_BL_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadultexp_sumLP_BL_model),sd_session_sumLP,mean_session_sumLP)
RFadultexp_sumLP_BL_model_fittable <- mixedmodel_fittable(RFadultexp_sumLP_BL_model, RFadultexp_sumLP_BL_nullmodel)

write.csv(RFadultexp_sumLP_BL_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_sumLP_Session_BL_model_paramtable.csv",sep=""))
write.csv(RFadultexp_sumLP_BL_model_fittable, paste(path_analysis_RFtest,"RFadultexp_sumLP_Session_BL_model_fittable.csv",sep=""))

##### Injection
RFadultexp_sumLP_INJ_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadultexp_sumLP_INJ_model),sd_session_sumLP,mean_session_sumLP)
RFadultexp_sumLP_INJ_model_fittable <- mixedmodel_fittable(RFadultexp_sumLP_INJ_model, RFadultexp_sumLP_INJ_nullmodel)
RFadultexp_sumLP_INJ_within_emmeans_contraststable <- mixedmodel_rescalecontrasts(mixedmodel_emmeanstable_Group(RFadultexp_sumLP_INJ_within_emmeans,RFadultexp_sumLP_INJ_model),sd_session_sumLP)
RFadultexp_sumLP_INJ_between_emmeans_contraststable <- mixedmodel_rescalecontrasts(mixedmodel_emmeanstable_Session(RFadultexp_sumLP_INJ_between_emmeans,RFadultexp_sumLP_INJ_model),sd_session_sumLP)

write.csv(RFadultexp_sumLP_INJ_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_sumLP_Session_INJ_model_paramtable.csv",sep=""))
write.csv(RFadultexp_sumLP_INJ_model_fittable, paste(path_analysis_RFtest,"RFadultexp_sumLP_Session_INJ_model_fittable.csv",sep=""))
write.csv(RFadultexp_sumLP_INJ_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_sumLP_Session_INJ_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_sumLP_INJ_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_sumLP_Session_INJ_emmeans_between_contraststable.csv",sep=""))


##### Post-Injection Early
RFadultexp_sumLP_POST_model_paramtable <- mixedmodel_rescaleparamtable(mixedmodel_paramtable(RFadultexp_sumLP_POST_model),sd_session_sumLP,mean_session_sumLP)
RFadultexp_sumLP_POST_model_fittable <- mixedmodel_fittable(RFadultexp_sumLP_POST_model, RFadultexp_sumLP_POST_nullmodel)
RFadultexp_sumLP_POST_within_emmeans_contraststable <- mixedmodel_rescalecontrasts(mixedmodel_emmeanstable_Group(RFadultexp_sumLP_POST_within_emmeans,RFadultexp_sumLP_POST_model),sd_session_sumLP)
RFadultexp_sumLP_POST_between_emmeans_contraststable <- mixedmodel_rescalecontrasts(mixedmodel_emmeanstable_Session(RFadultexp_sumLP_POST_between_emmeans,RFadultexp_sumLP_POST_model),sd_session_sumLP)

write.csv(RFadultexp_sumLP_POST_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_sumLP_Session_POST_model_paramtable.csv",sep=""))
write.csv(RFadultexp_sumLP_POST_model_fittable, paste(path_analysis_RFtest,"RFadultexp_sumLP_Session_POST_model_fittable.csv",sep=""))
write.csv(RFadultexp_sumLP_POST_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_sumLP_Session_POST_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_sumLP_POST_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_sumLP_Session_POST_emmeans_between_contraststable.csv",sep=""))


#### Percent Change Group Tables
RFadultexp_sumLP_PC_model_paramtable <- mixedmodel_paramtable(RFadultexp_sumLP_PC_model)
RFadultexp_sumLP_PC_model_fittable <- mixedmodel_fittable(RFadultexp_sumLP_PC_model, RFadultexp_sumLP_PC_nullmodel)
RFadultexp_sumLP_PC_within_emmeans_contraststable <- mixedmodel_emmeanstable_Group(RFadultexp_sumLP_PC_within_emmeans,RFadultexp_sumLP_PC_model)
RFadultexp_sumLP_PC_between_emmeans_contraststable <- mixedmodel_emmeanstable_Phase(RFadultexp_sumLP_PC_between_emmeans,RFadultexp_sumLP_PC_model)

write.csv(RFadultexp_sumLP_PC_model_paramtable, paste(path_analysis_RFtest,"RFadultexp_sumLP_PC_model_paramtable.csv",sep=""))
write.csv(RFadultexp_sumLP_PC_model_fittable, paste(path_analysis_RFtest,"RFadultexp_sumLP_PC_model_fittable.csv",sep=""))
write.csv(RFadultexp_sumLP_PC_within_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_sumLP_PC_emmeans_within_contraststable.csv",sep=""))
write.csv(RFadultexp_sumLP_PC_between_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultexp_sumLP_PC_emmeans_between_contraststable.csv",sep=""))




## Figure Panels -----
### Max Lever Presses -----
ymax_maxLPC <- 4
ymin_maxLPC <- -24

ymax_maxLP_bar <- 48
ymin_maxLP_bar <- 0

injmarker_maxLP <- c(rep(ymin_maxLPC,7))
injticks_maxLP <- data.frame(injsession,injmarker_maxLP)

# maxLP Change Line
maxLPC_line <- ggplot(subset(RFadultexpmeans_loco,Phase!="BL"), aes(x=sessionlabel, y=pass_maxLPC_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadultexpmeans_loco,Phase!="BL"), aes(ymin=pass_maxLPC_mean-pass_maxLPC_se, ymax=pass_maxLPC_mean+pass_maxLPC_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name ="Session", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLPC, ymax_maxLPC), breaks = seq(-24, ymax_maxLPC, 8))+ 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("Max Lever Presses", sep = "")) +
  geom_point(data=injticks_maxLP, aes(x=injsession, y=injmarker_maxLP, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE)+ # Add injection ticks
  annotate("text", x = 1, y = -17.5, label = "***", size = mystarsize) + # Add asterisks for morphine group session significance
  annotate("text", x = 2, y = -13, label = "***", size = mystarsize)+
  annotate("text", x = 3, y = -15, label = "***", size = mystarsize)+
  annotate("text", x = 4, y = -13.5, label = "***", size = mystarsize)+
  annotate("text", x = 5, y = -11.5, label = "***", size = mystarsize)+
  annotate("text", x = 6, y = -9.5, label = "*", size = mystarsize)+
  annotate("text", x = 7, y = -12.5, label = "*", size = mystarsize)+
  annotate("text", x = 1, y = -18.5, label = "###", size = 2)+ # Add hashtags for control group session significance
  annotate("text", x = 2, y = -14, label = "###", size = 2)+
  annotate("text", x = 3, y = -16, label = "##", size = 2)+
  annotate("text", x = 4, y = -14.5, label = "###", size = 2)+
  annotate("text", x = 5, y = -12.5, label = "###", size = 2)+
  annotate("text", x = 6, y = -10.5, label = "#", size = 2)+
  annotate("text", x = 7, y = -13.5, label = "###", size = 2)+
  annotate("text", x = 7, y = -15.5, label = "&", size = 2) # Add & sign for significance between control and morphine groups

maxLPC_line

# Adolescent Saline Group maxLP by Phase Bar
maxLP_AdolS_bar <- ggplot(subset(RFadultexpmeans_Phase_loco,Group=="AdolS"), aes(x=Phase, y=pass_maxLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadultexpmeans_SubjectPhase_loco, Group=="AdolS"), mapping = aes(x = Phase, y = pass_maxLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadultexpmeans_Phase_loco,Group=="AdolS"), aes(ymin=pass_maxLP_mean-pass_maxLP_se, ymax=pass_maxLP_mean+pass_maxLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLP_bar, ymax_maxLP_bar), breaks = seq(ymin_maxLP_bar, ymax_maxLP_bar, 15)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = AdolSMphasecolors) +
  scale_fill_manual("legend", values = AdolSMphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Max Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Adol Control", sep = ""))+
  geom_bracket(xmin=c(1), xmax = c(1.9), y.position=c(20) ,label.size = 2, label = "###",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2) +
  geom_bracket(xmin=c(2.1), xmax = c(3), y.position=c(20) ,label.size = 2, label = "###",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)+
  geom_bracket(xmin=c(2.1), xmax = c(4), y.position=c(27) ,label.size = 2, label = "###",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

maxLP_AdolS_bar

# Adolescent Morphine Group maxLP by Phase Bar
maxLP_AdolM_bar <- ggplot(subset(RFadultexpmeans_Phase_loco,Group=="AdolM"), aes(x=Phase, y=pass_maxLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadultexpmeans_SubjectPhase_loco, Group=="AdolM"), mapping = aes(x = Phase, y = pass_maxLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadultexpmeans_Phase_loco,Group=="AdolM"), aes(ymin=pass_maxLP_mean-pass_maxLP_se, ymax=pass_maxLP_mean+pass_maxLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLP_bar, ymax_maxLP_bar), breaks = seq(ymin_maxLP_bar, ymax_maxLP_bar, 15)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = AdolMMphasecolors) +
  scale_fill_manual("legend", values = AdolMMphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Max Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Adol Morphine", sep = ""))+
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(45) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(1.9), y.position=c(30) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.1), xmax = c(3), y.position=c(30) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.1), xmax = c(4), y.position=c(38) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

maxLP_AdolM_bar
# Create full panel
maxLP <- ggarrange(maxLP_AdolS_bar,maxLP_AdolM_bar,maxLPC_line, ncol=3, nrow=1, widths=c(1.8,1.8,4), labels = c("A", "B", "C"))
maxLP

### Total Session Lever Presses -----
ymax_sumLPC <- 260
ymin_sumLPC <- -270

ymax_sumLP_bar <- 1200
ymin_sumLP_bar <- 0

injmarker_sumLP <- c(rep(ymin_sumLPC,7))
injticks_sumLP <- data.frame(injsession,injmarker_sumLP)

# Total Session Lever Presses Change Line
sumLPC_line <- ggplot(subset(RFadultexpmeans_loco,Phase!="BL"), aes(x=sessionlabel, y=session_sumLPC_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadultexpmeans_loco,Phase!="BL"), aes(ymin=session_sumLPC_mean-session_sumLPC_se, ymax=session_sumLPC_mean+session_sumLPC_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name ="Session", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sumLPC, ymax_sumLPC), breaks = seq(-250, ymax_sumLPC, 125)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("Total Session Lever Presses", sep = "")) +
  geom_point(data=injticks_sumLP, aes(x=injsession, y=injmarker_sumLP, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE)+ # Add injection ticks
  annotate("text", x = 1, y = -200, label = "*", size = mystarsize) + # Add asterisks for morphine group session significance
  annotate("text", x = 6, y = 225, label = "*", size = mystarsize)+
  annotate("text", x = 8, y = -190, label = "***", size = mystarsize) + # Add asterisks for morphine group session significance
  annotate("text", x = 9, y = -145, label = "**", size = mystarsize)+
  annotate("text", x = 11, y = -135, label = "*", size = mystarsize)+
  annotate("text", x = 1, y = -215, label = "##", size = 2)+ # Add hashtags for control group session significance
  annotate("text", x = 2, y = -185, label = "##", size = 2)+
  annotate("text", x = 13, y = -90, label = "##", size = 2)

sumLPC_line

# Adolescent Saline Group sumLP by Phase Bar
sumLP_AdolS_bar <- ggplot(subset(RFadultexpmeans_Phase_loco,Group=="AdolS"), aes(x=Phase, y=session_sumLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadultexpmeans_SubjectPhase_loco, Group=="AdolS"), mapping = aes(x = Phase, y = session_sumLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadultexpmeans_Phase_loco,Group=="AdolS"), aes(ymin=session_sumLP_mean-session_sumLP_se, ymax=session_sumLP_mean+session_sumLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sumLP_bar, ymax_sumLP_bar), breaks = seq(ymin_sumLP_bar, ymax_sumLP_bar, 400)) +
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = AdolSMphasecolors) +
  scale_fill_manual("legend", values = AdolSMphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Session Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Adol Control", sep = ""))+
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(650) ,label.size = 2, label = "##",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

sumLP_AdolS_bar

# Adolescent Morphine Group sumLP by Phase Bar
sumLP_AdolM_bar <- ggplot(subset(RFadultexpmeans_Phase_loco,Group=="AdolM"), aes(x=Phase, y=session_sumLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(subset(RFadultexpmeans_SubjectPhase_loco, Group=="AdolM"), mapping = aes(x = Phase, y = session_sumLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=subset(RFadultexpmeans_Phase_loco,Group=="AdolM"), aes(ymin=session_sumLP_mean-session_sumLP_se, ymax=session_sumLP_mean+session_sumLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sumLP_bar, ymax_sumLP_bar), breaks = seq(ymin_sumLP_bar, ymax_sumLP_bar, 400)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = AdolMMphasecolors) +
  scale_fill_manual("legend", values = AdolMMphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Session Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Adol Morphine", sep = ""))+
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(950) ,label.size = mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2), xmax = c(2.9), y.position=c(800) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3.1), xmax = c(4), y.position=c(800) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

sumLP_AdolM_bar

# Create full panel
sumLP <- ggarrange(sumLP_AdolS_bar,sumLP_AdolM_bar,sumLPC_line, ncol=3, nrow=1, widths=c(1.8,1.8,4), labels = c("D", "E", "F"))
sumLP


### Max Lever Presses Percent Change -----
ymax_maxLPPC <- 111
ymin_maxLPPC <- -100

maxLP_PC_bar <- ggplot(subset(RFadultexpmeans_Phase_loco,Phase!="BL"), aes(x=Phase, y=pass_maxLPPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadultexpmeans_Phase_loco,Phase!="BL"), aes(ymin=pass_maxLPPC_mean-pass_maxLPPC_se, ymax=pass_maxLPPC_mean+pass_maxLPPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  #geom_point(subset(PCmeans_lm_subject, Phase!="BL"), mapping = aes(x = Phase, y = maxLPPC_mean), position = position_dodge(width = 0.95), size = mysubjectpointsize, stroke = mypointstroke, shape = 1, color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLPPC, ymax_maxLPPC), breaks = seq(-100, ymax_maxLPPC, 50)) + 
  scale_x_discrete(expand = c(0.25, 0), name ="Phase", labels=c('INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change from BL") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Max Lever Presses", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)+
  geom_bracket(xmin=c(1.83), xmax = c(2.17), y.position=c(102) ,label.size = 2, label = "&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

maxLP_PC_bar

### Total Session Lever Presses Percent Change -----
ymax_sumLPPC <- 32
ymin_sumLPPC <- -60

sumLP_PC_bar <- ggplot(subset(RFadultexpmeans_Phase_loco,Phase!="BL"), aes(x=Phase, y=session_sumLPPC_mean, Group=Group, fill=PCfill, color=PCfill)) +
  geom_col(position = position_dodge(width = 0.8), width=0.7, linetype="solid", linewidth = 0) +
  geom_errorbar(data=subset(RFadultexpmeans_Phase_loco,Phase!="BL"), aes(ymin=session_sumLPPC_mean-session_sumLPPC_se, ymax=session_sumLPPC_mean+session_sumLPPC_se), position = position_dodge(width = 0.8), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  #geom_point(subset(PCmeans_lm_subject, Phase!="BL"), mapping = aes(x = Phase, y = sumLPPC_mean), position = position_dodge(width = 0.95), size = mysubjectpointsize, stroke = mypointstroke, shape = 1, color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sumLPPC, ymax_sumLPPC), breaks = seq(ymin_sumLPPC, ymax_sumLPPC, 30)) + 
  scale_x_discrete(expand = c(0.25, 0), name ="Phase", labels=c('INJ','POST (Day 8-10)','POST (Day 11-16)')) + 
  scale_color_manual("legend", values = PCphasecolors) +
  scale_fill_manual("legend", values = PCphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Percent Change from BL") + xlab('Phase') +
  mytheme +
  ggtitle(paste("Total Session Lever Presses", sep = "")) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth = 0.5)+
  geom_bracket(xmin=c(0.83), xmax = c(1.17), y.position=c(28) ,label.size = 2, label = "&&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

sumLP_PC_bar

# Create PC Panel
sumLPmaxLPPC <- ggarrange(maxLP_PC_bar,sumLP_PC_bar, ncol=2, nrow=1, labels = c("G", "H"))
sumLPmaxLPPC

## Create Supplemental Figure -----
# Create full panel
suppfigure5 <- ggarrange(maxLP,sumLP,sumLPmaxLPPC, ncol=1, nrow=3)
suppfigure5
ggsave(plot=suppfigure5, path=path_figures, device="pdf", filename="SupplementaryF5_AllPanels.pdf", width=RFpanelwidth, height=RFpanelheight*3.5)


