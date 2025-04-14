############### NAIVE ADULT MORPHINE ADMINISTRATION ANALYSIS ###################
# This program takes in data processed by the script "RF Theta Calculation"    
# and run statistics to analyze findings. Data from each experiment are loaded,
# analyzed, and visualized.                                                    

# PREPARE R -----
## Set up paths
computerusergithubpath <- c('C:/Users/rmdon/OneDrive/Desktop/GitHub_Repositories/') # Path for Rachel's laptop to github repositories
computeruserboxpath <- c('C:/Users/rmdon/Box/') # Path for Rachel's laptop to Box drive

## Load required libraries
source(paste(computerusergithubpath,"Analysis-PBAdolescentMorphineReward/Functions/Functions_LoadPackages.R",sep=''))

## Load all custom functions
source.all(paste(computerusergithubpath,'JRoitmanICSS/Functions/Standard RF ICSS/',sep=''))
source.all(paste(computerusergithubpath,'Analysis-PBAdolescentMorphineReward/Functions/',sep=''))


## Set working directory to the raw data folder
setwd(paste(computeruserboxpath,'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Raw Data',sep=''))

## Set up paths to export figures and analysis to
### Overall paths
path_analysis <- c(paste(computeruserboxpath, 'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Analysis/',sep=''))
path_figures <- c(paste(computeruserboxpath, 'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Figures/Figure Panels/',sep=''))

## Specific paths
path_analysis_RFtest <- paste(path_analysis,"Adult Naive Rate Frequency Morphine Test/",sep="")
path_analysis_RFtest_means <- paste(path_analysis,"Adult Naive Rate Frequency Morphine Test/Means Tables/",sep="")
path_figures_RFtest <- paste(path_analysis_RFtest,"Figures/",sep="")
path_subjectfigures_RFtest <- paste(path_figures_RFtest,"Subject Figures/",sep="") # Subject shaping figures path

## Read in subject key
subjectkey <- read_csv("SubjectKey_AdultNaive.csv")
subjectexc <- c()

# PREPARE DATA -----
## Read in data and prepare variables -----
### Prep col names
freqcolnames <- c('141', '126', '112', '100', '89', '79', '71', '63', '56', '50', '45', '40','35','32','28') # Set column names of lever presses for each frequency

### Read in data and fit model
RFadultnaivedata_raw <- RFICSS_fitmodel("RawData_RFICSSMorphineTest_AdultNaive.csv")

### Add Sex and Group to raw
RFadultnaivedata_raw$Sex <- subjectkey$Sex[match(RFadultnaivedata_raw$SubjectID, subjectkey$SubjectID)]
RFadultnaivedata_raw$Group <- subjectkey$ExperimentalGroup[match(RFadultnaivedata_raw$SubjectID, subjectkey$SubjectID)]
RFadultnaivedata_raw$ShapingInclude <- subjectkey$Include[match(RFadultnaivedata_raw$SubjectID, subjectkey$SubjectID)]

### Prepare categorical variables
RFadultnaivedata_raw$Sex <- factor(RFadultnaivedata_raw$Sex, levels = c("M","F"))
RFadultnaivedata_raw$Group <- factor(RFadultnaivedata_raw$Group)

RFadultnaivedata_raw$Session <- factor(RFadultnaivedata_raw$Day) # Factor day by phase
RFadultnaivedata_raw$Session_collapsed <- ifelse(RFadultnaivedata_raw$Phase == "BL", 0, RFadultnaivedata_raw$Day) # Create collapsed Day variable where Baseline=0
RFadultnaivedata_raw$Session_collapsed <- factor(RFadultnaivedata_raw$Session_collapsed, levels = unique(RFadultnaivedata_raw$Session_collapsed)) # Factor collapsed Day variable



RFadultnaivedata_raw$Phase_raw <- factor(RFadultnaivedata_raw$Phase, levels = c("BL","INJ","POST"))
RFadultnaivedata_raw$Phase <- factor(if_else(RFadultnaivedata_raw$Phase=="POST" & RFadultnaivedata_raw$Day <=3,"POST_EARLY", 
                                                      if_else(RFadultnaivedata_raw$Phase=="POST" & RFadultnaivedata_raw$Day >=4,"POST_LATE", 
                                                              RFadultnaivedata_raw$Phase)), levels = c("BL","INJ","POST_EARLY","POST_LATE"))


### Add max lever presses, and total lever presses per pass
RFadultnaivedata_raw$maxLP <- rowMaxs(as.matrix(RFadultnaivedata_raw[freqcolnames])) # Max lever pressing per pass
RFadultnaivedata_raw$sumLP <- rowSums(as.matrix(RFadultnaivedata_raw[freqcolnames])) # Max lever pressing per pass

### Check Variables
unique(RFadultnaivedata_raw$SubjectID)
unique(RFadultnaivedata_raw$Sex)
unique(RFadultnaivedata_raw$Group)
unique(RFadultnaivedata_raw$ShapingInclude)

unique(RFadultnaivedata_raw$Day)
unique(RFadultnaivedata_raw$Session)
unique(RFadultnaivedata_raw$Session_collapsed)

unique(RFadultnaivedata_raw$Phase) 
unique(RFadultnaivedata_raw$Phase_raw) 
#view(RFadultnaivedata_raw)

## Prepare subject averages by day -----
### Find subject averages by day
RFadultnaivedata_lm <- RFadultnaivedata_raw %>% filter(FinalInclusion == "INCLUDE", Pass > 1, !SubjectID %in% subjectexc, ShapingInclude=="INCLUDE") %>%
  group_by(SubjectID, Sex, Group, ShapingInclude, UniqueSessionID, Phase, Phase_raw, Day, Session, Session_collapsed, Box, Amplitude) %>%
  dplyr::summarize(npass = n_distinct(Pass),
                   mean_theta = mean(theta, na.rm=TRUE), # find the means for theta, M50, alpha, sumLP, and maxLP
                   mean_M50 = mean(M50, na.rm=TRUE),
                   mean_alpha = mean(alpha, na.rm=TRUE),
                   mean_alpha95 = mean(alpha95, na.rm=TRUE))

RFadultnaivedata_loco <- RFadultnaivedata_raw %>% filter(Pass > 1, !SubjectID %in% subjectexc, ShapingInclude=="INCLUDE") %>%
  group_by(SubjectID, Sex, Group, ShapingInclude, UniqueSessionID, Phase, Phase_raw, Day, Session, Session_collapsed, Box, Amplitude) %>%
  dplyr::summarize(npass = n_distinct(Pass),
                   pass_sumLP = mean(TotalLeverPresses, na.rm=TRUE),
                   session_sumLP = sum(TotalLeverPresses, na.rm=TRUE),
                   pass_maxLP = mean(maxLP, na.rm=TRUE),
                   finalblock_loco = mean(BeamBreaks, na.rm=TRUE))

### Add percent change from baseline
RFadultnaivesubject_lm <- RFICSS_LM_findChange(RFadultnaivedata_lm,"BL")
RFadultnaivesubject_loco <- RFICSS_loco_findChange(RFadultnaivedata_loco,"BL")


# PREPARE MEANS -----
## Individual Session Means -----
sessionlabels <- (c(1:15))

RFadultnaivemeans_lm <- RFadultnaivesubject_lm %>%
  group_by(UniqueSessionID,Phase,Phase_raw,Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   theta_mean = mean(mean_theta, na.rm=TRUE), theta_sd = sd(mean_theta, na.rm=TRUE), theta_se = theta_sd / sqrt(n),
                   thetaC_mean = mean(mean_theta_C, na.rm=TRUE), thetaC_sd = sd(mean_theta_C, na.rm=TRUE), thetaC_se = thetaC_sd / sqrt(n),
                   
                   M50_mean = mean(mean_M50, na.rm=TRUE), M50_sd = sd(mean_M50, na.rm=TRUE), M50_se = M50_sd / sqrt(n),
                   M50C_mean = mean(mean_M50_C, na.rm=TRUE), M50C_sd = sd(mean_M50_C, na.rm=TRUE), M50C_se = M50C_sd / sqrt(n),
                   
                   alpha_mean = mean(mean_alpha, na.rm=TRUE), alpha_sd = sd(mean_alpha, na.rm=TRUE), alpha_se = alpha_sd / sqrt(n),
                   alphaC_mean = mean(mean_alpha_C, na.rm=TRUE), alphaC_sd = sd(mean_alpha_C, na.rm=TRUE), alphaC_se = alphaC_sd / sqrt(n))


RFadultnaivemeans_loco <- RFadultnaivesubject_loco %>%
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



RFadultnaivemeans_lm$sessionlabel = factor(RFadultnaivemeans_lm$UniqueSessionID, labels = c(1:19))
RFadultnaivemeans_lm

RFadultnaivemeans_loco$sessionlabel = factor(RFadultnaivemeans_loco$UniqueSessionID, labels = c(1:19))
RFadultnaivemeans_loco

## Overall Phase Means -----
RFadultnaivemeans_Phase_lm <- RFadultnaivesubject_lm %>%
  group_by(Phase,Phase_raw,Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   theta_mean = mean(mean_theta, na.rm=TRUE), theta_sd = sd(mean_theta, na.rm=TRUE), theta_se = theta_sd / sqrt(n),
                   thetaC_mean = mean(mean_theta_C, na.rm=TRUE), thetaC_sd = sd(mean_theta_C, na.rm=TRUE), thetaC_se = thetaC_sd / sqrt(n),
                   
                   M50_mean = mean(mean_M50, na.rm=TRUE), M50_sd = sd(mean_M50, na.rm=TRUE), M50_se = M50_sd / sqrt(n),
                   M50C_mean = mean(mean_M50_C, na.rm=TRUE), M50C_sd = sd(mean_M50_C, na.rm=TRUE), M50C_se = M50C_sd / sqrt(n),
                   
                   alpha_mean = mean(mean_alpha, na.rm=TRUE), alpha_sd = sd(mean_alpha, na.rm=TRUE), alpha_se = alpha_sd / sqrt(n),
                   alphaC_mean = mean(mean_alpha_C, na.rm=TRUE), alphaC_sd = sd(mean_alpha_C, na.rm=TRUE), alphaC_se = alphaC_sd / sqrt(n))


RFadultnaivemeans_Phase_loco <- RFadultnaivesubject_loco %>%
  group_by(Phase,Phase_raw,Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   
                   pass_maxLP_mean = mean(pass_maxLP, na.rm=TRUE), pass_maxLP_sd = sd(pass_maxLP,na.rm=TRUE), pass_maxLP_se = pass_maxLP_sd / sqrt(n),
                   pass_maxLPC_mean = mean(pass_maxLP_C, na.rm=TRUE), pass_maxLPC_sd = sd(pass_maxLP_C,na.rm=TRUE), pass_maxLPC_se = pass_maxLPC_sd / sqrt(n),
                   
                   session_sumLP_mean = mean(session_sumLP, na.rm=TRUE), session_sumLP_sd = sd(session_sumLP,na.rm=TRUE), session_sumLP_se = session_sumLP_sd / sqrt(n),
                   session_sumLPC_mean = mean(session_sumLP_C, na.rm=TRUE), session_sumLPC_sd = sd(session_sumLP_C,na.rm=TRUE), session_sumLPC_se = session_sumLPC_sd / sqrt(n),
                   
                   pass_sumLP_mean = mean(pass_sumLP, na.rm=TRUE), pass_sumLP_sd = sd(pass_sumLP,na.rm=TRUE), pass_sumLP_se = pass_sumLP_sd / sqrt(n),
                   pass_sumLPC_mean = mean(pass_sumLP_C, na.rm=TRUE), pass_sumLPC_sd = sd(pass_sumLP_C,na.rm=TRUE), pass_sumLPC_se = pass_sumLPC_sd / sqrt(n),
                   
                   finalblock_loco_mean = mean(finalblock_loco, na.rm=TRUE), finalblock_loco_sd = sd(finalblock_loco,na.rm=TRUE), finalblock_loco_se = finalblock_loco_sd / sqrt(n),
                   finalblock_locoC_mean = mean(finalblock_loco_C, na.rm=TRUE), finalblock_locoC_sd = sd(finalblock_loco_C,na.rm=TRUE), finalblock_locoC_se = finalblock_locoC_sd / sqrt(n))


## Subject Phase Means -----
RFadultnaivemeans_SubjectPhase_lm <- RFadultnaivesubject_lm %>%
  group_by(SubjectID,Sex,Phase,Phase_raw,Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   
                   theta_mean = mean(mean_theta, na.rm=TRUE), theta_sd = sd(mean_theta, na.rm=TRUE), theta_se = theta_sd / sqrt(n),
                   thetaC_mean = mean(mean_theta_C, na.rm=TRUE), thetaC_sd = sd(mean_theta_C, na.rm=TRUE), thetaC_se = thetaC_sd / sqrt(n),
                   
                   M50_mean = mean(mean_M50, na.rm=TRUE), M50_sd = sd(mean_M50, na.rm=TRUE), M50_se = M50_sd / sqrt(n),
                   M50C_mean = mean(mean_M50_C, na.rm=TRUE), M50C_sd = sd(mean_M50_C, na.rm=TRUE), M50C_se = M50C_sd / sqrt(n),
                   
                   alpha_mean = mean(mean_alpha, na.rm=TRUE), alpha_sd = sd(mean_alpha, na.rm=TRUE), alpha_se = alpha_sd / sqrt(n),
                   alphaC_mean = mean(mean_alpha_C, na.rm=TRUE), alphaC_sd = sd(mean_alpha_C, na.rm=TRUE), alphaC_se = alphaC_sd / sqrt(n))


RFadultnaivemeans_SubjectPhase_loco <- RFadultnaivesubject_loco %>%
  group_by(SubjectID,Sex,Phase,Phase_raw,Group) %>%
  dplyr::summarize(n = n_distinct(SubjectID),
                   
                   pass_maxLP_mean = mean(pass_maxLP, na.rm=TRUE), pass_maxLP_sd = sd(pass_maxLP,na.rm=TRUE), pass_maxLP_se = pass_maxLP_sd / sqrt(n),
                   pass_maxLPC_mean = mean(pass_maxLP_C, na.rm=TRUE), pass_maxLPC_sd = sd(pass_maxLP_C,na.rm=TRUE), pass_maxLPC_se = pass_maxLPC_sd / sqrt(n),
                   
                   session_sumLP_mean = mean(session_sumLP, na.rm=TRUE), session_sumLP_sd = sd(session_sumLP,na.rm=TRUE), session_sumLP_se = session_sumLP_sd / sqrt(n),
                   session_sumLPC_mean = mean(session_sumLP_C, na.rm=TRUE), session_sumLPC_sd = sd(session_sumLP_C,na.rm=TRUE), session_sumLPC_se = session_sumLPC_sd / sqrt(n),
                   
                   pass_sumLP_mean = mean(pass_sumLP, na.rm=TRUE), pass_sumLP_sd = sd(pass_sumLP,na.rm=TRUE), pass_sumLP_se = pass_sumLP_sd / sqrt(n),
                   pass_sumLPC_mean = mean(pass_sumLP_C, na.rm=TRUE), pass_sumLPC_sd = sd(pass_sumLP_C,na.rm=TRUE), pass_sumLPC_se = pass_sumLPC_sd / sqrt(n),
                   
                   finalblock_loco_mean = mean(finalblock_loco, na.rm=TRUE), finalblock_loco_sd = sd(finalblock_loco,na.rm=TRUE), finalblock_loco_se = finalblock_loco_sd / sqrt(n),
                   finalblock_locoC_mean = mean(finalblock_loco_C, na.rm=TRUE), finalblock_locoC_sd = sd(finalblock_loco_C,na.rm=TRUE), finalblock_locoC_se = finalblock_locoC_sd / sqrt(n))


### Export Means -----
write.csv(RFadultnaivemeans_Phase_lm, paste(path_analysis_RFtest_means,"RFadultnaive_Means_Phase_RFvariables.csv",sep=""))
write.csv(RFadultnaivemeans_Phase_loco, paste(path_analysis_RFtest_means,"RFadultnaive_Means_Phase_Motorvariables.csv",sep=""))

write.csv(RFadultnaivemeans_lm, paste(path_analysis_RFtest_means,"RFadultnaive_Means_PhaseandSession_RFvariables.csv",sep=""))
write.csv(RFadultnaivemeans_loco, paste(path_analysis_RFtest_means,"RFadultnaive_Means_PhaseandSession_Motorvariables.csv",sep=""))


# STATISTICS -----
## Descriptive Statistics -----

### Subject ns -----
RFadultnaive_ns <- RFadultnaivesubject_lm %>% group_by(Sex) %>% dplyr::summarize(n = n_distinct(SubjectID))
RFadultnaive_ns # Included subject ns

RFadultnaive_excluded_ns <- RFadultnaivedata_raw %>% filter(ShapingInclude!="INCLUDE") %>% group_by(Sex) %>% dplyr::summarize(n = n_distinct(SubjectID))
RFadultnaive_excluded_ns # Excluded subject ns

amps <- RFadultnaivedata_lm %>% group_by(SubjectID, Amplitude) %>% dplyr::summarize() # Included subject amplitudes

min(amps$Amplitude) # Overall group minimum amplitude
max(amps$Amplitude) # Overall group maximum amplitude

## M50 -----
### Exploratory Analysis: Sex -----
#### Mixed effects model with sex interaction and random intercept by subject
RFadultnaive_M50_sexInt_model <- glmmTMB(mean_M50 ~ Phase*Sex + (1|SubjectID), data = RFadultnaivedata_lm, family=stats::gaussian)
summary(RFadultnaive_M50_sexInt_model)
Anova(RFadultnaive_M50_sexInt_model, type=2)

#### Mixed effects model with sex fixed effect and random intercept by subject
RFadultnaive_M50_sex_model <- glmmTMB(mean_M50 ~ Phase + Sex + (1|SubjectID), data = RFadultnaivedata_lm, family=stats::gaussian)
summary(RFadultnaive_M50_sex_model)
Anova(RFadultnaive_M50_sex_model, type=2)

##### Overall, there is a significant main effect of sex (females have lower M50), but there is no interaction with phase.
##### As the interaction is lacking and there was not an a priori hypothesis of sex differences, sex is included as a
##### fixed effect in final models.

### Between Phase Analysis -----
#### Mixed effects model by phase and sex with random intercept by subject
RFadultnaive_M50_model <- glmmTMB(mean_M50 ~ Phase + Sex + (1|SubjectID), data = RFadultnaivedata_lm, family=stats::gaussian)

summary(RFadultnaive_M50_model)
Anova(RFadultnaive_M50_model, type=2)

#### Compare full model to the null model
RFadultnaive_M50_nullmodel <- glmmTMB(mean_M50 ~ 1 + Sex + (1|SubjectID), data = RFadultnaivedata_lm, family = stats::gaussian()) # Make the null model
anova(RFadultnaive_M50_nullmodel, RFadultnaive_M50_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultnaive_M50_emmeans <-emmeans(RFadultnaive_M50_model, pairwise ~ Phase, adjust = 'none')
summary(RFadultnaive_M50_emmeans)

### Within Phase Analysis -----
#### Baseline
##### Mixed effects model by day with random intercept by subject
RFadultnaive_M50_BL_model <- glmmTMB(mean_M50 ~ Session + Sex + (1|SubjectID), data = subset(RFadultnaivedata_lm, Phase=="BL"),
  family = stats::gaussian())

summary(RFadultnaive_M50_BL_model)
Anova(RFadultnaive_M50_BL_model, type=2)

##### Compare full model to the null model
RFadultnaive_M50_BL_nullmodel <- glmmTMB(mean_M50 ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_lm, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_M50_BL_nullmodel, RFadultnaive_M50_BL_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultnaive_M50_BL_emmeans <- emmeans(RFadultnaive_M50_BL_model, pairwise ~ Session, adjust = 'none')
summary(RFadultnaive_M50_BL_emmeans)

###### Overall, there is no significant effect of day in the baseline phase.


#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultnaive_M50_INJ_model <- glmmTMB(mean_M50 ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_lm, Phase %in% c("BL","INJ")),
                                     family = stats::gaussian())

summary(RFadultnaive_M50_INJ_model)
Anova(RFadultnaive_M50_INJ_model, type=2)

##### Compare full model to the null model
RFadultnaive_M50_INJ_nullmodel <- glmmTMB(mean_M50 ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_lm, Phase %in% c("BL","INJ")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_M50_INJ_nullmodel, RFadultnaive_M50_INJ_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultnaive_M50_INJ_emmeans <- emmeans(RFadultnaive_M50_INJ_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_M50_INJ_emmeans)
###### Overall, injections days 1, 4, 5, 6, and 7 are significantly different from baseline.
###### Comparing across only the injection days, there are no consistent significant differences.

#### Post-Injection - Early Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultnaive_M50_POSTE_model <- glmmTMB(mean_M50 ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_lm, Phase %in% c("BL","POST_EARLY")),
                                      family = stats::gaussian())

summary(RFadultnaive_M50_POSTE_model)
Anova(RFadultnaive_M50_POSTE_model, type=2)

##### Compare full model to the null model
RFadultnaive_M50_POSTE_nullmodel <- glmmTMB(mean_M50 ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_lm, Phase %in% c("BL","POST_EARLY")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_M50_POSTE_nullmodel, RFadultnaive_M50_POSTE_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultnaive_M50_POSTE_emmeans <- emmeans(RFadultnaive_M50_POSTE_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_M50_POSTE_emmeans)
##### Overall, early phase post injection days 1 and 3 are significantly different from baseline.
##### Comparing across only the post-injection days, there are no consistent significant differences.


#### Post-Injection - Late Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultnaive_M50_POSTL_model <- glmmTMB(mean_M50 ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_lm, Phase %in% c("BL","POST_LATE")),
                                        family = stats::gaussian())

summary(RFadultnaive_M50_POSTL_model)
Anova(RFadultnaive_M50_POSTL_model, type=2)

##### Compare full model to the null model
RFadultnaive_M50_POSTL_nullmodel <- glmmTMB(mean_M50 ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_lm, Phase %in% c("BL","POST_LATE")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_M50_POSTL_nullmodel, RFadultnaive_M50_POSTL_model, test = "LRT") # Full model is a better fit than the null mode

##### Follow up with emmeans
RFadultnaive_M50_POSTL_emmeans <- emmeans(RFadultnaive_M50_POSTL_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_M50_POSTL_emmeans)
##### Overall, late phase post injection days are not significantly different from baseline.
##### Only day 8 is significantly higher than baseline.
##### Comparing across only the post-injection days, there are no consistent significant differences.

### Export analysis tables -----
#### Between Phase Tables
RFadultnaive_M50_model_paramtable <- mixedmodel_paramtable(RFadultnaive_M50_model)
RFadultnaive_M50_model_fittable <- mixedmodel_fittable(RFadultnaive_M50_model, RFadultnaive_M50_nullmodel)
RFadultnaive_M50_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_M50_emmeans,RFadultnaive_M50_model)

write.csv(RFadultnaive_M50_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_M50_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_M50_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_M50_model_fittable.csv",sep=""))
write.csv(RFadultnaive_M50_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_M50_emmeans_contraststable.csv",sep=""))

#### Within Phase Tables
##### Baseline
RFadultnaive_M50_BL_model_paramtable <- mixedmodel_paramtable(RFadultnaive_M50_BL_model)
RFadultnaive_M50_BL_model_fittable <- mixedmodel_fittable(RFadultnaive_M50_BL_model, RFadultnaive_M50_BL_nullmodel)
RFadultnaive_M50_BL_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_M50_BL_emmeans,RFadultnaive_M50_BL_model)

write.csv(RFadultnaive_M50_BL_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_M50_BL_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_M50_BL_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_M50_BL_model_fittable.csv",sep=""))
write.csv(RFadultnaive_M50_BL_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_M50_BL_emmeans_contraststable.csv",sep=""))

##### Injection
RFadultnaive_M50_INJ_model_paramtable <- mixedmodel_paramtable(RFadultnaive_M50_INJ_model)
RFadultnaive_M50_INJ_model_fittable <- mixedmodel_fittable(RFadultnaive_M50_INJ_model, RFadultnaive_M50_INJ_nullmodel)
RFadultnaive_M50_INJ_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_M50_INJ_emmeans,RFadultnaive_M50_INJ_model)

write.csv(RFadultnaive_M50_INJ_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_M50_INJ_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_M50_INJ_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_M50_INJ_model_fittable.csv",sep=""))
write.csv(RFadultnaive_M50_INJ_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_M50_INJ_emmeans_contraststable.csv",sep=""))

##### Post-Injection Early
RFadultnaive_M50_POSTE_model_paramtable <- mixedmodel_paramtable(RFadultnaive_M50_POSTE_model)
RFadultnaive_M50_POSTE_model_fittable <- mixedmodel_fittable(RFadultnaive_M50_POSTE_model, RFadultnaive_M50_POSTE_nullmodel)
RFadultnaive_M50_POSTE_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_M50_POSTE_emmeans,RFadultnaive_M50_POSTE_model)

write.csv(RFadultnaive_M50_POSTE_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_M50_POSTE_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_M50_POSTE_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_M50_POSTE_model_fittable.csv",sep=""))
write.csv(RFadultnaive_M50_POSTE_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_M50_POSTE_emmeans_contraststable.csv",sep=""))

##### Post-Injection Late
RFadultnaive_M50_POSTL_model_paramtable <- mixedmodel_paramtable(RFadultnaive_M50_POSTL_model)
RFadultnaive_M50_POSTL_model_fittable <- mixedmodel_fittable(RFadultnaive_M50_POSTL_model, RFadultnaive_M50_POSTL_nullmodel)
RFadultnaive_M50_POSTL_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_M50_POSTL_emmeans,RFadultnaive_M50_POSTL_model)

write.csv(RFadultnaive_M50_POSTL_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_M50_POSTL_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_M50_POSTL_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_M50_POSTL_model_fittable.csv",sep=""))
write.csv(RFadultnaive_M50_POSTL_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_M50_POSTL_emmeans_contraststable.csv",sep=""))



## Block Locomotion -----
### Between Phase Analysis -----
#### Mixed effects model with random intercepts by subject and sex
RFadultnaive_blockloco_model <- glmmTMB(finalblock_loco ~ Phase + Sex + (1|SubjectID), data = RFadultnaivedata_loco, family=stats::gaussian)

summary(RFadultnaive_blockloco_model)
Anova(RFadultnaive_blockloco_model, type=2)

#### Compare full model to the null model
RFadultnaive_blockloco_nullmodel <- glmmTMB(finalblock_loco ~ 1 + Sex + (1|SubjectID), data = RFadultnaivedata_loco, family = stats::gaussian()) # Make the null model
anova(RFadultnaive_blockloco_nullmodel, RFadultnaive_blockloco_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultnaive_blockloco_emmeans <-emmeans(RFadultnaive_blockloco_model,pairwise~Phase, adjust = 'none')
summary(RFadultnaive_blockloco_emmeans)


### Within Phase Analysis -----
#### Baseline
##### Mixed effects model by day with random intercept by subject
RFadultnaive_blockloco_BL_model <- glmmTMB(finalblock_loco ~ Session + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase=="BL"),
                                       family = stats::gaussian())

summary(RFadultnaive_blockloco_BL_model)
Anova(RFadultnaive_blockloco_BL_model, type=2)

##### Compare full model to the null model
RFadultnaive_blockloco_BL_nullmodel <- glmmTMB(finalblock_loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_blockloco_BL_nullmodel, RFadultnaive_blockloco_BL_model, test = "LRT") # Full model is not a better fit than the null model

##### Follow up with emmeans
RFadultnaive_blockloco_BL_emmeans <- emmeans(RFadultnaive_blockloco_BL_model, pairwise ~ Session, adjust = 'none')
summary(RFadultnaive_blockloco_BL_emmeans)
##### Overall, there is no significant effect of day in the baseline phase.

#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultnaive_blockloco_INJ_model <- glmmTMB(finalblock_loco ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","INJ")),
                                        family = stats::gaussian())

summary(RFadultnaive_blockloco_INJ_model)
Anova(RFadultnaive_blockloco_INJ_model, type=2)

##### Compare full model to the null model
RFadultnaive_blockloco_INJ_nullmodel <- glmmTMB(finalblock_loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","INJ")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_blockloco_INJ_nullmodel, RFadultnaive_blockloco_INJ_model, test = "LRT") # Full model is not a better fit than the null model

##### Follow up with emmeans
RFadultnaive_blockloco_INJ_emmeans <- emmeans(RFadultnaive_blockloco_INJ_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_blockloco_INJ_emmeans)
##### Overall, injections days 1 through 4 are significantly lower than baseline.
##### Comparing across only the injection days, days 1 through 3 are significant or close to significantly lower than days 6 and 7


#### Post-Injection - Early Phase - relative to baseline
RFadultnaive_blockloco_POSTE_model <- glmmTMB(finalblock_loco ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_EARLY")),
                                          family = stats::gaussian())

summary(RFadultnaive_blockloco_POSTE_model)
Anova(RFadultnaive_blockloco_POSTE_model, type=2)

##### Compare full model to the null model
RFadultnaive_blockloco_POSTE_nullmodel <- glmmTMB(finalblock_loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_EARLY")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_blockloco_POSTE_nullmodel, RFadultnaive_blockloco_POSTE_model, test = "LRT") # Full model is not a better fit than the null model

##### Follow up with emmeans
RFadultnaive_blockloco_POSTE_emmeans <- emmeans(RFadultnaive_blockloco_POSTE_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_blockloco_POSTE_emmeans)
##### Overall, only early phase post injection days 1 and 2 are significantly lower than baseline.
##### Comparing across only the post-injection days, day 1 is significantly lower than day 3.


#### Post-Injection - Late Phase - relative to baseline
RFadultnaive_blockloco_POSTL_model <- glmmTMB(finalblock_loco ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_LATE")),
                                          family = stats::gaussian())

summary(RFadultnaive_blockloco_POSTL_model)
Anova(RFadultnaive_blockloco_POSTL_model, type=2)

##### Compare full model to the null model
RFadultnaive_blockloco_POSTL_nullmodel <- glmmTMB(finalblock_loco ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_LATE")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_blockloco_POSTL_nullmodel, RFadultnaive_blockloco_POSTL_model, test = "LRT") # Full model is not a better fit than the null model

##### Follow up with emmeans
RFadultnaive_blockloco_POSTL_emmeans <- emmeans(RFadultnaive_blockloco_POSTL_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_blockloco_POSTL_emmeans)
##### Overall, late phase post injection days are not significantly different from baseline.
##### Comparing across only the post-injection days, there are no consistent significant differences.


### Export analysis tables -----
#### Between Phase Tables
RFadultnaive_blockloco_model_paramtable <- mixedmodel_paramtable(RFadultnaive_blockloco_model)
RFadultnaive_blockloco_model_fittable <- mixedmodel_fittable(RFadultnaive_blockloco_model, RFadultnaive_blockloco_nullmodel)
RFadultnaive_blockloco_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_blockloco_emmeans,RFadultnaive_blockloco_model)

write.csv(RFadultnaive_blockloco_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_blockloco_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_model_fittable.csv",sep=""))
write.csv(RFadultnaive_blockloco_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_emmeans_contraststable.csv",sep=""))

#### Within Phase Tables
##### Baseline
RFadultnaive_blockloco_BL_model_paramtable <- mixedmodel_paramtable(RFadultnaive_blockloco_BL_model)
RFadultnaive_blockloco_BL_model_fittable <- mixedmodel_fittable(RFadultnaive_blockloco_BL_model, RFadultnaive_blockloco_BL_nullmodel)
RFadultnaive_blockloco_BL_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_blockloco_BL_emmeans,RFadultnaive_blockloco_BL_model)

write.csv(RFadultnaive_blockloco_BL_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_BL_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_blockloco_BL_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_BL_model_fittable.csv",sep=""))
write.csv(RFadultnaive_blockloco_BL_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_BL_emmeans_contraststable.csv",sep=""))

##### Injection
RFadultnaive_blockloco_INJ_model_paramtable <- mixedmodel_paramtable(RFadultnaive_blockloco_INJ_model)
RFadultnaive_blockloco_INJ_model_fittable <- mixedmodel_fittable(RFadultnaive_blockloco_INJ_model, RFadultnaive_blockloco_INJ_nullmodel)
RFadultnaive_blockloco_INJ_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_blockloco_INJ_emmeans,RFadultnaive_blockloco_INJ_model)

write.csv(RFadultnaive_blockloco_INJ_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_INJ_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_blockloco_INJ_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_INJ_model_fittable.csv",sep=""))
write.csv(RFadultnaive_blockloco_INJ_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_INJ_emmeans_contraststable.csv",sep=""))

##### Post-Injection Early
RFadultnaive_blockloco_POSTE_model_paramtable <- mixedmodel_paramtable(RFadultnaive_blockloco_POSTE_model)
RFadultnaive_blockloco_POSTE_model_fittable <- mixedmodel_fittable(RFadultnaive_blockloco_POSTE_model, RFadultnaive_blockloco_POSTE_nullmodel)
RFadultnaive_blockloco_POSTE_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_blockloco_POSTE_emmeans,RFadultnaive_blockloco_POSTE_model)

write.csv(RFadultnaive_blockloco_POSTE_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_POSTE_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_blockloco_POSTE_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_POSTE_model_fittable.csv",sep=""))
write.csv(RFadultnaive_blockloco_POSTE_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_POSTE_emmeans_contraststable.csv",sep=""))

##### Post-Injection Late
RFadultnaive_blockloco_POSTL_model_paramtable <- mixedmodel_paramtable(RFadultnaive_blockloco_POSTL_model)
RFadultnaive_blockloco_POSTL_model_fittable <- mixedmodel_fittable(RFadultnaive_blockloco_POSTL_model, RFadultnaive_blockloco_POSTL_nullmodel)
RFadultnaive_blockloco_POSTL_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_blockloco_POSTL_emmeans,RFadultnaive_blockloco_POSTL_model)

write.csv(RFadultnaive_blockloco_POSTL_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_POSTL_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_blockloco_POSTL_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_POSTL_model_fittable.csv",sep=""))
write.csv(RFadultnaive_blockloco_POSTL_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_blockloco_POSTL_emmeans_contraststable.csv",sep=""))

# FIGURES -----
## Example Subject Plots -----
### Prep data
examplesubjectdata <- RFadultnaivedata_raw %>% dplyr::filter(SubjectID == 141) %>% 
  pivot_longer(cols = all_of(freqcolnames), names_to = 'Frequency', values_to='LeverPresses')

examplesubjectdata$Pass <- factor(examplesubjectdata$Pass)

### Prep plotting variables
xorder <- c('28','32','35','40','45','50','56','63','71','79','89','100','112','126','141')
xlabels <- c('28','','35','','45','','56','','71','','89','','112','','141')


RFexample_ymax <- 65
RFexample_ymin <- 0

### Plot Reduced RF Example - Baseline
MBLplot <- ggplot(subset(examplesubjectdata, UniqueSessionID == "BL_3"), aes(x = Frequency, y = LeverPresses, group = Pass)) +
  geom_line(linewidth = mylinewidth, aes(color = Pass), show.legend = FALSE) + scale_x_discrete(limits = xorder, labels = xlabels) + 
  geom_point(aes(x = Frequency, y = LeverPresses, color = Pass), size = mypointsize, show.legend = FALSE) +
  scale_color_manual(values=c("#8d8d8d", "#0000FF", "#950EF2", "#FF0080", "#0FC127")) +
  scale_y_continuous(limits = c(RFexample_ymin, RFexample_ymax), breaks=seq(RFexample_ymin,RFexample_ymax,20), expand=c(0,0)) +
  mytheme +
  ylab("Lever Presses") +
  ggtitle("Baseline")

MBLplot

### Plot Reduced RF Example - Morphine Administration
MAplot <- ggplot(subset(examplesubjectdata, UniqueSessionID == "INJ_7"), aes(x = Frequency, y = LeverPresses, group = Pass)) +
  geom_line(linewidth = mylinewidth, aes(color = Pass), show.legend = FALSE) + scale_x_discrete(limits = xorder, labels = xlabels) + 
  geom_point(aes(x = Frequency, y = LeverPresses, color = Pass), size = mypointsize, show.legend = FALSE) +
  scale_color_manual(values=c("#8d8d8d", "#0000FF", "#950EF2", "#FF0080", "#0FC127")) +
  scale_y_continuous(limits = c(RFexample_ymin, RFexample_ymax), breaks=seq(RFexample_ymin,RFexample_ymax,20), expand=c(0,0)) +
  mytheme +
  ylab("Lever Presses") +
  ggtitle("Morphine (Day 7)")

MAplot

### Plot Reduced RF Example - Morphine Withdrawal
MWplot <- ggplot(subset(examplesubjectdata, UniqueSessionID == "POST_2"), aes(x = Frequency, y = LeverPresses, group = Pass)) +
  geom_line(linewidth = mylinewidth, aes(color = Pass), show.legend = FALSE) + scale_x_discrete(limits = xorder, labels = xlabels) + 
  geom_point(aes(x = Frequency, y = LeverPresses, color = Pass), size = mypointsize, show.legend = FALSE) +
  scale_color_manual(values=c("#8d8d8d", "#0000FF", "#950EF2", "#FF0080", "#0FC127")) + 
  scale_y_continuous(limits = c(RFexample_ymin, RFexample_ymax), breaks=seq(RFexample_ymin,RFexample_ymax,20), expand=c(0,0)) +
  mytheme +
  ylab("Lever Presses") +
  ggtitle("Withdrawal (Day 9)")

MWplot

### Make empty tile for key
keyplot <- ggplot() + mytheme

### Create full panel
exampleRFplots <- ggarrange(MBLplot,MAplot, MWplot,keyplot, ncol=4, nrow=1, widths=c(1,1,1,.8), labels = c("A", "B", "C", ""))
exampleRFplots
ggsave(plot=exampleRFplots, path=path_figures, device="pdf", filename="F3_NaiveAdultRF_MorphineTest_ExampleRFPlots.pdf", width=RFpanelwidth, height=RFpanelheight)


## Figure Panels -----
### Prepare plot colors and variables -----
groupcolors <- c("AdultNaiveM" = AdultNaiveMcolor)
Mphasecolors <- c("BL" = ICSSBLcolor, "INJ"=AdultNaiveMcolor,"POST_EARLY"=AdultNaiveEPcolor,"POST_LATE"=AdultNaiveLPcolor)

sexshapes <- c("M" = 0, "F" = 1)


# For line plots - Set up x axis ticks and labels - label day starting after baseline
sessionticks <- c(5:19)
sessionlabels <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15')

injsession <- c(1,2,3,4,5,6,7) # Sessions to label with injection triangle marker

### M50 -----
ymax_M50C <- 17
ymin_M50C <- -35

ymax_M50_bar <- 155
ymin_M50_bar <- 0

injmarker_M50 <- c(rep(ymin_M50C,7))
injticks_M50 <- data.frame(injsession,injmarker_M50)

#### M50 Change Line
M50C_line <- ggplot(subset(RFadultnaivemeans_lm,Phase!="BL"), aes(x=sessionlabel, y=M50C_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadultnaivemeans_lm,Phase!="BL"), aes(ymin=M50C_mean-M50C_se, ymax=M50C_mean+M50C_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name="Session", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_M50C, ymax_M50C), breaks = seq(-30, ymax_M50C, 15)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("M50", sep = "")) +
  geom_point(data=injticks_M50, aes(x=injsession, y=injmarker_M50, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE)+ # Add injection ticks
  annotate("text", x = 1, y = -27.5, label = "***", size = mystarsize)+
  annotate("text", x = 4, y = -31, label = "***", size = mystarsize)+
  annotate("text", x = 5, y = -25.5, label = "***", size = mystarsize)+
  annotate("text", x = 6, y = -19.5, label = "*", size = mystarsize)+
  annotate("text", x = 7, y = -26.5, label = "***", size = mystarsize)+
  annotate("text", x = 8, y = 14, label = "**", size = mystarsize)+
  annotate("text", x = 10, y = 13, label = "*", size = mystarsize)

M50C_line

#### M50 by Phase Bar
M50_bar <- ggplot(RFadultnaivemeans_Phase_lm, aes(x=Phase, y=M50_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(RFadultnaivemeans_SubjectPhase_lm, mapping = aes(x = Phase, y = M50_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=RFadultnaivemeans_Phase_lm, aes(ymin=M50_mean-M50_se, ymax=M50_mean+M50_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_M50_bar, ymax_M50_bar), breaks = seq(ymin_M50_bar, ymax_M50_bar,50)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-15)')) + 
  scale_color_manual("legend", values = Mphasecolors) + 
  scale_fill_manual("legend", values = Mphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("M50") + xlab('Phase') +
  mytheme +
  ggtitle(paste(" ", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(1.95), y.position=c(118) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.05), xmax = c(3), y.position=c(118) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(148) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.05), xmax = c(4), y.position=c(132) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

M50_bar

#### Create full panel
M50 <- ggarrange(M50_bar,M50C_line, ncol=2, nrow=1, widths=c(2.2,4), labels = c("D", "E"))
M50


### Final Block Locomotion -----
ymax_blocklocoC <- 10
ymin_blocklocoC <- -20

ymax_blockloco_bar <- 79
ymin_blockloco_bar <- 0

injmarker_blockloco <- c(rep(ymin_blocklocoC,7))
injticks_blockloco <- data.frame(injsession,injmarker_blockloco)

#### blockloco Change Line
blocklocoC_line <- ggplot(subset(RFadultnaivemeans_loco,Phase!="BL"), aes(x=sessionlabel, y=finalblock_locoC_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadultnaivemeans_loco,Phase!="BL"), aes(ymin=finalblock_locoC_mean-finalblock_locoC_se, ymax=finalblock_locoC_mean+finalblock_locoC_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name="Session", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blocklocoC, ymax_blocklocoC), breaks = seq(ymin_blocklocoC, ymax_blocklocoC, 10)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("Last Block Locomotion", sep = "")) +
  geom_point(data=injticks_blockloco, aes(x=injsession, y=injmarker_blockloco, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE)+ # Add injection ticks
  annotate("text", x = 1, y = -16.5, label = "***", size = mystarsize)+
  annotate("text", x = 2, y = -15, label = "**", size = mystarsize)+
  annotate("text", x = 3, y = -14.5, label = "**", size = mystarsize)+
  annotate("text", x = 4, y = -13, label = "*", size = mystarsize)+
  annotate("text", x = 8, y = -16.5, label = "***", size = mystarsize)+
  annotate("text", x = 9, y = -11, label = "*", size = mystarsize)
    
blocklocoC_line

#### blockloco by Phase Bar
blockloco_bar <- ggplot(RFadultnaivemeans_Phase_loco, aes(x=Phase, y=finalblock_loco_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(RFadultnaivemeans_SubjectPhase_loco, mapping = aes(x = Phase, y = finalblock_loco_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=RFadultnaivemeans_Phase_loco, aes(ymin=finalblock_loco_mean-finalblock_loco_se, ymax=finalblock_loco_mean+finalblock_loco_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_blockloco_bar, ymax_blockloco_bar), breaks = seq(ymin_blockloco_bar, ymax_blockloco_bar,25)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-15)')) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Block Locomotion") + xlab('Phase') +
  mytheme +
  ggtitle(paste(" ", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(55) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(64) ,label.size = mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(73) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

blockloco_bar

#### Create full panel
blockloco <- ggarrange(blockloco_bar, blocklocoC_line, ncol=2, nrow=1, widths=c(2.2,4), labels = c("F", "G"))
blockloco


## Figure 3 -----
### Create full panel
figure3 <- ggarrange(exampleRFplots,M50,blockloco, ncol=1, nrow=3)
figure3
ggsave(plot=figure3, path=path_figures, device="pdf", filename="F3_AllPanels.pdf", width=RFpanelwidth, height=RFpanelheight*3)


# SUPPLEMENTARY -----
## STATISTICS -----
## Max Lever Presses -----
### Between Phase Analysis -----
#### Mixed effects model with random intercepts by subject and sex
RFadultnaive_maxLP_model <- glmmTMB(pass_maxLP ~ Phase + Sex + (1|SubjectID), data = RFadultnaivedata_loco, family=stats::gaussian)

summary(RFadultnaive_maxLP_model)
Anova(RFadultnaive_maxLP_model, type=2)

#### Compare full model to the null model
RFadultnaive_maxLP_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = RFadultnaivedata_loco, family = stats::gaussian()) # Make the null model
anova(RFadultnaive_maxLP_nullmodel, RFadultnaive_maxLP_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultnaive_maxLP_emmeans <-emmeans(RFadultnaive_maxLP_model, pairwise ~ Phase, adjust = 'none')
summary(RFadultnaive_maxLP_emmeans)

### Within Phase Analysis -----
#### Baseline
##### Mixed effects model by day with random intercept by subject
RFadultnaive_maxLP_BL_model <- glmmTMB(pass_maxLP ~ Session + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase=="BL"),
                                       family = stats::gaussian())

summary(RFadultnaive_maxLP_BL_model)
Anova(RFadultnaive_maxLP_BL_model, type=2)

##### Compare full model to the null model
RFadultnaive_maxLP_BL_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_maxLP_BL_nullmodel, RFadultnaive_maxLP_BL_model, test = "LRT") # Full model is not a better fit than the null model

##### Follow up with emmeans
RFadultnaive_maxLP_BL_emmeans <- emmeans(RFadultnaive_maxLP_BL_model, pairwise ~ Session, adjust = 'none')
summary(RFadultnaive_maxLP_BL_emmeans)
##### Overall, there is no significant effect of day in the baseline phase.


#### Injection - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultnaive_maxLP_INJ_model <- glmmTMB(pass_maxLP ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","INJ")),
                                        family = stats::gaussian())

summary(RFadultnaive_maxLP_INJ_model)
Anova(RFadultnaive_maxLP_INJ_model, type=2)

##### Compare full model to the null model
RFadultnaive_maxLP_INJ_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","INJ")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_maxLP_INJ_nullmodel, RFadultnaive_maxLP_INJ_model, test = "LRT") # Full model is a better fit than the null model

##### Follow up with emmeans
RFadultnaive_maxLP_INJ_emmeans <- emmeans(RFadultnaive_maxLP_INJ_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_maxLP_INJ_emmeans)
##### Overall, all injections days are significantly lower than from baseline.
##### Comparing across only the injection days, there are no significant differences.


#### Post-Injection - Early Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultnaive_maxLP_POSTE_model <- glmmTMB(pass_maxLP ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_EARLY")),
                                          family = stats::gaussian())

summary(RFadultnaive_maxLP_POSTE_model)
Anova(RFadultnaive_maxLP_POSTE_model, type=2)

##### Compare full model to the null model
RFadultnaive_maxLP_POSTE_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_EARLY")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_maxLP_POSTE_nullmodel, RFadultnaive_maxLP_POSTE_model, test = "LRT") # Full model is a better fit than the null model

##### Follow up with emmeans
RFadultnaive_maxLP_POSTE_emmeans <- emmeans(RFadultnaive_maxLP_POSTE_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_maxLP_POSTE_emmeans)
##### Overall, only early phase post injection day 1 is significantly different from baseline.
##### Comparing across only the post-injection days, there are no consistent significant differences.


#### Post-Injection - Late Phase - relative to baseline
##### Mixed effects model by day with random intercept by subject
RFadultnaive_maxLP_POSTL_model <- glmmTMB(pass_maxLP ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_LATE")),
                                          family = stats::gaussian())

summary(RFadultnaive_maxLP_POSTL_model)
Anova(RFadultnaive_maxLP_POSTL_model, type=2)

##### Compare full model to the null model
RFadultnaive_maxLP_POSTL_nullmodel <- glmmTMB(pass_maxLP ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_LATE")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_maxLP_POSTL_nullmodel, RFadultnaive_maxLP_POSTL_model, test = "LRT") # Full model is not a better fit than the null model

##### Follow up with emmeans
RFadultnaive_maxLP_POSTL_emmeans <- emmeans(RFadultnaive_maxLP_POSTL_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_maxLP_POSTL_emmeans)
##### Overall, late phase post injection days 4 and 5 are significantly lower than baseline.
##### Comparing across only the post-injection days, there are no consistent significant differences.

### Export analysis tables -----
#### Between Phase Tables
RFadultnaive_maxLP_model_paramtable <- mixedmodel_paramtable(RFadultnaive_maxLP_model)
RFadultnaive_maxLP_model_fittable <- mixedmodel_fittable(RFadultnaive_maxLP_model, RFadultnaive_maxLP_nullmodel)
RFadultnaive_maxLP_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_maxLP_emmeans,RFadultnaive_maxLP_model)

write.csv(RFadultnaive_maxLP_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_maxLP_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_model_fittable.csv",sep=""))
write.csv(RFadultnaive_maxLP_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_emmeans_contraststable.csv",sep=""))

#### Within Phase Tables
##### Baseline
RFadultnaive_maxLP_BL_model_paramtable <- mixedmodel_paramtable(RFadultnaive_maxLP_BL_model)
RFadultnaive_maxLP_BL_model_fittable <- mixedmodel_fittable(RFadultnaive_maxLP_BL_model, RFadultnaive_maxLP_BL_nullmodel)
RFadultnaive_maxLP_BL_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_maxLP_BL_emmeans,RFadultnaive_maxLP_BL_model)

write.csv(RFadultnaive_maxLP_BL_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_BL_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_maxLP_BL_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_BL_model_fittable.csv",sep=""))
write.csv(RFadultnaive_maxLP_BL_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_BL_emmeans_contraststable.csv",sep=""))

##### Injection
RFadultnaive_maxLP_INJ_model_paramtable <- mixedmodel_paramtable(RFadultnaive_maxLP_INJ_model)
RFadultnaive_maxLP_INJ_model_fittable <- mixedmodel_fittable(RFadultnaive_maxLP_INJ_model, RFadultnaive_maxLP_INJ_nullmodel)
RFadultnaive_maxLP_INJ_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_maxLP_INJ_emmeans,RFadultnaive_maxLP_INJ_model)

write.csv(RFadultnaive_maxLP_INJ_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_INJ_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_maxLP_INJ_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_INJ_model_fittable.csv",sep=""))
write.csv(RFadultnaive_maxLP_INJ_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_INJ_emmeans_contraststable.csv",sep=""))

##### Post-Injection Early
RFadultnaive_maxLP_POSTE_model_paramtable <- mixedmodel_paramtable(RFadultnaive_maxLP_POSTE_model)
RFadultnaive_maxLP_POSTE_model_fittable <- mixedmodel_fittable(RFadultnaive_maxLP_POSTE_model, RFadultnaive_maxLP_POSTE_nullmodel)
RFadultnaive_maxLP_POSTE_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_maxLP_POSTE_emmeans,RFadultnaive_maxLP_POSTE_model)

write.csv(RFadultnaive_maxLP_POSTE_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_POSTE_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_maxLP_POSTE_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_POSTE_model_fittable.csv",sep=""))
write.csv(RFadultnaive_maxLP_POSTE_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_POSTE_emmeans_contraststable.csv",sep=""))

##### Post-Injection Late
RFadultnaive_maxLP_POSTL_model_paramtable <- mixedmodel_paramtable(RFadultnaive_maxLP_POSTL_model)
RFadultnaive_maxLP_POSTL_model_fittable <- mixedmodel_fittable(RFadultnaive_maxLP_POSTL_model, RFadultnaive_maxLP_POSTL_nullmodel)
RFadultnaive_maxLP_POSTL_emmeans_contraststable <- mixedmodel_emmeanstable(RFadultnaive_maxLP_POSTL_emmeans,RFadultnaive_maxLP_POSTL_model)

write.csv(RFadultnaive_maxLP_POSTL_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_POSTL_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_maxLP_POSTL_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_POSTL_model_fittable.csv",sep=""))
write.csv(RFadultnaive_maxLP_POSTL_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_maxLP_POSTL_emmeans_contraststable.csv",sep=""))



## Total Session Lever Presses -----
### Scale the data -----
#### NOTE: session_sumLP must be scaled to prevent convergence issues in the later within phase day models. 
#### All models are fit to the scaled and centered sumLP data, and tables will be back-transformed for interpretation at the end.
RFadultnaivedata_loco$session_sumLP_scaled <- scale(RFadultnaivedata_loco$session_sumLP)

# Save scaling parameters for back-transformation
mean_session_sumLP <- mean(RFadultnaivedata_loco$session_sumLP)
sd_session_sumLP <- sd(RFadultnaivedata_loco$session_sumLP)

### Between Phase Analysis -----
#### Mixed effects model with random intercepts by subject and sex
RFadultnaive_sumLP_model <- glmmTMB(session_sumLP_scaled ~ Phase + Sex + (1|SubjectID), data = RFadultnaivedata_loco, family=stats::gaussian)

summary(RFadultnaive_sumLP_model)
Anova(RFadultnaive_sumLP_model, type=2)

# Convert random-effect variance and SD back to raw scale
sumLP_var_scaled <- 0.5244
sumLP_sd_scaled  <- sqrt(sumLP_var_scaled)

sumLP_var_raw <- sumLP_var_scaled * (sd_session_sumLP^2)
sumLP_sd_raw  <- sumLP_sd_scaled  * sd_session_sumLP

sumLP_var_raw
sumLP_sd_raw

#### Compare full model to the null model
RFadultnaive_sumLP_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = RFadultnaivedata_loco, family = stats::gaussian()) # Make the null model
anova(RFadultnaive_sumLP_nullmodel, RFadultnaive_sumLP_model, test = "LRT") # Full model is a better fit than the null model

#### Post hoc tests with emmeans and pairs
RFadultnaive_sumLP_emmeans <-emmeans(RFadultnaive_sumLP_model, pairwise ~ Phase, adjust = 'none')
summary(RFadultnaive_sumLP_emmeans)

### Within Phase Analysis -----
#### Baseline
##### Find SD and means for later re-scaling
sumLP_BLdata <- subset(RFadultnaivedata_loco, Phase=="BL")
mean_session_sumLP_BL <- mean(sumLP_BLdata$session_sumLP)
sd_session_sumLP_BL <- sd(sumLP_BLdata$session_sumLP)

##### Mixed effects model by day with random intercept by subject
RFadultnaive_sumLP_BL_model <- glmmTMB(session_sumLP_scaled ~ Session + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase=="BL"),
                                       family = stats::gaussian())

summary(RFadultnaive_sumLP_BL_model)
Anova(RFadultnaive_sumLP_BL_model, type=2)

# Convert random-effect variance and SD back to raw scale
sumLP_BL_var_scaled <- 0.7087
sumLP_BL_sd_scaled  <- sqrt(sumLP_BL_var_scaled)

sumLP_BL_var_raw <- sumLP_BL_var_scaled * (sd_session_sumLP_BL^2)
sumLP_BL_sd_raw  <- sumLP_BL_sd_scaled  * sd_session_sumLP_BL

sumLP_BL_var_raw
sumLP_BL_sd_raw


##### Compare full model to the null model
RFadultnaive_sumLP_BL_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1  + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase=="BL"), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_sumLP_BL_nullmodel, RFadultnaive_sumLP_BL_model, test = "LRT") # Full model is not a better fit than the null model

##### Follow up with emmeans
RFadultnaive_sumLP_BL_emmeans <- emmeans(RFadultnaive_sumLP_BL_model, pairwise ~ Session, adjust = 'none')
summary(RFadultnaive_sumLP_BL_emmeans)
##### Overall, there is no significant effect of day in the baseline phase.


#### Injection - relative to baseline
##### Find SD and means for later re-scaling
sumLP_INJdata <- subset(RFadultnaivedata_loco, Phase %in% c("BL","INJ"))
mean_session_sumLP_INJ <- mean(sumLP_INJdata$session_sumLP)
sd_session_sumLP_INJ <- sd(sumLP_INJdata$session_sumLP)

##### Mixed effects model by day with random intercept by subject
RFadultnaive_sumLP_INJ_model <- glmmTMB(session_sumLP_scaled ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","INJ")),
                                        family = stats::gaussian())

summary(RFadultnaive_sumLP_INJ_model)
Anova(RFadultnaive_sumLP_INJ_model, type=2)

# Convert random-effect variance and SD back to raw scale
sumLP_INJ_var_scaled <- 0.7087
sumLP_INJ_sd_scaled  <- sqrt(sumLP_INJ_var_scaled)

sumLP_INJ_var_raw <- sumLP_INJ_var_scaled * (sd_session_sumLP_INJ^2)
sumLP_INJ_sd_raw  <- sumLP_INJ_sd_scaled  * sd_session_sumLP_INJ

sumLP_INJ_var_raw
sumLP_INJ_sd_raw


##### Compare full model to the null model
RFadultnaive_sumLP_INJ_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","INJ")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_sumLP_INJ_nullmodel, RFadultnaive_sumLP_INJ_model, test = "LRT") # Full model is not a better fit than the null model

##### Follow up with emmeans
RFadultnaive_sumLP_INJ_emmeans <- emmeans(RFadultnaive_sumLP_INJ_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_sumLP_INJ_emmeans)
##### Overall, injections days 1, 2, and 3 are significantly lower than baseline.
##### Comparing across only the injection days, there are no significant differences.

#### Post-Injection - Early Phase - relative to baseline
##### Find SD and means for later re-scaling
sumLP_POSTEdata <- subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_EARLY"))
mean_session_sumLP_POSTE <- mean(sumLP_POSTEdata$session_sumLP)
sd_session_sumLP_POSTE <- sd(sumLP_POSTEdata$session_sumLP)

##### Mixed effects model by day with random intercept by subject
RFadultnaive_sumLP_POSTE_model <- glmmTMB(session_sumLP_scaled ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_EARLY")),
                                          family = stats::gaussian())

summary(RFadultnaive_sumLP_POSTE_model)
Anova(RFadultnaive_sumLP_POSTE_model, type=2)

# Convert random-effect variance and SD back to raw scale
sumLP_POSTE_var_scaled <- 0.7076
sumLP_POSTE_sd_scaled  <- sqrt(sumLP_POSTE_var_scaled)

sumLP_POSTE_var_raw <- sumLP_POSTE_var_scaled * (sd_session_sumLP_POSTE^2)
sumLP_POSTE_sd_raw  <- sumLP_POSTE_sd_scaled  * sd_session_sumLP_POSTE

sumLP_POSTE_var_raw
sumLP_POSTE_sd_raw


##### Compare full model to the null model
RFadultnaive_sumLP_POSTE_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_EARLY")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_sumLP_POSTE_nullmodel, RFadultnaive_sumLP_POSTE_model, test = "LRT") # Full model is not a better fit than the null model

##### Follow up with emmeans
RFadultnaive_sumLP_POSTE_emmeans <- emmeans(RFadultnaive_sumLP_POSTE_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_sumLP_POSTE_emmeans)
##### Overall, all early phase post injection days are significantly lower than baseline.


#### Post-Injection - Late Phase - relative to baseline
##### Find SD and means for later re-scaling
sumLP_POSTLdata <- subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_LATE"))
mean_session_sumLP_POSTL <- mean(sumLP_POSTLdata$session_sumLP)
sd_session_sumLP_POSTL <- sd(sumLP_POSTLdata$session_sumLP)

##### Mixed effects model by day with random intercept by subject
RFadultnaive_sumLP_POSTL_model <- glmmTMB(session_sumLP_scaled ~ Session_collapsed + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_LATE")),
                                          family = stats::gaussian())

summary(RFadultnaive_sumLP_POSTL_model)
Anova(RFadultnaive_sumLP_POSTL_model, type=2)

# Convert random-effect variance and SD back to raw scale
sumLP_POSTL_var_scaled <- 0.8076
sumLP_POSTL_sd_scaled  <- sqrt(sumLP_POSTL_var_scaled)

sumLP_POSTL_var_raw <- sumLP_POSTL_var_scaled * (sd_session_sumLP_POSTL^2)
sumLP_POSTL_sd_raw  <- sumLP_POSTL_sd_scaled  * sd_session_sumLP_POSTL

sumLP_POSTL_var_raw
sumLP_POSTL_sd_raw

##### Compare full model to the null model
RFadultnaive_sumLP_POSTL_nullmodel <- glmmTMB(session_sumLP_scaled ~ 1 + Sex + (1|SubjectID), data = subset(RFadultnaivedata_loco, Phase %in% c("BL","POST_LATE")), family = stats::gaussian()) # Make the null model
anova(RFadultnaive_sumLP_POSTL_nullmodel, RFadultnaive_sumLP_POSTL_model, test = "LRT") # Full model is not a better fit than the null model

##### Follow up with emmeans
RFadultnaive_sumLP_POSTL_emmeans <- emmeans(RFadultnaive_sumLP_POSTL_model, pairwise ~ Session_collapsed, adjust = 'none')
summary(RFadultnaive_sumLP_POSTL_emmeans)
##### Overall, late phase post injection days 4, 5, and 7 are significantly lower than baseline.
##### Comparing across only the post-injection days, there are no consistent significant differences.


### Export analysis tables -----
#### Between Phase Tables
RFadultnaive_sumLP_model_paramtable <- mixedmodel_rescaleparamtable((mixedmodel_paramtable(RFadultnaive_sumLP_model)),sd_session_sumLP,mean_session_sumLP)
RFadultnaive_sumLP_model_fittable <- mixedmodel_fittable(RFadultnaive_sumLP_model, RFadultnaive_sumLP_nullmodel)
RFadultnaive_sumLP_emmeans_contraststable <- mixedmodel_rescalecontrasts((mixedmodel_emmeanstable(RFadultnaive_sumLP_emmeans,RFadultnaive_sumLP_model)),sd_session_sumLP)

write.csv(RFadultnaive_sumLP_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_sumLP_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_model_fittable.csv",sep=""))
write.csv(RFadultnaive_sumLP_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_emmeans_contraststable.csv",sep=""))

#### Within Phase Tables
##### Baseline
RFadultnaive_sumLP_BL_model_paramtable <- mixedmodel_rescaleparamtable((mixedmodel_paramtable(RFadultnaive_sumLP_BL_model)),sd_session_sumLP,mean_session_sumLP)
RFadultnaive_sumLP_BL_model_fittable <- mixedmodel_fittable(RFadultnaive_sumLP_BL_model, RFadultnaive_sumLP_BL_nullmodel)
RFadultnaive_sumLP_BL_emmeans_contraststable <- mixedmodel_rescalecontrasts((mixedmodel_emmeanstable(RFadultnaive_sumLP_BL_emmeans,RFadultnaive_sumLP_BL_model)),sd_session_sumLP)

write.csv(RFadultnaive_sumLP_BL_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_BL_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_sumLP_BL_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_BL_model_fittable.csv",sep=""))
write.csv(RFadultnaive_sumLP_BL_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_BL_emmeans_contraststable.csv",sep=""))

##### Injection
RFadultnaive_sumLP_INJ_model_paramtable <- mixedmodel_rescaleparamtable((mixedmodel_paramtable(RFadultnaive_sumLP_INJ_model)),sd_session_sumLP,mean_session_sumLP)
RFadultnaive_sumLP_INJ_model_fittable <- mixedmodel_fittable(RFadultnaive_sumLP_INJ_model, RFadultnaive_sumLP_INJ_nullmodel)
RFadultnaive_sumLP_INJ_emmeans_contraststable <- mixedmodel_rescalecontrasts((mixedmodel_emmeanstable(RFadultnaive_sumLP_INJ_emmeans,RFadultnaive_sumLP_INJ_model)),sd_session_sumLP)

write.csv(RFadultnaive_sumLP_INJ_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_INJ_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_sumLP_INJ_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_INJ_model_fittable.csv",sep=""))
write.csv(RFadultnaive_sumLP_INJ_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_INJ_emmeans_contraststable.csv",sep=""))

##### Post-Injection Early
RFadultnaive_sumLP_POSTE_model_paramtable <- mixedmodel_rescaleparamtable((mixedmodel_paramtable(RFadultnaive_sumLP_POSTE_model)),sd_session_sumLP,mean_session_sumLP)
RFadultnaive_sumLP_POSTE_model_fittable <- mixedmodel_fittable(RFadultnaive_sumLP_POSTE_model, RFadultnaive_sumLP_POSTE_nullmodel)
RFadultnaive_sumLP_POSTE_emmeans_contraststable <- mixedmodel_rescalecontrasts((mixedmodel_emmeanstable(RFadultnaive_sumLP_POSTE_emmeans,RFadultnaive_sumLP_POSTE_model)),sd_session_sumLP)

write.csv(RFadultnaive_sumLP_POSTE_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_POSTE_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_sumLP_POSTE_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_POSTE_model_fittable.csv",sep=""))
write.csv(RFadultnaive_sumLP_POSTE_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_POSTE_emmeans_contraststable.csv",sep=""))

##### Post-Injection Late
RFadultnaive_sumLP_POSTL_model_paramtable <- mixedmodel_rescaleparamtable((mixedmodel_paramtable(RFadultnaive_sumLP_POSTL_model)),sd_session_sumLP,mean_session_sumLP)
RFadultnaive_sumLP_POSTL_model_fittable <- mixedmodel_fittable(RFadultnaive_sumLP_POSTL_model, RFadultnaive_sumLP_POSTL_nullmodel)
RFadultnaive_sumLP_POSTL_emmeans_contraststable <- mixedmodel_rescalecontrasts((mixedmodel_emmeanstable(RFadultnaive_sumLP_POSTL_emmeans,RFadultnaive_sumLP_POSTL_model)),sd_session_sumLP)

write.csv(RFadultnaive_sumLP_POSTL_model_paramtable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_POSTL_model_paramtable.csv",sep=""))
write.csv(RFadultnaive_sumLP_POSTL_model_fittable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_POSTL_model_fittable.csv",sep=""))
write.csv(RFadultnaive_sumLP_POSTL_emmeans_contraststable, paste(path_analysis_RFtest,"RFadultnaive_sumLP_POSTL_emmeans_contraststable.csv",sep=""))

## FIGURES -----
### Figure panels -----


##### Max Lever Presses -----
ymax_maxLPC <- 4
ymin_maxLPC <- -32

ymax_maxLP_bar <- 110
ymin_maxLP_bar <- 0
injmarker_maxLP <- c(rep(ymin_maxLPC,7))
injticks_maxLP <- data.frame(injsession,injmarker_maxLP)

### maxLP Change Line
maxLPC_line <- ggplot(subset(RFadultnaivemeans_loco,Phase!="BL"), aes(x=sessionlabel, y=pass_maxLPC_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadultnaivemeans_loco,Phase!="BL"), aes(ymin=pass_maxLPC_mean-pass_maxLPC_se, ymax=pass_maxLPC_mean+pass_maxLPC_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name="Session", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLPC, ymax_maxLPC), breaks = seq(-30, ymax_maxLPC, 10)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("Max Lever Presses", sep = "")) +
  geom_point(data=injticks_maxLP, aes(x=injsession, y=injmarker_maxLP, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE)+ # Add injection ticks
  annotate("text", x = 1, y = -25, label = "***", size = mystarsize) +
  annotate("text", x = 2, y = -28, label = "***", size = mystarsize)+
  annotate("text", x = 3, y = -27, label = "***", size = mystarsize)+
  annotate("text", x = 4, y = -25.5, label = "***", size = mystarsize)+
  annotate("text", x = 5, y = -23.5, label = "***", size = mystarsize)+
  annotate("text", x = 6, y = -24, label = "***", size = mystarsize)+
  annotate("text", x = 7, y = -23, label = "***", size = mystarsize)+
  annotate("text", x = 8, y = -17, label = "**", size = mystarsize)

maxLPC_line

#### maxLP by Phase Bar
maxLP_bar <- ggplot(RFadultnaivemeans_Phase_loco, aes(x=Phase, y=pass_maxLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(RFadultnaivemeans_SubjectPhase_loco, mapping = aes(x = Phase, y = pass_maxLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=RFadultnaivemeans_Phase_loco, aes(ymin=pass_maxLP_mean-pass_maxLP_se, ymax=pass_maxLP_mean+pass_maxLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_maxLP_bar, ymax_maxLP_bar), breaks = seq(ymin_maxLP_bar, ymax_maxLP_bar,30)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-15)')) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Max Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste(" ", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(1.95), y.position=c(80) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(94) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(106) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.05), xmax = c(4), y.position=c(80) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.05), xmax = c(3), y.position=c(67) ,label.size = mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

maxLP_bar

#### Create full panel
maxLP <- ggarrange(maxLP_bar, maxLPC_line, ncol=2, nrow=1, widths=c(2.2,4), labels = c("A", "B"))
maxLP


#### Total Session Lever Presses -----
ymax_sumLPC <- 100
ymin_sumLPC <- -575

ymax_sumLP_bar <- 3400
ymin_sumLP_bar <- 0

injmarker_sumLP <- c(rep(ymin_sumLPC,7))
injticks_sumLP <- data.frame(injsession,injmarker_sumLP)

#### sumLP Change Line
sumLPC_line <- ggplot(subset(RFadultnaivemeans_loco,Phase!="BL"), aes(x=sessionlabel, y=session_sumLPC_mean, group=Group, fill=Group, color=Group)) +
  geom_hline(yintercept=0, linetype='dotted', col = 'grey', linewidth = mylinewidth) +
  geom_errorbar(data=subset(RFadultnaivemeans_loco,Phase!="BL"), aes(ymin=session_sumLPC_mean-session_sumLPC_se, ymax=session_sumLPC_mean+session_sumLPC_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black") +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_discrete(name="Session", breaks=sessionticks, labels=sessionlabels) + 
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sumLPC, ymax_sumLPC), breaks = seq(-450, ymax_sumLPC, 150)) + 
  scale_color_manual("legend", values = groupcolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Change from Baseline") + 
  mytheme +
  ggtitle(paste("Total Session Lever Presses", sep = "")) +
  geom_point(data=injticks_sumLP, aes(x=injsession, y=injmarker_sumLP, group=NULL, fill=NULL, color=NULL), size = myinjtrianglesize, shape=17, show.legend = FALSE) + # Add injection ticks
  annotate("text", x = 1, y = -460, label = "*", size = mystarsize)+
  annotate("text", x = 2, y = -490, label = "**", size = mystarsize)+
  annotate("text", x = 3, y = -515, label = "**", size = mystarsize)+
  annotate("text", x = 8, y = -475, label = "***", size = mystarsize)+
  annotate("text", x = 9, y = -335, label = "*", size = mystarsize)+
  annotate("text", x = 10, y = -340, label = "**", size = mystarsize)+
  annotate("text", x = 11, y = -320, label = "*", size = mystarsize)+
  annotate("text", x = 12, y = -300, label = "*", size = mystarsize)+
  annotate("text", x = 14, y = -330, label = "*", size = mystarsize)+
  annotate("text", x = 15, y = -310, label = "*", size = mystarsize)

sumLPC_line

#### sumLP by Phase Bar
sumLP_bar <- ggplot(RFadultnaivemeans_Phase_loco, aes(x=Phase, y=session_sumLP_mean, fill=Phase, color=Phase)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_point(RFadultnaivemeans_SubjectPhase_loco, mapping = aes(x = Phase, y = session_sumLP_mean, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") +
  scale_shape_manual(values=sexshapes)+
  geom_errorbar(data=RFadultnaivemeans_Phase_loco, aes(ymin=session_sumLP_mean-session_sumLP_se, ymax=session_sumLP_mean+session_sumLP_se), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin_sumLP_bar, ymax_sumLP_bar), breaks = seq(ymin_sumLP_bar, ymax_sumLP_bar,1000)) + 
  scale_x_discrete(name ="Phase", labels=c('BL','INJ','POST (Day 8-10)','POST (Day 11-15)')) + 
  scale_color_manual("legend", values = Mphasecolors) +
  scale_fill_manual("legend", values = Mphasecolors) +
  theme_prism(base_size = mybasesize, base_fontface = "plain") +
  ylab("Session Lever Presses") + xlab('Phase') +
  mytheme +
  ggtitle(paste(" ", sep = "")) +
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(2450) ,label.size = mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(2800) ,label.size = mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(3150) ,label.size = mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = .5)

sumLP_bar

#### Create full panel
sumLP <- ggarrange(sumLP_bar, sumLPC_line, ncol=2, nrow=1, widths=c(2.2,4), labels = c("C", "D"))
sumLP


## Supplementary Figure 2 -----
### Create full panel
sfigure2 <- ggarrange(maxLP,sumLP, ncol=1, nrow=2)
sfigure2
ggsave(plot=sfigure2, path=path_figures, device="pdf", filename="SupplementaryF2_AllPanels.pdf", width=RFpanelwidth, height=RFpanelheight*2)
