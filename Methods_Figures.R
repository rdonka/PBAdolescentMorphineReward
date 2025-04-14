############### NAIVE ADULT MORPHINE ADMINISTRATION ANALYSIS ###################
# This program takes in data processed by the script "RF Theta Calculation"    
# and run statistics to analyze findings. Data from each experiment are loaded,
# analyzed, and visualized.                                                    

# PREPARE R -----
## Set up paths
computerusergithubpath <- c('C:/Users/rmdon/OneDrive/Desktop/GitHub_Repositories/') # Path for Rachel's laptop to github repositories
computeruserboxpath <- c('C:/Users/rmdon/Box/') # Path for Rachel's laptop to Box drive

## Load required libraries
source(paste(computerusergithubpath,"JRoitmanICSS/Functions/Functions_LoadPackages.R",sep=''))

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
path_figures_RFtest <- paste(path_analysis_RFtest,"Figures/",sep="")
path_subjectfigures_RFtest <- paste(path_figures_RFtest,"Subject Figures/",sep="") # Subject shaping figures path




# STANDARD RF FIGURES -----
RFblocks <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)

RFfreqs_char <- c("141","126","112","100","89", "79","71", "63","56", "50","45", "40","35", "32","28")
RFfreqs_num <- as.numeric(RFfreqs_char)

RFfreqtrial <- data.frame(Block = RFblocks, Frequency = RFfreqs_num)

## Figure Panels -----

### Prepare plotting variables -----
FreqTriallabels <- c(1,5,10,15)

### Frequency by Trial -----
RF_FreqTrial <- ggplot(RFfreqtrial, aes(x=Block, y=Frequency)) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize*2, show.legend = FALSE) +
  scale_x_reverse(breaks = FreqTriallabels) + 
  scale_y_continuous(limits = c(0,150), breaks = c(0,35,70,105,141)) + 
  mytheme +
  ggtitle(paste("Standard RF Pass Design", sep = ""))

RF_FreqTrial


# ABRIDGED RF FIGURES -----
aRFblocks <- c(1,2,3,4,5,6,7)
aRFfreqs_char <- c("141","100","89", "79","71", "63","50")
aRFfreqs_num <- as.numeric(aRFfreqs_char)

aRFfreqtrial <- data.frame(Block = aRFblocks, Frequency = aRFfreqs_num)

## Figure Panels -----

### Prepare plotting variables -----
aFreqTriallabels <- c(1,3,5,7)

### Frequency by Trial -----
aRF_FreqTrial <- ggplot(aRFfreqtrial, aes(x=Block, y=Frequency)) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize*2, show.legend = FALSE) +
  scale_x_reverse(breaks = aFreqTriallabels) + 
  scale_y_continuous(limits = c(0,150), breaks = c(0,35,70,105,141)) + 
  mytheme +
  ggtitle(paste("Abridged RF Pass Design", sep = ""))

aRF_FreqTrial

## Figure 2 -----
# Create full panel
figure2 <- ggarrange(RF_FreqTrial,aRF_FreqTrial,ncol=1, nrow=2, labels = c("B", "D"))
figure2
ggsave(plot=figure2, path=path_figures, device="pdf", filename="F2_AllPanels.pdf", width=3.5, height=RFpanelheight*3)




