############### ADOLESCENT MORPHINE PATCH EXCITABILITY ANALYSIS ################

# Prepare R for Analysis -----
## Set up paths
computerusergithubpath <- c('C:/Users/rmdon/OneDrive/Desktop/GitHub_Repositories/') # Path for Rachel's laptop to github repositories
computeruserboxpath <- c('C:/Users/rmdon/Box/') # Path for Rachel's laptop to Box drive
rawdatapath <- c('JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Raw Data/')

## Load functions
source(paste(computerusergithubpath, '/Analysis-PBAdolescentMorphineReward/Functions/Functions_LoadPackages.R', sep=''))
source(paste(computerusergithubpath,"/Analysis-PBAdolescentMorphineReward/Functions/Functions_SetPlotVariables.R",sep=''))
library(effsize)
##### Set up paths to export figures and analysis to
# Overall paths
path_analysis <- c(paste(computeruserboxpath, 'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Analysis/',sep=''))
path_figures <- c(paste(computeruserboxpath, 'JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Figures/Figure Panels/',sep=''))


# PREPARE DATA -----
## Read in data and prepare variables -----
subjectkey <- read_csv(paste(computeruserboxpath,rawdatapath,'SubjectKey_AdolescentMorphineWithdrawalPatch.csv',sep=''))
nonspikedata_raw <- read_csv(paste(computeruserboxpath,rawdatapath,'RawData_AdolescentMorphineWithdrawalPatch_NonSpikeExcitabilityMeasurements.csv',sep=''))
spikedata_raw <- read_csv(paste(computeruserboxpath,rawdatapath,'RawData_AdolescentMorphineWithdrawalPatch_Spikes.csv',sep=''))
rebounddata_raw <- read_csv(paste(computeruserboxpath,rawdatapath,'RawData_AdolescentMorphineWithdrawalPatch_ReboundSpike.csv',sep=''))

# Combine subject key
nonspikedata <- left_join(nonspikedata_raw, subjectkey, by=c('ID_Subject'))
spikedata <- left_join(spikedata_raw, subjectkey, by=c('ID_Subject'))
rebounddata <- left_join(rebounddata_raw, subjectkey, by=c('ID_Subject'))

# Factor group
nonspikedata$Group <- factor(nonspikedata$Group, levels=c('Saline','Morphine'))
spikedata$Group <- factor(spikedata$Group, levels=c('Saline','Morphine'))
rebounddata$Group <- factor(rebounddata$Group, levels=c('Saline','Morphine'))

nonspikedata$UniqueID <- paste(nonspikedata$ID_Subject, nonspikedata$ID_Cell,sep='_')
spikedata$UniqueID <- paste(spikedata$ID_Subject, spikedata$ID_Cell,sep='_')
rebounddata$UniqueID <- paste(rebounddata$ID_Subject, rebounddata$ID_Cell,sep='_')

# PREPARE MEANS -----
nonspikemeans <- nonspikedata %>%
  group_by(Group, Variable) %>% 
  summarise(n = n_distinct(ID_Subject),
            valuemean = mean(Value, na.rm=TRUE),
            valuesd = sd(Value,na.rm=TRUE),             
            valuese = valuesd / sqrt(n))

spikemeans <- spikedata %>% filter(Step<15) %>%
  group_by(Group, Step, StepSize) %>% 
  summarise(n = n_distinct(ID_Subject),
            SpikeCountmean = mean(SpikeCount, na.rm=TRUE),
            SpikeCountsd = sd(SpikeCount,na.rm=TRUE),             
            SpikeCountse = SpikeCountsd / sqrt(n))

reboundmeans <- rebounddata %>% group_by(Group) %>% 
  summarise(n = n_distinct(ID_Subject),
            ReboundCountmean = mean(ReboundCount, na.rm=TRUE),
            ReboundCountsd = sd(ReboundCount,na.rm=TRUE),             
            ReboundCountse = ReboundCountsd / sqrt(n),
            
            ReboundBinarymean = mean(ReboundBinary, na.rm=TRUE),
            ReboundBinarysd = sd(ReboundBinary,na.rm=TRUE),             
            ReboundBinaryse = ReboundBinarysd / sqrt(n))

view(nonspikemeans)

# STATISTICS -----
RMPdata <- nonspikedata %>% subset(Variable=='RMP')
t.test(Value ~ Group, data=RMPdata)
cohen.d(Value ~ Group, data=RMPdata)

t.test(ReboundCount ~ Group, data=rebounddata)
cohen.d(ReboundCount ~ Group, data=rebounddata)

Thresholddata <- nonspikedata %>% subset(Variable=='Threshold')
t.test(Value ~ Group, data=Thresholddata)
cohen.d(Value ~ Group, data=Thresholddata)

InputResistancedata <- nonspikedata %>% subset(Variable=='InputResistance')
t.test(Value ~ Group, data=InputResistancedata)
cohen.d(Value ~ Group, data=InputResistancedata)

Rheodata <- nonspikedata %>% subset(Variable=='Rheo')
t.test(Value ~ Group, data=Rheodata)
cohen.d(Value ~ Group, data=Rheodata)



spikeaov <- aov(SpikeCount ~ Step*Group, data = subset(spikedata, Step<15), )
summary(spikeaov)
effectsize::eta_squared(spikeaov, partial = TRUE)

cohens_f(spikeaov, partial = TRUE)

# FIGURES -----
## Figure 6 Panels -----
### Prepare plot colors and variables -----
groupcolors <- c('Saline'=AdolScolor, 'Morphine'=AdolEPcolor)
grouplabels <- c('Saline'='Control', 'Morphine'='Morphine')
sexshapes <- c("M" = 0, "F" = 1)

# Override set default column width
mycolwidth <- 0.8

# Make empty tile for spacing
emptyplot <- ggplot() + mytheme

# Make empty tile for spacing
holdfortraceplot <- ggplot() + mytheme

### RMP -----
ymax_RMP <- 0
ymin_RMP <- -80

RMP_bar <- ggplot(subset(nonspikemeans, Variable=='RMP'), aes(x=Group, y=valuemean, fill=Group)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_errorbar(aes(ymin=valuemean-valuese, ymax=valuemean+valuese), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  geom_point(subset(nonspikedata, Variable=='RMP'), mapping = aes(x=Group, y=Value, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") + # Add subject points
  scale_shape_manual(values=sexshapes)+
  scale_y_continuous(expand = c(0,0), limits = c(ymin_RMP, ymax_RMP), breaks = seq(-75, ymax_RMP, 25)) + 
  mytheme +
  scale_fill_manual(values=groupcolors) +
  scale_x_discrete(name ="Group",labels=grouplabels) +
  theme(axis.title.x = element_blank())+
  ylab("RMP (mV)") +
  ggtitle("")+
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(-72),label.size = 2, label = "&&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2) +
  geom_hline(yintercept = 0, col = 'black', linewidth = 1)

RMP_bar

RMP <- ggarrange(RMP_bar,holdfortraceplot,holdfortraceplot, ncol=3, nrow=1, widths=c(1,1,1.8), labels = c("A", "B", ""))
RMP

### Spikes During Release from Hold -----
ymax_Reb <- 12
ymin_Reb <- -0.25

ReboundSpikeCount_bar <- ggplot(subset(reboundmeans), aes(x=Group, y=ReboundCountmean, fill=Group)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_errorbar(aes(ymin=ReboundCountmean-ReboundCountse, ymax=ReboundCountmean+ReboundCountse), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  geom_point(subset(rebounddata), mapping = aes(x=Group, y=ReboundCount, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") + # Add subject points
  scale_shape_manual(values=sexshapes)+
  scale_y_continuous(expand = c(0,0), limits = c(ymin_Reb, ymax_Reb), breaks = seq(0, ymax_Reb, 4)) + 
  mytheme +
  scale_fill_manual(values=groupcolors) +
  scale_x_discrete(name ="Group",labels=grouplabels) +
  theme(axis.title.x = element_blank())+
  ylab("Spikes During Release from Hold")+ 
  ggtitle("") +
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(9),label.size = 2, label = "&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2)

ReboundSpikeCount_bar

### Threshold -----
ymax_Thresh <- 0
ymin_Thresh <- -60

Threshold_bar <- ggplot(subset(nonspikemeans, Variable=='Threshold'), aes(x=Group, y=valuemean, fill=Group)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_errorbar(aes(ymin=valuemean-valuese, ymax=valuemean+valuese), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  geom_point(subset(nonspikedata, Variable=='Threshold'), mapping = aes(x=Group, y=Value, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") + # Add subject points
  scale_shape_manual(values=sexshapes)+
  scale_y_continuous(expand = c(0,0), limits = c(ymin_Thresh, ymax_Thresh), breaks = seq(ymin_Thresh, ymax_Thresh, 20)) + 
  mytheme +
  scale_fill_manual(values=groupcolors) +
  scale_x_discrete(name ="Group",labels=grouplabels) +
  theme(axis.title.x = element_blank())+
  ylab("Threshold (mV)")+
  ggtitle("")+
  geom_bracket(xmin=c(1), xmax = c(2), y.position=c(-53),label.size = 2, label = "&&",inherit.aes = FALSE, tip.length = c(.04, .04), size = .4, vjust = -0.2) +
  geom_hline(yintercept = 0, col = 'black', linewidth = 1)

Threshold_bar

ReboundThreshold <- ggarrange(ReboundSpikeCount_bar,Threshold_bar,holdfortraceplot, ncol=3, nrow=1, widths=c(1,1,1.8), labels = c("C", "D", "E"))
ReboundThreshold

### Input Resistance -----
ymax_IR <- 900
ymin_IR <- 0

InputResistance_bar <- ggplot(subset(nonspikemeans, Variable=='InputResistance'), aes(x=Group, y=valuemean, fill=Group)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_errorbar(aes(ymin=valuemean-valuese, ymax=valuemean+valuese), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  geom_point(subset(nonspikedata, Variable=='InputResistance'), mapping = aes(x=Group, y=Value, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") + # Add subject points
  scale_shape_manual(values=sexshapes)+
  scale_y_continuous(expand = c(0,0), limits = c(ymin_IR, ymax_IR), breaks = seq(ymin_IR, ymax_IR, 300)) + 
  mytheme +
  scale_fill_manual(values=groupcolors) +
  scale_x_discrete(name ="Group",labels=grouplabels) +
  theme(axis.title.x = element_blank())+
  ylab("Input Resistance (mOhms)")+ 
  ggtitle("")+
  geom_hline(yintercept = 0, col = 'black', linewidth = 1)

InputResistance_bar

### Rheobase -----
ymax_Rheo <- 180
ymin_Rheo <- 0

Rheo_bar <- ggplot(subset(nonspikemeans, Variable=='Rheo'), aes(x=Group, y=valuemean, fill=Group)) +
  geom_col(width=mycolwidth, linetype="solid", linewidth = mylinewidth) +
  geom_errorbar(aes(ymin=valuemean-valuese, ymax=valuemean+valuese), width=myerrorbarwidth, linewidth = myerrorbarlinewidth, color="black", position=position_dodge()) +
  geom_point(subset(nonspikedata, Variable=='Rheo'), mapping = aes(x=Group, y=Value, shape=Sex), size = mysubjectpointsize, stroke = mypointstroke, color = "black") + # Add subject points
  scale_shape_manual(values=sexshapes)+
  scale_y_continuous(expand = c(0,0), limits = c(ymin_Rheo, ymax_Rheo), breaks = seq(ymin_Rheo, ymax_Rheo, 60)) + 
  mytheme +
  scale_fill_manual(values=groupcolors) +
  scale_x_discrete(name ="Group",labels=grouplabels) +
  theme(axis.title.x = element_blank())+
  ylab("Rheobase (pA)")+
  ggtitle("")

Rheo_bar

IRRheo <- ggarrange(InputResistance_bar,Rheo_bar,holdfortraceplot, ncol=3, nrow=1, widths=c(1,1,1.8), labels = c("F", "G", "H"))
IRRheo


### Spikes Per Step -----
ymax_Spike <- 2.15
ymin_Spike <- 0

spikecount_line <- ggplot(spikemeans, aes(x=StepSize, y=SpikeCountmean, group=Group, fill=Group, color=Group)) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  mytheme +
  scale_y_continuous(limits = c(ymin_Spike, ymax_Spike), breaks = seq(ymin_Spike, ymax_Spike, 0.5)) + 
  scale_x_continuous(name="Step", limits = c(-30, 105), breaks = seq(-20, 105, 20)) +
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values=groupcolors) +
  scale_color_manual(values=groupcolors) +
  ggtitle("")+
  ylab("Spikes per Current Step")

spikecount_line

Spike <- ggarrange(spikecount_line,holdfortraceplot, ncol=2, nrow=1, widths=c(2,1), labels = c("I", "", ""))
Spike


## Create Figure 6 -----
figure6_1 <- ggarrange(RMP_bar, ReboundSpikeCount_bar, emptyplot, emptyplot, ncol=4, nrow=1, labels = c("A", "B", "C", ""), widths=c(1,1,1,1), heights=c(1,1,1,1))
figure6_2 <- ggarrange(InputResistance_bar, Threshold_bar, Rheo_bar, emptyplot, ncol=4, nrow=1, labels = c("D", "E", "F", "G"), widths=c(1,1,1,1), heights=c(1,1,1,1))
figure6_3 <- ggarrange(spikecount_line, emptyplot, ncol=2, nrow=1, labels = c("H", "I"), widths=c(2,2), heights=c(1,1))


# Create full panel
figure6 <- ggarrange(figure6_1, figure6_2, figure6_3, ncol=1, nrow=3, widths=c(1,1,1), heights=c(1,1,1))
figure6
ggsave(plot=figure6, path=path_figures, device="pdf", filename="F6_AllPanels.pdf", width=RFpanelwidth, height=RFpanelheight*3)

