### Set plot variables
# Font sizes
mybasesize <- 10
mytitlesize <- 10
myaxislabelsize <- 8
mystarsize <- 5

# General plot area
myaxislinewidth <- 0.5
mylinewidth <- 1

# Plot variables
mylinewidth <- .8
mycolwidth <- .6
mycollinewidth <- 0.8
myerrorbarwidth <- .3
myerrorbarlinewidth <- .28
mypointsize <- 1.5
mysubjectpointsize <- 1
mypointstroke <- .5

# Inj label size
myinjtrianglesize <- 4

### Make plot theme
mytheme <-  theme_prism(base_size = mybasesize, base_fontface = "plain") + 
            theme(legend.position = 'none', axis.line.x=element_line(linewidth=myaxislinewidth), axis.line.y=element_line(linewidth=myaxislinewidth),
                  axis.ticks=element_line(linewidth=myaxislinewidth), axis.ticks.length=unit(1, "mm"),
                  axis.text.x=element_text(size=myaxislabelsize), axis.text.y=element_text(size=myaxislabelsize), 
                  axis.title=element_text(size=mybasesize, face='plain'), title=element_text(size=mytitlesize, face='plain'), 
                  axis.title.y=element_text(margin=margin(r=10)))



### Load plot colors
colorkey <- read_csv("C:/Users/rmdon/Box/JRoit Lab/Publications/In Process/PB Adolescent Morphine Reward/Figures/Adolescent Morphine ICSS Color Key.csv") # read in raw data

ICSSBLcolor <- colorkey$HexCode[which(colorkey$VariableName == 'ICSSBLcolor')]
AdolScolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdolScolor')]
AdolMcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdolMcolor')]
AdolEPcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdolEPcolor')]
AdolLPcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdolLPcolor')]
Injcolor <- colorkey$HexCode[which(colorkey$VariableName == 'Injcolor')]

AdolSAdultcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdolSAdultcolor')]
AdolMAdultcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdolMAdultcolor')]

AdolSAdultEPcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdolSAdultEPcolor')]
AdolSAdultLPcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdolSAdultLPcolor')]

AdolMAdultEPcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdolMAdultEPcolor')]
AdolMAdultLPcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdolMAdultLPcolor')]


AdultNaiveMcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdultNaiveMcolor')]
AdultNaiveEPcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdultNaiveEPcolor')]
AdultNaiveLPcolor <- colorkey$HexCode[which(colorkey$VariableName == 'AdultNaiveLPcolor')]

### Set plot output widths
RFpanelwidth <- 7.75
RFpanelheight <- 2.15

RFplotwidth_line <- 4.5
RFplotheight_line <- 1.75

RFplotwidth_bar <- 1.5
RFplotheight_bar <- 1.75
