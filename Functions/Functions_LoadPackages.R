#### Load packages I always use
library(rstudioapi)
library(miceadds)
library(conflicted)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(rstatix)
library(stringr)
library(ggprism)
library(drc)

# Mixed effects modeling
library(glmmTMB)
library(emmeans)
library(parameters)
library(performance)
library(effectsize)

# Set up preferences
conflicted::conflicts_prefer(dplyr::summarize, dplyr::mutate, dplyr::filter, dplyr::select, lmerTest::lmer, )
