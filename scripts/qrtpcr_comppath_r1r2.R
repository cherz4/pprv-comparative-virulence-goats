
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: June 23, 2024
# Code purpose: analyze and visualize qRT-PCR antigen data for comparative pathology goat experiment both rounds

#############################################################################################################


######################
#  Load libraries ----
######################

library(ggplot2)   # for plotting figures
library(gridExtra) # for grid.arrange function in doing multiple plots for jpeg
library(tidyverse) # for piping %>%


###############################################
#  Molecular Data Import, Data Management  ----
###############################################

# Round 1
# Isolates: 24475/21 in Barn 2, 2014 Isolate in Barn 3, 38920/19 in Barn 6, controls in Barn 5
# Isolates: 26173/19 in Barn 1, 35368/19 in Barn 4, 25870/21 in Barn 6, Controls in Barn 3


# Load the data by pointing to the location within the project that this data resides (do not need to setwd)
# Note: parenthesis are backward compared to Windows Explorer paths
# Round 1
molecdat_r1 <- read.csv("data/comp_path/mol/molecdat_clean_comppath_r1_v2.csv", 
                        header = TRUE, stringsAsFactors = FALSE)

# Round 2
molecdat_r2 <- read.csv("data/comp_path/mol/molecdat_clean_comppath_r2_v2.csv", 
                        header = TRUE, stringsAsFactors = FALSE) %>%
               mutate(barn = case_when(barn == 6 ~ 7, 
                                       barn == 3 ~ 8,
                                       TRUE ~ barn))
# Change barn 6 to read barn 7 to have a unique set of barn numbers for combining files

#Bind
molecdat <- rbind(molecdat_r1, molecdat_r2)

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
molecdat$eartag <- as.factor(molecdat$eartag)
molecdat$expt <- as.factor(molecdat$expt)
molecdat$species <- as.factor(molecdat$species)
molecdat$barn <- as.factor(molecdat$barn)
molecdat$status <- as.factor(molecdat$status)
molecdat$status_f <- as.factor(molecdat$status_f)

# Setting the line and symbol colors for each barn and status

barncol <- c("1" = "yellow", "2" = "cyan1", "3" = "cyan4", "4" = "aquamarine4", "5" = "black",
             "6" = "red", "7" = "seagreen4", "8" = "black")


#####################
#  All Barns ----
#####################
# Smooth: LOESS 

# ocular
pcr_o <- ggplot(molecdat, aes(dpi, 35-omean, color = barn)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Ocular", x = "", y = "") +
  scale_color_manual(name = "Isolates", 
                     #labels = c("Goats", "Goats"),
                     values = barncol) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
pcr_n <- ggplot(molecdat, aes(dpi, 35-nmean, color = barn)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Nasal", x = "", y = "") +
  scale_color_manual(name = "Isolates", 
                     #labels = c("Goats", "Goats"),
                     values = barncol) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
pcr_r<- ggplot(molecdat, aes(dpi, 35-rmean, color = barn)) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  geom_point(size = 1, alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  scale_y_continuous(limits = c(0,20)) +
  labs(title = "Rectal", x = "", y = "") +
  scale_color_manual(name = "Isolates", 
                     #labels = c("Goats", "Goats"),
                     values = barncol) +
  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

########################################
#  Saving figure to file (as jpeg) ----
########################################

jpeg("output/comp_path/pprv_comppath_r1r2_qRTPCR_all_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(pcr_o, pcr_n, pcr_r, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())