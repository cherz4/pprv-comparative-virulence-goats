
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: March 24, 2023
# Code purpose: analyze and visualize serology data for compartive pathology goat experiment ROUND 1

#############################################################################################################


######################
#  Load libraries ----
######################

library(ggplot2)   # for plotting figures
#library(gridExtra) # for grid.arrange function in doing multiple plots for jpeg
#library(tidyverse) # for piping %>%


###############################################
#  Serology Data Import, Data Management  ----
###############################################

# Load the data by pointing to the location within the project that this data resides (do not need to setwd)
# Note: parenthesis are backward compared to Windows Explorer paths
serodat <- read.csv("data/comp_path/sero/serodat_clean_comppath_r1_final.csv", header = TRUE)

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
serodat$eartag <- as.factor(serodat$eartag)
serodat$barn <- as.factor(serodat$barn)

# Setting the line and symbol colors for each barn and status
barncol <- c("2" = "cyan1", "3" = "cyan3", "5" = "grey30", "6" = "cyan4")


# Split data into each barn, and have combined barns
barn2_sero <- subset(serodat, barn == 2) 
barn3_sero <- subset(serodat, barn == 3)
barn5_sero <- subset(serodat, barn == 5) # controls
barn6_sero <- subset(serodat, barn == 6)
barns2356_sero <- subset(serodat, barn == 2 | barn == 3 | barn == 5 | barn == 6)


########################################
#  Saving figures to file (as jpeg) ----
########################################

#####################
#  Barn 2 ----
#####################
# Barn 2: 3 Experimental Goats

jpeg("output/comp_path/pprv_comppath_r1_sero_barn2_todpi12.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn2_sero, aes(dpi, sn, color ="cyan3")) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12, by=1)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn 2", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan3") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
invisible(dev.off())


#####################
#  Barn 3 ----
#####################
# Barn 3: 3 Experimental Goats

jpeg("output/comp_path/pprv_comppath_r1_sero_barn3_todpi12.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn3_sero, aes(dpi, sn, color ="cyan3")) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12, by=1)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn 3", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan3") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
invisible(dev.off())

#####################
#  Barn 5 ----
#####################
# Barn 5: 3 Negative Control Goats

jpeg("output/comp_path/pprv_comppath_r1_sero_barn5controls_todpi12.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn5_sero, aes(dpi, sn, color ="grey30")) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12, by=1)) +
  scale_y_reverse(limits = c(110, 0), breaks = c(0, 25, 50, 75, 100)) +  # modified this line to force specific breaks and get 107 SN in
  scale_color_manual(name = "Barn 5", 
                     labels = c("Goats (Control)"),
                     values = "grey30") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
invisible(dev.off())

#####################
#  Barn 6 ----
#####################
# Barn 6: 3 Experimental Goats

jpeg("output/comp_path/pprv_comppath_r1_sero_barn6_todpi12.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barn6_sero, aes(dpi, sn, color="cyan3")) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12, by=1)) +
  scale_y_reverse(limits = c(105, 0)) +
  scale_color_manual(name = "Barn 6", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan3") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
invisible(dev.off())


######################
#  Combined Barns ----
######################
# All experimental barns together, total 16 experimental Goats

jpeg("output/comp_path/pprv_comppath_r1_sero_barns2356_todpi12.jpeg", width = 6, height = 5, units = "in", quality = 100, res = 600)
ggplot(barns2356_sero, aes(dpi, sn, color=barn)) +
  geom_line(aes(group = eartag)) +
  geom_point() +
  labs(title = "", x = "Days Post Infection", y = "S/N % Sample Competition") +
  scale_x_continuous(limits = c(0,12), breaks = seq(0,12, by=1)) +
  scale_y_reverse(limits = c(105.1, 0)) +    # made this 105.1 instead of 105, to accommodate 1 point
  scale_color_manual(name = "", 
                     labels = c("Barn 2", "Barn 3", "Barn 5", "Barn 6"),
                     values = barncol) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 47, label = "Positive", colour = "red") +
  annotate("text", x = 1, y = 53, label = "Negative", colour = "red") +
  theme_minimal() + 
  theme(legend.position = "none", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
invisible(dev.off())


# Need to put the three plots above each into a named object and then put names in here in line 162 to plot a panel
# jpeg("output/comp_path/pprv_comppath_r1_clin_3panel_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
# grid.arrange(barn2_sero, barn3_sero, barn5_sero, barn6_sero, barns2356_sero, ncol=2, nrow = 3)
# invisible(dev.off())


