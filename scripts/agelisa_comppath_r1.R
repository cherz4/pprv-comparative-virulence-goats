
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: March 24, 2023
# Code purpose: analyze and visualize AgELISA antigen data for comparative pathology goat experiment ROUND 1

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

# Load the data by pointing to the location within the project that this data resides (do not need to setwd)
# Note: parenthesis are backward compared to Windows Explorer paths
molecdat <- read.csv("data/comp_path/mol/molecdat_clean_comppath_r1.csv", header = TRUE, stringsAsFactors = FALSE)

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
molecdat_noc <- molecdat %>% select(eartag, expt, species, dpi, barn, od, sp, status, 
                                    od_f, sp_f, status_f, notes1, notes2)

molecdat_noc$eartag <- as.factor(molecdat_noc$eartag)
molecdat_noc$expt <- as.factor(molecdat_noc$expt)
molecdat_noc$species <- as.factor(molecdat_noc$species)
molecdat_noc$barn <- as.factor(molecdat_noc$barn)
molecdat_noc$status <- as.factor(molecdat_noc$status)
molecdat_noc$status_f <- as.factor(molecdat_noc$status_f)

# Setting the line and symbol colors for each barn and status
barncol <- c("2" = "cyan1", "3" = "cyan3", "5" = "grey30", "6" = "cyan4")


# Split data into each barn, and have combined barns
barn2_ag <- subset(molecdat_noc, barn == 2) 
barn3_ag <- subset(molecdat_noc, barn == 3)
#barn5_ag <- subset(molecdat_noc, barn == 5) # controls - no AgELISA samples taken
barn6_ag <- subset(molecdat_noc, barn == 6)
#barns2356_ag <- subset(molecdat_noc, barn == 2 | barn == 3 | barn == 5 | barn == 6) # controls barn 5, no samples
barns236_ag <- subset(molecdat_noc, barn == 2 | barn == 3 | barn == 6)


#####################
#  Barn 2 ----
#####################
# Barn 2: 3 Experimental Goats
# Smooth: LOESS 

ag_n_barn2 <- ggplot(barn2_ag, aes(dpi, sp, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Nasal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 2: Nasal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barn2 <- ggplot(barn2_ag, aes(dpi, sp_f, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Rectal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 2: Rectal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


#####################
#  Barn 3 ----
#####################
# Barn 3: 3 Experimental Goats
# Smooth: LOESS 

ag_n_barn3 <- ggplot(barn3_ag, aes(dpi, sp, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Nasal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 3: Nasal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barn3 <- ggplot(barn3_ag, aes(dpi, sp_f, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Rectal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 3: Rectal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


#####################
#  Barn 5 ----
#####################
# Barn 5: 3 Negative Control Goats
# Smooth: LOESS 

# NO SAMPLES TAKEN
# 
# ag_n_barn5 <- ggplot(barn5_ag, aes(dpi, sp, color = "grey30")) +
#               geom_line(aes(group = eartag), alpha = 0.4) +
#               geom_point(size = 1, alpha = 0.4) +
#               stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
#               scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
#               scale_y_continuous(limits = c(-5, 515)) +
#               labs(title = "Nasal", x = "Days Post Infection", y = "") +
#               scale_color_manual(name = "Barn 5: Nasal Status", 
#                                  #labels = c("Goats", "Goats"),
#                                  values = "grey30") +
#               geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
#               theme_minimal() + 
#               theme(legend.position = "none", axis.title.y = element_text(size = 8), 
#                     axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
# 
# ag_r_barn5 <- ggplot(barn5_ag, aes(dpi, sp_f, color = "grey30")) +
#               geom_line(aes(group = eartag), alpha = 0.4) +
#               geom_point(size = 1, alpha = 0.4) +
#               stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
#               scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
#               scale_y_continuous(limits = c(-5, 515)) +
#               labs(title = "Rectal", x = "Days Post Infection", y = "") +
#               scale_color_manual(name = "Barn 5: Rectal Status", 
#                                  #labels = c("Goats", "Goats"),
#                                  values = "grey30") +
#               geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
#               theme_minimal() + 
#               theme(legend.position = "none", axis.title.y = element_text(size = 8), 
#                     axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
# 

#####################
#  Barn 6 ----
#####################
# Barn 6: 3 Experimental Goats
# Smooth: LOESS 

ag_n_barn6 <- ggplot(barn6_ag, aes(dpi, sp, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Nasal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 6: Nasal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") + 
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


ag_r_barn6 <- ggplot(barn6_ag, aes(dpi, sp_f, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Rectal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 6: Rectal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() +
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))




######################
#  Combined Barns ----
######################
# All experimental barns together, total 16 experimental Goats
# Smooth: LOESS 

ag_n_barns236 <- ggplot(barns236_ag, aes(dpi, sp, color = barn, linetype = barn)) +
                  geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
                  geom_point(size = 1, alpha = 0.4) +
                  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                  scale_y_continuous(limits = c(-5, 515)) +
                  labs(title = "Nasal", x = "Days Post Infection", y = "") +
                  scale_color_manual(#name = "Nasal Status", 
                                     #labels = c("Goats", "Goats"),
                                     values = barncol) +  
                  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
                  theme_minimal() +
                  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barns236 <- ggplot(barns236_ag, aes(dpi, sp_f, color = barn, linetype = barn)) +
                  geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
                  geom_point(size = 1, alpha = 0.4) +
                  stat_smooth(aes(group=barn), method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                  scale_y_continuous(limits = c(-5, 515)) +
                  labs(title = "Rectal", x = "Days Post Infection", y = "") +
                  scale_color_manual(#name = "Rectal Status", 
                                     #labels = c("Goats", "Goats"),
                                     values = barncol) +  
                  geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
                  theme_minimal() +
                  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


########################################
#  Saving figures to file (as jpeg) ----
########################################

jpeg("output/comp_path/pprv_comppath_r1_agelisa_barn2_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barn2, ag_r_barn2, ncol=2, left = "S/P Ratio")
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r1_agelisa_barn3_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barn3, ag_r_barn3, ncol=2, left = "S/P Ratio")
invisible(dev.off())

# jpeg("output/comp_path/pprv_comppath_r1_agelisa_barn5_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
# grid.arrange(ag_n_barn5, ag_r_barn5, ncol=2, left = "S/P Ratio")
# invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r1_agelisa_barn6_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barn6, ag_r_barn6, ncol=2, left = "S/P Ratio")
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r1_agelisa_barns236_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barns236, ag_r_barns236, ncol=2, left = "S/P Ratio")
invisible(dev.off())
