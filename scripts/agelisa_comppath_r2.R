
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: March 24, 2023
# Code purpose: analyze and visualize AgELISA antigen data for comparative pathology goat experiment ROUND 2

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
molecdat <- read.csv("data/comp_path/mol/molecdat_clean_comppath_r2.csv", header = TRUE, stringsAsFactors = FALSE)

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
barncol <- c("1" = "cyan1", "3" = "grey30", "4" = "cyan3", "6" = "cyan4")


# Split data into each barn, and have combined barns
barn1_ag <- subset(molecdat_noc, barn == 1) 
# barn3_ag <- subset(molecdat_noc, barn == 3) # controls, no AgELISA samples taken
barn4_ag <- subset(molecdat_noc, barn == 4) 
barn6_ag <- subset(molecdat_noc, barn == 6)
# barns1346_ag <- subset(molecdat_noc, barn == 1 | barn == 3 | barn == 4 | barn == 6)
barns146_ag <- subset(molecdat_noc, barn == 1 | barn == 4 | barn == 6)


#####################
#  Barn 1 ----
#####################
# Barn 1: 3 Experimental Goats
# Smooth: LOESS 

ag_n_barn1 <- ggplot(barn1_ag, aes(dpi, sp, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Nasal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 1: Nasal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

ag_r_barn1 <- ggplot(barn1_ag, aes(dpi, sp_f, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Rectal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 1: Rectal Status", 
                                 #labels = c("Goats", "Goats"),
                                 values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() + 
              theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# #####################
# #  Barn 3 ----
# #####################
# # Barn 3: 3 Negative Control Goats
# # Smooth: LOESS 
# 
# ag_n_barn3 <- ggplot(barn3_ag, aes(dpi, sp, color = "grey30")) +
#               geom_line(aes(group = eartag), alpha = 0.4) +
#               geom_point(size = 1, alpha = 0.4) +
#               stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
#               scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
#               scale_y_continuous(limits = c(-5, 515)) +
#               labs(title = "Nasal", x = "Days Post Infection", y = "") +
#               scale_color_manual(name = "Barn 3: Nasal Status", 
#                                  #labels = c("Goats", "Goats"),
#                                  values = "grey30") +
#               geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
#               theme_minimal() + 
#               theme(legend.position = "none", axis.title.y = element_text(size = 8), 
#                     axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
# 
# ag_r_barn3 <- ggplot(barn3_ag, aes(dpi, sp_f, color = "grey30")) +
#               geom_line(aes(group = eartag), alpha = 0.4) +
#               geom_point(size = 1, alpha = 0.4) +
#               stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
#               scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
#               scale_y_continuous(limits = c(-5, 515)) +
#               labs(title = "Rectal", x = "Days Post Infection", y = "") +
#               scale_color_manual(name = "Barn 3: Rectal Status", 
#                                  #labels = c("Goats", "Goats"),
#                                  values = "grey30") +
#               geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
#               theme_minimal() + 
#               theme(legend.position = "none", axis.title.y = element_text(size = 8), 
#                     axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


#####################
#  Barn 4 ----
#####################
# Barn 4: 3 Experimental Goats
# Smooth: LOESS

# NO SAMPLES TAKEN

ag_n_barn4 <- ggplot(barn4_ag, aes(dpi, sp, color = "cyan3")) +
              geom_line(aes(group = eartag), alpha = 0.4) +
              geom_point(size = 1, alpha = 0.4) +
              stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
              scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
              scale_y_continuous(limits = c(-5, 515)) +
              labs(title = "Nasal", x = "Days Post Infection", y = "") +
              scale_color_manual(name = "Barn 4: Nasal Status",
                                #labels = c("Goats", "Goats"),
                                values = "cyan3") +
              geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
              theme_minimal() +
              theme(legend.position = "none", axis.title.y = element_text(size = 8),
                   axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

 ag_r_barn4 <- ggplot(barn4_ag, aes(dpi, sp_f, color = "cyan3")) +
               geom_line(aes(group = eartag), alpha = 0.4) +
               geom_point(size = 1, alpha = 0.4) +
               stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
               scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
               scale_y_continuous(limits = c(-5, 515)) +
               labs(title = "Rectal", x = "Days Post Infection", y = "") +
               scale_color_manual(name = "Barn 4: Rectal Status",
                                  #labels = c("Goats", "Goats"),
                                  values = "cyan3") +
               geom_hline(yintercept = 20, linetype ="dashed", color = "red") +
               theme_minimal() +
               theme(legend.position = "none", axis.title.y = element_text(size = 8),
                     axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


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

ag_n_barns146 <- ggplot(barns146_ag, aes(dpi, sp, color = barn, linetype = barn)) +
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

ag_r_barns146 <- ggplot(barns146_ag, aes(dpi, sp_f, color = barn, linetype = barn)) +
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

jpeg("output/comp_path/pprv_comppath_r2_agelisa_barn1_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barn1, ag_r_barn1, ncol=2, left = "S/P Ratio")
invisible(dev.off())

# jpeg("output/comp_path/pprv_comppath_r2_agelisa_barn3_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
# grid.arrange(ag_n_barn3, ag_r_barn3, ncol=2, left = "S/P Ratio")
# invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r2_agelisa_barn4_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barn4, ag_r_barn4, ncol=2, left = "S/P Ratio")
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r2_agelisa_barn6_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barn6, ag_r_barn6, ncol=2, left = "S/P Ratio")
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r2_agelisa_barns146_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(ag_n_barns146, ag_r_barns146, ncol=2, left = "S/P Ratio")
invisible(dev.off())
