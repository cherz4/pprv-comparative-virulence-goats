
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: June 22, 2024
# Code purpose: analyze and visualize qRT-PCR antigen data for comparative pathology goat experiment ROUND 2

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
molecdat <- read.csv("data/comp_path/mol/molecdat_clean_comppath_r2_v2.csv", header = TRUE, stringsAsFactors = FALSE)

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
molecdat$eartag <- as.factor(molecdat$eartag)
molecdat$expt <- as.factor(molecdat$expt)
molecdat$species <- as.factor(molecdat$species)
molecdat$barn <- as.factor(molecdat$barn)
molecdat$status <- as.factor(molecdat$status)
molecdat$status_f <- as.factor(molecdat$status_f)

# Setting the line and symbol colors for each barn and status
barncol <- c("1" = "cyan1", "3" = "grey30", "4" = "cyan3", "6" = "cyan4")


# Split data into each barn, and have combined barns
barn1_pcr <- subset(molecdat, barn == 1) 
# barn3_pcr <- subset(molecdat, barn == 3) # controls, no samples taken
barn4_pcr <- subset(molecdat, barn == 4) 
barn6_pcr <- subset(molecdat, barn == 6)
# barns1346_pcr <- subset(molecdat, barn == 1 | barn == 3 | barn == 4 | barn == 6)
barns146_pcr <- subset(molecdat, barn == 1 | barn == 4 | barn == 6)


#####################
#  Barn 1 ----
#####################
# Barn 1: 3 Experimental Goats
# Smooth: LOESS 

# ocular
pcr_o_barn1 <- ggplot(barn1_pcr, aes(dpi, 35-omean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Ocular", x = "", y = "") +
                scale_color_manual(name = "Barn 1", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
pcr_n_barn1 <- ggplot(barn1_pcr, aes(dpi, 35-nmean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Nasal", x = "", y = "") +
                scale_color_manual(name = "Barn 1", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
pcr_r_barn1 <- ggplot(barn1_pcr, aes(dpi, 35-rmean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Rectal", x = "", y = "") +
                scale_color_manual(name = "Barn 1", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))



# #####################
# #  Barn 3 ----
# #####################
# # Barn 3: 3 Negative Control Goats
# # Smooth: LOESS 

# # ocular
# pcr_o_barn3 <- ggplot(barn3_pcr, aes(dpi, 35-omean, color = "grey30")) +
#                 geom_line(aes(group = eartag), alpha = 0.4) +
#                 geom_point(size = 1, alpha = 0.4) +
#                 stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
#                 scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
#                 scale_y_continuous(limits = c(0,20)) +
#                 labs(title = "Ocular", x = "", y = "") +
#                 scale_color_manual(name = "Barn 3", 
#                                    #labels = c("Goats", "Goats"),
#                                    values = "grey30") +
#                 geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
#                 theme_minimal() + 
#                 theme(legend.position = "none", axis.title.y = element_text(size = 8), 
#                       axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
# 
# # nasal
# pcr_n_barn3 <- ggplot(barn3_pcr, aes(dpi, 35-nmean, color = "grey30")) +
#                 geom_line(aes(group = eartag), alpha = 0.4) +
#                 geom_point(size = 1, alpha = 0.4) +
#                 stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
#                 scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
#                 scale_y_continuous(limits = c(0,20)) +
#                 labs(title = "Nasal", x = "", y = "") +
#                 scale_color_manual(name = "Barn 3", 
#                                    #labels = c("Goats", "Goats"),
#                                    values = "grey30") +
#                 geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
#                 theme_minimal() + 
#                 theme(legend.position = "none", axis.title.y = element_text(size = 8), 
#                       axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
# 
# # rectal
# pcr_r_barn3 <- ggplot(barn3_pcr, aes(dpi, 35-rmean, color = "grey30")) +
#                 geom_line(aes(group = eartag), alpha = 0.4) +
#                 geom_point(size = 1, alpha = 0.4) +
#                 stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
#                 scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
#                 scale_y_continuous(limits = c(0,20)) +
#                 labs(title = "Rectal", x = "", y = "") +
#                 scale_color_manual(name = "Barn 3", 
#                                    #labels = c("Goats", "Goats"),
#                                    values = "grey30") +
#                 geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
#                 theme_minimal() + 
#                 theme(legend.position = "none", axis.title.y = element_text(size = 8), 
#                       axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

#####################
#  Barn 4 ----
#####################
# Barn 4: 3 Experimental Goats
# Smooth: LOESS

# ocular
pcr_o_barn4 <- ggplot(barn4_pcr, aes(dpi, 35-omean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Ocular", x = "", y = "") +
                scale_color_manual(name = "Barn 4",
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() +
                theme(legend.position = "none", axis.title.y = element_text(size = 8),
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
pcr_n_barn4 <- ggplot(barn4_pcr, aes(dpi, 35-nmean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Nasal", x = "", y = "") +
                scale_color_manual(name = "Barn 4",
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() +
                theme(legend.position = "none", axis.title.y = element_text(size = 8),
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
pcr_r_barn4 <- ggplot(barn4_pcr, aes(dpi, 35-rmean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Rectal", x = "", y = "") +
                scale_color_manual(name = "Barn 4",
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() +
                theme(legend.position = "none", axis.title.y = element_text(size = 8),
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


#####################
#  Barn 6 ----
#####################
# Barn 6: 3 Experimental Goats
# Smooth: LOESS 

# ocular
pcr_o_barn6 <- ggplot(barn6_pcr, aes(dpi, 35-omean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Ocular", x = "", y = "") +
                scale_color_manual(name = "Barn 6", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
pcr_n_barn6 <- ggplot(barn6_pcr, aes(dpi, 35-nmean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Nasal", x = "", y = "") +
                scale_color_manual(name = "Barn 6", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
pcr_r_barn6 <- ggplot(barn6_pcr, aes(dpi, 35-rmean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Rectal", x = "", y = "") +
                scale_color_manual(name = "Barn 6", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))



######################
#  Combined Barns ----
######################
# All experimental barns together, total 16 experimental Goats
# Smooth: LOESS 

#ocular
pcr_o_barns146 <- ggplot(barns146_pcr, aes(dpi, 35-omean, color = barn)) +
                  geom_line(aes(group = eartag), alpha = 0.4) +
                  geom_point(size = 1, alpha = 0.4) +
                  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                  scale_y_continuous(limits = c(0,20)) +
                  labs(title = "Ocular", x = "", y = "") +
                  scale_color_manual(name = "Barn", 
                                     #labels = c("Goats", "Goats"),
                                     values = barncol) +
                  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                  theme_minimal() + 
                  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

## nasal
pcr_n_barns146 <- ggplot(barns146_pcr, aes(dpi, 35-nmean, color = barn)) +
                  geom_line(aes(group = eartag), alpha = 0.4) +
                  geom_point(size = 1, alpha = 0.4) +
                  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                  scale_y_continuous(limits = c(0,20)) +
                  labs(title = "Nasal", x = "", y = "") +
                  scale_color_manual(name = "Barn", 
                                     #labels = c("Goats", "Goats"),
                                     values = barncol) +
                  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                  theme_minimal() + 
                  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
pcr_r_barns146 <- ggplot(barns146_pcr, aes(dpi, 35-rmean, color = barn)) +
                  geom_line(aes(group = eartag), alpha = 0.4) +
                  geom_point(size = 1, alpha = 0.4) +
                  stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                  scale_y_continuous(limits = c(0,20)) +
                  labs(title = "Rectal", x = "", y = "") +
                  scale_color_manual(name = "Barn", 
                                     #labels = c("Goats", "Goats"),
                                     values = barncol) +
                  geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                  theme_minimal() + 
                  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
                        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


########################################
#  Saving figures to file (as jpeg) ----
########################################

jpeg("output/comp_path/pprv_comppath_r2_qRTPCR_barn1_todpi12_v2.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(pcr_o_barn1, pcr_n_barn1, pcr_r_barn1, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

# jpeg("output/comp_path/pprv_comppath_r2_qRTPCR_barn3_todpi12_v2.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
# grid.arrange(pcr_o_barn3, pcr_n_barn3, pcr_r_barn3, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
# invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r2_qRTPCR_barn4_todpi12_v2.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(pcr_o_barn4, pcr_n_barn4, pcr_r_barn4, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r2_qRTPCR_barn6_todpi12_v2.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(pcr_o_barn6, pcr_n_barn6, pcr_r_barn6, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r2_qRTPCR_barns146_todpi12_v2.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(pcr_o_barns146, pcr_n_barns146, pcr_r_barns146, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())
