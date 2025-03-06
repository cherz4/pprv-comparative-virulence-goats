
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: June 22, 2024
# Code purpose: analyze and visualize qRT-PCR antigen data for comparative pathology goat experiment ROUND 1

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
molecdat <- read.csv("data/comp_path/mol/molecdat_clean_comppath_r1_v2.csv", header = TRUE, stringsAsFactors = FALSE)

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
molecdat$eartag <- as.factor(molecdat$eartag)
molecdat$expt <- as.factor(molecdat$expt)
molecdat$species <- as.factor(molecdat$species)
molecdat$barn <- as.factor(molecdat$barn)
molecdat$status <- as.factor(molecdat$status)
molecdat$status_f <- as.factor(molecdat$status_f)

# Setting the line and symbol colors for each barn and status
barncol <- c("2" = "cyan1", "3" = "cyan3", "5" = "grey30", "6" = "cyan4")


# Split data into each barn, and have combined barns
barn2_pcr <- subset(molecdat, barn == 2) 
barn3_pcr <- subset(molecdat, barn == 3)
#barn5_pcr <- subset(molecdat, barn == 5) # controls - no qRT-PCR samples taken
barn6_pcr <- subset(molecdat, barn == 6)
#barns2356_pcr <- subset(molecdat, barn == 2 | barn == 3 | barn == 5 | barn == 6) # controls barn 5, no samples
barns236_pcr <- subset(molecdat, barn == 2 | barn == 3 | barn == 6)


#####################
#  Barn 2 ----
#####################
# Barn 2: 3 Experimental Goats
# Smooth: LOESS 

# ocular
pcr_o_barn2 <- ggplot(barn2_pcr, aes(dpi, 35-omean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Ocular", x = "", y = "") +
                scale_color_manual(name = "Barn 2", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
pcr_n_barn2 <- ggplot(barn2_pcr, aes(dpi, 35-nmean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Nasal", x = "", y = "") +
                scale_color_manual(name = "Barn 2", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
pcr_r_barn2 <- ggplot(barn2_pcr, aes(dpi, 35-rmean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Rectal", x = "", y = "") +
                scale_color_manual(name = "Barn 2", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))



#####################
#  Barn 3 ----
#####################
# Barn 3: 3 Experimental Goats
# Smooth: LOESS 

# ocular
pcr_o_barn3 <- ggplot(barn3_pcr, aes(dpi, 35-omean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Ocular", x = "", y = "") +
                scale_color_manual(name = "Barn 3", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# nasal
pcr_n_barn3 <- ggplot(barn3_pcr, aes(dpi, 35-nmean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Nasal", x = "", y = "") +
                scale_color_manual(name = "Barn 3", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
                theme_minimal() + 
                theme(legend.position = "none", axis.title.y = element_text(size = 8), 
                      axis.title.x = element_text(size = 8), axis.text=element_text(size=6))

# rectal
pcr_r_barn3 <- ggplot(barn3_pcr, aes(dpi, 35-rmean, color = "cyan3")) +
                geom_line(aes(group = eartag), alpha = 0.4) +
                geom_point(size = 1, alpha = 0.4) +
                stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
                scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
                scale_y_continuous(limits = c(0,20)) +
                labs(title = "Rectal", x = "", y = "") +
                scale_color_manual(name = "Barn 3", 
                                   #labels = c("Goats", "Goats"),
                                   values = "cyan3") +
                geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
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
# ocular
# pcr_o_barn5 <- ggplot(barn5_pcr, aes(dpi, 35-omean, color = "cyan3")) +
#                 geom_line(aes(group = eartag), alpha = 0.4) +
#                 geom_point(size = 1, alpha = 0.4) +
#                 stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
#                 scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
#                 scale_y_continuous(limits = c(0,20)) +
#                 labs(title = "Ocular", x = "", y = "") +
#                 scale_color_manual(name = "Barn 5", 
#                                    #labels = c("Goats", "Goats"),
#                                    values = "cyan3") +
#                 geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
#                 theme_minimal() + 
#                 theme(legend.position = "none", axis.title.y = element_text(size = 8), 
#                       axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
# 
# # nasal
# pcr_n_barn5 <- ggplot(barn5_pcr, aes(dpi, 35-nmean, color = "cyan3")) +
#                 geom_line(aes(group = eartag), alpha = 0.4) +
#                 geom_point(size = 1, alpha = 0.4) +
#                 stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
#                 scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
#                 scale_y_continuous(limits = c(0,20)) +
#                 labs(title = "Nasal", x = "", y = "") +
#                 scale_color_manual(name = "Barn 5", 
#                                    #labels = c("Goats", "Goats"),
#                                    values = "cyan3") +
#                 geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
#                 theme_minimal() + 
#                 theme(legend.position = "none", axis.title.y = element_text(size = 8), 
#                       axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
# 
# # rectal
# pcr_r_barn5 <- ggplot(barn5_pcr, aes(dpi, 35-rmean, color = "cyan3")) +
#                 geom_line(aes(group = eartag), alpha = 0.4) +
#                 geom_point(size = 1, alpha = 0.4) +
#                 stat_smooth(method = "loess", formula = y ~ x, size =1, alpha = 0.2) +
#                 scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
#                 scale_y_continuous(limits = c(0,20)) +
#                 labs(title = "Rectal", x = "", y = "") +
#                 scale_color_manual(name = "Barn 5", 
#                                    #labels = c("Goats", "Goats"),
#                                    values = "cyan3") +
#                 geom_hline(yintercept = 35, linetype ="dashed", color = "red") +
#                 theme_minimal() + 
#                 theme(legend.position = "none", axis.title.y = element_text(size = 8), 
#                       axis.title.x = element_text(size = 8), axis.text=element_text(size=6))
# 

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
pcr_o_barns236 <- ggplot(barns236_pcr, aes(dpi, 35-omean, color = barn)) +
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
pcr_n_barns236 <- ggplot(barns236_pcr, aes(dpi, 35-nmean, color = barn)) +
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
pcr_r_barns236 <- ggplot(barns236_pcr, aes(dpi, 35-rmean, color = barn)) +
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

jpeg("output/comp_path/pprv_comppath_r1_qRTPCR_barn2_todpi12_v2.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(pcr_o_barn2, pcr_n_barn2, pcr_r_barn2, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r1_qRTPCR_barn3_todpi12_v2.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(pcr_o_barn3, pcr_n_barn3, pcr_r_barn3, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

# jpeg("output/comp_path/pprv_comppath_r1_qRTPCR_barn5_todpi12_v2.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
# grid.arrange(pcr_o_barn5, pcr_n_barn5, pcr_r_barn5, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
# invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r1_qRTPCR_barn6_todpi12_v2.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(pcr_o_barn6, pcr_n_barn6, pcr_r_barn6, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r1_qRTPCR_barns236_todpi12_v2.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(pcr_o_barns236, pcr_n_barns236, pcr_r_barns236, ncol=3, left = "qRT-PCR 35 - mean Ct", bottom = "Days Post Infection")
invisible(dev.off())
