
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: March 23, 2023
# Code purpose: analyze and visualize clinical data for comparative pathology goat experiment ROUND 2

#############################################################################################################


######################
#  Load libraries ----
######################

library(ggplot2)   # for plotting figures
# library(mgcv)    # need if doing a GAM (generalized additive model) for the smooth
library(gridExtra) # for grid.arrange function in doing multiple plots for jpeg
library(tidyverse) # for piping %>%


###############################################
#  Clinical Data Import, Data Management  ----
###############################################

# Load the data by pointing to the location within the project that this data resides (do not need to setwd)
# Note: parenthesis are backward compared to Windows Explorer paths
clindat <- read.csv("data/comp_path/clin/clindat_clean_comppath_r2.csv", header = TRUE) # use csv file name on your computer here

# Making sure each variable (column in data frame) is of the correct data type (factor, numeric, character, etc)
clindat$eartag <- as.factor(clindat$eartag)
clindat$barn <- as.factor(clindat$barn)
clindat$tempc <- as.numeric(clindat$tempc)
clindat <- subset(clindat, !dpi <0) # making sure there are no dpi prior to experiment start in the data
#clindat$dpi <- as.integer(clindat$dpi)

# Create clinical score without temperature - this is used throughout the file
clindat <- clindat %>% mutate(clinscore_mod = clindat$general + clindat$feces + clindat$odnd + clindat$headmucos +
                                                  clindat$resp)

# Setting the line and symbol colors for each barn and status
barncol <- c("1" = "cyan1", "3" = "grey30", "4" = "cyan3", "6" = "cyan4")


# Split data into each barn, and have combined barns
barn1 <- subset(clindat, barn == 1) 
barn3 <- subset(clindat, barn == 3) # controls
barn4 <- subset(clindat, barn == 4) 
barn6 <- subset(clindat, barn == 6)
barns1346 <- subset(clindat, barn == 1 | barn == 3 | barn == 4 | barn == 6)


#####################
#  Barn 1 ----
#####################
# Barn 1: 3 Experimental Goats
# Smooth: LOESS 

# PEAKS
# Getting information about where there is a peak each clinical variable is
# rectal temp
barn1i1 <- barn1[barn1$tempc>=36 & barn1$tempc<=41,]
barn1i1_fit <- loess(tempc ~ dpi, barn1i1)
barn1i1_nd <- data.frame(dpi=seq(min(barn1i1$dpi), max(barn1i1$dpi), length=100))
barn1i1_nd$fit <- predict(barn1i1_fit, newdata=barn1i1_nd)
barn1i1tmax <- barn1i1_nd$dpi[which.max(barn1i1_nd$fit)]
#plot(barn1i1$tempc~barn1i1$dpi)
#lines(barn1i1_nd$fit~barn1i1_nd$dpi, col="grey")

# clinical score
barn1i1c <- barn1
barn1i1c_fit <- loess(clinscore_mod ~ dpi, barn1i1c)
barn1i1c_nd <- data.frame(dpi=seq(min(barn1i1c$dpi), max(barn1i1c$dpi), length=100))
barn1i1c_nd$fit <- predict(barn1i1c_fit, newdata=barn1i1c_nd)
barn1i1cmax <- barn1i1c_nd$dpi[which.max(barn1i1c_nd$fit)]
# plot(barn1i1c$clinscore_mod~barn1i1c$dpi)
# lines(barn1i1c_nd$fit~barn1i1c_nd$dpi, col="grey")


# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barn1 <- ggplot(barn1, aes(dpi, tempc, color="cyan1")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan1") +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  scale_y_continuous(limits = c(37,41)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barn1i1tmax), color = "cyan1") +
  annotate(geom = "text", x = (barn1i1tmax + 1), y = 41, label = paste(round(barn1i1tmax,1)), color = "cyan1") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barn1 <- ggplot(barn1, aes(dpi, clinscore_mod, color="cyan1")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan1") +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barn1i1cmax), color = "cyan1") +
  annotate(geom = "text", x = (barn1i1cmax + 1), y = 10, label = paste(round(barn1i1cmax,1)), color = "cyan1") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


#####################
#  Barn 3 ----
#####################
# Barn 3: 3 Negative Control Goats
# Smooth: LOESS 

# PEAKS
# Getting information about where there is a peak each clinical variable is
# rectal temp
barn3i1 <- barn3[barn3$tempc>=36 & barn3$tempc<=41,]
barn3i1_fit <- loess(tempc ~ dpi, barn3i1)
barn3i1_nd <- data.frame(dpi=seq(min(barn3i1$dpi), max(barn3i1$dpi), length=100))
barn3i1_nd$fit <- predict(barn3i1_fit, newdata=barn3i1_nd)
barn3i1tmax <- barn3i1_nd$dpi[which.max(barn3i1_nd$fit)]
#plot(barn3i1$tempc~barn3i1$dpi)
#lines(barn3i1_nd$fit~barn3i1_nd$dpi, col="grey")

# clinical score
barn3i1c <- barn3
barn3i1c_fit <- loess(clinscore_mod ~ dpi, barn3i1c)
barn3i1c_nd <- data.frame(dpi=seq(min(barn3i1c$dpi), max(barn3i1c$dpi), length=100))
barn3i1c_nd$fit <- predict(barn3i1c_fit, newdata=barn3i1c_nd)
barn3i1cmax <- barn3i1c_nd$dpi[which.max(barn3i1c_nd$fit)]
# plot(barn3i1c$clinscore_mod~barn3i1c$dpi)
# lines(barn3i1c_nd$fit~barn3i1c_nd$dpi, col="grey")


# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barn3 <- ggplot(barn3, aes(dpi, tempc, color="grey30")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goats (Control)"),
                     values = "grey30") +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  scale_y_continuous(limits = c(37,41)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barn3i1tmax), color = "grey30") +
  annotate(geom = "text", x = (barn3i1tmax + 1), y = 41, label = paste(round(barn3i1tmax,1)), color = "grey30") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barn3 <- ggplot(barn3, aes(dpi, clinscore_mod, color="grey30")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats (Control)"),
                     values = "grey30") +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barn3i1cmax), color = "grey30") +
  annotate(geom = "text", x = (barn3i1cmax + 1), y = 10, label = paste(round(barn3i1cmax,1)), color = "grey30") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


#####################
#  Barn 4 ----
#####################
# Barn 4: 3 Experimental Goats
# Smooth: LOESS 

# PEAKS
# Getting information about where there is a peak each clinical variable is
# rectal temp
barn4i1 <- barn4[barn4$tempc>=36 & barn4$tempc<=41,]
barn4i1_fit <- loess(tempc ~ dpi, barn4i1)
barn4i1_nd <- data.frame(dpi=seq(min(barn4i1$dpi), max(barn4i1$dpi), length=100))
barn4i1_nd$fit <- predict(barn4i1_fit, newdata=barn4i1_nd)
barn4i1tmax <- barn4i1_nd$dpi[which.max(barn4i1_nd$fit)]
#plot(barn4i1$tempc~barn4i1$dpi)
#lines(barn4i1_nd$fit~barn4i1_nd$dpi, col="grey")

# clinical score
barn4i1c <- barn4
barn4i1c_fit <- loess(clinscore_mod ~ dpi, barn4i1c)
barn4i1c_nd <- data.frame(dpi=seq(min(barn4i1c$dpi), max(barn4i1c$dpi), length=100))
barn4i1c_nd$fit <- predict(barn4i1c_fit, newdata=barn4i1c_nd)
barn4i1cmax <- barn4i1c_nd$dpi[which.max(barn4i1c_nd$fit)]
# plot(barn4i1c$clinscore_mod~barn4i1c$dpi)
# lines(barn4i1c_nd$fit~barn4i1c_nd$dpi, col="grey")


# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barn4 <- ggplot(barn4, aes(dpi, tempc, color="cyan3")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan3") +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  scale_y_continuous(limits = c(37,41)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barn4i1tmax), color = "cyan3") +
  annotate(geom = "text", x = (barn4i1tmax - 1), y = 41, label = paste(round(barn4i1tmax,1)), color = "cyan3") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barn4 <- ggplot(barn4, aes(dpi, clinscore_mod, color="cyan3")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan3") +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barn4i1cmax), color = "cyan3") +
  annotate(geom = "text", x = (barn4i1cmax + 1), y = 10, label = paste(round(barn4i1cmax,1)), color = "cyan3") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))



#####################
#  Barn 6 ----
#####################
# Barn 6: 3 Experimental Goats
# Smooth: LOESS 

# PEAKS
# Getting information about where there is a peak each clinical variable is
# rectal temp
barn6i1 <- barn6[barn6$tempc>=36 & barn6$tempc<=41,]
barn6i1_fit <- loess(tempc ~ dpi, barn6i1)
barn6i1_nd <- data.frame(dpi=seq(min(barn6i1$dpi), max(barn6i1$dpi), length=100))
barn6i1_nd$fit <- predict(barn6i1_fit, newdata=barn6i1_nd)
barn6i1tmax <- barn6i1_nd$dpi[which.max(barn6i1_nd$fit)]
# plot(barn6i1$tempc~barn6i1$dpi)
# lines(barn6i1_nd$fit~barn6i1_nd$dpi, col="grey")

#clinical score
barn6i1c <- barn6
barn6i1c_fit <- loess(clinscore_mod ~ dpi, barn6i1c)
barn6i1c_nd <- data.frame(dpi=seq(min(barn6i1c$dpi), max(barn6i1c$dpi), length=100))
barn6i1c_nd$fit <- predict(barn6i1c_fit, newdata=barn6i1c_nd)
barn6i1cmax <- barn6i1c_nd$dpi[which.max(barn6i1c_nd$fit)]
# plot(barn6i1c$clinscore_mod~barn6i1c$dpi)
# lines(barn6i1c_nd$fit~barn6i1c_nd$dpi, col="grey")

# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barn6 <- ggplot(barn6, aes(dpi, tempc, color="cyan4")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan4") +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  scale_y_continuous(limits = c(37,41)) +
  # Next two lines are for plotting peaks, if desired  NOTE: had to shift backward to see peak label in second line
  geom_vline(aes(xintercept = barn6i1tmax), color = "cyan4") +
  annotate(geom = "text", x = (barn6i1tmax - 1), y = 41, label = paste(round(barn6i1tmax,1)), color = "cyan4") + 
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barn6 <- ggplot(barn6, aes(dpi, clinscore_mod, color="cyan4")) +
  geom_line(aes(group = eartag), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_color_manual(name = "", 
                     labels = c("Goats (Inoculated)"),
                     values = "cyan4") +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barn6i1cmax), color = "cyan4") +
  annotate(geom = "text", x = (barn6i1cmax + 1), y = 10, label = paste(round(barn6i1cmax,1)), color = "cyan4") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))




######################
#  Combined Barns ----
######################
# All experimental barns together, total 16 experimental Goats
# Smooth: LOESS 

# PEAKS
# Getting information about where there is a peak each clinical variable is
# rectal temp
barns1346i1 <- barns1346[barns1346$tempc>=36 & barns1346$tempc<=41,]
barns1346i1_fit <- loess(tempc ~ dpi, barns1346i1)
barns1346i1_nd <- data.frame(dpi=seq(min(barns1346i1$dpi), max(barns1346i1$dpi), length=100))
barns1346i1_nd$fit <- predict(barns1346i1_fit, newdata=barns1346i1_nd)
barns1346i1tmax <- barns1346i1_nd$dpi[which.max(barns1346i1_nd$fit)]
# plot(barns1346i1$tempc~barns1346i1$dpi)
# lines(barns1346i1_nd$fit~barns1346i1_nd$dpi, col="grey")

#clinical score
barns1346i1c <- barns1346
barns1346i1c_fit <- loess(clinscore_mod ~ dpi, barns1346i1c)
barns1346i1c_nd <- data.frame(dpi=seq(min(barns1346i1c$dpi), max(barns1346i1c$dpi), length=100))
barns1346i1c_nd$fit <- predict(barns1346i1c_fit, newdata=barns1346i1c_nd)
barns1346i1cmax <- barns1346i1c_nd$dpi[which.max(barns1346i1c_nd$fit)]
# plot(barns1346i1c$clinscore_mod~barns1346i1c$dpi)
# lines(barns1346i1c_nd$fit~barns1346i1c_nd$dpi, col="grey")


# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barns1346 <- ggplot(barns1346, aes(dpi, tempc, color=barn, linetype = barn)) +
  geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_linetype_manual(name = "", values = c("1" = "solid", "3" = "dotted", "4" = "dashed", "6" = "twodash")) +
  scale_color_manual(name = "", 
                     #labels = c("Goats", "Goats"),
                     values = barncol) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  scale_y_continuous(limits = c(37,41)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barns1346i1tmax), color = "grey50") +
  annotate(geom = "text", x = (barns1346i1tmax + 1), y = 41, label = paste(round(barns1346i1tmax,1)), color = "grey50") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barns1346 <- ggplot(barns1346, aes(dpi, clinscore_mod, color=barn, linetype = barn)) +
  geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_linetype_manual(name = "", values = c("1" = "solid", "3" = "dotted", "4" = "dashed", "6" = "twodash")) +
  scale_color_manual(name = "", 
                     #labels = c("Goats", "Goats"),
                     values = barncol) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barns1346i1cmax), color = "grey50") +
  annotate(geom = "text", x = (barns1346i1cmax + 1), y = 10, label = paste(round(barns1346i1cmax,1)), color = "grey50") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


########################################
#  Saving figures to file (as jpeg) ----
########################################

jpeg("output/comp_path/pprv_comppath_r2_clin_barn1_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barn1, cs_barn1, ncol=2)
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r2_clin_barn3_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barn3, cs_barn3, ncol=2)
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r2_clin_barn4controls_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barn4, cs_barn4, ncol=2)
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r2_clin_barn6_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barn6, cs_barn6, ncol=2)
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r2_clin_barns1346_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barns1346, cs_barns1346, ncol=2)
invisible(dev.off())

# If need symbols by unique eartag, sample code for this is in PPRV Trial Stage 1.1A code
