
#############################################################################################################
# Code written by: Catherine M. Herzog, PhD MPH
# Code last modified on: March 23, 2023
# Code purpose: analyze and visualize clinical data for comparative pathology goat experiment ROUND 1

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
clindat <- read.csv("data/comp_path/clin/clindat_clean_comppath_r1.csv", header = TRUE) # use csv file name on your computer here

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
barncol <- c("2" = "cyan1", "3" = "cyan3", "5" = "grey30", "6" = "cyan4")


# Split data into each barn, and have combined barns
barn2 <- subset(clindat, barn == 2) 
barn3 <- subset(clindat, barn == 3)
barn5 <- subset(clindat, barn == 5) # controls
barn6 <- subset(clindat, barn == 6)
barns2356 <- subset(clindat, barn == 2 | barn == 3 | barn == 5 | barn == 6)


#####################
#  Barn 2 ----
#####################
# Barn 2: 3 Experimental Goats
# Smooth: LOESS 

# PEAKS
# Getting information about where there is a peak each clinical variable is
# rectal temp
barn2i1 <- barn2[barn2$tempc>=36 & barn2$tempc<=41,]
barn2i1_fit <- loess(tempc ~ dpi, barn2i1)
barn2i1_nd <- data.frame(dpi=seq(min(barn2i1$dpi), max(barn2i1$dpi), length=100))
barn2i1_nd$fit <- predict(barn2i1_fit, newdata=barn2i1_nd)
barn2i1tmax <- barn2i1_nd$dpi[which.max(barn2i1_nd$fit)]
#plot(barn2i1$tempc~barn2i1$dpi)
#lines(barn2i1_nd$fit~barn2i1_nd$dpi, col="grey")

# clinical score
barn2i1c <- barn2
barn2i1c_fit <- loess(clinscore_mod ~ dpi, barn2i1c)
barn2i1c_nd <- data.frame(dpi=seq(min(barn2i1c$dpi), max(barn2i1c$dpi), length=100))
barn2i1c_nd$fit <- predict(barn2i1c_fit, newdata=barn2i1c_nd)
barn2i1cmax <- barn2i1c_nd$dpi[which.max(barn2i1c_nd$fit)]
# plot(barn2i1c$clinscore_mod~barn2i1c$dpi)
# lines(barn2i1c_nd$fit~barn2i1c_nd$dpi, col="grey")


# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barn2 <- ggplot(barn2, aes(dpi, tempc, color="cyan1")) +
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
  geom_vline(aes(xintercept = barn2i1tmax), color = "cyan1") +
  annotate(geom = "text", x = (barn2i1tmax + 1), y = 41, label = paste(round(barn2i1tmax,1)), color = "cyan1") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barn2 <- ggplot(barn2, aes(dpi, clinscore_mod, color="cyan1")) +
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
  geom_vline(aes(xintercept = barn2i1cmax), color = "cyan1") +
  annotate(geom = "text", x = (barn2i1cmax + 1), y = 10, label = paste(round(barn2i1cmax,1)), color = "cyan1") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


#####################
#  Barn 3 ----
#####################
# Barn 3: 3 Experimental Goats
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
rt_barn3 <- ggplot(barn3, aes(dpi, tempc, color="cyan3")) +
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
  geom_vline(aes(xintercept = barn3i1tmax), color = "cyan3") +
  annotate(geom = "text", x = (barn3i1tmax + 1), y = 41, label = paste(round(barn3i1tmax,1)), color = "cyan3") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barn3 <- ggplot(barn3, aes(dpi, clinscore_mod, color="cyan3")) +
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
  geom_vline(aes(xintercept = barn3i1cmax), color = "cyan3") +
  annotate(geom = "text", x = (barn3i1cmax + 1), y = 10, label = paste(round(barn3i1cmax,1)), color = "cyan3") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


#####################
#  Barn 5 ----
#####################
# Barn 5: 3 Negative Control Goats
# Smooth: LOESS 

# PEAKS
# Getting information about where there is a peak each clinical variable is
# rectal temp
barn5i1 <- barn5[barn5$tempc>=36 & barn5$tempc<=41,]
barn5i1_fit <- loess(tempc ~ dpi, barn5i1)
barn5i1_nd <- data.frame(dpi=seq(min(barn5i1$dpi), max(barn5i1$dpi), length=100))
barn5i1_nd$fit <- predict(barn5i1_fit, newdata=barn5i1_nd)
barn5i1tmax <- barn5i1_nd$dpi[which.max(barn5i1_nd$fit)]
#plot(barn5i1$tempc~barn5i1$dpi)
#lines(barn5i1_nd$fit~barn5i1_nd$dpi, col="grey")

# clinical score
barn5i1c <- barn5
barn5i1c_fit <- loess(clinscore_mod ~ dpi, barn5i1c)
barn5i1c_nd <- data.frame(dpi=seq(min(barn5i1c$dpi), max(barn5i1c$dpi), length=100))
barn5i1c_nd$fit <- predict(barn5i1c_fit, newdata=barn5i1c_nd)
barn5i1cmax <- barn5i1c_nd$dpi[which.max(barn5i1c_nd$fit)]
# plot(barn5i1c$clinscore_mod~barn5i1c$dpi)
# lines(barn5i1c_nd$fit~barn5i1c_nd$dpi, col="grey")


# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barn5 <- ggplot(barn5, aes(dpi, tempc, color="grey30")) +
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
  geom_vline(aes(xintercept = barn5i1tmax), color = "grey30") +
  annotate(geom = "text", x = (barn5i1tmax - 1), y = 41, label = paste(round(barn5i1tmax,1)), color = "grey30") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barn5 <- ggplot(barn5, aes(dpi, clinscore_mod, color="grey30")) +
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
  geom_vline(aes(xintercept = barn5i1cmax), color = "grey30") +
  annotate(geom = "text", x = (barn5i1cmax + 1), y = 10, label = paste(round(barn5i1cmax,1)), color = "grey30") +
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
barns2356i1 <- barns2356[barns2356$tempc>=36 & barns2356$tempc<=41,]
barns2356i1_fit <- loess(tempc ~ dpi, barns2356i1)
barns2356i1_nd <- data.frame(dpi=seq(min(barns2356i1$dpi), max(barns2356i1$dpi), length=100))
barns2356i1_nd$fit <- predict(barns2356i1_fit, newdata=barns2356i1_nd)
barns2356i1tmax <- barns2356i1_nd$dpi[which.max(barns2356i1_nd$fit)]
# plot(barns2356i1$tempc~barns2356i1$dpi)
# lines(barns2356i1_nd$fit~barns2356i1_nd$dpi, col="grey")

#clinical score
barns2356i1c <- barns2356
barns2356i1c_fit <- loess(clinscore_mod ~ dpi, barns2356i1c)
barns2356i1c_nd <- data.frame(dpi=seq(min(barns2356i1c$dpi), max(barns2356i1c$dpi), length=100))
barns2356i1c_nd$fit <- predict(barns2356i1c_fit, newdata=barns2356i1c_nd)
barns2356i1cmax <- barns2356i1c_nd$dpi[which.max(barns2356i1c_nd$fit)]
# plot(barns2356i1c$clinscore_mod~barns2356i1c$dpi)
# lines(barns2356i1c_nd$fit~barns2356i1c_nd$dpi, col="grey")


# PLOTTING
# Developing figure for rectal temperature and saving to an object
rt_barns2356 <- ggplot(barns2356, aes(dpi, tempc, color=barn, linetype = barn)) +
  geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x, , k=7), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = expression("Rectal Temperature " ( degree*C))) +
  scale_linetype_manual(name = "", values = c("2" = "solid", "3" = "dashed", "5" = "dotted", "6" = "twodash")) +
  scale_color_manual(name = "", 
                     #labels = c("Goats", "Goats"),
                     values = barncol) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  scale_y_continuous(limits = c(37,41)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barns2356i1tmax), color = "grey50") +
  annotate(geom = "text", x = (barns2356i1tmax + 1), y = 41, label = paste(round(barns2356i1tmax,1)), color = "grey50") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


# Developing figure for clinical score and saving to an object
cs_barns2356 <- ggplot(barns2356, aes(dpi, clinscore_mod, color=barn, linetype = barn)) +
  geom_line(aes(group = eartag, linetype = barn), alpha = 0.4) +
  stat_smooth(method = "loess", formula = y ~ x, size =1) +
  #stat_smooth(method = "gam", formula = y ~ s(x), size =1) +  # code for GAM instead of LOESS (not used)
  labs(x = "Days Post Infection", y = "Modified Clinical Score") +
  scale_linetype_manual(name = "", values = c("2" = "solid", "3" = "dashed", "5" = "dotted", "6" = "twodash")) +
  scale_color_manual(name = "", 
                     #labels = c("Goats", "Goats"),
                     values = barncol) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, 1)) +
  scale_x_continuous(limits = c(0,12), breaks = seq(0:12)) +
  # Next two lines are for plotting peaks, if desired
  geom_vline(aes(xintercept = barns2356i1cmax), color = "grey50") +
  annotate(geom = "text", x = (barns2356i1cmax + 1), y = 10, label = paste(round(barns2356i1cmax,1)), color = "grey50") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), axis.text=element_text(size=6))


########################################
#  Saving figures to file (as jpeg) ----
########################################

jpeg("output/comp_path/pprv_comppath_r1_clin_barn2_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barn2, cs_barn2, ncol=2)
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r1_clin_barn3_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barn3, cs_barn3, ncol=2)
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r1_clin_barn5controls_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barn5, cs_barn5, ncol=2)
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r1_clin_barn6_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barn6, cs_barn6, ncol=2)
invisible(dev.off())

jpeg("output/comp_path/pprv_comppath_r1_clin_barns2356_todpi12.jpeg", width = 6, height = 3, units = "in", quality = 100, res = 600)
grid.arrange(rt_barns2356, cs_barns2356, ncol=2)
invisible(dev.off())

# If need symbols by unique eartag, sample code for this is in PPRV Trial Stage 1.1A code
