############################
######    PACKAGES    ######
############################

library(data.table)
library(lattice)
library(reshape2)
library(dplyr)
library(ggsci)


##################################################################################################
## Plotting spectrum for the BASE FREQUENCIES
##################################################################################################

################################
#######    DATA FILES    #######
################################
rm(list = setdiff(ls(), c("make_list")))	

data_fund <- data.frame(t(fread(".\\z-score g_avg fft base.txt", header = F, sep ="\t")))


#################################
######    CHANNEL LABELS   ######
#################################

channel_labels <- read.csv2("labels.csv", header = F)
channel_labels <- as.vector(channel_labels[,1])

############################
######    DATASETS    ######
############################


colnames(data_fund) <- c(channel_labels,"avg")
rownames(data_fund) <- 1:nrow(data_fund)
data_fund$timepoint <- 1:nrow(data_fund)
head(data_fund)
tail(data_fund)


### Checking for harmonics

idxs = data_fund$avg > 2.32

freqs =c(1: nrow(data_fund))/45.601
freqs[idxs]

base1 <- c(1:10)*2.5
base2 <- c(1:10)*3.16

# for frequency 2.5, the first 10 harmonics are significant
# for frequency 3.16, the first 8 harmonics are significant



####################################################
#######    Plotting spectrum up to 10 Hz    #######
####################################################

last_freq = 10
data_fund_reduced <- data_fund[1:(last_freq*45.601), ] #data_fund[1:456, ] # freq x 45.601
data_fund_reduced_plot <- data.frame(mutate(rowwise(data_fund_reduced), cond = "data_fund", ROI = mean(c(PO4, PO8, P6, P8, O2, PO3, PO7, P5, P7, O1))))

## NOTE: I am using a bilateral ROI because in this data I have 3.16Hz on left electrodes


data_reduced_plot <- rbind(data_fund_reduced_plot)
data_reduced_plot <- dcast(timepoint + cond ~ "ROI", data = data_reduced_plot, value.var = "ROI", mean)
data_reduced_plot$cond <- factor(data_reduced_plot$cond, levels = c("data_fund"))
unique(data_reduced_plot$cond)
#scales::show_col(pal_nejm(alpha = 0.8)(8))
my_cols <- '#000000' #pal_nejm(alpha = 1)(8)

mytheme = trellis.par.get()
mytheme$superpose.line$col = my_cols
mytheme$superpose.line$alpha = 1
mytheme$superpose.symbol$col = my_cols
mytheme$superpose.symbol$alpha = 0.8
mytheme$superpose.polygon$col = my_cols
mytheme$superpose.polygon$alpha = 1
mytheme$superpose.line$lwd = 5
mytheme$fontsize$text = 22
mytheme$axis.text$cex = 1
mytheme$axis.line = list(col = 0, lwd = 3)

x_labels <- as.character(0:(last_freq+1))
x_labels_loc <- c(0.4,ceiling(c(1:last_freq)*45.601))
y_max = max(data_reduced_plot$ROI)
y_labels <- as.character(format(seq(-0.25,1,1),digits = 2))
y_labels_loc <- c(-0.247,seq(0,0.75,0.25),0.997)


# my_plot <- xyplot(ROI ~ timepoint, groups = cond, type = "l", data = data_reduced_plot)
# print(my_plot)

my_plot1 <- xyplot(ROI ~ timepoint, groups = cond, type = "l", data = data_reduced_plot,
	# xlim = c(0.00, 456), ylim = c(-0.25, 1),
	# scales = list(col = 1, tck = c(1,0), x = list(labels = x_labels, at = x_labels_loc), y = list(labels = y_labels, at = y_labels_loc)),
	scales = list(col = 1, tck = c(1,0), x = list(labels = x_labels, at = x_labels_loc)),
	xlab = list("Hz", cex = 1.1), ylab = list("z-scored Amplitudes", cex = 1),
	par.settings = mytheme,
	panel = function(x, y, ...){
		panel.xyplot(x, y, ...)
		# panel.arrows(x0 = c(117,233), x1 = c(117,233), y0 = c(0.15+0.05,0.05+0.05), y1 = c(0.15+0.15,0.05+0.15), code = 1, fill = "black", type = "closed", angle = 30, length = 0.3, lwd = 3, unit = "cm")
		# panel.text(x = c(117,233), y = c(0.15+0.2,0.05+0.2), labels = c("2.4 Hz","4.8 Hz"), font = 2, cex = 0.75)
		lims <- current.panel.limits()
		panel.abline(h = 2.32, lwd = 2, col = 'gray', lty = 2)
		panel.abline(v = lims$xlim[1], lwd = 6)
		panel.abline(h = lims$ylim[1], lwd = 6)
	}
)

ppi <- 300

print(my_plot1)
tiff(filename = "Spectrum_base_10.tiff", compression = "zip", width = 7*2*ppi, height = 7*ppi, pointsize = 7*(9/5), res = ppi)
print(my_plot1)

dev.off()





##################################################################################################
## Plotting spectrum for the IM RESPONSE
##################################################################################################

rm(list = setdiff(ls(), c("make_list", "my_plot1")))	


data_im <- data.frame(t(fread(".\\z-score IM.txt", header = F, sep ="\t")))

#################################
######    CHANNEL LABELS   ######
#################################

channel_labels <- read.csv2("labels.csv", header = F)
channel_labels <- as.vector(channel_labels[,1])

############################
######    DATASETS    ######
############################


colnames(data_im) <- c(channel_labels,"avg")
rownames(data_im) <- 1:nrow(data_im)
data_im$timepoint <- 1:nrow(data_im)
head(data_im)
tail(data_im)

####################################################
#######    Plotting spectrum up to 10 Hz    #######
####################################################

last_freq = 10
data_im_reduced <- data_im[1:(last_freq*45.601), ]
data_im_reduced_plot <- data.frame(mutate(rowwise(data_im_reduced), cond = "data_im", ROI = mean(c(PO4, PO8, O2, Oz, OI1h, P6, P8, FCz, FC1, FC2, F1, Fz, F2))))


data_reduced_plot <- rbind(data_im_reduced_plot)
data_reduced_plot <- dcast(timepoint + cond ~ "ROI", data = data_reduced_plot, value.var = "ROI", mean)
data_reduced_plot$cond <- factor(data_reduced_plot$cond, levels = c("data_im"))
unique(data_reduced_plot$cond)
scales::show_col(pal_nejm(alpha = 0.8)(8))
my_cols <- '#000000' #pal_nejm(alpha = 1)(8)

mytheme = trellis.par.get()
mytheme$superpose.line$col = my_cols
mytheme$superpose.line$alpha = 1
mytheme$superpose.symbol$col = my_cols
mytheme$superpose.symbol$alpha = 0.8
mytheme$superpose.polygon$col = my_cols
mytheme$superpose.polygon$alpha = 1
mytheme$superpose.line$lwd = 5
mytheme$fontsize$text = 22
mytheme$axis.text$cex = 1
mytheme$axis.line = list(col = 0, lwd = 3)

x_labels <- as.character(0:(last_freq+1))
x_labels_loc <- c(0.4,ceiling(c(1:last_freq)*45.601))
y_max = max(data_reduced_plot$ROI)
y_labels <- as.character(format(seq(-0.25,1,1),digits = 2))
y_labels_loc <- c(-0.247,seq(0,0.75,0.25),0.997)


my_plot2 <- xyplot(ROI ~ timepoint, groups = cond, type = "l", data = data_reduced_plot,
                  # xlim = c(0.00, 456), ylim = c(-0.25, 1),
                  # scales = list(col = 1, tck = c(1,0), x = list(labels = x_labels, at = x_labels_loc), y = list(labels = y_labels, at = y_labels_loc)),
                  scales = list(col = 1, tck = c(1,0), x = list(labels = x_labels, at = x_labels_loc)),
                  xlab = list("Hz", cex = 1.1), ylab = list("z-scored Amplitudes", cex = 1),
                  par.settings = mytheme,
                  panel = function(x, y, ...){
                    panel.xyplot(x, y, ...)
                    # panel.arrows(x0 = c(117,233), x1 = c(117,233), y0 = c(0.15+0.05,0.05+0.05), y1 = c(0.15+0.15,0.05+0.15), code = 1, fill = "black", type = "closed", angle = 30, length = 0.3, lwd = 3, unit = "cm")
                    # panel.text(x = c(117,233), y = c(0.15+0.2,0.05+0.2), labels = c("2.4 Hz","4.8 Hz"), font = 2, cex = 0.75)
                    lims <- current.panel.limits()
                    panel.abline(h = 2.32, lwd = 2, col = 'gray', lty = 2)
                    panel.abline(v = lims$xlim[1], lwd = 6)
                    panel.abline(h = lims$ylim[1], lwd = 6)
                  }
)

ppi <- 300
print(my_plot)

tiff(filename = "Spectrum_IM_10.tiff", compression = "zip", width = 7*2*ppi, height = 7*ppi, pointsize = 7*(9/5), res = ppi)

print(my_plot2)
dev.off()



##################################################################################################
## Repeating everything but with the full spectrum up to 30
##################################################################################################

####################################
######    DATA PREPARATION    ######
####################################

rm(list = ls())

################################
#######    DATA FILES    #######
################################

data_fund <- data.frame(t(fread(".\\z-score g_avg fft base.txt", header = F, sep ="\t")))

#################################
######    CHANNEL LABELS   ######
#################################

channel_labels <- read.csv2("labels.csv", header = F)
channel_labels <- as.vector(channel_labels[,1])

############################
######    DATASETS    ######
############################


colnames(data_fund) <- c(channel_labels,"avg")
rownames(data_fund) <- 1:nrow(data_fund)
data_fund$timepoint <- 1:nrow(data_fund)
head(data_fund)
tail(data_fund)

last_freq = 30
data_fund_reduced <- data_fund[1:(last_freq*45.601), ] #data_fund[1:456, ] # freq x 45.601
data_fund_reduced_plot <- data.frame(mutate(rowwise(data_fund_reduced), cond = "data_fund", ROI = mean(c(PO4, PO8, P6, P8, O2, PO3, PO7, P5, P7, O1))))


data_reduced_plot <- rbind(data_fund_reduced_plot)
data_reduced_plot <- dcast(timepoint + cond ~ "ROI", data = data_reduced_plot, value.var = "ROI", mean)
data_reduced_plot$cond <- factor(data_reduced_plot$cond, levels = c("data_fund"))
unique(data_reduced_plot$cond)
#scales::show_col(pal_nejm(alpha = 0.8)(8))
my_cols <- '#000000' #pal_nejm(alpha = 1)(8)

mytheme = trellis.par.get()
mytheme$superpose.line$col = my_cols
mytheme$superpose.line$alpha = 1
mytheme$superpose.symbol$col = my_cols
mytheme$superpose.symbol$alpha = 0.8
mytheme$superpose.polygon$col = my_cols
mytheme$superpose.polygon$alpha = 1
mytheme$superpose.line$lwd = 5
mytheme$fontsize$text = 22
mytheme$axis.text$cex = 1
mytheme$axis.line = list(col = 0, lwd = 3)

x_labels <- as.character(0:(last_freq+1))
x_labels_loc <- c(0.4,ceiling(c(1:last_freq)*45.601))
y_max = max(data_reduced_plot$ROI)
y_labels <- as.character(format(seq(-0.25,1,1),digits = 2))
y_labels_loc <- c(-0.247,seq(0,0.75,0.25),0.997)


# my_plot <- xyplot(ROI ~ timepoint, groups = cond, type = "l", data = data_reduced_plot)
# print(my_plot)

my_plot1 <- xyplot(ROI ~ timepoint, groups = cond, type = "l", data = data_reduced_plot,
                  # xlim = c(0.00, 456), ylim = c(-0.25, 1),
                  # scales = list(col = 1, tck = c(1,0), x = list(labels = x_labels, at = x_labels_loc), y = list(labels = y_labels, at = y_labels_loc)),
                  scales = list(col = 1, tck = c(1,0), x = list(labels = x_labels, at = x_labels_loc)),
                  xlab = list("Hz", cex = 1.1), ylab = list("z-scored Amplitudes", cex = 1),
                  par.settings = mytheme,
                  panel = function(x, y, ...){
                    panel.xyplot(x, y, ...)
                    # panel.arrows(x0 = c(117,233), x1 = c(117,233), y0 = c(0.15+0.05,0.05+0.05), y1 = c(0.15+0.15,0.05+0.15), code = 1, fill = "black", type = "closed", angle = 30, length = 0.3, lwd = 3, unit = "cm")
                    # panel.text(x = c(117,233), y = c(0.15+0.2,0.05+0.2), labels = c("2.4 Hz","4.8 Hz"), font = 2, cex = 0.75)
                    lims <- current.panel.limits()
                    panel.abline(h = 2.32, lwd = 2, col = 'gray', lty = 2)
                    panel.abline(v = lims$xlim[1], lwd = 6)
                    panel.abline(h = lims$ylim[1], lwd = 6)
                  }
)

ppi <- 300
print(my_plot1)

tiff(filename = "Spectrum_base_30.tiff", compression = "zip", width = 7*2*ppi, height = 7*ppi, pointsize = 7*(9/5), res = ppi)

print(my_plot1)

dev.off()





##################################################################################################
## IM response
##################################################################################################
rm(list = setdiff(ls(), c("make_list", "my_plot1")))	


data_im <- data.frame(t(fread(".\\z-score IM.txt", header = F, sep ="\t")))

#################################
######    CHANNEL LABELS   ######
#################################

channel_labels <- read.csv2("labels.csv", header = F)
channel_labels <- as.vector(channel_labels[,1])

############################
######    DATASETS    ######
############################


colnames(data_im) <- c(channel_labels,"avg")
rownames(data_im) <- 1:nrow(data_im)
data_im$timepoint <- 1:nrow(data_im)
head(data_im)
tail(data_im)


last_freq = 30
data_im_reduced <- data_im[1:(last_freq*45.601), ]
data_im_reduced_plot <- data.frame(mutate(rowwise(data_im_reduced), cond = "data_im", ROI = mean(c(PO4, PO8, O2, Oz, OI1h, P6, P8, FCz, FC1, FC2, F1, Fz, F2))))

data_reduced_plot <- rbind(data_im_reduced_plot)
data_reduced_plot <- dcast(timepoint + cond ~ "ROI", data = data_reduced_plot, value.var = "ROI", mean)
data_reduced_plot$cond <- factor(data_reduced_plot$cond, levels = c("data_im"))
unique(data_reduced_plot$cond)
#scales::show_col(pal_nejm(alpha = 0.8)(8))
my_cols <- '#000000' #pal_nejm(alpha = 1)(8)

mytheme = trellis.par.get()
mytheme$superpose.line$col = my_cols
mytheme$superpose.line$alpha = 1
mytheme$superpose.symbol$col = my_cols
mytheme$superpose.symbol$alpha = 0.8
mytheme$superpose.polygon$col = my_cols
mytheme$superpose.polygon$alpha = 1
mytheme$superpose.line$lwd = 5
mytheme$fontsize$text = 22
mytheme$axis.text$cex = 1
mytheme$axis.line = list(col = 0, lwd = 3)

x_labels <- as.character(0:(last_freq+1))
x_labels_loc <- c(0.4,ceiling(c(1:last_freq)*45.601))
y_max = max(data_reduced_plot$ROI)
y_labels <- as.character(format(seq(-0.25,1,1),digits = 2))
y_labels_loc <- c(-0.247,seq(0,0.75,0.25),0.997)

# my_plot <- xyplot(ROI ~ timepoint, groups = cond, type = "l", data = data_reduced_plot)
# print(my_plot)

my_plot2 <- xyplot(ROI ~ timepoint, groups = cond, type = "l", data = data_reduced_plot,
                  # xlim = c(0.00, 456), ylim = c(-0.25, 1),
                  # scales = list(col = 1, tck = c(1,0), x = list(labels = x_labels, at = x_labels_loc), y = list(labels = y_labels, at = y_labels_loc)),
                  scales = list(col = 1, tck = c(1,0), x = list(labels = x_labels, at = x_labels_loc)),
                  xlab = list("Hz", cex = 1.1), ylab = list("z-scored Amplitudes", cex = 1),
                  par.settings = mytheme,
                  panel = function(x, y, ...){
                    panel.xyplot(x, y, ...)
                    # panel.arrows(x0 = c(117,233), x1 = c(117,233), y0 = c(0.15+0.05,0.05+0.05), y1 = c(0.15+0.15,0.05+0.15), code = 1, fill = "black", type = "closed", angle = 30, length = 0.3, lwd = 3, unit = "cm")
                    # panel.text(x = c(117,233), y = c(0.15+0.2,0.05+0.2), labels = c("2.4 Hz","4.8 Hz"), font = 2, cex = 0.75)
                    lims <- current.panel.limits()
                    panel.abline(h = 2.32, lwd = 2, col = 'gray', lty = 2)
                    panel.abline(v = lims$xlim[1], lwd = 6)
                    panel.abline(h = lims$ylim[1], lwd = 6)
                  }
)

ppi <- 300
tiff(filename = "Spectrum_IM_30.tiff", compression = "zip", width = 7*2*ppi, height = 7*ppi, pointsize = 7*(9/5), res = ppi)

print(my_plot2)

dev.off()

