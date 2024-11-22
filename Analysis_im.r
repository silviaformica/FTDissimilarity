############################
######    PACKAGES    ######
############################

library(lattice)
library(latticeExtra)
library(reshape2)
library(car)
library(plyr)
library(dplyr)
library(magrittr)
library(heplots)
library(eegUtils)
library(BayesFactor)
library(data.table)
library(ggplot2)
library(MBESS)
library(lsr)
library(scales)
library(colormap)
library(ggsci)

############################
######    SETTINGS    ######
############################

options(contrasts = c("contr.sum", "contr.poly"))

#############################
######    FUNCTIONS    ######
#############################

######################################
######    READ DATA	FUNCTION    ######
######################################

make_list <- function(ppn){

	###### READ DATA ######
	
	if(ppn < 10){
		file.name <- sprintf(".\\data\\subject%02d_im_sns.txt", ppn)
	}else{
		file.name <- sprintf(".\\data\\subject%d_im_sns.txt", ppn)

	}

	data.file <- read.table(file.name, header = T, sep = ",")

	###### ADD FLUENCY AND SYNCHRONY COLUMNS ######

	data.file <- mutate(data.file, movs = factor(case_when(epoch %in% c(1,2) ~ "diff", epoch %in% c(3,4) ~ "same"), levels = c("diff","same")))
	data.file <- mutate(data.file, agents = factor(case_when(epoch %in% c(1,3) ~ "one", epoch %in% c(2,4) ~ "two"), levels = c("one","two")))
			
	###### CHANGE DATASET TO PPN ######
	
	data.file$dataset <- rep(ppn, times = 65*4)
	names(data.file)[1] <- "ppn"	

	###### RETURN DATA ######
			
	return(data.file)
}

	####################################
	######    DATA PREPARATION    ######
	####################################

#############################
######    READ DATA    ######
#############################

	###### CLEAR WORKING SPACE ######
	
rm(list = setdiff(ls(), c("make_list")))	

	###### MAKE LIST ######
	
my_ppn <- c(1:3,5:14,16:26,28,30:33)
data.list <- lapply(my_ppn, make_list)
sapply(data.list, nrow)
sapply(data.list, function(x) mean(x$H1, na.rm = T))

	###### MAKE DATA FRAME ######
	
data.full <- do.call(rbind, data.list)
head(data.full)
tail(data.full)
summary(data.full)
str(data.full)


######################################
######    READ CHANNEL LABELS   ######
######################################

channel_labels <- read.csv2("labels.csv", header = F)
channel_labels <- as.vector(channel_labels[,1])

data.full$channel <- factor(data.full$channel, levels = 1:65, labels = channel_labels)
str(data.full)	

##############################################
######    REMOVE IRRELEVANT CHANNELS    ######
##############################################

data.full <- filter(data.full, !(channel %in% c("65")))

##########################################
######    CHECK HIGHEST HARMONIC    ######
##########################################

colMeans(select(data.full, "H1","H2","H3"))

###############################################
######    CHECK SIGNIFICANT HARMONICS    ######
###############################################

data.check.sig <- data.full %>%
	filter(channel == "PO8") %>%
	reshape2::dcast(ppn + movs + agents ~ "my_h", value.var = "H1", mean) %>%
	reshape2::dcast(ppn + agents ~ "my_h", value.var = "my_h", mean) %>%
	reshape2::dcast(ppn ~ "my_h", value.var = "my_h", mean)
mean(data.check.sig$my_h)

t.test(data.check.sig$my_h, mu = 0, alternative = "greater")

#############################
######    TOPOPLOTS    ######
#############################

################################
######    COLOR SCALES    ######
################################

jet.colors <- colorRampPalette(c('#000052', '#0C44AC', '#FFFFFF', '#ED0101', '#970005')) # original


###############################################################
######    LOCALIZER across frequencies and conditions   ######
###############################################################

###### MAKE DATAFRAME ######

data.topo.avg <- data.full %>%
  select(channel, movs, agents, ppn, H1) %>%
  data.table %>%
  dcast(channel + movs + agents ~ 'H1', value.var = "H1", mean) %>%
  dcast(channel + movs ~ 'H1', value.var = "H1", mean) %>%
  dcast(channel ~ 'H1', value.var = "H1", mean) %>%
  setDF

###### PLOT ######	

###### DATA ######	

data.topo.avg.plot <- data.topo.avg %>%
  mutate(H_avg = H1) %>%
  select(channel, H_avg) %>%	
  filter(channel != "FT9", channel != "FT10")
names(data.topo.avg.plot) <- c("electrode","amplitude")

my_max <- print(max(data.topo.avg.plot$amplitude))

highlight = c("PO4","PO8","O2","OZ","OI2H","OI1H","P6","P8","FCZ","FC1","FC2","F1","FZ","F2")
my_topo <- topoplot(data.topo.avg.plot,chan_marker = "none", interp_limit = "head", method = "Biharmonic", contour = T)+
  #stat_contour(aes(z = fill,linetype = ggplot2::after_stat(level) < 0),bins = 6,colour = "black",size = rel(1.1 * 1),show.legend = FALSE)+   
  annotate("point", x = filter(electrode_locations(data.topo.avg.plot),!(electrode %in% highlight))$x, y = filter(electrode_locations(data.topo.avg.plot),!(electrode %in% highlight))$y, colour = "black", size = 1.5)+
  annotate("point", x = filter(electrode_locations(data.topo.avg.plot),electrode %in% highlight)$x, y = filter(electrode_locations(data.topo.avg.plot),electrode %in% highlight)$y, colour = "black", fill = "white", pch = 21, stroke = 3, size = 10)+#rel(2 * 1.5))+
  #annotate("point", x = filter(electrode_locations(data.topo.avg.plot),electrode %in% highlight)$x, y = filter(electrode_locations(data.topo.avg.plot),electrode %in% highlight)$y, colour = "black", pch = 3, stroke = 1.5, size = rel(2 * 1.5))+
  scale_fill_gradientn(colours = jet.colors(256), labels = c(0, "max"), breaks = c(0.00, my_max), limits = c(0.00, my_max), oob = squish) + 
  #guides(fill = guide_colorbar(barwidth = 1.5, barheight = 15, label = T, label.theme = element_text(angle = 0, size = 34), ticks = F, title = "", title.position = "right", title.hjust = 0.5, title.theme = element_text(angle = 270)))
  guides(fill = guide_colorbar(barwidth = 0, barheight = 0, label = F, title = ""))
suppressWarnings(print(my_topo))

###### FIND HIGHEST ######

data.topo.avg.plot %>% arrange(desc(amplitude)) %>% slice(1:10)

###### PRINT PLOT ###### 

ppi <- 300	
tiff(filename = "./localizer_im.tiff", compression = "zip", width = 7*1*ppi, height = 7*1*ppi, pointsize = 7*(9/5), res = ppi)
my_topo #+ labs(title = "Base Rate") + theme(plot.title = element_text(size = 50, hjust = 0.5))
dev.off()	


#########################################
######    CONDITION DIFFERENCES    ######
#########################################

main_movs <- reshape2::dcast(ppn + channel + movs + agents ~ "H1", data = data.full, value.var = "H1", mean) %>%
	reshape2::dcast(ppn + channel + movs ~ "H1", value.var = "H1", mean) %>%
	reshape2::dcast(channel ~ movs, value.var = "H1", mean) %>%
	filter(channel != "FT9", channel != "FT10") %>% 
	mutate(amplitude = diff - same) %>% select(-diff,-same) %>% setnames("channel", "electrode")
main_agents <- reshape2::dcast(ppn + channel + movs + agents ~ "H1", data = data.full, value.var = "H1", mean) %>%
	reshape2::dcast(ppn + channel + agents ~ "H1", value.var = "H1", mean) %>%
	reshape2::dcast(channel ~ agents, value.var = "H1", mean) %>%
	filter(channel != "FT9", channel != "FT10") %>% 
	mutate(amplitude = two - one) %>% select(-one,-two) %>% setnames("channel", "electrode")
inter <- reshape2::dcast(ppn + channel + movs + agents ~ "H1", data = data.full, value.var = "H1", mean) %>%
	reshape2::dcast(channel ~ movs + agents, value.var = "H1", mean) %>%
	filter(channel != "FT9", channel != "FT10") %>% 
	mutate(amplitude = (diff_two - diff_one) - (same_two - same_one)) %>% select(-diff_two,-diff_one,-same_two,-same_one) %>% setnames("channel", "electrode")	


#################################################
######    Main Effect: MOVEMENT TYPE      ######
#################################################

my_plot_data <- get("main_movs")	
my_max <- print(max(my_plot_data$amplitude))
my_min <- print(min(my_plot_data$amplitude))
my_limit <- print(max(c(abs(my_max),abs(my_min))))
highlight = c("PO4","PO8","O2","OZ","OI2H","OI1H","P6","P8","FCZ","FC1","FC2","F1","FZ","F2")
my_topo <- topoplot(my_plot_data, chan_marker = "none", interp_limit = "head", contour = T)+
    #stat_contour(aes(z = fill,linetype = ggplot2::after_stat(level) < 0),bins = 6,colour = "black",size = rel(1.1 * 1),show.legend = FALSE)+   
  annotate("point", x = filter(electrode_locations(my_plot_data),!(electrode %in% highlight))$x, y = filter(electrode_locations(my_plot_data),!(electrode %in% highlight))$y, colour = "black", size = 1.5)+
  annotate("point", x = filter(electrode_locations(my_plot_data),electrode %in% highlight)$x, y = filter(electrode_locations(my_plot_data),electrode %in% highlight)$y, colour = "black", fill = "white", pch = 21, stroke = 3, size = 10)+#rel(2 * 1.5))+
  #annotate("point", x = filter(electrode_locations(data.topo.cond.plot),electrode %in% highlight)$x, y = filter(electrode_locations(data.topo.cond.plot),electrode %in% highlight)$y, colour = "black", pch = 3, stroke = 1.5, size = rel(2 * 1.5))+
	scale_fill_gradientn(colours = jet.colors(256), labels = c(round(-my_limit, 3), round(my_limit, 3)), breaks = c(-my_limit, my_limit), limits = c(-my_limit, my_limit), oob = squish) + 
	guides(fill = guide_colorbar(barwidth = 1.5, barheight = 15, label = T, label.theme = element_text(angle = 0, size = 34), ticks = F, title = "", title.position = "right", title.hjust = 0.5, title.theme = element_text(angle = 270)))
	guides(fill = guide_colorbar(barwidth = 0, barheight = 0, label = F, title = "")) 
suppressWarnings(print(my_topo))


ppi <- 300	
tiff(filename = "./EffectMovs_im.tiff", compression = "zip", width = 7*1*ppi, height = 7*1*ppi, pointsize = 7*(9/5), res = ppi)
my_topo #+ labs(title = "Base Rate") + theme(plot.title = element_text(size = 50, hjust = 0.5))
dev.off()	


#########################################
######    Main Effect: AGENTS      ######
#########################################

my_plot_data <- get("main_agents")	
my_max <- print(max(my_plot_data$amplitude))
my_min <- print(min(my_plot_data$amplitude))
my_limit <- print(max(c(abs(my_max),abs(my_min))))
my_topo <- topoplot(my_plot_data, chan_marker = "none", interp_limit = "head", contour = T)+
  #stat_contour(aes(z = fill,linetype = ggplot2::after_stat(level) < 0),bins = 6,colour = "black",size = rel(1.1 * 1),show.legend = FALSE)+   
  annotate("point", x = filter(electrode_locations(my_plot_data),!(electrode %in% highlight))$x, y = filter(electrode_locations(my_plot_data),!(electrode %in% highlight))$y, colour = "black", size = 1.5)+
  annotate("point", x = filter(electrode_locations(my_plot_data),electrode %in% highlight)$x, y = filter(electrode_locations(my_plot_data),electrode %in% highlight)$y, colour = "black", fill = "white", pch = 21, stroke = 3, size = 10)+#rel(2 * 1.5))+
  #annotate("point", x = filter(electrode_locations(data.topo.cond.plot),electrode %in% highlight)$x, y = filter(electrode_locations(data.topo.cond.plot),electrode %in% highlight)$y, colour = "black", pch = 3, stroke = 1.5, size = rel(2 * 1.5))+
  scale_fill_gradientn(colours = jet.colors(256), labels = c(round(-my_limit, 3), round(my_limit, 3)), breaks = c(-my_limit, my_limit), limits = c(-my_limit, my_limit), oob = squish) + 
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 15, label = T, label.theme = element_text(angle = 0, size = 34), ticks = F, title = "", title.position = "right", title.hjust = 0.5, title.theme = element_text(angle = 270)))
  guides(fill = guide_colorbar(barwidth = 0, barheight = 0, label = F, title = "")) 
suppressWarnings(print(my_topo))


ppi <- 300	
tiff(filename = "./EffectAgents_im.tiff", compression = "zip", width = 7*1*ppi, height = 7*1*ppi, pointsize = 7*(9/5), res = ppi)
my_topo #+ labs(title = "Base Rate") + theme(plot.title = element_text(size = 50, hjust = 0.5))
dev.off()	


	
#########################
######    BOXPLOTS    ######
#########################

#########################
######    THEME    ######
#########################

mytheme = trellis.par.get() 
#scales::show_col(ggsci::pal_nejm(alpha = 0.7)(8))	
mytheme$superpose.symbol = list(fill = c("#5F5FCF", '#E2E2F6'),col= "black") # inside boxplot	
mytheme$superpose.polygon$col = c("#5F5FCF", '#E2E2F6') # legend	
mytheme$box.umbrella = list(col = "black", lwd = 2) # whiskers	
mytheme$box.rectangle = list(col = "black", lwd = 2) # box outline
mytheme$plot.symbol = list(col = "black", pch = 16, cex = 1) # outliers
mytheme$add.line$lwd = 1 # background lines
mytheme$fontsize$text = 20
mytheme$axis.text$cex = 1
mytheme$axis.line$lwd = 2

########################
######    PLOT    ######
########################

data.plot.H1.ppn <- data.full %>%
	#filter(channel %in% c("PO4","PO8","O2","Oz","OI2h","OI1h","P6","P8")) %>%
	#filter(channel %in% c("FCz","FC1","FC2","F1","Fz","F2")) %>%
	filter(channel %in% c("PO4","PO8","O2","Oz","OI2h","OI1h","P6","P8","FCz","FC1","FC2","F1","Fz","F2")) %>%
	reshape2::dcast(ppn + movs + agents ~ "H1", value.var = "H1", mean)
data.plot.H1.avg <- reshape2::dcast(movs + agents ~ "H1", data = data.plot.H1.ppn, value.var = "H1", mean)

boxplot.H1 <- bwplot(H1 ~ agents, groups = movs, data = data.plot.H1.ppn,
	ylim = c(-0.15,0.3), scales = list(tck = c(1,0), x = list(labels = c("1P", "2P")), y = list(at = round(seq(-0.15, 0.25, 0.05),2))),
	ylab = "Baseline-Subtracted Amplitudes",
	auto.key = list(corner = c(0.5,.975), text = c("Different","Same"), points = F, rectangles = T, columns = 2, between.columns = 0.75, height = 0.85, size = 3.5, padding.text = 2.5),
	panel = panel.superpose, par.settings = mytheme,
	panel.groups = function(x, y, pch, ... , group.number) {
		panel.abline(h = seq(-0.15,0.25, 0.025), col="lightgrey", lty = 3, alpha = 0.5)
	})+
bwplot(H1 ~ agents, groups = movs, data = data.plot.H1.ppn, 
    box.width = 1/4, do.out = T,
    panel = panel.superpose, par.settings = mytheme,
    panel.groups = function(x, y, pch, ... , group.number, cex) {		
    	panel.bwplot(x + (group.number-1.5)/3.5, y, pch = "", ...) #"x": 1 to 2, "group.number": 1 to 2; "(group.number-1.5)": middle group has to be 0; "/ 3.5": "(box.width^-1)-0.5"
    })+
xyplot(H1 ~ agents, groups = movs, data = data.plot.H1.avg, type = c("p"), pch = 16, cex = 1.25, lwd = 2,
	panel = panel.superpose, 
	panel.groups = function(x, y, ... , group.number){
		panel.xyplot(x + (group.number-1.5)/3.5, y, ...) #x: 1 or 2, group.number: 1 or 2
		panel.segments(x0 = (x + (group.number-1.5)/3.5) - 1/8, x1 = (x + (group.number-1.5)/3.5) + 1/8, y0 = y, y1 = y, ...)
	})
print(boxplot.H1)

ppi <- 300	
tiff(filename = "./BoxPlot_im.tiff", compression = "zip", width = 7*1*ppi, height = 7*1*ppi, pointsize = 7*(9/5), res = ppi)
print(boxplot.H1)
dev.off()	

# dev.off()

	############################
	######    ANALYSIS    ######
	############################

###############################
######    WITH REGION    ######
###############################

	###### ANOVA ######

data.analysis.H1 <- data.full %>%
	filter(channel %in% c("PO4","PO8","O2","Oz","OI2h","OI1h","P6","P8","FCz","FC1","FC2","F1","Fz","F2")) %>%
	mutate(region = factor(ifelse(channel %in% c("PO4","PO8","O2","Oz","OI2h","OI1h","P6","P8"), "post","front"))) %>%
	reshape2::dcast(ppn ~ region + movs + agents, value.var = "H1", mean)
colMeans(data.analysis.H1)

region <- factor(rep(1:2,each = 4))
movs <- factor(rep(1:2, each = 2,times=2))
agents <- factor(rep(1:2, times = 4))
idata <- data.frame(region,movs,agents)

fit.H1 <- lm(cbind(front_diff_one, front_diff_two, front_same_one, front_same_two, post_diff_one, post_diff_two, post_same_one, post_same_two) ~ 1, data = data.analysis.H1)
etasq(Anova(fit.H1, type = "III", test = "Wilks", idata = idata, idesign = ~region*movs*agents), anova = T, partial = T)
etasq(Anova(fit.H1, type = "III", test = "Wilks", idata = idata, idesign = ~region*movs*agents), anova = T, partial = T)$"approx F"


# cohens d for main effects of agents
sqrt(6.8816)/sqrt(29)

# cohens d for main effects of movs
sqrt(6.1990)/sqrt(29)


## computing means
two_ag = c(data.analysis.H1$post_same_two, data.analysis.H1$post_diff_two, data.analysis.H1$front_same_two, data.analysis.H1$front_diff_two)
mean(two_ag)
sd(two_ag)

one_ag = c(data.analysis.H1$post_diff_one, data.analysis.H1$post_same_one, data.analysis.H1$front_same_one, data.analysis.H1$front_diff_one)
mean(one_ag)
sd(one_ag)

## computing means
same = c(data.analysis.H1$post_same_two, data.analysis.H1$post_same_one, data.analysis.H1$front_same_two, data.analysis.H1$front_same_one)
mean(same)
sd(same)

diff = c(data.analysis.H1$post_diff_one, data.analysis.H1$post_diff_two, data.analysis.H1$front_diff_one, data.analysis.H1$front_diff_two)
mean(diff)
sd(diff)


#####################################
#### EXTRA 
#####################################

## Plotting the boxplots for the two regions separately

data.plot.H1.ppn <- data.full %>%
  filter(channel %in% c("PO4","PO8","O2","Oz","OI2h","OI1h","P6","P8","FCz","FC1","FC2","F1","Fz","F2")) %>%
  mutate(region = factor(ifelse(channel %in% c("PO4","PO8","O2","Oz","OI2h","OI1h","P6","P8"), "post","front"))) %>%
  reshape2::dcast(ppn + region + movs + agents ~ "H1", value.var = "H1", mean)


## plotting only the frontal

boxplot.H1 <- bwplot(H1 ~ agents, groups = movs, data = data.plot.H1.ppn[data.plot.H1.ppn$region == 'front',],
                     ylim = c(-0.15,0.3), scales = list(tck = c(1,0), x = list(labels = c("1P", "2P")), y = list(at = round(seq(-0.15, 0.25, 0.05),2))),
                     ylab = "Baseline-Subtracted Amplitudes",
                     main = 'Frontal ROI',
                     auto.key = list(corner = c(0.5,.975), text = c("Different","Same"), points = F, rectangles = T, columns = 2, between.columns = 0.75, height = 0.85, size = 3.5, padding.text = 2.5),
                     panel = panel.superpose, par.settings = mytheme,
                     panel.groups = function(x, y, pch, ... , group.number) {
                       panel.abline(h = seq(-0.15,0.25, 0.025), col="lightgrey", lty = 3, alpha = 0.5)
                     })+
  bwplot(H1 ~ agents, groups = movs, data = data.plot.H1.ppn[data.plot.H1.ppn$region == 'front',], 
         box.width = 1/4, do.out = T,
         panel = panel.superpose, par.settings = mytheme,
         panel.groups = function(x, y, pch, ... , group.number, cex) {		
           panel.bwplot(x + (group.number-1.5)/3.5, y, pch = "", ...) #"x": 1 to 2, "group.number": 1 to 2; "(group.number-1.5)": middle group has to be 0; "/ 3.5": "(box.width^-1)-0.5"
         })+
  xyplot(H1 ~ agents, groups = movs, data = data.plot.H1.avg, type = c("p"), pch = 16, cex = 1.25, lwd = 2,
         panel = panel.superpose, 
         panel.groups = function(x, y, ... , group.number){
           panel.xyplot(x + (group.number-1.5)/3.5, y, ...) #x: 1 or 2, group.number: 1 or 2
           panel.segments(x0 = (x + (group.number-1.5)/3.5) - 1/8, x1 = (x + (group.number-1.5)/3.5) + 1/8, y0 = y, y1 = y, ...)
         })
# print(boxplot.H1)


## plotting only the posterior

boxplot.H2 <- bwplot(H1 ~ agents, groups = movs, data = data.plot.H1.ppn[data.plot.H1.ppn$region == 'post',],
                     ylim = c(-0.15,0.3), scales = list(tck = c(1,0), x = list(labels = c("1P", "2P")), y = list(at = round(seq(-0.15, 0.25, 0.05),2))),
                     ylab = "Baseline-Subtracted Amplitudes",
                     main = 'Posterior ROI',
                     auto.key = list(corner = c(0.5,.975), text = c("Different","Same"), points = F, rectangles = T, columns = 2, between.columns = 0.75, height = 0.85, size = 3.5, padding.text = 2.5),
                     panel = panel.superpose, par.settings = mytheme,
                     panel.groups = function(x, y, pch, ... , group.number) {
                       panel.abline(h = seq(-0.15,0.25, 0.025), col="lightgrey", lty = 3, alpha = 0.5)
                     })+
  bwplot(H1 ~ agents, groups = movs, data = data.plot.H1.ppn[data.plot.H1.ppn$region == 'post',], 
         box.width = 1/4, do.out = T,
         panel = panel.superpose, par.settings = mytheme,
         panel.groups = function(x, y, pch, ... , group.number, cex) {		
           panel.bwplot(x + (group.number-1.5)/3.5, y, pch = "", ...) #"x": 1 to 2, "group.number": 1 to 2; "(group.number-1.5)": middle group has to be 0; "/ 3.5": "(box.width^-1)-0.5"
         })+
  xyplot(H1 ~ agents, groups = movs, data = data.plot.H1.avg, type = c("p"), pch = 16, cex = 1.25, lwd = 2,
         panel = panel.superpose, 
         panel.groups = function(x, y, ... , group.number){
           panel.xyplot(x + (group.number-1.5)/3.5, y, ...) #x: 1 or 2, group.number: 1 or 2
           panel.segments(x0 = (x + (group.number-1.5)/3.5) - 1/8, x1 = (x + (group.number-1.5)/3.5) + 1/8, y0 = y, y1 = y, ...)
         })
# print(boxplot.H1)


if (!require(gridExtra)) {
  install.packages("gridExtra")
}
grid.arrange(boxplot.H1, boxplot.H2, nrow = 1)
ppi <- 300	

tiff(filename = "./boxPlot_IM_byROI.tiff", compression = "zip", width = 14*1*ppi, height = 7*1*ppi, pointsize = 7*(9/5), res = ppi)
grid.arrange(boxplot.H1, boxplot.H2, nrow = 1)

dev.off()	


