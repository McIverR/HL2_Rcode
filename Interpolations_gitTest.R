# Find the Right interpolation method for the functional groups of Algae at HL-2 1999-2016 time series - mean cell abundance grouped by month. ###

## Update : September 26, 2017: Tested interpolations include linear, spline, Stienman, and multiple weighting types for Moving Weighted Mean interpolation (K-2, larger is Not valid for our data - monthly grouped phytoplankton)
# --> RM choosing MWM interpolation to go forward with 

#Missing data interpolation using package imputeTS in R --> multiple algorithms for data filling for univariate time series
#Description Imputation (replacement) of missing values in univariate time series. https://cran.r-project.org/web/packages/imputeTS/imputeTS.pdf
# REF: Moritz, Steffen, et al. "Comparison of different Methods for Univariate Time Series Imputation in R." arXiv preprint arXiv:1510.03924 (2015).


#Offers several imputation functions and missing data plots. Available imputation algorithms include: 'Mean', 'LOCF', 'Interpolation','Moving Average', 'Seasonal Decomposition', 'Kalman Smoothing on Structural Time Series models', 'Kalman Smoothing on ARIMA models'.

##To impute (fill all missing values) in a time series x, run:
#na.interpolation(x)
##To plot missing data statistics for a time series x, run: 
#plotNA.distribution(x)
##To print missing data statistics for a time series x, run:
#statsNA(x)

#Every other imputation function (starting with na.’algorithm name’) and plotting function (starting with plotNA.’plot name’) work the same way as in this example.

# begin: visualize distribution of missing values in the series (This function visualizes the distribution of missing values within a time series. Therefore, the time series is plotted and whenever a value is NA the background is colored differently. This gives a nice overview, where in the time series most of the missing values occur)

#1 NA interpolation - options = linear, splines, Stineman

#2 Moving Mean Window interpolation - options = linear, linearm exponential, and can specify the width of window

#3 Seasonally decomposed missing value imputation - removes seasonal component, and imputation by chosen algorithm


## Create the data frame, subsetted as wanted
rm(list=ls())
setwd("~/phytoplanktontaxonomy/spreadsheets")
#rawdat = read.csv("Phyto_HL2-1999-2015.csv") #### This will bring up the CSV without the column with genus groupings
rawdat = read.csv("Phyto_HL2-1999-2016_spGroupedit.csv") ### This will bring up the edited CSV with the correct column and genus grouping. Other abnormalities with the dataset (misspellings etc) have also been fixed.
rawdat$Date <- (as.Date(substr(as.character(rawdat$Date),1,10),format="%d/%m/%Y")) # remove time stamp from the date and time
rawdat$Month <- month(ymd(rawdat$Date), label = TRUE, abbr = FALSE) # turn month integers into month names (full or abbreviated)
rawdat$Month_integer <- (substr(as.character(rawdat$Date),6,7)) # Add month column to the table
head(rawdat)

#### Create the Year Month Interval for the frequent and most abundant subset dataframes ####
rawdat$year_month <- paste(rawdat$Year, rawdat$Month_integer, sep = "")

mon <- (as.character(01:12))
mon[mon == "1"] <- "01"
mon[mon == "2"] <- "02"
mon[mon == "3"] <- "03"
mon[mon == "4"] <- "04"
mon[mon == "5"] <- "05"
mon[mon == "6"] <- "06"
mon[mon == "7"] <- "07"
mon[mon == "8"] <- "08"
mon[mon == "9"] <- "09"

Year <- as.character(rep(1999:2016, each = 52))

#### matrix with mean cells/L/month from most frequent occuring genus', monthly resolution ####
code = paste(Year,mon, sep="")
ids = sort(unique(code))
nbs = length(ids)

abun_sort <- arrange(rawdat, year_month) ## The dataframe with all observations sorted in ascending year_month
head(abun_sort)

top_abun = matrix(NaN, nbs, 10) # for each of the 10 'functional' groups
j = 10
for (j in 1:nbs)
{
  ind = abun_sort$CATEGORY == 1 & abun_sort$year_month == ids[j]
  top_abun[j,1] = mean(abun_sort$CELLS_LITRE[ind])
  ind = abun_sort$CATEGORY == 2 & abun_sort$year_month == ids[j]
  top_abun[j,2] = mean(abun_sort$CELLS_LITRE[ind])
  ind = abun_sort$CATEGORY == 3 & abun_sort$year_month == ids[j]
  top_abun[j,3] = mean(abun_sort$CELLS_LITRE[ind])
  ind = abun_sort$CATEGORY == 4 & abun_sort$year_month == ids[j]
  top_abun[j,4] = mean(abun_sort$CELLS_LITRE[ind])
  ind = abun_sort$CATEGORY == 5 & abun_sort$year_month == ids[j]
  top_abun[j,5] = mean(abun_sort$CELLS_LITRE[ind])
  ind = abun_sort$CATEGORY == 6 & abun_sort$year_month == ids[j]
  top_abun[j,6] = mean(abun_sort$CELLS_LITRE[ind])
  ind = abun_sort$CATEGORY == 7 & abun_sort$year_month == ids[j]
  top_abun[j,7] = mean(abun_sort$CELLS_LITRE[ind])
  ind = abun_sort$CATEGORY == 8 & abun_sort$year_month == ids[j]
  top_abun[j,8] = mean(abun_sort$CELLS_LITRE[ind])
  ind = abun_sort$CATEGORY == 9 & abun_sort$year_month == ids[j]
  top_abun[j,9] = mean(abun_sort$CELLS_LITRE[ind])
  ind = abun_sort$CATEGORY == 10 & abun_sort$year_month == ids[j]
  top_abun[j,10] = mean(abun_sort$CELLS_LITRE[ind])
  
}

head(top_abun)

abun_df <- as.data.frame(top_abun)
mon_vec <- sort(as.character(ids))
abun_df <- cbind(abun_df, mon_vec)
colnames(abun_df) <- c("Centric Diatom", "Pennate Diatom", "HAB Diatom", "HAB Flagellate", "Thecate/Non-Theca", "Silicoflagellates", "Flagellates & Chryptophores", "Ciliates", "Euglena", "Terres. Fragments", "Year_month")
head(abun_df) 

# Want to get index of the rows that are all = zero (no sampling so NA) vs those that just have no values (zero values) --> distingushes between zero abundance and zero sampling effort (NA)
abun_df[is.na(abun_df)] <- 0

abun_df$sum_all <- rowSums(abun_df[,c("Centric Diatom", "Pennate Diatom", "HAB Diatom", "HAB Flagellate", "Thecate/Non-Theca", "Silicoflagellates", "Flagellates & Chryptophores", "Ciliates", "Euglena", "Terres. Fragments")])
abun_df$sum_all[abun_df$sum_all == 0] <- NA
abun_df$sum_all

abun_df$`Centric Diatom` <- ifelse(is.na(abun_df$sum_all), NA, abun_df$`Centric Diatom`) 
abun_df$`Pennate Diatom` <- ifelse(is.na(abun_df$sum_all), NA, abun_df$`Pennate Diatom`) 
abun_df$`HAB Diatom`<- ifelse(is.na(abun_df$sum_all), NA, abun_df$`HAB Diatom`) 
abun_df$`HAB Flagellate`<- ifelse(is.na(abun_df$sum_all), NA, abun_df$`HAB Flagellate`) 
abun_df$`Thecate/Non-Theca`<- ifelse(is.na(abun_df$sum_all), NA, abun_df$`Thecate/Non-Theca`) 
abun_df$Silicoflagellates<- ifelse(is.na(abun_df$sum_all), NA, abun_df$Silicoflagellates)
abun_df$`Flagellates & Chryptophores`<- ifelse(is.na(abun_df$sum_all), NA, abun_df$`Flagellates & Chryptophores` ) 
abun_df$Ciliates<- ifelse(is.na(abun_df$sum_all), NA, abun_df$Ciliates ) 
abun_df$Euglena<- ifelse(is.na(abun_df$sum_all), NA, abun_df$Euglena) 
abun_df$`Terres. Fragments`<- ifelse(is.na(abun_df$sum_all), NA, abun_df$`Terres. Fragments`) 


#### plot the data - what do they look like? With transformations?? All exported and plotted in the "HL-2 Data and timeseries" keynote presentation ####
setwd("~/Desktop")
cat <- c("CentricDiatom", "PennateDiatom", "HABDiatom", "HABFlagellate", "ThecateNonTheca", "Silicoflagellates", "FlagellatesChryptophores", "Ciliates", "Euglena", "TerresFragments")

for(ic in c(1:10)){
  png(paste(cat[ic],"_notrans.png", sep = ","))
  func_group <- abun_df[,ic]   
  hist(func_group)
  dev.off()
  png(paste(cat[ic],"_log10.png", sep = ","))
  hist(log10(func_group))
  dev.off()
  png(paste(cat[ic],"_log.png", sep = ","))
  hist(log(func_group))
  dev.off()
  png(paste(cat[ic],"_sqrt.png", sep = ","))
  hist(sqrt(func_group))
  dev.off()
}



library(imputeTS)
### visualize NAs
plotNA.distribution(log10(abun_df$`Thecate/Non-Theca`), colPoints = "white",
                    colBackgroundMV = "indianred2", main = "Distribution of NAs",
                    xlab = "Time", ylab = "Value", pch = 20, cexPoints = 0.8,
                    col = "white")

### run stats on missing data
statsNA(log10(abun_df$`Thecate/Non-Theca`)) # shows the percentage of NA, the number of NAs in a row, etc

#### loop through if you want to write out results to table, or plot for each functional group ####
ic = 2
for (ic in c(1:10)){
  
  x <- abun_df[,ic]
  x <- x + 1
  x <- log10(x)
  
  ##### 1 #####
  ##### Linear Interpolation and LOESS + smoothing #####
  ## if want to use the log10 of the data frame, need to add 1 to all values to remove zeros
  #x <- abun_df$`Centric Diatom`
  #x <- x + 1
  #x = log10(x)
  linear_interp <- na.interpolation((x), option = "linear") #for linear interpolation using 'approx'
  ## LOESS for smoothing - span = 10
  #len <- 1:length(abun_df$`Centric Diatom`)
  #loessMod10 <- loess(linear_interp ~ len, span=0.10) # 10% smoothing span
  # predict line
  #smoothed10 <- predict(loessMod10) 
  ## LOESS for smoothing - smoothing span = 25
  #loessMod25 <- loess(linear_interp ~ len, span=0.25) # 25% smoothing span
  # predict line
  #smoothed25 <- predict(loessMod25) 
  ## LOESS for smoothing - smoothing span = .50
  #loessMod50 <- loess(linear_interp ~ len, span=0.50) # 50% smoothing span, every 4 years
  # predict line
  #smoothed50 <- predict(loessMod50) 
  
  
  # plot the value (log10(x)) with the smoothing LOESS results
  #plot((x), x=abun_df$Year_month, type="l", main="Loess Smoothing and Prediction 10, 25, 50 span- linear interp", xlab="Year/Month", ylab="Mean cell abundance `HAB Diatom`")
  #lines(smoothed10, x=1:length(smoothed10), col="red")
  #lines(smoothed25, x=1:length(smoothed25), col="blue")
  #lines(smoothed50, x=1:length(smoothed50), col="green")
  
  
  #par(mfrow= c(2,1))
  #plot(loessMod10)
  #plot(smoothed10, col = "red")
  #par(mfrow= c(2,1))
  #plot(loessMod25)
  #plot(smoothed25, col = "green")
  #par(mfrow= c(2,1))
  #plot(loessMod50)
  #plot(smoothed50, col = "blue")
  
  ##### Spline interpolation and LOESS + smoothing ####
  #x <- abun_df$`Centric Diatom`
  #x <- x + 1
  spline_interp <- na.interpolation(x, option = "spline") # for spling interpolation using 'spline'
  
  ## LOESS for smoothing - span = 10
  #len <- 1:length(abun_df$`Centric Diatom`)
  #loessMod10 <- loess(spline_interp ~ len, span=0.10) # 10% smoothing span
  # predict line
  #smoothed10 <- predict(loessMod10) 
  ## LOESS for smoothing - smoothing span = 25
  #loessMod25 <- loess(spline_interp ~ len, span=0.25) # 10% smoothing span
  # predict line
  #smoothed25 <- predict(loessMod25) 
  ## LOESS for smoothing - smoothing span = .50
  #loessMod50 <- loess(spline_interp ~ len, span=0.50) # 10% smoothing span
  # predict line
  #smoothed50 <- predict(loessMod50) 
  
  # plot the value (log10(x)) with the smoothing LOESS results
  #plot(log10(x), x=abun_df$Year_month, type="l", main="Loess Smoothing and Prediction 10, 25, 50 span- spline interp", xlab="Year/Month", ylab="Mean cell abundance")
  #lines(smoothed10, x=1:length(smoothed10), col="red")
  #lines(smoothed25, x=1:length(smoothed25), col="blue")
  #lines(smoothed50, x=1:length(smoothed50), col="green")
  ##### Stine interpolation and LOESS + smoothing ####
  #x <- abun_df$`Centric Diatom`
  #x <- x + 1
  #x = log10(x)
  stine_interp <- na.interpolation(x, option = "stine") # for spling interpolation using 'stine'
  
  ## LOESS for smoothing - span = 10
  #len <- 1:length(abun_df$`Centric Diatom`)
  #loessMod10 <- loess(stine_interp ~ len, span=0.10) # 10% smoothing span
  # predict line
  #smoothed10 <- predict(loessMod10) 
  ## LOESS for smoothing - smoothing span = 25
  #loessMod25 <- loess(stine_interp ~ len, span=0.25) # 10% smoothing span
  # predict line
  #smoothed25 <- predict(loessMod25) 
  ## LOESS for smoothing - smoothing span = .50
  #loessMod50 <- loess(stine_interp ~ len, span=0.50) # 10% smoothing span
  # predict line
  #smoothed50 <- predict(loessMod50) 
  
  # plot the value (log10(x)) with the smoothing LOESS results
  #plot(log10(x), x=abun_df$Year_month, type="l", main="Loess Smoothing and Prediction 10, 25, 50 span- Stineman interp", xlab="Year/Month", ylab="Mean cell abundance")
  #lines(smoothed10, x=1:length(smoothed10), col="red")
  #lines(smoothed25, x=1:length(smoothed25), col="blue")
  #lines(smoothed50, x=1:length(smoothed50), col="green")
  
  #par(mfrow=c(4,1))
  #plotNA.distribution(x, colPoints = "steelblue", colBackgroundMV = "indianred2", main = paste(cat[ic]," Distribution of NAs filled with Linear Interp", sep = ","),xlab = "Time", ylab = "Value", pch = 20, cexPoints = 0.8,col = "black")
  #plot(linear_interp, col = "red")
  #plot(spline_interp, col = "blue")
  #plot(stine_interp, col = "green")
  
}



#### Write out the table of the missing values to compare #####
missingTempDataIndices=which(is.na(abun_df[,ic]))
date_na <- abun_df$Year_month[missingTempDataIndices]
lin_filled <- linear_interp[missingTempDataIndices]
spline_filled <- spline_interp[missingTempDataIndices]
stine_filled <- stine_interp[missingTempDataIndices]

interp_df <- cbind(lin_filled, spline_filled, stine_filled)

setwd("~/phytoplanktontaxonomy/spreadsheets/Interpolation Data")
write.table(interp_df, paste("InterpolatedVals_multiMeth_HL2_FuncGroup.csv", sep = ""), append = TRUE, sep = ",", col.names = TRUE, row.names = date_na) 

}


###### 2 Moving Weighted Moving Average #####
#This means for an NA value at position i of a time series, the observations i-1,i+1 and i+1, i+2 (assuming a window size of k=2) are used to calculate the mean.
# I should not use a window greater than k=2 as I am already grouping the cell counts by month, very broad in terms of phytoplankton temporal trends. 
#Whenever there are less than 2 non-NA values in the complete window available, the window size is incrementally increased, till at least 2 non-NA values are there.

##### Loop if running a table out or multiple images to export ####
ic =1
for (ic in c(2:10)) {
  ##### Simple weighting: SMA: all observations in the window are equally weighted for calculating the mean. ####
  x <- abun_df[,ic]
  x <- x + 1
  x <- log10(x)
  
  # k = 2
  k2_s <- na.ma(x, k = 2, weighting = "simple")
  # k = 4
  #k4_s <- na.ma(x, k = 4, weighting = "simple")
  # k = 6
  #k6_s <- na.ma(x, k = 6, weighting = "simple")
  # k = 8
  #k8_s <- na.ma(x, k = 10, weighting = "simple")
  # k = 10
  #k10_s <- na.ma(x, k = 10, weighting = "simple")
  
  
  #smoothed2 <- predict(k2) 
  #smoothed4 <- predict(k4) 
  #smoothed6 <- predict(k6) 
  #smoothed8 <- predict(k8) 
  #smoothed10 <- predict(k10) 
  
  ##### Linear weighting: LWMA: weights decrease in arithmetical progression. The observations directly next to a central value i, have weight 1/2, the observations one further away (i-2,i+2) have weight 1/3, the next (i-3,i+3) have weight 1/4, ... #####
  x <- abun_df[,ic]
  x <- x + 1
  x <- log10(x)
  
  # k = 2
  k2_l <- na.ma(x, k = 2, weighting = "linear")
  # k = 4
  #k4_l <- na.ma(x, k = 4, weighting = "linear")
  # k = 6
  #k6_l <- na.ma(x, k = 6, weighting = "linear")
  # k = 8
  #k8_l <- na.ma(x, k = 8, weighting = "linear")
  # k = 10
  #k10_l <- na.ma(x, k = 10, weighting = "linear")
  
  
  #smoothed2 <- predict(k2) 
  #smoothed4 <- predict(k4) 
  #smoothed6 <- predict(k6) 
  #smoothed8 <- predict(k8) 
  #smoothed10 <- predict(k10) 
  
  ##### Exponential weighting: EWMA: uses weighting factors which decrease exponentially. The observations directly next to a central value i, have weight 1/2^1, the observations one further away (i-2,i+2) have weight 1/2^2, the next (i-3,i+3) have weight 1/2^3, ...---> USING 26/09/2017 RM ####
  x <- abun_df[,ic]
  x <- x + 1
  x <- log10(x)
  
  # k = 2
  k2_e <- na.ma(x, k = 2, weighting = "exponential")
  # k = 4
  #k4_e <- na.ma(x, k = 4, weighting = "exponential")
  # k = 6
  #k6_e <- na.ma(x, k = 6, weighting = "exponential")
  # k = 8
  #k8_e <- na.ma(x, k = 8, weighting = "exponential")
  # k = 10
  #k10_e <- na.ma(x, k = 10, weighting = "exponential")
  
  #smoothed2 <- predict(k2) 
  #smoothed4 <- predict(k4) 
  #smoothed6 <- predict(k6) 
  #smoothed8 <- predict(k8) 
  #smoothed10 <- predict(k10) 
  ###### plot and table out the MWM ####
  missingTempDataIndices=which(is.na(abun_df[,ic]))
  date_na <- abun_df$Year_month[missingTempDataIndices]
  k2s_na <- k2_s[missingTempDataIndices]
  k2l_na <- k2_l[missingTempDataIndices]
  k2e_na <- k2_e[missingTempDataIndices]
  
  interp_df <- cbind(k2s_na, k2l_na, k2e_na)
  
  setwd("~/phytoplanktontaxonomy/spreadsheets/Interpolation Data")
  write.table(interp_df, paste("MWMk2_multimod_HLFuncGroup_.csv", sep = ""), append = TRUE, sep = ",", col.names = TRUE, row.names = date_na) 
  
  
  
  par(mfrow=c(4,1))
  plotNA.distribution(x, colPoints = "steelblue", colBackgroundMV = "indianred2", main = paste(cat[ic]," Distribution of NAs filled with MWA: Simple, Linear, Expon. weight, K=2", sep = ","),xlab = "Time", ylab = "Value", pch = 20, cexPoints = 0.8,col = "black")
  plot(k2_s, col = "red")
  plot(k2s_na, add = TRUE)
  plot(k2_l, col = "blue")
  plot(k2_e, col = "green")
  
}

smoothed2 <- predict(k2e_na) 



#### De-seasoning time series, in R --> come back to TOMORRROW (27/09/2017) ####


#create time series object
## Want to minimize the residuals = transformations of each phytoplankton functional group at HL-2:

###############################
#LOESS --> input needed is imputed Time series (aka missing values filled in)
###############################



for (ic in c(1:10)) {
  #create the time series object
  TS <- ts((k10), start = c(1999, 1), frequency = 12) # can transform the data within the equation
  
  # Decomposition using LOESS method, stl() package in R
  stl_ts <- stl(TS, s.window = "periodic", t.window = 12, robust = TRUE) # robust = TRUE deals with strong outliers, reduces residual noise, might be good for the TS where there are values of 1 that obstruct the pattern
  
  setwd("~/Desktop")
  png(paste("TS_LOESS_weightedmeank10impute_",cat[ic],"_HL2.png",sep =","), width = 650)
  
  plot(stl_ts, col = c(1:4), lwd = 2, main = paste("Decomposed Time series ",cat[ic]), cex.axis = 20, cex.main = 20)
  
  dev.off()
  
}
# ts object of deseasoned trend
seasonal   <- stl_ts$time.series[,1]
trend	   <- stl_ts$time.series[,2]
residual  <- stl_ts$time.series[,3]
# check the distributuon of the residual
par(mfrow=c(1,2))
hist(residual)
qqnorm(residual)
sh <- shapiro.test(residual) 
# get the normality test object for printing
test <- sh$method
W <- format(sh$statistic, digits = 3)
P.val <- format(sh$p.value, digits = 3)
text(0,0.4, test)
text(0, 0.35, paste("W = ",W, sep = ""))
text(0, 0.3, paste("P-value = ", P.val, sep = ""))


## check the distribution of the residuals + overal trend for the deseasoned time series trend
hist(trend + residual)
qqnorm(trend + residual)
sh <- shapiro.test(trend + residual) 

# get the normality test object for printing
test <- sh$method
W <- format(sh$statistic, digits = 3)
P.val <- format(sh$p.value, digits = 3)
text(0, 0.2, test)
text(0, 0.15, paste("W = ",W, sep = ""))
text(0, 0.1, paste("P-value = ", P.val, sep = ""))


## Create the de-seasoned trend and plot 
TS_ds <- ts(c(trend + residual), start = c(1999, 1), frequency = 12)
fit <- tslm((TS_ds ~ abun_df$Year_month)) # y ~ x 
summary(fit)
# make the plot
png("deseasoned_DynamicSSEBSA_AnnComp_log10CHL_loess.png", width = 650)

par(mfrow = c(1,1))
plot(trend+residual, lwd = 2, col = "darkblue", 
     main="CHL over Time, Seasonally Adjusted",
     ylab = "Chla (log10 ug/L)",
     xlab = "Year")

abline(fit, col = "darkred", lwd = 2)

# get the r2 and the p value of the line to paste on the figure
modsum <- summary(fit)
r2 = modsum$adj.r.squared
my.p =modsum$coefficients[2,4]
#get rounded intercept and slope estimates from the linear model - round to 2 significant digits
cf <- round(coef(fit), 2) 
eq <- paste("Y =",cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x")

## create legend of all the line info to paste
rp = vector('expression',3)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=2)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
rp[3] = eq
#rp

# add legend, top right 
legend("bottom", legend = rp, bty = 'n', cex = 1.2, x.intersp = 1, y.intersp = 1.2)

setwd("~/EBSAs/Time_series_andAnomalies/decomposed and deseasoned TS Annual Composite Dynamic SS EBSA")
dev.off()


oddZero<- (abun_df[c(45,50),])
dim(oddZero)
oddZero


##### correlation table for each of the phytoplankton groups - using log10 and non transformed monthly mean cell abundances ####

abun_df$year <- (substr(as.character(abun_df$Year_month),1,4))
abun_df$month <- substr(as.character(abun_df$Year_month), 5,6)

cormat <- as.matrix(abun_df[,1:10])

# get average of each per month
##### Cut data frame into autotrophs only ####
rawdat_cut <- filter(rawdat, CATEGORY <11) # exclude categories 9 and 10 

##### Specify year and month character vectors  ####
#month <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
year <- c("1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016")
cater = c("Centric Diatoms", "Pennate Diatoms", "HAB Diatoms", "HAB Flagellates", "Thecate Non-Thecates", "Silicates", "Flagellates", "Cilliates", "Euglenoids", "Terrestrial")  #### 

#### LINE GRAPHS --> Get general monthly stats (# samplding days, sum cells, mean cells/cast) ALL cells of all groups for each year --> for just catergories 1-8 (autotrophs) ####



#### THis is a test - what happens when I need to update the file on Git Hub #####


