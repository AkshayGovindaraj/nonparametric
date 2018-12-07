install.packages("wavelets")
install.packages("wmtsa")
install.packages("waveslim")
install.packages("wavethresh")
install.packages("Metrics")
install.packages("tseries")

library(wavelets)
library(wmtsa)
library(waveslim)
library(wavethresh)
library(Metrics)
library(tseries)

# Part 1 & Part 2
# Plot the AAPL stock price data time series for 1024 days
AMZN_data <- read.csv('AMZN.csv')
AMZN_data_matrix <-as.matrix(AMZN_data)
plot.ts(AMZN_data$Close,pch=10,cex=0.2, xlab="Days", ylab="Stock price",main="AMZN Stock Price Close")

# Haar transform
AMZN_haar <- mra(AMZN_data$Close, "haar", 9, "dwt")
names(AMZN_haar) <- c("d1", "d2", "d3", "d4","d5","d6","d7","d8","d9", "s")
mar = c(0,0,0,0)
par(mfrow=c(2,5))
# plot.ts(AMZN_data$Close,pch=10,cex=0.2, xlab="Days", ylab="Stock price",main="Haar Transform Coefficients")
for(i in 1:10)
  plot.ts(AMZN_haar[[i]],  ylab=names(AMZN_haar)[i])


# Daubechies 8 transform
AMZN_la8 <- mra(AMZN_data$Close, "la8", 9, "dwt")
names(AMZN_la8) <- c("d1", "d2", "d3", "d4","d5","d6","d7","d8","d9","s")
mar = c(0,0,0,0)
par(mfrow=c(2,5))
# plot.ts(AMZN_data$Close,pch=10,cex=0.2, xlab="Days", ylab="Stock price",main="La8 Transform Coefficients")
for(i in 1:10)
  plot.ts(AMZN_la8[[i]], axes=TRUE, ylab=names(AMZN_la8)[i])


# Part - 3
# Thresholding and inverse with Daubachies

AMZN.daub.thres = lapply(c("universal","adaptive"),function(k,x) wavShrink(x,wavelet = "s8",thresh.fun = k, xform = "dwt"),x =AMZN_data$Close)
par(mfrow=c(2,1))
ifultools::stackPlot(x=AMZN_data$Days,y=data.frame(AMZN_data$Close,AMZN.daub.thres), ylab=c("true price", "universal","adaptive"), xlab="Days")
plot(AMZN_data$Days,AMZN_data$Close,type="l",col="black")
lines(AMZN_data$Days,matrix(unlist(AMZN.daub.thres[1])),col="blue")
lines(AMZN_data$Days,matrix(unlist(AMZN.daub.thres[2])),col="red")

# Thresholding and inverse with Haar
AMZN.haar.thres = lapply(c("universal","adaptive"),function(k,x) wavShrink(x,wavelet = "haar",thresh.fun = k, xform = "dwt"),x =AMZN_data$Close)
ifultools::stackPlot(x=AMZN_data$Day,y=data.frame(AMZN_data$Close,AMZN.haar.thres), ylab=c("true price", "universal","adaptive"), xlab="Days")
plot(AMZN_data$Days,AMZN_data$Close,type="l",col="black")
lines(AMZN_data$Days,matrix(unlist(AMZN.haar.thres[1])),col="blue")
lines(AMZN_data$Days,matrix(unlist(AMZN.haar.thres[2])),col="red")

Error_haar_universal<- mse(AMZN_data$Close,AMZN.haar.thres[[1]]) # MSE of haar wavelet using universal thresholding = 1.870186
Error_haar_adaptive<- mse(AMZN_data$Close,AMZN.haar.thres[[2]]) # MSE of haar wavelet using adaptive thresholding = 357.904
Error_daub_universal<- mse(AMZN_data$Close,AMZN.daub.thres[[1]]) # MSE of Daubachies wavelet using universal threshold = 1.074561
Error_daub_adaptive<- mse(AMZN_data$Close,AMZN.daub.thres[[2]]) # MSE of Daubachies wavelet using adaptive threshold =  55.86791

Error_haar_universal <- data.frame(Error_haar_universal)
Error_haar_adaptive <- data.frame(Error_haar_adaptive)
Error_daub_universal <- data.frame(Error_daub_universal)
Error_daub_adaptive <- data.frame(Error_daub_adaptive)
Error_comparison <- cbind(Error_haar_universal,Error_haar_adaptive,Error_daub_universal,Error_daub_adaptive)
Error_comparison

# Calculating the Haar Data reduction ratio Universal
AMZN_haar_dwt <- dwt(AMZN_data$Close, wf='haar', n.levels=9, boundary = 'periodic')
AMZN_haar_thres_dwt <- dwt(as.vector(AMZN.haar.thres[[1]]), wf='haar', n.levels=9, boundary = 'periodic')

zero_haar <- sum(AMZN_haar_dwt$d1==0) + sum(AMZN_haar_dwt$d2==0) + sum(AMZN_haar_dwt$d3==0) + sum(AMZN_haar_dwt$d4==0) + sum(AMZN_haar_dwt$d5==0)
zero_haar <- zero_haar + sum(AMZN_haar_dwt$d6==0) + sum(AMZN_haar_dwt$d7==0) + sum(AMZN_haar_dwt$d8==0) + sum(AMZN_haar_dwt$d9==0) + sum(AMZN_haar_dwt$s9==0)
zero_haar_thres <- sum(AMZN_haar_thres_dwt$d1==0) + sum(AMZN_haar_thres_dwt$d2==0) + sum(AMZN_haar_thres_dwt$d3==0) + sum(AMZN_haar_thres_dwt$d4==0) + sum(AMZN_haar_thres_dwt$d5==0)
zero_haar_thres <- zero_haar_thres + sum(AMZN_haar_thres_dwt$d6==0) + sum(AMZN_haar_thres_dwt$d7==0) + sum(AMZN_haar_thres_dwt$d8==0) + sum(AMZN_haar_thres_dwt$d9==0) + sum(AMZN_haar_thres_dwt$s9==0)

DRR_haar_univ = (zero_haar_thres - zero_haar)*100/1024

# Daub Data Reduction Ratio Universal
AMZN_daub_dwt <- dwt(AMZN_data$Close, wf='la8', n.levels=9, boundary = 'periodic')
AMZN_daub_thres_dwt <- dwt(as.vector(AMZN.daub.thres[[1]]), wf='la8', n.levels=9, boundary = 'periodic')

zero_daub <- sum(AMZN_daub_dwt$d1==0) + sum(AMZN_daub_dwt$d2==0) + sum(AMZN_daub_dwt$d3==0) + sum(AMZN_daub_dwt$d4==0) + sum(AMZN_daub_dwt$d5==0)
zero_daub <- zero_daub + sum(AMZN_daub_dwt$d6==0) + sum(AMZN_daub_dwt$d7==0) + sum(AMZN_daub_dwt$d8==0) + sum(AMZN_daub_dwt$d9==0) + sum(AMZN_daub_dwt$s9==0)
zero_daub_thres <- sum(AMZN_daub_thres_dwt$d1==0) + sum(AMZN_daub_thres_dwt$d2==0) + sum(AMZN_daub_thres_dwt$d3==0) + sum(AMZN_daub_thres_dwt$d4==0) + sum(AMZN_daub_thres_dwt$d5==0)
zero_daub_thres <- zero_daub_thres + sum(AMZN_daub_thres_dwt$d6==0) + sum(AMZN_daub_thres_dwt$d7==0) + sum(AMZN_daub_thres_dwt$d8==0) + sum(AMZN_daub_thres_dwt$d9==0) + sum(AMZN_daub_thres_dwt$s9==0)

DRR_daub_univ = (zero_daub_thres - zero_daub)*100/1024

# Calculating the Haar Data reduction ratio Adaptive
AMZN_haar_dwt <- dwt(AMZN_data$Close, wf='haar', n.levels=9, boundary = 'periodic')
AMZN_haar_thres_dwt <- dwt(as.vector(AMZN.haar.thres[[2]]), wf='haar', n.levels=9, boundary = 'periodic')

zero_haar <- sum(AMZN_haar_dwt$d1==0) + sum(AMZN_haar_dwt$d2==0) + sum(AMZN_haar_dwt$d3==0) + sum(AMZN_haar_dwt$d4==0) + sum(AMZN_haar_dwt$d5==0)
zero_haar <- zero_haar + sum(AMZN_haar_dwt$d6==0) + sum(AMZN_haar_dwt$d7==0) + sum(AMZN_haar_dwt$d8==0) + sum(AMZN_haar_dwt$d9==0) + sum(AMZN_haar_dwt$s9==0)
zero_haar_thres <- sum(AMZN_haar_thres_dwt$d1==0) + sum(AMZN_haar_thres_dwt$d2==0) + sum(AMZN_haar_thres_dwt$d3==0) + sum(AMZN_haar_thres_dwt$d4==0) + sum(AMZN_haar_thres_dwt$d5==0)
zero_haar_thres <- zero_haar_thres + sum(AMZN_haar_thres_dwt$d6==0) + sum(AMZN_haar_thres_dwt$d7==0) + sum(AMZN_haar_thres_dwt$d8==0) + sum(AMZN_haar_thres_dwt$d9==0) + sum(AMZN_haar_thres_dwt$s9==0)

DRR_haar_adap = (zero_haar_thres - zero_haar)*100/1024

# Daub Data Reduction Ratio Adaptive
AMZN_daub_dwt <- dwt(AMZN_data$Close, wf='la8', n.levels=9, boundary = 'periodic')
AMZN_daub_thres_dwt <- dwt(as.vector(AMZN.daub.thres[[2]]), wf='la8', n.levels=9, boundary = 'periodic')

zero_daub <- sum(AMZN_daub_dwt$d1==0) + sum(AMZN_daub_dwt$d2==0) + sum(AMZN_daub_dwt$d3==0) + sum(AMZN_daub_dwt$d4==0) + sum(AMZN_daub_dwt$d5==0)
zero_daub <- zero_daub + sum(AMZN_daub_dwt$d6==0) + sum(AMZN_daub_dwt$d7==0) + sum(AMZN_daub_dwt$d8==0) + sum(AMZN_daub_dwt$d9==0) + sum(AMZN_daub_dwt$s9==0)
zero_daub_thres <- sum(AMZN_daub_thres_dwt$d1==0) + sum(AMZN_daub_thres_dwt$d2==0) + sum(AMZN_daub_thres_dwt$d3==0) + sum(AMZN_daub_thres_dwt$d4==0) + sum(AMZN_daub_thres_dwt$d5==0)
zero_daub_thres <- zero_daub_thres + sum(AMZN_daub_thres_dwt$d6==0) + sum(AMZN_daub_thres_dwt$d7==0) + sum(AMZN_daub_thres_dwt$d8==0) + sum(AMZN_daub_thres_dwt$d9==0) + sum(AMZN_daub_thres_dwt$s9==0)

DRR_daub_adap = (zero_daub_thres - zero_daub)*100/1024

## Part 4 : Local and Global Changes

AMZN_data_c1 <- AMZN_data
AMZN_data_c2 <- AMZN_data
AMZN_data_c3 <- AMZN_data

# Local change

AMZN_data_c1$Close[16] = AMZN_data_c1$Close[16] + 100
AMZN_data_c1$Close[17] = AMZN_data_c1$Close[17] - 100
AMZN_data_c1$Close[300] = AMZN_data_c1$Close[300] + 100
AMZN_data_c1$Close[301] = AMZN_data_c1$Close[301] - 100
AMZN_data_c1$Close[485] = AMZN_data_c1$Close[485] + 100
AMZN_data_c1$Close[486] = AMZN_data_c1$Close[486] - 100
AMZN_data_c1$Close[813] = AMZN_data_c1$Close[813] + 100
AMZN_data_c1$Close[814] = AMZN_data_c1$Close[814] - 100

# Plotting Haar MRA
AMZN_haar_c1 <- mra(AMZN_data_c1$Close, "haar", 9, "dwt")
names(AMZN_haar_c1) <- c("d1", "d2", "d3", "d4","d5","d6","d7","d8","d9", "s")
mar = c(0,0,0,0)
par(mfrow=c(2,5))
# plot.ts(AMZN_data_c1$Close,pch=10,cex=0.2, xlab="Days", ylab="Stock price",main="Local - Haar MRA")
for(i in 1:10)
  plot.ts(AMZN_haar_c1[[i]], axes=TRUE, ylab=names(AMZN_haar_c1)[i])

# Plotting La8 MRA
AMZN_daub_c1 <- mra(AMZN_data_c1$Close, "la8", 9, "dwt")
names(AMZN_daub_c1) <- c("d1", "d2", "d3", "d4","d5","d6","d7","d8","d9", "s")
mar = c(0,0,0,0)
par(mfrow=c(2,5))
# plot.ts(AMZN_data_c1$Close,pch=10,cex=0.2, xlab="Days", ylab="Stock price",main="Local - La8 MRA")
for(i in 1:10)
  plot.ts(AMZN_daub_c1[[i]], axes=TRUE, ylab=names(AMZN_daub_c1)[i])

# Global Change

for(i in 501:650)
  AMZN_data_c2$Close[i] = AMZN_data_c2$Close[i] + 50
  
# Plotting Haar MRA
AMZN_haar_c2 <- mra(AMZN_data_c2$Close, "haar", 9, "dwt")
names(AMZN_haar_c2) <- c("d1", "d2", "d3", "d4","d5","d6","d7","d8","d9", "s")
mar = c(0,0,0,0)
par(mfrow=c(2,5))
# plot.ts(AMZN_data_c2$Close,pch=10,cex=0.2, xlab="Days", ylab="Stock price",main="Global- Haar MRA")
for(i in 1:10)
  plot.ts(AMZN_haar_c2[[i]], axes=TRUE, ylab=names(AMZN_haar_c2)[i])

# Plotting La8 MRA
AMZN_daub_c2 <- mra(AMZN_data_c2$Close, "la8", 9, "dwt")
names(AMZN_daub_c2) <- c("d1", "d2", "d3", "d4","d5","d6","d7","d8","d9", "s")
mar = c(0,0,0,0)
par(mfrow=c(2,5))
# plot.ts(AMZN_data_c2$Close,pch=10,cex=0.2, xlab="Days", ylab="Stock price",main="Global- La8 MRA")
for(i in 1:10)
  plot.ts(AMZN_daub_c2[[i]], axes=TRUE, ylab=names(AMZN_daub_c2)[i])
  
# Local and Global changes

AMZN_data_c3$Close[16] = AMZN_data_c3$Close[16] + 100
AMZN_data_c3$Close[17] = AMZN_data_c3$Close[17] - 100
AMZN_data_c3$Close[300] = AMZN_data_c3$Close[300] + 100
AMZN_data_c3$Close[301] = AMZN_data_c3$Close[301] - 100
AMZN_data_c3$Close[485] = AMZN_data_c3$Close[485] + 100
AMZN_data_c3$Close[486] = AMZN_data_c3$Close[486] - 100
AMZN_data_c3$Close[813] = AMZN_data_c3$Close[813] + 100
AMZN_data_c3$Close[814] = AMZN_data_c3$Close[814] - 100

for(i in 501:650)
  AMZN_data_c3$Close[i] = AMZN_data_c3$Close[i] + 50

# Plotting Haar MRA
AMZN_haar_c3 <- mra(AMZN_data_c3$Close, "haar", 9, "dwt")
names(AMZN_haar_c3) <- c("d1", "d2", "d3", "d4","d5","d6","d7","d8","d9", "s")
mar = c(0,0,0,0)
par(mfrow=c(2,5))
# plot.ts(AMZN_data_c3$Close,pch=10,cex=0.2, xlab="Days", ylab="Stock price",main="Local & Global- Haar MRA")
for(i in 1:10)
  plot.ts(AMZN_haar_c3[[i]], axes=TRUE, ylab=names(AMZN_haar_c3)[i])

# Plotting La8 MRA
AMZN_daub_c3 <- mra(AMZN_data_c3$Close, "la8", 9, "dwt")
names(AMZN_daub_c3) <- c("d1", "d2", "d3", "d4","d5","d6","d7","d8","d9", "s")
mar = c(0,0,0,0)
par(mfrow=c(2,5))
# plot.ts(AMZN_data_c3$Close,pch=10,cex=0.2, xlab="Days", ylab="Stock price",main="Local & Global- La8 MRA")
for(i in 1:10)
  plot.ts(AMZN_daub_c3[[i]], axes=TRUE, ylab=names(AMZN_daub_c3)[i])
