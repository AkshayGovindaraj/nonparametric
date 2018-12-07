# install.packages('locpol','splines','sm')
# setwd('/Users/visswanathvenkataraman/Desktop/ISYE_6404_Exam3')
library(locpol) 
library(splines)
data<- read.csv("motor.csv",header=TRUE)

y <- data$accel
x1 <- data$t1
x2 <- data$t2
x3 <- data$t3
nrows <- length(y)
plot(x1,y)
plot(x2,y)
plot(x3,y)

# Question 3
# Part 1
# LOOCV Smoothing Spline
Q3_P1_sspline <- function(variable){
  ss_res_sspline <- c(0)
  for(i in 1:nrows){
    y_pred <- predict(smooth.spline(x=variable[-i],y=y[-i],spar=0.6),variable[i])
    ss_res_sspline<- ss_res_sspline + (y[i]-y_pred$y)^2
  }
  ss_res_sspline <- ss_res_sspline/nrows
  return(ss_res_sspline)
}

Q3_P1_sspline(x1)
Q3_P1_sspline(x2)
Q3_P1_sspline(x3)

sspline_fit_x1 <- smooth.spline(x1,y)
fit.ss_x1 <- fitted(sspline_fit_x1)

sspline_fit_x2 <- smooth.spline(x2,y)
fit.ss_x2 <- fitted(sspline_fit_x2)

sspline_fit_x3 <- smooth.spline(x3,y)
fit.ss_x3 <- fitted(sspline_fit_x3)



# Bandwidth Selection for Kernel Regression (using different x-variables)
bandwidth_selection <- function(x_variable){
  n = length(x_variable)
  # n: sample size
  h_seq = seq(from=1,to=5, by=0.1)
  # smoothing bandwidths we are using
  CV_err_h = rep(NA,length(h_seq))
  for(j in 1:length(h_seq)){
    h_using = h_seq[j]
    CV_err = rep(NA, n)
    for(i in 1:n){
      X_val = x_variable[i]
      Y_val = y[i]
      # validation set
      X_tr = x_variable[-i]
      Y_tr = y[-i]
      # training set
      Y_val_predict = ksmooth(x=X_tr,y=Y_tr,kernel = "normal",bandwidth=h_using,x.points = X_val)
      CV_err[i] = (Y_val - Y_val_predict$y)^2
      # we measure the error in terms of difference square
    }
    CV_err_h[j] = mean(CV_err)
  }
  x=min(CV_err_h,na.rm=TRUE)
  plot(x=h_seq,y=CV_err_h,type="b",lwd=3,col="blue",xlab="Smoothing bandwidth",ylab="LOOCV prediction error")
  return(h_seq[which(CV_err_h==x)])
}

bandwidth_selection(x1) # Kernel Band-width for x1 = 3.7
bandwidth_selection(x2) # Kernel Band-width for x2 = 3.4
bandwidth_selection(x3) # Kernel Band-width for x3 = 3.9


# LOOCV Kernel
Q3_P1_kernel <- function(variable,BW){
  ss_res_kreg <- c(0)
  for(i in 1:nrows){
    y_pred <- ksmooth(x=variable[-i],y=y[-i],bandwidth = BW,kernel="normal",x.points=variable[i])
    ss_res_kreg <- ss_res_kreg + (y[i]-y_pred$y)^2
  }
  ss_res_kreg <- ss_res_kreg/nrows
  return( ss_res_kreg)
  
}

Q3_P1_kernel(x1,3.7)
Q3_P1_kernel(x2,3.4)
Q3_P1_kernel(x3,3.9)

kernel_fit_x1 <- ksmooth(x1,y,kernel="normal",bandwidth=3.7)
kernel_fit_x2 <- ksmooth(x2,y,kernel="normal",bandwidth=3.4)
kernel_fit_x3 <- ksmooth(x3,y,kernel="normal",bandwidth=3.9)

# Comparing Fit Values
compare_x1_both_methods <- data.frame(Q3_P1_sspline(x1),Q3_P1_kernel(x1,68))
compare1 <- data.frame(LOOCV_Value_x1=t(compare_x1_both_methods))

compare_x2_both_methods <- data.frame(Q3_P1_sspline(x2),Q3_P1_kernel(x2,34))
compare2 <- data.frame(LOOCV_Value_x2=t(compare_x2_both_methods))

compare_x3_both_methods <- data.frame(Q3_P1_sspline(x3),Q3_P1_kernel(x3,2975))
compare3 <- data.frame(LOOCV_Value_x3=t(compare_x3_both_methods))

compare <- cbind(compare1,compare2,compare3)

# The LOOCV values are lower for Smoothing Spline Method as compared to the Kernel Regression 
# for all three cases.

par(mfrow=c(3,1))
plot(x1,y,ylim=c(-160,120),main="Spline vs Kernel Fit for variable x1")
lines(kernel_fit_x1,lwd=2,type='l',col="limegreen")
lines(fit.ss_x1,lwd=2,type='l',col="red")
legend(x=2,y=100,legend=c("Kernel Fit - x1","Spline Fit - x1"),lty=1,col=c("limegreen","red"))

plot(x2,y,ylim=c(-160,120),main="Spline vs Kernel Fit for variable x2")
lines(kernel_fit_x1,lwd=2,type='l',col="blue")
lines(fit.ss_x2,lwd=2,type='l',col="red")
legend(x=2,y=100,legend=c("Kernel Fit - x2","Spline Fit - x2"),lty=1,col=c("blue","red"))


plot(x2,y,ylim=c(-160,120),main="Spline vs Kernel Fit for variable x3")
lines(kernel_fit_x1,lwd=2,type='l',col="purple")
lines(fit.ss_x2,lwd=2,type='l',col="yellow")
legend(x=2,y=100,legend=c("Kernel Fit - x3","Spline Fit - x2"),lty=1,col=c("purple","yellow"))



# Part 2
# Data from Middle rows
x11 = 85
x12 = 52
x13 = 2035
x_middle = cbind(x1,x2,x3)

# Data from the Edge
x21 = 119
x22 = 82
x23 = 2720
x_edge = cbind(x21,x22,x23)

# Prediction of Y at x_middle data points using the Spline Regression