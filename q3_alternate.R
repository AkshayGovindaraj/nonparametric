#### Spline and Kernel 

library(locpol) 
library(splines)
library(boot)
library(gam)
library(mgcv)
data<- read.csv("motor.csv",header=TRUE)

y <- data$accel
x1 <- data$t1
x2 <- data$t2
x3 <- data$t3
nrows <- length(y)
plot(x1,y)
plot(x2,y)
plot(x3,y)

## Spline regression

sspline_fit_x1 <- smooth.spline(x1,y)
fit.ss_x1 <- fitted(sspline_fit_x1)

plot(y-fit.ss_x1, fit.ss_x1)

sspline_fit_x2 <- smooth.spline(x2,y-fit.ss_x1)
fit.ss_x2 <- fitted(sspline_fit_x2)

plot(fit.ss_x1 + fit.ss_x2)      # Combined predicions

sspline_fit_x3 <- smooth.spline(x3,y-fit.ss_x1-fit.ss_x2)
fit.ss_x3 <- fitted(sspline_fit_x3)

plot(fit.ss_x1 + fit.ss_x2 + fit.ss_x3)      # Combined predicions

plot(fit.ss_x1 + fit.ss_x2 + fit.ss_x3 - y)      # Final Residuals

# LOOCV Spline

Q3_P1_sspline <- function(variable){
  ss_res_sspline <- c(0)
  x_1 = variable$t1
  x_2 = variable$t2
  x_3 = variable$t3
  y_1 = variable$accel
  for(i in 1:nrows){
    # y_pred <- predict(smooth.spline(x=x1[-i],y=y[-i],spar=0.6),variable[i])
    x1 = x_1[-i]
    x2 = x_2[-i]
    x3 = x_3[-i]
    y1 = y_1[-i]
    sspline_fit_x1 <- smooth.spline(x1,y1)
    fit.ss_x1 <- fitted(sspline_fit_x1)
    sspline_fit_x2 <- smooth.spline(x2,y1-fit.ss_x1)
    fit.ss_x2 <- fitted(sspline_fit_x2)
    sspline_fit_x3 <- smooth.spline(x3,y1-fit.ss_x1-fit.ss_x2)
    fit.ss_x3 <- fitted(sspline_fit_x3)
    y_pred <- predict(sspline_fit_x1,x1_1)
    ss_res_sspline<- ss_res_sspline + (y_1[i]-y_pred$y)^2
  }
  ss_res_sspline <- ss_res_sspline/nrows
  return(ss_res_sspline)
}

## Kernel regression

kernel_fit_x1 <- ksmooth(x1,y,kernel="normal",bandwidth=3.7)

plot(y-kernel_fit_x1$y)   # Residual plot

kernel_fit_x2 <- ksmooth(x2,y-kernel_fit_x1$y,kernel="normal",bandwidth=3.4)

plot(y-kernel_fit_x1$y-kernel_fit_x2$y)   # Residual plot

kernel_fit_x3 <- ksmooth(x3,y-kernel_fit_x1$y-kernel_fit_x2$y,kernel="normal",bandwidth=3.9)

plot(y-kernel_fit_x1$y-kernel_fit_x2$y-kernel_fit_x3$y)   # Residual plot

## Combined plots

plot(x3,y,ylim=c(-160,120),main="Spline vs Kernel Fit for variable x3")
lines(kernel_fit_x1$y+kernel_fit_x2$y+kernel_fit_x3$y,lwd=2,type='l',col="purple")
lines(fit.ss_x1+fit.ss_x2+fit.ss_x3,lwd=2,type='l',col="yellow")
legend(x=2,y=100,legend=c("Kernel Fit - x3","Spline Fit - x2"),lty=1,col=c("purple","yellow"))

### Part 2

# Choosing first point

x1_1 = 2.6
x2_1 = 1.0
x3_1 = 2.0

# Choosing second point

x1_2 = 24.0
x2_2 = 22.9
x3_2 = 23.5

## Predictions based on these points

y1_pred_spline = predict(sspline_fit_x1,x1_1)$y + predict(sspline_fit_x2,x2_1)$y + predict(sspline_fit_x3,x3_1)$y
y2_pred_spline = predict(sspline_fit_x1,x1_2)$y + predict(sspline_fit_x2,x2_2)$y + predict(sspline_fit_x3,x3_2)$y

kernel_fit <- gam(accel~s(t1,bs='gp'),data=data)
data1 <- data
data1$residuals <- kernel_fit$residuals
kernel_fit2 <- gam(residuals~s(t2,bs='gp'),data=data1)
data1$residuals1 <- kernel_fit2$residuals
kernel_fit3 <- gam(residuals1~s(t3,bs='gp'),data=data1)

x1_middle <- data.frame(t1=c(24.6))#,t2=c(23.5),t3=c(23.2))
x2_middle <- data.frame(t2=c(23.5))
x3_middle <- data.frame(t3=c(23.2))

prediction <- predict.gam(kernel_fit,x1_middle) + predict.gam(kernel_fit2,x2_middle) + predict.gam(kernel_fit3,x3_middle)

# y1_pred_kernel = predict(sspline_fit_x1,x1_1)$y + predict(sspline_fit_x2,x2_1)$y + predict(sspline_fit_x3,x3_1)$y
# y2_pred_kernel = predict(sspline_fit_x1,x1_2)$y + predict(sspline_fit_x2,x2_2)$y + predict(sspline_fit_x3,x3_2)$y

## Part 3 : Bootstrap

bootstrap_fn <- function(d,i,x){
  data_source <- d[i,]
  y <- data_source$accel
  x1 <- data_source$t1
  x2 <- data_source$t2
  x3 <- data_source$t3
  # bootstrap_fit <- gam(accel ~s(t1) +s(t2)+s(t3),data = data_source)
  # prediction <- predict.gam(bootstrap_fit,x)
  bsspline_fit_x1 <- smooth.spline(x1,y)
  fit.ss_x1 <- fitted(bsspline_fit_x1)
  bsspline_fit_x2 <- smooth.spline(x2,y-fit.ss_x1)
  fit.ss_x2 <- fitted(bsspline_fit_x2)
  bsspline_fit_x3 <- smooth.spline(x3,y-fit.ss_x1-fit.ss_x2)
  fit.ss_x3 <- fitted(bsspline_fit_x3)
  prediction <- predict(bsspline_fit_x1,24.6)$y + predict(bsspline_fit_x2,23.5)$y + predict(bsspline_fit_x3,23.2)$y
  return(prediction)
}

B <- 1e04
bs_CI <- boot(data, bootstrap_fn, R=B,parallel = 'multicore')
boot.ci(boot.out = bs_CI,conf = 0.9,type = 'bca')

bootstrap_kernel <- function(d,i,x){
  data_source <- d[i,]
  kernel_fit <- gam(accel~s(t1,bs='gp'),data=data_source)
  data1 <- data_source
  data1$residuals <- kernel_fit$residuals
  kernel_fit2 <- gam(residuals~s(t2,bs='gp'),data=data1)
  data1$residuals1 <- kernel_fit2$residuals
  kernel_fit3 <- gam(residuals1~s(t3,bs='gp'),data=data1)
  
  x1_middle <- data.frame(t1=c(24.6))#,t2=c(23.5),t3=c(23.2))
  x2_middle <- data.frame(t2=c(23.5))
  x3_middle <- data.frame(t3=c(23.2))
  
  prediction <- predict.gam(kernel_fit,x1_middle) + predict.gam(kernel_fit2,x2_middle) + predict.gam(kernel_fit3,x3_middle)
  return(prediction)
}

bs_kernel <- boot(data, bootstrap_kernel, R=B,parallel = 'multicore')
boot.ci(boot.out = bs_kernel,type = 'bca')
