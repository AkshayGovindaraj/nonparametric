#### Spline and Kernel Final

library(locpol) 
library(splines)
library(boot)
library(gam)
library(mgcv)
data<- read.csv("octopus.csv",header=TRUE)

y <- data$totWt
x1 <- data$upBeak
x2 <- data$loBeak
x3 <- data$latWall
nrows <- length(y)
plot(x1,y,xlab="Upper Beak length", ylab="Total Weight")
plot(x2,y,xlab="Lower Beak length", ylab="Total Weight")
plot(x3,y,xlab="Lateral Wall length", ylab="Total Weight")

## Spline regression

sspline_fit_x1 <- smooth.spline(x1,y)
fit.ss_x1 <- fitted(sspline_fit_x1)

plot(fit.ss_x1,y-fit.ss_x1)    # Residuals plot
plot(fit.ss_x1)                # Predictions

sspline_fit_x2 <- smooth.spline(x2,y-fit.ss_x1)
fit.ss_x2 <- fitted(sspline_fit_x2)

plot(fit.ss_x1 + fit.ss_x2,y)      # Combined predicions
plot(fit.ss_x2,y-fit.ss_x1 - fit.ss_x2)      # Combined Residuals

sspline_fit_x3 <- smooth.spline(x3,y-fit.ss_x1-fit.ss_x2)
fit.ss_x3 <- fitted(sspline_fit_x3)

plot(fit.ss_x1 + fit.ss_x2 + fit.ss_x3, y)      # Combined predicions

plot(fit.ss_x1 + fit.ss_x2 + fit.ss_x3 - y)      # Final Residuals
plot(fit.ss_x1 + fit.ss_x2  - y)   
plot(fit.ss_x1  - y)   
# LOOCV Spline

Q3_P1_sspline <- function(variable){
  ss_res_sspline <- c(0)
  x_1 = variable$upBeak
  x_2 = variable$loBeak
  x_3 = variable$latWall
  y_1 = variable$totWt
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
    y_pred <- predict(sspline_fit_x1,x_1[i])$y + predict(sspline_fit_x2,x_2[i])$y + predict(sspline_fit_x3,x_3[i])$y
    ss_res_sspline<- ss_res_sspline + (y_1[i]-y_pred)^2
  }
  ss_res_sspline <- ss_res_sspline/nrows
  return(sqrt(ss_res_sspline))
}

Q3_P1_sspline(data)     # Convert to root mean square error


## Kernel regression

# kernel_fit_x1 <- ksmooth(x1,y,kernel="normal",bandwidth=20)
kernel_fit <- gam(totWt~s(upBeak,bs='gp'),data=data)
data1 <- data
data1$residuals <- kernel_fit$residuals
# plot(y-kernel_fit_x1$y)   # Residual plot

kernel_fit2 <- gam(residuals~s(loBeak,bs='gp'),data=data1)
data1$residuals2 <- kernel_fit2$residuals
# plot(kernel_fit2$residuals)   # Residual plot

# kernel_fit_x3 <- ksmooth(x3,y-kernel_fit_x1$y-kernel_fit_x2$y,kernel="normal",bandwidth=20)
kernel_fit3 <- gam(residuals2~s(latWall,bs='gp'),data=data1)

# plot(y-kernel_fit_x1$y-kernel_fit_x2$y-kernel_fit_x3$y)   # Residual plot

# LOOCV Kernel

Q3_P1_kernel <- function(variable){
  ss_res_kernel <- c(0)
  x_1 = variable$upBeak
  x_2 = variable$loBeak
  x_3 = variable$latWall
  y_1 = variable$totWt
  for(i in 1:nrows){
    x1 = x_1[-i]
    x2 = x_2[-i]
    x3 = x_3[-i]
    y1 = y_1[-i]
    data_kn = data.frame(x1,x2,x3,y1)
    kernel_fit <- gam(y1~s(x1,bs='gp'),data=data_kn)
    data_kn$residuals <- kernel_fit$residuals
    kernel_fit2 <- gam(residuals~s(x2,bs='gp'),data=data_kn)
    data_kn$residuals2 <- kernel_fit2$residuals
    kernel_fit3<- gam(residuals2~s(x3,bs='gp'),data=data_kn)
    data_kn$residuals2 <- kernel_fit2$residuals
    y_pred <- predict.gam(kernel_fit, data.frame(x1=c(x_1[i]))) + predict.gam(kernel_fit2, data.frame(x2=c(x_2[i]))) + predict.gam(kernel_fit3, data.frame(x3=c(x_3[i])))
    ss_res_kernel<- ss_res_kernel + (y_1[i]-y_pred)^2
  }
  ss_res_kernel <- ss_res_kernel/nrows
  ss_res_kernel <- data.frame(ss_res_kernel)
  return(sqrt(ss_res_kernel))
}

Q3_P1_kernel(data)
## Combined plots

plot(y,ylim=c(0,6000),main="Spline vs Kernel Fit for variable x3")
lines(kernel_fit$fitted.values+kernel_fit2$fitted.values+kernel_fit3$fitted.values,lwd=2,type='l',col="purple")
lines(fit.ss_x1+fit.ss_x2+fit.ss_x3,lwd=2,type='l',col="red")
legend(x=2,y=5000,legend=c("Kernel Fit","Spline Fit "),lty=1,col=c("purple","red"))

### Part 2

# Choosing first point

x1_1 = 137
x2_1 = 172
x3_1 = 153

# Choosing second point

x1_2 = 136
x2_2 = 166
x3_2 = 227

## Predictions based on these points

y1_pred_spline = predict(sspline_fit_x1,x1_1)$y + predict(sspline_fit_x2,x2_1)$y + predict(sspline_fit_x3,x3_1)$y
y2_pred_spline = predict(sspline_fit_x1,x1_2)$y + predict(sspline_fit_x2,x2_2)$y + predict(sspline_fit_x3,x3_2)$y

kernel_fit <- gam(totWt~s(upBeak,bs='gp'),data=data)
data1 <- data
data1$residuals <- kernel_fit$residuals
kernel_fit2 <- gam(residuals~s(loBeak,bs='gp'),data=data1)
data1$residuals1 <- kernel_fit2$residuals
kernel_fit3 <- gam(residuals1~s(latWall,bs='gp'),data=data1)

x1_middle <- data.frame(upBeak=c(x1_1))
x2_middle <- data.frame(loBeak=c(x2_1))
x3_middle <- data.frame(latWall=c(x3_1))

prediction_mid <- predict.gam(kernel_fit,x1_middle) + predict.gam(kernel_fit2,x2_middle) + predict.gam(kernel_fit3,x3_middle)

x1_edge <- data.frame(upBeak=c(x1_2))#
x2_edge <- data.frame(loBeak=c(x2_2))
x3_edge <- data.frame(latWall=c(x3_2))

prediction_edge <- predict.gam(kernel_fit,x1_edge) + predict.gam(kernel_fit2,x2_edge) + predict.gam(kernel_fit3,x3_edge)


## Part 3 : Bootstrap
# Spline
#Middle
bootstrap_spline1 <- function(d,i,x){
  data_source <- d[i,]
  y <- data_source$totWt
  x1 <- data_source$upBeak
  x2 <- data_source$loBeak
  x3 <- data_source$latWall
  bsspline_fit_x1 <- smooth.spline(x1,y)
  fit.ss_x1 <- fitted(bsspline_fit_x1)
  bsspline_fit_x2 <- smooth.spline(x2,y-fit.ss_x1)
  fit.ss_x2 <- fitted(bsspline_fit_x2)
  bsspline_fit_x3 <- smooth.spline(x3,y-fit.ss_x1-fit.ss_x2)
  fit.ss_x3 <- fitted(bsspline_fit_x3)
  prediction <- predict(bsspline_fit_x1,x1_1)$y + predict(bsspline_fit_x2,x2_1)$y + predict(bsspline_fit_x3,x3_1)$y
  return(prediction)
}

B <- 1e02
bs_spline <- boot(data, bootstrap_spline1, R=B,parallel = 'multicore')
CI_mid_spline <- boot.ci(boot.out = bs_spline,conf = 0.9,type = 'bca')

#Edge

bootstrap_spline2 <- function(d,i,x){
  data_source <- d[i,]
  y <- data_source$totWt
  x1 <- data_source$upBeak
  x2 <- data_source$loBeak
  x3 <- data_source$latWall
  bsspline_fit_x1 <- smooth.spline(x1,y)
  fit.ss_x1 <- fitted(bsspline_fit_x1)
  bsspline_fit_x2 <- smooth.spline(x2,y-fit.ss_x1)
  fit.ss_x2 <- fitted(bsspline_fit_x2)
  bsspline_fit_x3 <- smooth.spline(x3,y-fit.ss_x1-fit.ss_x2)
  fit.ss_x3 <- fitted(bsspline_fit_x3)
  prediction <- predict(bsspline_fit_x1,x1_2)$y + predict(bsspline_fit_x2,x2_2)$y + predict(bsspline_fit_x3,x3_2)$y
  return(prediction)
}

bs_spline2 <- boot(data, bootstrap_spline2, R=B,parallel = 'multicore')
CI_edge_spline <- boot.ci(boot.out = bs_spline2,conf = 0.9,type = 'bca')

# Kernel
# Mid

bootstrap_kernel1 <- function(d,i,x){
  data_source <- d[i,]
  kernel_fit <- gam(totWt~s(upBeak,bs='gp'),data=data_source)
  data1 <- data_source
  data1$residuals <- kernel_fit$residuals
  kernel_fit2 <- gam(residuals~s(loBeak,bs='gp'),data=data1)
  data1$residuals1 <- kernel_fit2$residuals
  kernel_fit3 <- gam(residuals1~s(latWall,bs='gp'),data=data1)
  
  x1_middle <- data.frame(upBeak=c(x1_1))
  x2_middle <- data.frame(loBeak=c(x2_1))
  x3_middle <- data.frame(latWall=c(x3_1))
  
  prediction <- predict.gam(kernel_fit,x1_middle) + predict.gam(kernel_fit2,x2_middle) + predict.gam(kernel_fit3,x3_middle)
  return(prediction)
}

bs_kernel1 <- boot(data, bootstrap_kernel1, R=B,parallel = 'multicore')
CI_mid_kernel1 <- boot.ci(boot.out = bs_kernel1,type = 'bca')

#Edge
bootstrap_kernel2 <- function(d,i,x){
  data_source <- d[i,]
  kernel_fit <- gam(totWt~s(upBeak,bs='gp'),data=data_source)
  data1 <- data_source
  data1$residuals <- kernel_fit$residuals
  kernel_fit2 <- gam(residuals~s(loBeak,bs='gp'),data=data1)
  data1$residuals1 <- kernel_fit2$residuals
  kernel_fit3 <- gam(residuals1~s(latWall,bs='gp'),data=data1)
  
  x1_middle <- data.frame(upBeak=c(x1_2))
  x2_middle <- data.frame(loBeak=c(x2_2))
  x3_middle <- data.frame(latWall=c(x3_2))
  
  prediction <- predict.gam(kernel_fit,x1_middle) + predict.gam(kernel_fit2,x2_middle) + predict.gam(kernel_fit3,x3_middle)
  return(prediction)
}

bs_kernel2 <- boot(data, bootstrap_kernel2, R=B,parallel = 'multicore')
CI_edge_kernel2 <- boot.ci(boot.out = bs_kernel2,type = 'bca')