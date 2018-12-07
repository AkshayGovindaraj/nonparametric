#### Script for wavelet transform

library(waveslim)
library(wmtsa)

x = rnorm(32)

## Part 1
#  x.haar <- mra(x,"haar",4,"dwt")
#  names(x.haar)<- c("d1","d2","d3","d4","s4")
# 
#  
#  x.la8 <- mra(x,"la8",4,"dwt")
#  names(x.la8)<- c("d1","d2","d3","d4","s4")
# 
# ## Part 2
# par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
# plot.ts(x, axes=TRUE, ylab="", main="(a)")
# 
# for(i in 1:5)
# plot.ts(x.haar[[i]], axes=TRUE, ylab=names(x.haar)[i])
# axis(side=1, at=seq(0,368,by=23),
#      labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
# 
# 
# par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
# plot.ts(x, axes=FALSE, ylab="", main="(b)")
# for(i in 1:5)
#   plot.ts(x.la8[[i]], axes=FALSE, ylab=names(x.la8)[i])
# axis(side=1, at=seq(0,368,by=23),
#      labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))

## Using dwt

x_haar = dwt(x, wf='haar', n.levels=4, boundary = 'periodic')
x_la8 = dwt(x, wf='la8', n.levels=4, boundary = 'periodic')

#### Part 3
# Shrinkage Methods
# x.haar.uni <- wavShrink(x_haar,wavelet="haar",shrink.fun="hard",thresh.fun="universal",thresh.scale=1)
# x.la8.uni <- wavShrink(x_la8,wavelet="la8",shrink.fun="hard",thresh.fun="adaptive",thresh.scale=1)

#Haar thresholding
x_sure_haar = sure.thresh(x_haar, max.level = 4, hard = TRUE)
x_univ_haar = universal.thresh(x_haar, max.level = 4, hard = TRUE)

x_sure_inv_haar = idwt(x_sure_haar) 
x_univ_inv_haar = idwt(x_univ_haar)

#la8 thresholding
x_sure_la8 = sure.thresh(x_la8, max.level = 4, hard = TRUE)
x_univ_la8 = universal.thresh(x_la8, max.level = 4, hard = TRUE)

x_sure_inv_la8 = idwt(x_sure_la8)
x_univ_inv_la8 = idwt(x_univ_la8)

#### Part 4
## Local Transformation

fc1 = x_haar
fc2 = x_haar
fc3 = x_haar

#Local changes
fc1$d1[1] = 0
fc2$d1[3] = 0
fc3$d1[9] = 0

#Global changes
fc1$d4[1] = -2
fc2$d4[2] = 2
fc3$d4[1] = 0

#Sure thresholding
fc1_sure = sure.thresh(fc1, max.level = 4, hard = TRUE)
fc2_sure = sure.thresh(fc2, max.level = 4, hard = TRUE)
fc3_sure = sure.thresh(fc3, max.level = 4, hard = TRUE)

### To Add: Multiresolution plot for x_haar, fc1, fc2, fc3 together

x_haar_mra <- mra(idwt(x_haar),"haar",4,"dwt")
fc1_sure_mra <- mra(idwt(fc1_sure),"haar",4,"dwt")
fc2_sure_mra <- mra(idwt(fc2_sure),"haar",4,"dwt")
fc3_sure_mra <- mra(idwt(fc3_sure),"haar",4,"dwt")

par(mfrow=c(5,1))
 for(i in 1:5)
 plot.ts(x_haar_mra[[i]], axes=TRUE, ylab=names(x_haar_mra)[i])
 axis(side=1, at=seq(0,368,by=23),
      labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))


 par(mfrow=c(5,1))
 for(i in 1:5)
   plot.ts(fc1_sure_mra[[i]], axes=TRUE, ylab=names(fc1_sure_mra)[i])
 axis(side=1, at=seq(0,368,by=23),
      labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
 
 par(mfrow=c(5,1))
 for(i in 1:5)
   plot.ts(fc2_sure_mra[[i]], axes=TRUE, ylab=names(fc2_sure_mra)[i])
 axis(side=1, at=seq(0,368,by=23),
      labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
 
 par(mfrow=c(5,1))
 for(i in 1:5)
   plot.ts(fc3_sure_mra[[i]], axes=TRUE, ylab=names(fc3_sure_mra)[i])
 axis(side=1, at=seq(0,368,by=23),
      labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
 




# x_haar$d1[3] = 0 
# x_haarinv = idwt(x_haar)
# 
# x_la8$d1[3] = 0
# x_la8inv = idwt(x_la8)

#Plot


## Global Transformation
# x_haar = dwt(x, wf='haar', n.levels=4, boundary = 'periodic')
# x_la8 = dwt(x, wf='la8', n.levels=4, boundary = 'periodic')
# 
# x_haar$d4[2]=0
# x_haargl_inv = idwt(x_haar)
# 
# x_la8$d4[1] = 1
# x_la8gl_inv = idwt(x_la8)

#Plot


