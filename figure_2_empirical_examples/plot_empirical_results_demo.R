########################################################
# Important Note:
#  Dataset and documentation of all analytical procedures from a scientific paper, 
#  "Anticipating the occurrence and type of critical transitions", that is now under revision in Science Advances.

#   This R file demonstrates the DEV analysis applied in the following three datasets (Fig. 2D-L)
#   (1) Voice experiment data are provided by Dr. Isao Tokuda (Fig. 2D-F) 
#   (2) Cellular ATP experiment data are offered by Dr. Markus Schwarzlander (Fig. 2G-I)
#   (3) Greenhouse earth climate data is publicly available and were downloaded from the World Data Center for Paleoclimatology, 
#       National Geophysical Data Center, Boulder, Colorado (http:// www.ncdc.noaa.gov/paleo/data.html). 
#       This dataset has also been used for computing the other EWS in a previous study (Dakos et al. PNAS 2008)(Fig. 2J-L)

###############################################
## Terms of Use for the dataset (1) & (2)
# I.   You must acknowledge the use of content.
# II.  Monitoring data is made available for use in activities of a non-profit nature only.
# III. Users must contact the database administrator (1) Dr. Isao Tokuda (voice data; isao@fc.ritsumei.ac.jp) 
#      (2) Dr. Markus Schwarzlander (cellular ATP data; markus.schwarzlander@uni-muenster.de) before using the dataset for any publications, 
#      including conference presentations as well as handouts and presentation materials for meetings such as committees and councils. 
#      Co-authorship may be required for some publications depending on the way the time series data is to be used.
###############################################


rm(list = ls())
### Load SSR ###
library(rEDM)

###
save.plot.file <- F
if(save.plot.file){pdf("figure_2.pdf", width = 9, height = 9)}

par(mfrow=c(3,3),oma=c(0,0,2,0))

### voice: calculate dev

raw_time <- read.table("voice.txt", header=F)[,1]
raw_time_series <- read.table("voice.txt", header=F)[,2]

x_trend <- 0.3
x_bif <- 0.49

seg1 <- 0.38995
seg2 <- 0.41495
seg3 <- 0.43995
seg4 <- 0.46495

E <- 6
tau <- 6
theta <- seq(0,2.5,by=0.5)
window_size <- 200
step_size <- 50

### start algorithm
window_indices <- seq(window_size, NROW(raw_time_series), step_size)
matrix_result <- matrix(NaN, nrow = length(window_indices), ncol = 4)
index <- 0

for(j in window_indices)
{
    index <- index + 1
    rolling_window <- raw_time_series[(j-window_size+1):j]
    
    norm_rolling_window <- (rolling_window - mean(rolling_window, na.rm=TRUE))/sd(rolling_window, na.rm=TRUE)
    
    # calculate best theta
    smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta, silent=TRUE)
    best <- order(-smap$rho)[1]
    theta_best <- smap[best,]$theta
    
    # calculate eigenvalues for best theta
    smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta_best, silent=TRUE, save_smap_coefficients=TRUE)
    
    #smap_co <- smap[[1]]$smap_coefficients # In earlier versions of rEDM package
    smap_co <- smap$smap_coefficients[[1]]

    matrix_eigen <- matrix(NA, nrow = NROW(smap_co), ncol = 3)
    
    for(k in 1:NROW(smap_co))
    {
        if(!is.na(smap_co[k,1]))
        {
            M <- rbind(as.numeric(smap_co[k, 1:E]), cbind(diag(E - 1), rep(0, E - 1)))
            M_eigen <- eigen(M)$values
            lambda1 <- M_eigen[order(abs(M_eigen))[E]]
            
            matrix_eigen[k,1] <- abs(lambda1)
            matrix_eigen[k,2] <- Re(lambda1)
            matrix_eigen[k,3] <- Im(lambda1)
        }
    }
    
    # save results
    matrix_result[index,1] <- j
    matrix_result[index,2] <- mean(matrix_eigen[,1],na.rm=TRUE)
    matrix_result[index,3] <- mean(matrix_eigen[,2],na.rm=TRUE)
    matrix_result[index,4] <- mean(matrix_eigen[,3],na.rm=TRUE)
}

result <- matrix_result[which(matrix_result[,1]<min(which(raw_time>x_bif)) & matrix_result[,1]>min(which(raw_time>x_trend))),]

# plot ts

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0, xlim=c(min(raw_time),max(raw_time)), ylim=c(0.2,1.2), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.65, main="", xlab="Time (s)", ylab="Subglottal pressure, kPa")

polygon(c(raw_time[min(which(raw_time>x_bif))],max(raw_time),max(raw_time),raw_time[min(which(raw_time>x_bif))]),c(-100,-100,100,100), col="red",density=15, border="white")

points(c(seg1,seg1),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg2,seg2),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg3,seg3),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg4,seg4),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(x_bif,x_bif),c(-100,100), type="l", col="red", lwd=1.5, lty=2)

points(raw_time,raw_time_series, type="l", col="black")

mtext("A",side=3, at=min(raw_time), line=0.5,cex=2)

# plot dev

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0,xlim=c(min(raw_time),max(raw_time)), ylim=c(0.6,1.1), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.65, main="", xlab="Time (s)", ylab="| DEV |")

polygon(c(raw_time[min(which(raw_time>x_bif))],max(raw_time),max(raw_time),raw_time[min(which(raw_time>x_bif))]),c(-100,-100,100,100), col="red",density=15, border="white")

points(c(x_bif,x_bif),c(-100,100), type="l", col="red", lwd=1.5, lty=2)

points(c(0,12500),c(1,1), type="l", col="gray", lwd=3, lty=2)

mat = cbind(result[,1], result[,2], 1:nrow(result))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(raw_time[result[,1]], result[,2], type="p", cex=1, pch=20, col=cols)

lines(raw_time[result[,1]],predict(lm(result[,2]~result[,1])), col='red', lwd=2)

mtext("B",side=3, at=min(raw_time), line=0.5,cex=2)

## plot dev in complex

par(mar=c(4.1,2.5,1,4.5)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.65, cex.axis=1.75, type="l", xlab="Re(DEV)", ylab="", xlim=c(-0.2,1), ylim=c(-0.6,0.6), yaxt="n")

re <- result[,3]
im <- result[,4]
mat = cbind(re, im, 1:length(re))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

re_seg1 <- mean(result[which(result[,1]<min(which(raw_time>seg2)) & result[,1]>min(which(raw_time>seg1))),3])
im_seg1 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg2)) & result[,1]>min(which(raw_time>seg1))),4]))
re_seg2 <- mean(result[which(result[,1]<min(which(raw_time>seg3)) & result[,1]>min(which(raw_time>seg2))),3])
im_seg2 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg3)) & result[,1]>min(which(raw_time>seg2))),4]))
re_seg3 <- mean(result[which(result[,1]<min(which(raw_time>seg4)) & result[,1]>min(which(raw_time>seg3))),3])
im_seg3 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg4)) & result[,1]>min(which(raw_time>seg3))),4]))
re_seg4 <- mean(result[which(result[,1]<min(which(raw_time>x_bif)) & result[,1]>min(which(raw_time>seg4))),3])
im_seg4 <- -mean(abs(result[which(result[,1]<min(which(raw_time>x_bif)) & result[,1]>min(which(raw_time>seg4))),4]))

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(c(re_seg1, re_seg2, re_seg3, re_seg4), c(im_seg1, im_seg2, im_seg3, im_seg4), lwd=1.5, type="l", cex=4.75, pch=20, col="darkgray")
points(re_seg1, im_seg1, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[15])
points(re_seg2, im_seg2, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[25])
points(re_seg3, im_seg3, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[40])
points(re_seg4, im_seg4, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[60])

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DEV)",side=4,line=3,cex=1.15)

mtext("C",side=3, at=-0.2, line=0.5,cex=2)

### Mitochondria: calculate dev

raw_time <- read.table("mitochondria.txt", header=F)[,1]
raw_time_series <- read.table("mitochondria.txt", header=F)[,2]

x_trend <- 3.2
x_bif <- 4.7

seg1 <- 3.9
seg2 <- 4.1
seg3 <- 4.3
seg4 <- 4.5

E <- 3
tau <- 10
theta <- seq(0,2.5,by=0.5)
window_size <- 125
step_size <- 1

### start algorithm
window_indices <- seq(window_size, NROW(raw_time_series), step_size)
matrix_result <- matrix(NaN, nrow = length(window_indices), ncol = 4)
index <- 0

for(j in window_indices)
{
    index <- index + 1
    rolling_window <- raw_time_series[(j-window_size+1):j]
    
    norm_rolling_window <- (rolling_window - mean(rolling_window, na.rm=TRUE))/sd(rolling_window, na.rm=TRUE)
    
    # calculate best theta
    smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta, silent=TRUE)
    best <- order(-smap$rho)[1]
    theta_best <- smap[best,]$theta
    
    # calculate eigenvalues for best theta
    smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta_best, silent=TRUE, save_smap_coefficients=TRUE)
    
    #smap_co <- smap[[1]]$smap_coefficients # In earlier versions of rEDM package
    smap_co <- smap$smap_coefficients[[1]]
    
    matrix_eigen <- matrix(NA, nrow = NROW(smap_co), ncol = 3)
    
    for(k in 1:NROW(smap_co))
    {
        if(!is.na(smap_co[k,1]))
        {
            M <- rbind(as.numeric(smap_co[k, 1:E]), cbind(diag(E - 1), rep(0, E - 1)))
            M_eigen <- eigen(M)$values
            lambda1 <- M_eigen[order(abs(M_eigen))[E]]
            
            matrix_eigen[k,1] <- abs(lambda1)
            matrix_eigen[k,2] <- Re(lambda1)
            matrix_eigen[k,3] <- Im(lambda1)
        }
    }
    
    # save results
    matrix_result[index,1] <- j
    matrix_result[index,2] <- mean(matrix_eigen[,1],na.rm=TRUE)
    matrix_result[index,3] <- mean(matrix_eigen[,2],na.rm=TRUE)
    matrix_result[index,4] <- mean(matrix_eigen[,3],na.rm=TRUE)
}

result <- matrix_result[which(matrix_result[,1]<min(which(raw_time>x_bif)) & matrix_result[,1]>min(which(raw_time>x_trend))),]

# plot ts

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0, xlim=c(min(raw_time),max(raw_time)), ylim=c(70,90), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.65, main="", xlab="Time (s)", ylab="Relative fluorescence ratio, %")

polygon(c(raw_time[min(which(raw_time>x_bif))],max(raw_time),max(raw_time),raw_time[min(which(raw_time>x_bif))]),c(-100,-100,100,100), col="red",density=15, border="white")

points(c(seg1,seg1),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg2,seg2),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg3,seg3),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg4,seg4),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(x_bif,x_bif),c(-100,100), type="l", col="red", lwd=1.5, lty=2)

points(raw_time,raw_time_series*100, type="l", col="black")

mtext("D",side=3, at=min(raw_time), line=0.5,cex=2)

# plot dev

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0,xlim=c(min(raw_time),max(raw_time)), ylim=c(0.7,1.1), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.65, main="", xlab="Time (s)", ylab="| DEV |")

polygon(c(raw_time[min(which(raw_time>x_bif))],max(raw_time),max(raw_time),raw_time[min(which(raw_time>x_bif))]),c(-100,-100,100,100), col="red",density=15, border="white")

points(c(x_bif,x_bif),c(-100,100), type="l", col="red", lwd=1.5, lty=2)

points(c(0,12500),c(1,1), type="l", col="gray", lwd=3, lty=2)

mat = cbind(result[,1], result[,2], 1:nrow(result))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(raw_time[result[,1]], result[,2], type="p", cex=1, pch=20, col=cols)

lines(raw_time[result[,1]],predict(lm(result[,2]~result[,1])), col='red', lwd=2)

mtext("E",side=3, at=min(raw_time), line=0.5,cex=2)

## plot dev in complex

par(mar=c(4.1,2.5,1,4.5)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.65, cex.axis=1.75, type="l", xlab="Re(DEV)", ylab="", xlim=c(0.8,1), ylim=c(-0.2,0.2), yaxt="n")

re <- result[,3]
im <- result[,4]
mat = cbind(re, im, 1:length(re))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

re_seg1 <- mean(result[which(result[,1]<min(which(raw_time>seg2)) & result[,1]>min(which(raw_time>seg1))),3])
im_seg1 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg2)) & result[,1]>min(which(raw_time>seg1))),4]))
re_seg2 <- mean(result[which(result[,1]<min(which(raw_time>seg3)) & result[,1]>min(which(raw_time>seg2))),3])
im_seg2 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg3)) & result[,1]>min(which(raw_time>seg2))),4]))
re_seg3 <- mean(result[which(result[,1]<min(which(raw_time>seg4)) & result[,1]>min(which(raw_time>seg3))),3])
im_seg3 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg4)) & result[,1]>min(which(raw_time>seg3))),4]))
re_seg4 <- mean(result[which(result[,1]<min(which(raw_time>x_bif)) & result[,1]>min(which(raw_time>seg4))),3])
im_seg4 <- -mean(abs(result[which(result[,1]<min(which(raw_time>x_bif)) & result[,1]>min(which(raw_time>seg4))),4]))

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(c(re_seg1, re_seg2, re_seg3, re_seg4), c(im_seg1, im_seg2, im_seg3, im_seg4), lwd=1.5, type="l", cex=4.75, pch=20, col="darkgray")
points(re_seg1, im_seg1, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[15])
points(re_seg2, im_seg2, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[25])
points(re_seg3, im_seg3, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[40])
points(re_seg4, im_seg4, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[60])

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DEV)",side=4,line=3,cex=1.15)

mtext("F",side=3, at=0.8, line=0.5,cex=2)

### Climate: calculate dev

raw_time <- read.table("climate.txt", header=F)[,2]
raw_time_series <- read.table("climate.txt", header=F)[,1]

x_trend <- -38
x_bif <- -34

seg0 <- -37.831
seg1 <- -37.191
seg2 <- -36.551
seg3 <- -35.911
seg4 <- -35.271
seg5 <- -34.632

E <- 5
tau <- 2
theta <- seq(0,2.5,by=0.5)
window_size <- 150
step_size <- 5

### start algorithm
window_indices <- seq(window_size, NROW(raw_time_series), step_size)
matrix_result <- matrix(NaN, nrow = length(window_indices), ncol = 4)
index <- 0

for(j in window_indices)
{
    index <- index + 1
    rolling_window <- raw_time_series[(j-window_size+1):j]
    
    norm_rolling_window <- (rolling_window - mean(rolling_window, na.rm=TRUE))/sd(rolling_window, na.rm=TRUE)
    
    # calculate best theta
    smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta, silent=TRUE)
    best <- order(-smap$rho)[1]
    theta_best <- smap[best,]$theta
    
    # calculate eigenvalues for best theta
    smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta_best, silent=TRUE, save_smap_coefficients=TRUE)
    
    #smap_co <- smap[[1]]$smap_coefficients # In earlier versions of rEDM package
    smap_co <- smap$smap_coefficients[[1]]
    
    matrix_eigen <- matrix(NA, nrow = NROW(smap_co), ncol = 3)
    
    for(k in 1:NROW(smap_co))
    {
        if(!is.na(smap_co[k,1]))
        {
            M <- rbind(as.numeric(smap_co[k, 1:E]), cbind(diag(E - 1), rep(0, E - 1)))
            M_eigen <- eigen(M)$values
            lambda1 <- M_eigen[order(abs(M_eigen))[E]]
            
            matrix_eigen[k,1] <- abs(lambda1)
            matrix_eigen[k,2] <- Re(lambda1)
            matrix_eigen[k,3] <- Im(lambda1)
        }
    }
    
    # save results
    matrix_result[index,1] <- j
    matrix_result[index,2] <- mean(matrix_eigen[,1],na.rm=TRUE)
    matrix_result[index,3] <- mean(matrix_eigen[,2],na.rm=TRUE)
    matrix_result[index,4] <- mean(matrix_eigen[,3],na.rm=TRUE)
}

result <- matrix_result[which(matrix_result[,1]<min(which(raw_time>x_bif)) & matrix_result[,1]>min(which(raw_time>x_trend))),]

# plot ts

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0, xlim=c(min(raw_time),max(raw_time)), ylim=c(0,100), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.65, main="", xlab=expression(paste("Time (10"^"6"," years ago)")), ylab=expression(paste("CaCO"[3],", %")))

polygon(c(raw_time[min(which(raw_time>x_bif))],max(raw_time),max(raw_time),raw_time[min(which(raw_time>x_bif))]),c(-100,-100,100,100), col="red",density=15, border="white")

points(c(seg0,seg0),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg1,seg1),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg2,seg2),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg3,seg3),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg4,seg4),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg5,seg5),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(x_bif,x_bif),c(-100,100), type="l", col="red", lwd=1.5, lty=2)

points(raw_time,raw_time_series, type="l", col="black")

mtext("G",side=3, at=min(raw_time), line=0.5,cex=2)

# plot dev

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0,xlim=c(min(raw_time),max(raw_time)), ylim=c(0.7,1.1), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.65, main="", xlab=expression(paste("Time (10"^"6"," years ago)")), ylab="| DEV |")

polygon(c(raw_time[min(which(raw_time>x_bif))],max(raw_time),max(raw_time),raw_time[min(which(raw_time>x_bif))]),c(-100,-100,100,100), col="red",density=15, border="white")

points(c(x_bif,x_bif),c(-100,100), type="l", col="red", lwd=1.5, lty=2)

points(c(-1000,12500),c(1,1), type="l", col="gray", lwd=3, lty=2)

mat = cbind(result[,1], result[,2], 1:nrow(result))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(raw_time[result[,1]], result[,2], type="p", cex=1, pch=20, col=cols)

lines(raw_time[result[,1]],predict(lm(result[,2]~result[,1])), col='red', lwd=2)

mtext("H",side=3, at=min(raw_time), line=0.5,cex=2)

## plot dev in complex

par(mar=c(4.1,2.5,1,4.5)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.65, cex.axis=1.75, type="l", xlab="Re(DEV)", ylab="", xlim=c(0.6,1), ylim=c(-0.2,0.2), yaxt="n")

re <- result[,3]
im <- result[,4]
mat = cbind(re, im, 1:length(re))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

re_seg0 <- mean(result[which(result[,1]<min(which(raw_time>seg1)) & result[,1]>min(which(raw_time>seg0))),3])
im_seg0 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg1)) & result[,1]>min(which(raw_time>seg0))),4]))
re_seg1 <- mean(result[which(result[,1]<min(which(raw_time>seg2)) & result[,1]>min(which(raw_time>seg1))),3])
im_seg1 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg2)) & result[,1]>min(which(raw_time>seg1))),4]))
re_seg2 <- mean(result[which(result[,1]<min(which(raw_time>seg3)) & result[,1]>min(which(raw_time>seg2))),3])
im_seg2 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg3)) & result[,1]>min(which(raw_time>seg2))),4]))
re_seg3 <- mean(result[which(result[,1]<min(which(raw_time>seg4)) & result[,1]>min(which(raw_time>seg3))),3])
im_seg3 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg4)) & result[,1]>min(which(raw_time>seg3))),4]))
re_seg4 <- mean(result[which(result[,1]<min(which(raw_time>seg5)) & result[,1]>min(which(raw_time>seg4))),3])
im_seg4 <- -mean(abs(result[which(result[,1]<min(which(raw_time>seg5)) & result[,1]>min(which(raw_time>seg4))),4]))
re_seg5 <- mean(result[which(result[,1]<min(which(raw_time>x_bif)) & result[,1]>min(which(raw_time>seg5))),3])
im_seg5 <- -mean(abs(result[which(result[,1]<min(which(raw_time>x_bif)) & result[,1]>min(which(raw_time>seg5))),4]))

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(c(re_seg0, re_seg1, re_seg2, re_seg3, re_seg4, re_seg5), c(im_seg0, im_seg1, im_seg2, im_seg3, im_seg4, im_seg5), lwd=1.5, type="l", cex=4.75, pch=20, col="darkgray")
points(re_seg0, im_seg0, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[10])
points(re_seg1, im_seg1, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[15])
points(re_seg2, im_seg2, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[25])
points(re_seg3, im_seg3, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[40])
points(re_seg4, im_seg4, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[60])
points(re_seg5, im_seg5, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[70])

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DEV)",side=4,line=3,cex=1.15)

mtext("I",side=3, at=0.6, line=0.5,cex=2)

#dev.off()
