########################################################
# Important Note:
#  Dataset and documentation of all analytical procedures from a scientific paper, 
#  "Anticipating the occurrence and type of critical transitions", that is now in revision pending with Science Advances.
#
#   This R file demonstrates the DEV analysis applied in the following five datasets (Fig. 2A-O)
#   (1) Microcosm experiment data provided by Dr. Egbert H. van Nes (Fig. 2A-C; microcosm.txt). 
#       The complete dataset analyzed in the original paper (Veraart et al. 2012 Nature) has been published on online DANS repository, doi:10.17026/dans-ztg-93aw 
#   (2) Voice experiment data are provided by Dr. Isao Tokuda (Fig. 2D-F; voice.txt) 
#   (3) Cellular ATP experiment data are offered by Dr. Markus Schwarzlander (Fig. 2G-I; mitochondria.txt)
#   (4) Greenhouse earth climate data is publicly available and were downloaded from the World Data Center for Paleoclimatology, 
#       National Geophysical Data Center, Boulder, Colorado (http:// www.ncdc.noaa.gov/paleo/data.html). 
#       This dataset has also been used for computing the other EWS in a previous study (Dakos et al. PNAS 2008)(Fig. 2J-L; climate.txt)
#   (5) Frequency data provided by Bonneville Power Administration (Fig. 2M-O; blackout.txt).

###############################################
## Terms of Use for the datasets
# (I)   You must acknowledge the use of content.
# (II)  Monitoring data is made available for use in activities of a non-profit nature only.
# (III) Users must contact the database administrator listed in Table S2 before using the dataset for any publications, 
#      including conference presentations as well as handouts and presentation materials for meetings such as committees and councils. 
#      Co-authorship may be required for some publications depending on the way the time series data is to be used.
###############################################


rm(list = ls())
### Load SSR ###
library(rEDM)

###
save.plot.file <- T
if(save.plot.file){pdf("figure_2_new.pdf", width = 9, height = 15)}else{win.graph(80,100)}
par(mfrow=c(5,3),oma=c(0,0,2,0))



############################################
### microcosm: calculate dev (Fig. 2A-C)

raw_time <- read.table("microcosm.txt", header=F)[,1]
raw_time_series <- read.table("microcosm.txt", header=F)[,2]

x_trend <- 0
x_bif <- 27

seg0 <- 0.57092
seg0_end <- 1.4539
seg1 <- 4.5496
seg1_end <- 5.4326
seg2 <- 10.642
seg2_end <- 11.525
seg3 <- 14.454
seg3_end <- 15.337
seg4 <- 19.499
seg4_end <- 20.382
seg5 <- 24.502
seg5_end <- 25.385

E <- 6
tau <- 2
theta <- seq(0,2.5,by=0.5)
window_size <- 200
step_size <- 5

### start algorithm
matrix_segments <- matrix(NaN, nrow = 250, ncol = 6)
matrix_segments[,1] <- intersect(which(raw_time>=seg0), which(raw_time<=seg0_end))
matrix_segments[,2] <- intersect(which(raw_time>=seg1), which(raw_time<=seg1_end))
matrix_segments[,3] <- intersect(which(raw_time>=seg2), which(raw_time<=seg2_end))
matrix_segments[,4] <- intersect(which(raw_time>=seg3), which(raw_time<=seg3_end))
matrix_segments[,5] <- intersect(which(raw_time>=seg4), which(raw_time<=seg4_end))
matrix_segments[,6] <- intersect(which(raw_time>=seg5), which(raw_time<=seg5_end))

result <- c()

for(i in 1:6)
{
  time_series_segment <- raw_time_series[matrix_segments[,i]]
  
  window_indices <- seq(window_size, NROW(time_series_segment), step_size)
  matrix_result <- matrix(NaN, nrow = length(window_indices), ncol = 4)
  index <- 0
  
  for(j in window_indices)
  {
    index <- index + 1
    rolling_window <- time_series_segment[(j-window_size+1):j]
    
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
    matrix_result[index,1] <- matrix_segments[,i][j]
    matrix_result[index,2] <- mean(matrix_eigen[,1],na.rm=TRUE)
    matrix_result[index,3] <- mean(matrix_eigen[,2],na.rm=TRUE)
    matrix_result[index,4] <- mean(matrix_eigen[,3],na.rm=TRUE)
  }
  result <- rbind(result,matrix_result)
}

# plot ts

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0, xlim=c(min(raw_time),max(raw_time)), ylim=c(80,170), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.65, main="", xlab="Time (day)", ylab=expression(paste("Light attenuation (m"^-1,")")))

polygon(c(raw_time[min(which(raw_time>x_bif))],max(raw_time),max(raw_time),raw_time[min(which(raw_time>x_bif))]),c(-100,-100,200,200), col="red",density=15, border="white")


polygon(c(seg0,seg0_end,seg0_end,seg0),c(-100,-100,200,200), col=adjustcolor("darkgray", alpha.f = 0.4), border="white")
polygon(c(seg1,seg1_end,seg1_end,seg1),c(-100,-100,200,200), col=adjustcolor("darkgray", alpha.f = 0.4), border="white")
polygon(c(seg2,seg2_end,seg2_end,seg2),c(-100,-100,200,200), col=adjustcolor("darkgray", alpha.f = 0.4), border="white")
polygon(c(seg3,seg3_end,seg3_end,seg3),c(-100,-100,200,200), col=adjustcolor("darkgray", alpha.f = 0.4), border="white")
polygon(c(seg4,seg4_end,seg4_end,seg4),c(-100,-100,200,200), col=adjustcolor("darkgray", alpha.f = 0.4), border="white")
polygon(c(seg5,seg5_end,seg5_end,seg5),c(-100,-100,200,200), col=adjustcolor("darkgray", alpha.f = 0.4), border="white")
points(c(x_bif,x_bif),c(-100,200), type="l", col="red", lwd=1.5, lty=2)

points(raw_time,raw_time_series, type="l", col="black")

text(seg0_end,147, "P1", col="red", cex=1.2)
arrows(seg0_end,151,seg0_end,157, col="red", length=0.1, angle=20, lwd=1.5)
text(seg1_end,146, "P2", col="red", cex=1.2)
arrows(seg1_end,150,seg1_end,156, col="red", length=0.1, angle=20, lwd=1.5)
text(seg2_end,146, "P3", col="red", cex=1.2)
arrows(seg2_end,150,seg2_end,156, col="red", length=0.1, angle=20, lwd=1.5)
text(seg3_end,142.5, "P4", col="red", cex=1.2)
arrows(seg3_end,146.5,seg3_end,152.5, col="red", length=0.1, angle=20, lwd=1.5)
text(seg4_end,139.5, "P5", col="red", cex=1.2)
arrows(seg4_end,143.5,seg4_end,149.5, col="red", length=0.1, angle=20, lwd=1.5)
text(seg5_end,128.5, "P6", col="red", cex=1.2)
arrows(seg5_end,132.5,seg5_end,138.5, col="red", length=0.1, angle=20, lwd=1.5)

mtext("A",side=3, at=min(raw_time), line=0.5,cex=1.5)

# plot dev

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0,xlim=c(min(raw_time),max(raw_time)), ylim=c(0.7,1.1), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.65, main="", xlab="Time (day)", ylab="| DEV |")

polygon(c(raw_time[min(which(raw_time>x_bif))],max(raw_time),max(raw_time),raw_time[min(which(raw_time>x_bif))]),c(-100,-100,200,200), col="red",density=15, border="white")

points(c(x_bif,x_bif),c(-100,200), type="l", col="red", lwd=1.5, lty=2)

points(c(0,12500),c(1,1), type="l", col="gray", lwd=3, lty=2)

mat = cbind(raw_time, raw_time, 1:length(raw_time))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(raw_time[result[1:11,1]], result[1:11,2], type="p", cex=1, pch=20, col=cols[1])
points(raw_time[result[12:22,1]], result[12:22,2], type="p", cex=1, pch=20, col=cols[1409])
points(raw_time[result[23:33,1]], result[23:33,2], type="p", cex=1, pch=20, col=cols[3221])
points(raw_time[result[34:44,1]], result[34:44,2], type="p", cex=1, pch=20, col=cols[4250])
points(raw_time[result[45:55,1]], result[45:55,2], type="p", cex=1, pch=20, col=cols[5300])
points(raw_time[result[56:66,1]], result[56:66,2], type="p", cex=1, pch=20, col=cols[6800])

lines(raw_time[result[,1]],predict(lm(result[,2]~raw_time[result[,1]])), col='red', lwd=2)

mtext("B",side=3, at=min(raw_time), line=0.5,cex=1.5)

## plot dev in complex

par(mar=c(4.1,2.5,1,4.5)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.65, cex.axis=1.75, type="l", xlab="Re(DEV)", ylab="", xlim=c(0.6,1), ylim=c(-0.2,0.2), yaxt="n")

n <- NROW(result)

re <- rep(NA, 6)
im <- rep(NA, 6)

for(i in 1:6)
{
  re[i] <- mean(result[(n/6*(i-1)+1):(n/6*i),3])
  im[i] <- mean(result[(n/6*(i-1)+1):(n/6*i),4])
}

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")
points(re,im, lwd=1.5, type="l", cex=4.75, pch=20, col="darkgray")
points(re[1], im[1], lty=1, lwd=1.5, type="p", cex=4.75, pch=20, col=cols[1])
points(re[2], im[2], lty=1, lwd=1.5, type="p", cex=4.75, pch=20, col=cols[1409])
points(re[3], im[3], lty=1, lwd=1.5, type="p", cex=4.75, pch=20, col=cols[3221])
points(re[4], im[4], lty=1, lwd=1.5, type="p", cex=4.75, pch=20, col=cols[4250])
points(re[5], im[5], lty=1, lwd=1.5, type="p", cex=4.75, pch=20, col=cols[5300])
points(re[6], im[6], lty=1, lwd=1.5, type="p", cex=4.75, pch=20, col=cols[6800])

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DEV)",side=4,line=3,cex=1.15)

mtext("C",side=3, at=0.6, line=0.5,cex=1.5)





############################################
### voice: calculate dev (Fig. 2D-F)

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

mtext("D",side=3, at=min(raw_time), line=0.5,cex=1.5)

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

mtext("E",side=3, at=min(raw_time), line=0.5,cex=1.5)

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
points(re_seg4, im_seg4, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[trunc(quantile(1:length(cols),0.9))])#60

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DEV)",side=4,line=3,cex=1.15)

mtext("F",side=3, at=-0.2, line=0.5,cex=1.5)





############################################
### Mitochondria: calculate dev (Fig. 2G-I)

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

mtext("G",side=3, at=min(raw_time), line=0.5,cex=1.5)

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

mtext("H",side=3, at=min(raw_time), line=0.5,cex=1.5)

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
points(re_seg4, im_seg4, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[trunc(quantile(1:length(cols),0.9))])#70

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DEV)",side=4,line=3,cex=1.15)

mtext("I",side=3, at=0.8, line=0.5,cex=1.5)





############################################
### Climate: calculate dev (Fig. 2J-L)

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

mtext("J",side=3, at=min(raw_time), line=0.5,cex=1.5)

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

mtext("K",side=3, at=min(raw_time), line=0.5,cex=1.5)

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

mtext("L",side=3, at=0.6, line=0.5,cex=1.5)



############################################
### blackout: calculate dev (Fig. 2M-O)

raw_time <- read.table("blackout.txt", header=F)[,1]
raw_time_series <- read.table("blackout.txt", header=F)[,2]

x_trend <- 200
x_bif <- 550

seg1 <- 290
seg2 <- 355
seg3 <- 420
seg4 <- 485

E <- 4
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 1000
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

plot(0,0, xlim=c(min(raw_time),max(raw_time)), ylim=c(-0.4,0.4), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.65, main="", xlab="Time (s)", ylab="Frequency, Hz")

polygon(c(raw_time[min(which(raw_time>x_bif))],max(raw_time),max(raw_time),raw_time[min(which(raw_time>x_bif))]),c(-100,-100,100,100), col="red",density=15, border="white")

points(c(seg1,seg1),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg2,seg2),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg3,seg3),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(seg4,seg4),c(-100,100), type="l", col="gray", lwd=1.5, lty=2)
points(c(x_bif,x_bif),c(-100,100), type="l", col="red", lwd=1.5, lty=2)

points(raw_time,raw_time_series, type="l", col="black")

mtext("M",side=3, at=min(raw_time), line=0.5,cex=1.5)

# plot dev

par(mar=c(4.1,5.3,1,1)+0.1)

plot(0,0,xlim=c(min(raw_time),max(raw_time)), ylim=c(0.8,1.1), type="p", col="white", las=0, cex.axis=1.75, cex.lab=1.65, main="", xlab="Time (s)", ylab="| DEV |")

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

mtext("N",side=3, at=min(raw_time), line=0.5,cex=1.5)

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
points(re_seg4, im_seg4, lty=1, lwd=1, type="p", cex=4.75, pch=20, col=cols[trunc(quantile(1:length(cols),0.9))])# 70

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DEV)",side=4,line=3,cex=1.15)

mtext("O",side=3, at=0.8, line=0.5,cex=1.5)

if(save.plot.file){dev.off()}
