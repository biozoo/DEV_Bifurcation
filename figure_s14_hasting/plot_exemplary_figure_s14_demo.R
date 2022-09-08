rm(list = ls())
### load packages and set seed
library(rEDM)
#set.seed(45278)
set.seed(2^13+1)

########################################################
# Important Note:
# To significantly reduce the computation time running the analysis including 100 replicates, all results presented in this R-file are based on one time series replicate per mathematical model. 
# Due to process noise (as well as observation error), the DEV results produced from analyzing one replicate look more variable and change less smoothly. 
# However, the findings based on analyzing one replicate still reveal the same dynamical behaviors qualitatively (e.g., increasing trend of |DEV|) 
# and quantitatively (e.g., |DEV|->1 at the bifurcation point).

# To obtain the results presented in the main text, the R-file can be easily modified in four steps:
#  1. Calculate all results of this R-file 100 times.
#  2. Save each set of the results.
#  3. Average across all results associated with the same time step.
#  4. Plot the averaged results instead of the results based on one replicate as shown in this file.
#########################################################

### Stochasticity-driven patch dynamics model (Hastings and Wysham Ecol. Lett. 2010): simulate time series
npatch=8
Tmax=250
a1=1.35;a2=2.37
Tcri=190

N=rep(0.5,npatch)
N1t=NULL
if(T){ # run the model for 1000 steps w.o. changing any parameters
    for(tind in 1:1000){
        zt=runif(1,-0.001,0.001)
        if(zt>0){At=matrix(runif(npatch^2,0.02,0.03),npatch,npatch)}else{At=matrix(runif(npatch^2,0.01,0.02),npatch,npatch)}
        diag(At)=0
        diag(At)=-apply(At,2,sum)
        a=a1
        At=matrix(0.02,npatch,npatch)
        diag(At)=0
        diag(At)=-apply(At,2,sum)
        Nt=N*exp(a*(1-N+zt))
        # dispersal
        N=Nt+c(At%*%Nt)
        N1t=rbind(N1t,N)
    }
}

N1=NULL
a.seq=seq(a1,a2,length.out =Tmax)
zt.seq=runif(Tmax,-0.001,0.001)
for(tind in 1:Tmax){
    # local population growth
    zt=zt.seq[tind] # stochaticity
    a=a.seq[tind] # gradually changed the parameter (growth rate)
    if(zt>0){At=matrix(runif(npatch^2,0.02,0.03),npatch,npatch)
    }else{At=matrix(runif(npatch^2,0.01,0.02),npatch,npatch)}
    
    diag(At)=0;diag(At)=-apply(At,2,sum)
    # local population growth
    Nt=N*exp(a*(1-N+zt))
    # dispersal
    N1=rbind(N1,Nt+c(At%*%Nt))
    N=N1[tind,]
}

real_ts <- N1[,4]
real_x <- c(1:length(real_ts))

### hasting: calculate univariate dev
raw_time_series <- real_ts[1:1:190]
E <- 8
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 150
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

uni_result <- matrix_result

### hasting: calculate multivariate dev
emp_time_series <- N1[1:190,]
window_indices <- seq(window_size, NROW(emp_time_series), step_size)
matrix_result <- matrix(NaN, nrow = length(window_indices), ncol = 4)
index <- 0

#create matrix to save data
for(j in window_indices)
{
    index <- index + 1
    rolling_window_1 <- emp_time_series[(j-window_size+1):j,1]
    rolling_window_2 <- emp_time_series[(j-window_size+1):j,2]
    rolling_window_3 <- emp_time_series[(j-window_size+1):j,3]
    rolling_window_4 <- emp_time_series[(j-window_size+1):j,4]
    rolling_window_5 <- emp_time_series[(j-window_size+1):j,5]
    rolling_window_6 <- emp_time_series[(j-window_size+1):j,6]
    rolling_window_7 <- emp_time_series[(j-window_size+1):j,7]
    rolling_window_8 <- emp_time_series[(j-window_size+1):j,8]
    
    n <- window_size
    vectors <- matrix(NaN, nrow = n, ncol = E+1)
    vectors[,1] <- c(1:n)
    vectors[,2] <- (rolling_window_1-mean(rolling_window_1))/sd(rolling_window_1)
    vectors[,3] <- (rolling_window_2-mean(rolling_window_2))/sd(rolling_window_2)
    vectors[,4] <- (rolling_window_3-mean(rolling_window_3))/sd(rolling_window_3)
    vectors[,5] <- (rolling_window_4-mean(rolling_window_4))/sd(rolling_window_4)
    vectors[,6] <- (rolling_window_5-mean(rolling_window_5))/sd(rolling_window_5)
    vectors[,7] <- (rolling_window_6-mean(rolling_window_6))/sd(rolling_window_6)
    vectors[,8] <- (rolling_window_7-mean(rolling_window_7))/sd(rolling_window_7)
    vectors[,9] <- (rolling_window_8-mean(rolling_window_8))/sd(rolling_window_8)
    colnames(vectors) <- c("time","N1","N2","N3","N4","N5","N6","N7","N8")
    
    
    if(T){ # rEDM package ver 1.2.3
      out_block1 <- block_lnlp(vectors, norm=2, method="s-map", theta=theta_best, target=1, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_1 <- out_block1$smap_coefficients[[1]]
      
      out_block2 <- block_lnlp(vectors, norm=2, method="s-map", theta=theta_best, target=2, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_2 <- out_block2$smap_coefficients[[1]]
      
      out_block3 <- block_lnlp(vectors, norm=2, method="s-map", theta=theta_best, target=3, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_3 <- out_block3$smap_coefficients[[1]]
      
      out_block4 <- block_lnlp(vectors, norm=2, method="s-map", theta=theta_best, target=4, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_4 <- out_block4$smap_coefficients[[1]]
      
      out_block5 <- block_lnlp(vectors, norm=2, method="s-map", theta=theta_best, target=5, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_5 <- out_block5$smap_coefficients[[1]]
      
      out_block6 <- block_lnlp(vectors, norm=2, method="s-map", theta=theta_best, target=6, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_6 <- out_block6$smap_coefficients[[1]]
      
      out_block7 <- block_lnlp(vectors, norm=2, method="s-map", theta=theta_best, target=7, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_7 <- out_block7$smap_coefficients[[1]]
      
      out_block8 <- block_lnlp(vectors, norm=2, method="s-map", theta=theta_best, target=8, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_8 <- out_block8$smap_coefficients[[1]]
    }
    if(F){ # # In earlier version of rEDM package 
      out_block1 <- block_lnlp(vectors, norm_type="L2 norm", method="s-map", theta=theta_best, target=1, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_1 <- out_block1[[1]]$smap_coefficients
      
      out_block2 <- block_lnlp(vectors, norm_type="L2 norm", method="s-map", theta=theta_best, target=2, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_2 <- out_block2[[1]]$smap_coefficients
      
      out_block3 <- block_lnlp(vectors, norm_type="L2 norm", method="s-map", theta=theta_best, target=3, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_3 <- out_block3[[1]]$smap_coefficients
      
      out_block4 <- block_lnlp(vectors, norm_type="L2 norm", method="s-map", theta=theta_best, target=4, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_4 <- out_block4[[1]]$smap_coefficients
      
      out_block5 <- block_lnlp(vectors, norm_type="L2 norm", method="s-map", theta=theta_best, target=5, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_5 <- out_block5[[1]]$smap_coefficients
      
      out_block6 <- block_lnlp(vectors, norm_type="L2 norm", method="s-map", theta=theta_best, target=6, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_6 <- out_block6[[1]]$smap_coefficients
      
      out_block7 <- block_lnlp(vectors, norm_type="L2 norm", method="s-map", theta=theta_best, target=7, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_7 <- out_block7[[1]]$smap_coefficients
      
      out_block8 <- block_lnlp(vectors, norm_type="L2 norm", method="s-map", theta=theta_best, target=8, stats_only=FALSE, exclusion_radius = 0, first_column_time = TRUE, columns = c("time","N1","N2","N3","N4","N5","N6","N7","N8"), save_smap_coefficients=TRUE, tp=1)
      smap_co_8 <- out_block8[[1]]$smap_coefficients
    }
    

    
    
    matrix_eigen <- matrix(NA, nrow = NROW(smap_co_1), ncol = 3)
    
    for(k in 1:NROW(smap_co_1))
    {
        if(!is.na(smap_co_1[k,1]))
        {
            M <- rbind(smap_co_1[k, 1:E], smap_co_2[k, 1:E], smap_co_3[k, 1:E], smap_co_4[k, 1:E], smap_co_5[k, 1:E], smap_co_6[k, 1:E], smap_co_7[k, 1:E], smap_co_8[k, 1:E])
            M_eigen <- eigen(M)$values
            lambda1 <- M_eigen[order(abs(M_eigen))[E]]
            
            matrix_eigen[k,1] <- abs(lambda1)
            matrix_eigen[k,2] <- Re(lambda1)
            matrix_eigen[k,3] <- Im(lambda1)
        }
    }
    
    # calculat save data
    matrix_result[index,1] <- j
    matrix_result[index,2] <- mean(matrix_eigen[,1],na.rm=TRUE)
    matrix_result[index,3] <- mean(matrix_eigen[,2],na.rm=TRUE)
    matrix_result[index,4] <- mean(matrix_eigen[,3],na.rm=TRUE)
}

multi_result <- matrix_result

### plot exemplary figure s14

save.plot.file <- F
if(save.plot.file){pdf("figure_s14.pdf", width = 6, height = 6)}

par(mfrow=c(2,2),oma=c(0,0,2,0))

## plot A
n <- length(real_ts)

par(mar=c(4.1,5.3,1,1)+0.1)

plot(real_x,real_ts, ylim=c(0,2), type="l", col="white", las=0, cex.axis=1.75, cex.lab=1.75, main="", xlab="Time", ylab="Exemplary abundance")

polygon(c(190,250,250,190),c(0,0,2,2), col="red",density=15, border="white")
points(c(190,190),c(min(0),max(2)),type="l", lty=2, col="red")

points(real_x, real_ts, lwd=1, type="l", col="black")

mtext("A",side=3, at=real_x[1], line=0.5,cex=2)

## plot B
n <- NROW(uni_result)

par(mar=c(4.1,5,1,1)+0.1)

plot(uni_result[,1], uni_result[,2], ylim=c(0.4,1.125), lty=1, lwd=0.5, type="p", cex=1, pch=20, col="white", las=0, cex.axis=1.75, cex.lab=1.75, xlim=c(min(real_x), max(real_x)), main="", xlab="Time", ylab="univariate |DEV|")

polygon(c(190,250,250,190),c(0,0,2,2), col="red",density=15, border="white")
points(c(190,190),c(min(0),max(2)),type="l", lty=2, col="red")

points(c(0,12500),c(1,1), type="l", col="gray", lwd=3, lty=2)

bluePal <- colorRampPalette(c('lightblue','blue'))
col <- bluePal(10)[as.numeric(cut(1:length(uni_result[,2]),breaks = 10))]

points(uni_result[,1], uni_result[,2], lty=1, lwd=0.5, type="p", cex=1, pch=20, col=col)

mtext("B",side=3, at=real_x[1], line=0.5, cex=2)

## plot C
n <- NROW(multi_result)

par(mar=c(4.1,5,1,1)+0.1)

plot(multi_result[,1], multi_result[,2], ylim=c(0.4,1.125), lty=1, lwd=0.5, type="p", cex=1, pch=20, col="white", las=0, cex.axis=1.75, cex.lab=1.75, xlim=c(min(real_x), max(real_x)), main="", xlab="Time", ylab="multivariate |DEV|")

polygon(c(190,250,250,190),c(0,0,2,2), col="red",density=15, border="white")
points(c(190,190),c(min(0),max(2)),type="l", lty=2, col="red")

points(c(0,12500),c(1,1), type="l", col="gray", lwd=3, lty=2)


bluePal <- colorRampPalette(c('lightblue','blue'))
col <- bluePal(10)[as.numeric(cut(1:length(multi_result[,2]),breaks = 10))]

points(multi_result[,1], multi_result[,2], lty=1, lwd=0.5, type="p", cex=1, pch=20, col=col)


mtext("C",side=3, at=real_x[1], line=0.5, cex=2)


## plot D

par(mar=c(4.1,5,1,1)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), cex.lab=1.75, cex.axis=1.75, type="l", xlab="Re(DEV)", ylab="Im(DEV)", xlim=c(-1,-0.5), ylim=c(-0.25,0.25))

## complex eigenvalue

n <- NROW(multi_result)

points(c(-1,1), c(0,0), lty=1, lwd=1.5, type="l", col="black")

bluePal <- colorRampPalette(c('lightblue','blue'))
col <- bluePal(10)[as.numeric(cut(1:length(multi_result[,3]),breaks = 10))]

points(multi_result[,3], multi_result[,4], lty=1, lwd=1.5, type="p", cex=2, pch=20, col=col)
points(multi_result[,3], -multi_result[,4], lty=1, lwd=1.5, type="p", cex=2, pch=20, col=col)

mtext("D",side=3, at=-1, line=0.5,cex=2)

#dev.off()
