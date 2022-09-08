rm(list = ls())
### load packages and set seed
library(rEDM)
library(pracma)
set.seed(2^14+1)

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

### noy meir: simulate time series
noy_meir <- function(r=0.75, b=0.1, h=0.75, p=2, F=1, t_start=1, t_end=10000, x0=7.5, sigma=0, noise_on=0)
{
    F <- if(length(F)==1){rep(F,t_end)}else{F}

    ricker <- function(N, r, b, h, p, F, noise_add=0)
    {
        out <- N*exp(r-b*N)-F*(N^p/(N^p+h^p)) + noise_add
        
        return(out)
    }
    
    L <- matrix(NaN, nrow = t_end+1-t_start, ncol = 1)
    L[1] <- x0
    
    for(i in 1:(NROW(L)-1))
    {
        L_noise_add <- if(noise_on!=0){sigma*L[i]*rnorm(1)}else{0}
        L[i+1] <- ricker(N=L[i], noise_add=L_noise_add, r=r, b=b, h=h, p=p, F=F[i])
    }
    
    return(L)
}

F_bifurcation <- 2/10000*c(0:(10000-1))
ts_noy_meir <- noy_meir(F=F_bifurcation, sigma=0.01, noise_on=1)

### henon: simulate time series
henon <- function(r=0.5, b=0.3, t_start=1, t_end=10000, x0=c(1.21,1.21), sigma=0, noise_on=0)
{
    r <- if(length(r)==1){rep(r,t_end)}else{r}
    
    hen_map <- function(N, r, b, noise_add=0)
    {
        N1 <- 1 - r*N[1]^2 + b*N[2] + noise_add
        N2 <- N[1]
        
        out <- c(N1, N2)
        
        return(out)
    }
    
    L <- matrix(NaN, nrow = t_end+1-t_start, ncol = 2)
    L[1,] <- x0
    
    for(i in 1:(NROW(L)-1))
    {
        L_noise_add <- if(noise_on!=0){sigma*L[i,1]*rnorm(1)}else{0}
        L[i+1,] <- hen_map(N=L[i,], noise_add=L_noise_add, r=r[i], b=b)
    }
    
    return(L)
}

r_bifurcation <- 0.1 + 0.3/10000*c(0:(10000-1))
ts_henon <- henon(r=r_bifurcation, sigma=0.01, noise_on=1)

### mac arthur: simulate time series
arthur <- function(l=3.4, t_start=1, t_end=10000, x0=c(0.62,1.4), sigma=0, noise_on=0)
{
  l <- if(length(l)==1){rep(l,t_end)}else{l}
  
  hopf <- function(N, l, noise_add=0)
  {
    a <- 4
    b <- 1
    c <- 6
    m <- 0.5
    mu <- -3
    
    N1 <- (1+l)*N[1] - a*N[1]^2 - b*N[1]*N[2]/(1+m*N[1]) + noise_add
    N2 <- (1+mu)*N[2] + c*N[1]*N[2]/(1+m*N[1])
    
    out <- c(N1, N2)
    
    return(out)
  }
  
  L <- matrix(NaN, nrow = t_end+1-t_start, ncol = 2)
  L[1,] <- x0
  
  for(i in 1:(NROW(L)-1))
  {
    L_noise_add <- if(noise_on!=0){rnorm(1)*sigma}else{0}
    
    L[i+1,] <- hopf(N=L[i,], l=l[i], noise_add=L_noise_add)
  }
  
  return(L)
}

l_bifurcation <- 3.48 + 0.28/10000*c(0:(10000-1))
ts_arthur <- arthur(l=l_bifurcation, sigma=0.001, noise_on=1)


### noy meir: calculate dev
raw_time_series <- ts_noy_meir[1:9000]
E <- 2
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 250
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

noy_meir_result <- matrix_result

### henon: calculate dev
raw_time_series <- ts_henon[1:8500,1]
E <- 3
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 100
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

henon_result <- matrix_result

### mac arthur: calculate dev
raw_time_series <- ts_arthur[1:8500,1]
E <- 6
tau <- 1
theta <- seq(0,2.5,by=0.5)
window_size <- 100
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

arthur_result <- matrix_result

### plot figure 1 for single exemplary time series
save.plot.file <- F
if(save.plot.file){pdf("figure_1.pdf", width = 9, height = 7.5)}

win.graph(70,60)
par(mfrow=c(2,2),oma=c(0,0,2,0))

# plot A
par(mar=c(4.1,5,2,4)+0.1)
plot(0,0, col="white", ylim=c(0, 9), xlim=c(0, 10000), las=0, cex.axis=1.65, cex.lab=1.75, main="", xlab="", ylab="Variable N", yaxt="n")

raw_time_series <- ts_noy_meir
mat = cbind(raw_time_series, raw_time_series, 1:length(raw_time_series))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T)]

for(i in 2:10000)
{
    points(c((i-1):i),raw_time_series[(i-1):i], type="l", lwd=1, pch=20, col=cols[i])
}

mtext(text = "Time", side = 1, line = 2.5, cex = 1.45)
mtext("A",side=3, at=0, line=0.5,cex=2)

#col of axis
axis(2, las=0, col="blue", cex.axis=1.65, lwd=2)

par(new=TRUE)
plot(0,0, lty=1, lwd=3, ylim=c(0,1.75), xlim=c(0, 10000), type="l", col="white", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

raw_time_series <- ts_henon[,1]
mat = cbind(raw_time_series, raw_time_series, 1:length(raw_time_series))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('tomato','red'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

for(i in 2:10000)
{
    points(c((i-1):i),raw_time_series[(i-1):i], type="l", lwd=1, pch=20, col=cols[i])
}

raw_time_series <- ts_arthur[,1]
mat = cbind(raw_time_series, raw_time_series, 1:length(raw_time_series))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('plum','purple'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

for(i in 2:10000)
{
    points(c((i-1):i),raw_time_series[(i-1):i], type="l", lwd=1, pch=20, col=cols[i])
}

axis(4, las=0, col="black", cex.axis=1.65, lwd=2)
mtext("Variable x",side=4,line=3,cex=1.45)

legend("bottomleft", lwd=4.5, cex=1.5, col=c("blue", "red", "purple"), legend=c("Noy-Meir model", "Henon map", "MacArthur model"),bty="n")

# plot B
par(mar=c(4.1,5,2,4)+0.1)
plot(-1, -1, xlim = c(0, 1), col="white", ylim=c(0, 9), las=0, cex.axis=1.65, cex.lab=1.75, main="", xlab="", ylab="Variable N", yaxt="n", xaxt="n")

rmax <- 2
F <- seq(0, rmax, by = 0.001)
h <- 0.75
b <- 0.1
r <- 0.75

mat = cbind(1:length(F), 1:length(F), 1:length(F))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

n <- 100

for (z in 1:length(F)) {
    xl <- vector()
    xl[1] <- 10
    for (i in 2:n) {
        
        xl[i] <- xl[i - 1] * exp(r - b * xl[i - 1]) - F[z] * (xl[i - 1]^2)/(xl[i - 1]^2 + h^2)
        
    }
    uval <- unique(xl[70:n])
    points(rep(F[z], length(uval))/rmax, uval, cex = 0.3, pch = 19, col=cols[z])
}

mtext(text = "Bifurcation parameter", side = 1, line = 2.5, cex = 1.45)
mtext("B",side=3, at=0, line=0.5,cex=2)

#col of axis
axis(2, las=0, col="blue", cex.axis=1.65, lwd=2)

par(new=TRUE)

plot(-1, -1, xlim = c(0, 1), col="white", ylim=c(0, 1.75), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

rmax <- 0.38
r <- seq(0.1, rmax, by = 0.0001)

mat = cbind(1:length(r), 1:length(r), 1:length(r))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('tomato','red'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

n <- 100

for (z in 1:length(r)) {
    xl <- vector()
    xy <- vector()
    xl[1] <- 0.1
    xy[1] <- 0.4
    for (i in 2:n) {
        
        xl[i] <- 1 - r[z] * xl[i - 1]^2 + 0.3*xy[i - 1]
        xy[i] <- xl[i - 1]
        
    }
    uval <- unique(xl[80:n])
    points((rep(r[z], length(uval))-0.1)/(rmax-0.1), uval, cex = 0.3, pch = 19, col=cols[z])
}

rmax <- 3.76
F <- seq(3.48, rmax, by = 0.0001)

mat = cbind(1:length(F), 1:length(F), 1:length(F))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('plum','purple'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

n <- 40

for (z in 1:length(F)) {
    
    xl <- arthur(l=F[z],t_end=1000)
    
    uval <- unique(xl[(1000-n):1000,1])
    points((rep(F[z], length(uval))-3.48)/(rmax-3.48), uval, cex = 0.3, pch = 19, col=cols[z])
}

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65, lwd=2)
mtext("Variable x",side=4,line=3,cex=1.45)

#plot C
par(mar=c(4.1,5,2,4)+0.1)

plot(0,0, col="white", ylim=c(0, 1.25), xlim=c(0, 10000), las=0, cex.axis=1.65, cex.lab=1.75, main="", xlab="", ylab="| DEV |")

points(c(0,10000), c(1,1), col="gray", type="l", lwd=3, lty=2)

mat = cbind(henon_result[1:169,1], henon_result[1:169,2], 1:169)
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('tomato','red'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(henon_result[1:169,1], henon_result[1:169,2], lty=1, lwd=5, type="p", cex=0.5, pch=20, col=cols)

mat = cbind(arthur_result[1:169,1], arthur_result[1:169,2], 1:169)
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('plum','purple'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(arthur_result[1:169,1], arthur_result[1:169,2], lty=1, lwd=5, type="p", cex=0.5, pch=20, col=cols)

mat = cbind(noy_meir_result[,1], noy_meir_result[,2], 1:nrow(noy_meir_result))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(noy_meir_result[,1], noy_meir_result[,2], lty=1, lwd=5, type="p", cex=0.5, pch=20, col=cols)

mtext(text = "Time", side = 1, line = 2.5, cex = 1.45)
mtext("C",side=3, at=-1, line=0.5,cex=2)

# plot D
noy_meir_jacob <- function(x, F=1)
{
    
    r <- 0.75
    b <- 0.1
    h <- 0.75
    
    out <- x * exp(r - b * x) - F * (x^2)/(x^2 + h^2)
    
    return(out)
}

F <- seq(0, 2, by = 0.0002)
L <- matrix(NaN, nrow = length(F), ncol = 3)
L[,1] <- F

for (z in 1:length(F))
{
    # calculate equilibrium
    equilibrium <- mean(noy_meir(x0=8, F=F[z],t_end=1000)[900:1000])
    
    eigenvalues <- fderiv(noy_meir_jacob, x=equilibrium, F=F[z])
    L[z,2] <- Re(eigenvalues)
    L[z,3] <- Im(eigenvalues)
}

L_noy <- L

r <- seq(0.1, 0.37, by = 0.1/10000)
L <- matrix(NaN, nrow = length(r), ncol = 2)
L[,1] <- r

for (z in 1:length(r))
{
    x <- 1/(2*r[z]) * (0.3 - 1 + sqrt((0.3-1)^2 + 4*r[z]) )
    L[z,2] <- -r[z]*x-sqrt(r[z]^2*x^2+0.3)
}

L_hen <- L

r <- seq(3.48, 3.73, by = 0.25/1000)
L <- matrix(NaN, nrow = length(r), ncol = 3)
L[,1] <- r

for (z in 1:length(r))
{
    a <- 4
    b <- 1
    c <- 6
    m <- 0.5
    mu <- -3
    
    eig <- eigen(matrix(c((c^2+c*(a+m-m*r[z])*mu-m*(a+m*r[z])*mu^2)/(c*(c+m*mu)),b*mu/c,(c*r[z]+a*mu+m*r[z]*mu)/b,1), nrow=2, ncol=2, byrow=T))$values[1]
    L[z,2] <- Re(eig)
    L[z,3] <- Im(eig)
}

L_hopf <- L

## unit-circle
par(mar=c(4.1,5,2,6)+0.1)

f <- function(x) exp(-(0+1i)*x)
x <- seq(0, 2*pi, by=0.01)

plot(f(x), type="l", cex.lab=1.75, cex.axis=1.65, xlab="", ylab="", yaxt="n")

re_hen <- henon_result[,3]
im_hen <- henon_result[,4]
mat = cbind(re_hen, im_hen, 1:length(re_hen))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('tomato','red'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(re_hen, im_hen, lty=1, lwd=5, type="p", cex=0.5, pch=20, col=cols)
points(re_hen, -im_hen, lty=1, lwd=5, type="p", cex=0.5, pch=20, col=cols)
points(Re(L_hen[,2]), Im(L_hen[,2]), col="black",type="l", lty=1, lwd=3.5)

re_noy <- noy_meir_result[,3]
im_noy <- noy_meir_result[,4]
mat = cbind(re_noy, im_noy, 1:length(re_noy))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('lightblue','blue'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(mat[, 1:2], lty=1, lwd=5, type="p", cex=0.5, pch=20, col=cols)
points(re_noy, -im_noy, lty=1, lwd=5, type="p", cex=0.5, pch=20, col=cols)
points(L_noy[,2],L_noy[,3], col="black",type="l", lty=1, lwd=3.5)

re_hopf <- arthur_result[,3]
im_hopf <- arthur_result[,4]
mat = cbind(re_hopf, im_hopf, 1:length(re_hopf))
n = 255
data_seq = seq(min(mat[,3]), max(mat[,3]), length=n)
col_pal = colorRampPalette(c('plum','purple'))(n+1)
cols = col_pal[ cut(mat[,3], data_seq, include.lowest=T) ]

points(re_hopf, im_hopf, lty=1, lwd=5, type="p", cex=0.5, pch=20, col=cols)
points(re_hopf, -im_hopf, lty=1, lwd=5, type="p", cex=0.5, pch=20, col=cols)
points(L_hopf[200:1001,2],L_hopf[200:1001,3], col="black",type="l", lty=1, lwd=3.5)
points(L_hopf[200:1001,2],-L_hopf[200:1001,3], col="black",type="l", lty=1, lwd=3.5)

#2nd y-axis
axis(4, las=0, col="black", cex.axis=1.65)
mtext("Im(DEV)",side=4,line=3,cex=1.5)

mtext(text = "Re(DEV)", side = 1, line = 2.5, cex = 1.45)
mtext("D",side=3, at=-1, line=0.5,cex=2)

#dev.off()
