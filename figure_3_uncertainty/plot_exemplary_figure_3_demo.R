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

## observation error
ts_obs_error <- function(TS, CV=0)
{
    TS_out <- TS
    
    for(i in 1:NROW(TS))
    {
        TS_out[i] <- TS_out[i] + rnorm(1)*CV*TS_out[i]
    }
    
    return(TS_out)
}

## noy meir: simulate time series
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
ts_noy_meir <- matrix(NaN, nrow = 10000, ncol = 9)

for(o.n in 1:3)
{
    for(p.n in 1:3)
    {
        o.noise <- c(0,0.05,0.1)[o.n]
        p.noise <- c(0.01,0.05,0.1)[p.n]
        ts_noy_meir[,3*(o.n-1)+p.n] <- ts_obs_error(noy_meir(F=F_bifurcation, sigma=p.noise, noise_on=1), CV=o.noise)
    }
}

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
ts_henon <- matrix(NaN, nrow = 10000, ncol = 9)

for(o.n in 1:3)
{
    for(p.n in 1:3)
    {
        o.noise <- c(0,0.05,0.1)[o.n]
        p.noise <- c(0.01,0.05,0.1)[p.n]
        ts_henon[,3*(o.n-1)+p.n] <- ts_obs_error(henon(r=r_bifurcation, sigma=p.noise, noise_on=1)[,1], CV=o.noise)
    }
}

### mac arthur: simulate time series
arthur <- function(l=3.4, t_start=1, t_end=10000, x0=c(0.665,1.085), sigma=0, noise_on=0)
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
ts_arthur <- matrix(NaN, nrow = 10000, ncol = 9)

for(o.n in 1:3)
{
    for(p.n in 1:3)
    {
        o.noise <- c(0,0.05,0.1)[o.n]
        p.noise <- c(0.001,0.0025,0.004)[p.n]
        ts_arthur[,3*(o.n-1)+p.n] <- ts_obs_error(arthur(l=l_bifurcation, sigma=p.noise, noise_on=1)[,1], CV=o.noise)
    }
}

### plot figure 3
save.plot.file <- F
if(save.plot.file){pdf("figure_3.pdf", width = 6, height = 5)}

par(mfrow=c(1,1),oma=c(0,0,2,0))

## dev vs rho
par(mar=c(4.1,5,1,1)+0.1)

plot(0, 0, ylim=c(0,1), lty=1, lwd=0.5, type="p", cex=1.75, pch=20, col="white", las=0, cex.axis=2, cex.lab=2.25, xlim=c(-0.1,1), main="", xlab=expression(rho), ylab="Eucledean deviation")

r <- seq(0.1, 0.37, by = 0.27/10000)
L <- matrix(NaN, nrow = length(r), ncol = 2)
L[,1] <- r

for (z in 1:length(r))
{
    x <- 1/(2*r[z]) * (0.3 - 1 + sqrt((0.3-1)^2 + 4*r[z]) )
    L[z,2] <- -r[z]*x-sqrt(r[z]^2*x^2+0.3)
}

ana.re <- Re(L[,2])
ana.im <- Im(L[,2])

for(i in 1:9)
{
    raw_time_series <- ts_henon[1:8500,i]
    E <- 3
    tau <- 1
    theta <- seq(0,2.5,by=0.5)
    window_size <- 100
    step_size <- 500
    
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
        rho_best <- smap[best,]$rho
        
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
        matrix_result[index,2] <- rho_best
        matrix_result[index,3] <- mean(sqrt((matrix_eigen[,2]-ana.re[j])^2+(matrix_eigen[,3]+ana.im[j])^2),na.rm=TRUE)
    }
    
    points(matrix_result[,2],matrix_result[,3], type="p", cex=1, pch=21, col="red")
}

r <- seq(3.48, 3.73, by = 0.25/10000)
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

ana.re <- L[,2]
ana.im <- L[,3]

for(i in 1:9)
{
    raw_time_series <- ts_arthur[1:8500,i]
    E <- 6
    tau <- 1
    theta <- seq(0,2.5,by=0.5)
    window_size <- 250
    step_size <- 500
    
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
        rho_best <- smap[best,]$rho
        
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
        matrix_result[index,2] <- rho_best
        matrix_result[index,3] <- mean(sqrt((matrix_eigen[,2]-ana.re[j])^2+(matrix_eigen[,3]+ana.im[j])^2),na.rm=TRUE)
    }
    
    points(matrix_result[,2],matrix_result[,3], type="p", cex=1, pch=21, col="purple")
}

noy_meir_jacob <- function(x, F=1)
{
    
    r <- 0.75
    b <- 0.1
    h <- 0.75
    
    out <- x * exp(r - b * x) - F * (x^2)/(x^2 + h^2)
    
    return(out)
}

F <- seq(0, 2, by = 2/10000)
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

ana.re <- L[,2]
ana.im <- L[,3]

for(i in 1:9)
{
    raw_time_series <- ts_noy_meir[1:9000,i]
    E <- 2
    tau <- 1
    theta <- seq(0,2.5,by=0.5)
    window_size <- 250
    step_size <- 500
    
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
        rho_best <- smap[best,]$rho
        
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
        matrix_result[index,2] <- rho_best
        matrix_result[index,3] <- mean(sqrt((matrix_eigen[,2]-ana.re[j])^2+(matrix_eigen[,3]+ana.im[j])^2),na.rm=TRUE)
    }
    
    points(matrix_result[,2],matrix_result[,3], type="p", cex=1, pch=21, col="blue")
}

#dev.off()
