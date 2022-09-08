rm(list = ls())
### Load SSR ###
library(rEDM)

### Note that S-Map calculations for large window sizes and/or time series need a lot of computational time

### noy meir: simulate time series
noy_meir <- function(r=0.75, b=0.1, h=0.75, p=2, F=1, t_start=1, t_end=1000, x0=7.5, sigma=0, noise_on=0)
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

F_bifurcation <- 2/1000*c(0:(1000-1))
ts_noy_meir <- noy_meir(F=F_bifurcation, sigma=0.01, noise_on=1)

### set up parameter
raw_time_series <- ts_noy_meir[1:820] #from noy_meir.R
E <- seq(1,5,by=1)
tau <- seq(1,5,by=1)
theta <- 0
window_size <- c(100,150,200,250,300,400,500,600)
step_size <- 25

### start algorithm
matrix_results <- data.frame()

for(i in 1:length(window_size))
{
    for(j in 1:length(E))
    {
        for(k in 1:length(tau))
        {
            for(l in 1:length(theta))
            {
                matrix_rho <- c()
                m <- 0
                
                while(m <= length(raw_time_series) - window_size[i] - step_size)
                {
                    raw_ts_part <- raw_time_series[(m + 1):(m + window_size[i])]
                    time_series <- (raw_ts_part - mean(raw_ts_part, na.rm=TRUE))/sd(raw_ts_part, na.rm=TRUE)
                    
                    smap <- s_map(time_series, E=E[j], tau=tau[k], theta=theta[l], silent=TRUE)
                    matrix_rho <- cbind(matrix_rho, smap$rho)
                    
                    m <- m + step_size
                }
                
                matrix_results <- rbind(matrix_results, data.frame(
                w = window_size[i],
                E = E[j],
                tau = tau[k],
                theta = theta[l],
                rho = mean(matrix_rho)
                ))
                
            }
        }
    }
}

results <- matrix_results

### sort results
matrix_rho <- c()
for(i in 1:length(window_size))
{
    best_w <- order(-results[results$w == window_size[i],]$rho)[1]
    matrix_rho <- cbind(matrix_rho, results[results$w == window_size[i],]$rho[best_w])
}


### plot rho vs window size
save.plot.file <- F
if(save.plot.file){pdf("figure_s2a.pdf", width = 4, height = 3)}

par(mar=c(4,4.5,2,1)+0.1)

plot(window_size, matrix_rho, lty=1, lwd=2, type="l", col="blue", las=0, cex.axis=1.5, cex.lab=1.5, xlim=c(min(window_size), max(window_size)), main="", xlab="Window size", ylab=expression("Predictability, "~rho))

#dev.off()
