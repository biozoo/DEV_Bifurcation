####################################################################
# The R code used to calculate DEV from the time series of each population in multi-speices Ricker models showing critical transition
# The synthetic time series data was derived from the R file, "21sp_Ricker_CT_demo2.R" 
# For simplicity, we only show the scenario that sequentially removes the species with the smallest growth rate (Fig. S16 a-e)
####################################################################

rm(list = ls())
# Change Working Directory
OPath <- 'D:\\data\\Cooperation\\Florian\\Code\\Code_test\\Code_demo\figure_s16_HighDimensionalModels'
setwd(OPath)
SaveFile <- F #T/F for saving files
seed=68365
set.seed(seed)
library(deSolve)
library(dplyr)
library(rEDM)

###Load datasets with various dimensionality 
da=list()
da[[1]]=read.csv('500RemoveSeq1_RK21.csv',header=T)
da[[2]]=read.csv('500RemoveSeq1_RK16.csv',header=T)
da[[3]]=read.csv('500RemoveSeq1_RK11.csv',header=T)
da[[4]]=read.csv('500RemoveSeq1_RK6.csv',header=T)
da[[5]]=read.csv('500RemoveSeq1_RK3.csv',header=T)

# Number of species in each dataset
(nvb=unlist(lapply(da,ncol))-1)
result.dev=NULL
for(ls.i in 1:length(da)){
  nvb.t=nvb[ls.i]
  da.t=da[[ls.i]]
  rg.ts=da.t[,1]-da.t[1,1]+1
  da.t=da.t[,-1]
  print(ls.i)
  #vb.i=2
  for(vb.i in 1:nvb.t){
    raw_time_series <- da.t[,vb.i]
    ###set up the parameter ranges for performing cross-validation
    E <- seq(1,12,by=1)
    tau <- seq(1,1,by=1)
    theta <- 0
    window_size <- seq(25, 250, by=25)
    step_size <- 25
    
    ### Note that S-Map calculations for large window sizes and/or time series need a lot of computational time
    
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
    
    
    d.t=diff(c(matrix_rho));
    ind.w1=which((d.t<median(d.t))&c(matrix_rho[-length(window_size)]>0))
    if(length(ind.w1)>0){ind.w2=ind.w1[1]}else{ind.w2=which.max(c(matrix_rho))}
    # Select the optimal moving window
    opt.w=window_size[ind.w2]

    results.w=filter(results,w==opt.w)
    # Select the optimal parameter sets
    param.opt=results.w[which.max(results.w[,'rho']),]
    
    
    ################################################################################
    E <- as.numeric(param.opt['E'])
    tau <- as.numeric(param.opt['tau'])
    theta <- seq(0,2.5,by=0.5)
    window_size <- as.numeric(param.opt['w'])
    step_size <- 25
    
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
      param.opt['theta'] <- theta_best
      
      # calculate eigenvalues for best theta
      smap <- s_map(norm_rolling_window, E=E, tau=tau, theta=theta_best, silent=TRUE, save_smap_coefficients=TRUE)
      smap_co <- as.matrix(smap$smap_coefficients[[1]])
      colnames(smap_co)=rownames(smap_co)=NULL
      
      matrix_eigen <- matrix(NA, nrow = NROW(smap_co), ncol = 3)
      
      for(k in 1:NROW(smap_co))
      {
        if(!is.na(smap_co[k,1]))
        {
          M <- rbind(c(smap_co[k, 1:E]), cbind(diag(E - 1), rep(0, E - 1)))
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
    
    ### plot dev vs time and in complex plain
    result.t=unlist(c(nvb.t,vb.i,param.opt,matrix_result[nrow(matrix_result),-1],matrix_result[nrow(matrix_result),-1]-matrix_result[1,-1]))
    names(result.t)=c('SpeciesNumber','SpeciesID',names(param.opt),paste('Term',c('Abs','Re','Im'),sep='_'),paste('Diff',c('Abs','Re','Im'),sep='_'))
    result.dev=rbind(result.dev,result.t)
    
  }# end of vb.i
}# end of ls.i

rownames(result.dev)=NULL
if(SaveFile){write.csv(result.dev,'DEV_Ricker_n500_RemoveSeq1.csv',row.names=F)}

# Panel a & e: Removing sequentially the species of the smallest r
win.graph(45,90);par(mfcol=c(2,1),mar=c(4,4,2,1))
boxplot(Term_Abs~SpeciesNumber,result.dev,xlab='',ylab='|DEV|', main='', outline=T)
stripchart(Term_Abs~SpeciesNumber, vertical = TRUE, data = result.dev,method = "jitter", add = TRUE, pch = 20, col = 'blue')
abline(h=1.0)
boxplot(Diff_Abs~SpeciesNumber,result.dev,xlab='Number of species',ylab='delta |DEV|', outline=T)
stripchart(Diff_Abs~SpeciesNumber, vertical = TRUE, data = result.dev,method = "jitter", add = TRUE, pch = 20, col = 'blue')
abline(h=0)

