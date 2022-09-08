####################################################################
# The R code used to obtain multi-speices Ricker models with various dimensions and critical transitions
# For simplicity, we only show the scenario that sequentially removes the species with the smallest grwoth rate (Fig. S16 a-e)
####################################################################
rm(list = ls())
# Change Working Directory
OPath <- 'D:\\data\\Cooperation\\Florian\\Code\\Code_test\\Code_demo\figure_s16_HighDimensionalModels'
setwd(OPath)
SaveFile <- F #T/F for saving files
#initialize random seed
seed.t1=1012
set.seed(seed.t1)
library(dplyr)

####moving window data
mvda=function(x,wd=12,move=1){
  da=NULL
  for(i in 1:(length(x)-wd+1)){
    da=cbind(da,x[i:(wd+i-1)])
  }
  dat=da
  if(move>1){
    ind.t=c(1+move*c(0:(ncol(da)%/%move)))
    ind.t=ind.t[ind.t<=ncol(da)]
    dat=da[,ind.t]
  }
  if((length(x)-wd)%%move>0){dat=cbind(dat,da[,ncol(da)])}
  return(dat)  
}

sp_no_i <-  97

sp_no <- sp_no_i
truncated_norm <- function(average, SD, thres)
{
  temp <- rnorm(1, mean=average, sd=SD)
  if (abs(temp) > thres) {
    temp <- 0
  }
  return(temp)
}

#define parameters for DE 
nv <- numeric(sp_no)  #true value
nv2 <- numeric(sp_no) #observation with error
initial_nv <- numeric(sp_no) #common initial values
r <- numeric(sp_no)

####code for generating growth rate r####
max_var_r <- 0.5

for(j in 1:sp_no) {
  r[j] = 1.5*(1.0 + runif(1, -max_var_r, max_var_r)) 
}

initial_nv <- numeric(sp_no) # initial values
for(j in 1: sp_no) initial_nv[j] <- 0.1*(1 + truncated_norm(average=0, SD=0.2, thres=0.9))

############generate interaction matrix##################
interaction_matrix <- matrix(0:0,nrow=sp_no, ncol=sp_no)
max_int<-0.5
sdd<-2.0

for(k in 1:(sp_no-1)) {
  for(j in (k+1):sp_no) {
    interaction_matrix[k, j] = truncated_norm(average=0, SD=sdd, thres=max_int)  #effect of sp. j on sp.k
  }
}

for(k in 1:sp_no) {
  interaction_matrix[k, k] = -1.0
}
interaction_matrix[1,3]
for(k in 2:sp_no) {
  for(j in 1:(k-1)) {
    if(interaction_matrix[j, k] < 0) {
      interaction_matrix[k, j] = interaction_matrix[j, k]*(1.0 + 0.9*runif(1, -1, 1))  #competition
    }
    else {
      interaction_matrix[k, j] = -1*interaction_matrix[j, k]*(1.0 + 0.9*runif(1, -1, 1)) #trophic interaction, sp.j is prey while sp.k is predator
    }
  }
}

Ricker_map <- function(invec, t, int_matrix, rr)
{
  growth_vec <- int_matrix%*%invec
  next_vec <- invec*exp(rr*(1.0 + growth_vec))
  return(next_vec)
}

#Function to calculate theoretical Jacobian matrix
Ricker_map_diff <- function(invec, t, int_matrix, rr)
{
  growth_vec <- int_matrix%*%invec
  temp_vec <- invec*exp(rr*(1.0 + growth_vec))
  
  #differential coefficient component (1) (all element)
  DF <- diag(unlist(c(rr*temp_vec)))%*%(int_matrix)
  
  #differential coefficient component (2) (diagnal only)
  for(k in 1: length(int_matrix[1,])) {
    DF[k, k] = DF[k, k] + exp(rr*(1.0 + growth_vec))[k,]
  }
  return(DF)
}

######################The first step analysis with the first interaction matrix and the first set of noise series##############
end_time <- 2000
start <- 1000
end <- 2000

max_time <- 2000
env<-numeric(max_time)
cor <- 0.8
max_red <- 0.5
env[1]<-0.1
for(i in 2:max_time) env[i] = cor*env[i-1] + (1-cor)*runif(1, -max_red, max_red)
plot(c(1:1000), env[1:1000], type="l")
max(env)

#Add noises on species growth
env_dep<-numeric(sp_no)
env_dep <- runif(sp_no, -1.0, 1.0)

r_env <- matrix(0, nrow=sp_no, ncol=end_time)
forcing_size <- 0.05*1
for(k in 1: sp_no) {
  for(t in 1: end_time) {
    r_env[k, t] <- r[k]*(1.0 + forcing_size*env_dep[k]*env[t])
  }
}

#Initial settings
time <- 0
nv <- initial_nv
result2<-t(rbind(time, as.data.frame(nv)))
for(time in 1:end_time)
{
  r_t=r_env[, time]
  nv <- Ricker_map(nv, time, interaction_matrix, r_t)
  result2_t<-t(rbind(time, as.data.frame(nv)))
  result2<-rbind(result2, result2_t)
}

colnames(result2) <- c('time',paste('V',1:sp_no,sep=''))
rownames(result2) <- NULL
result3 <- result2[start:end,]
ncol(result3)


########################################################
## Remove species with low abundance and stable dynamics
cri.m=10^-3; cri.sd=10^-2

dnsp=1
while(dnsp>0){
  # Delete sp without dynamics (1st round)
  exind=sort(union(which(is.na(result3[nrow(result3),-1])),
                   which(apply(result3[,-1]<cri.m,2,all)|(apply(result3[,-1],2,sd)<cri.sd))))
  sp_no_t=sp_no-length(exind);dnsp=sp_no-sp_no_t;(sp_no=sp_no_t)
  if(length(exind)==0){exind=10^8}
  if(sp_no_t<2){break}
  r=r[-exind];interaction_matrix=interaction_matrix[-exind,][,-exind];
  initial_nv=initial_nv[-exind];r_env=r_env[-exind,]
  #Initial settings
  time <- 0
  nv <- initial_nv
  result2<-t(rbind(time, as.data.frame(nv)))
  for(time in 1:end_time)
  {
    r_t=r_env[, time]
    nv <- Ricker_map(nv, time, interaction_matrix, r_t)
    result2_t<-t(rbind(time, as.data.frame(nv)))
    result2<-rbind(result2, result2_t)
  }
  
  if(sp_no_t>0 & all(!is.na(result2[nrow(result2),]))) {
    colnames(result2) <- c('time',paste('V',1:sp_no,sep=''))
    rownames(result2) <- NULL
    result3 <- result2[start:end,]
  }
}
print(paste(sp_no_i,sp_no,sep='->'))

(startsp=paste('RK_',sp_no,sep=''))

# Save a multi-species Ricker model with 21 species coexistence 
com.ls <- r.ls <- intM.ls <- list()
com.ls[[1]] <- result3
r.ls[[1]] <- as.matrix(r)
intM.ls[[1]] <- interaction_matrix

#####################################################################
# Species were sequentially removed depending on their growth rate, r
# 1: Removing sequentially the species of the smallest r 
# 2: Removing sequentially the species of the largest r 
# 3: Removing sequentially the species of the most extreme r (i.e., large deviations from the average); 
remove.seq <- 1 
# The number of species removed at each step 
rmsp.seq=c(5,5,5,3)

for(k in 2:5){
  result3 <- com.ls[[k-1]]
  r <- r.ls[[k-1]]
  interaction_matrix <- intM.ls[[k-1]]
  sp_no.o=ncol(interaction_matrix)
  rmsp=rmsp.seq[k-1] # The number of species removed from the community
  if(remove.seq==1){
    rv_sp=order(r)[1:rmsp]
  }else if(remove.seq==2){
    rv_sp=order(r,decreasing=T)[1:rmsp]
  }else if(remove.seq==3){
    rv_sp=order(abs(r-mean(r)),decreasing=T)[1:rmsp]
  }
  
  r.t=r[-rv_sp,];
  intM.t=interaction_matrix[-rv_sp,][,-rv_sp]
  initial_nv=result3[nrow(result3),-c(1,rv_sp+1)]
  sp_no=ncol(intM.t)
  
  colnames(intM.t)=rownames(intM.t)=NULL;intM.t=apply(intM.t,2,as.numeric)
  initial_nv=matrix(unlist(c(initial_nv)),ncol=1)
  
  #Initial settings
  time <- 0
  nv <- initial_nv
  result2<-c(time, nv)
  end_time=1000
  for(time in 1:end_time)
  {
    nv <- Ricker_map(invec=nv, t=time, int_matrix=intM.t, rr=r.t)
    result2_t<-t(rbind(time, as.data.frame(nv)))
    result2<-rbind(result2, result2_t)
  }
  
  colnames(result2) <- c('time',paste('V',1:sp_no,sep=''))
  rownames(result2) <- NULL
  
  com.ls[[k]] <- result2
  r.ls[[k]] <- as.matrix(r.t)
  intM.ls[[k]] <- intM.t
}

#########################################################################################
### Trigger the critical transition in multi-species Ricker model
tslength=500 # The length of time series data collected before reaching the tipping point
####################################################################
# Tuning the bifurcation parameters to trigger critical transitions
# Note that the tuning parameters for each community are different
initialF.seq=c(0.70,0.60,0.7,0.68,0.7)
terminalF.seq=c(1.3,1.45,1.51,1.52,1.5)
trend.seq=c(1.2,1.256,1.256,1.256,1.256)
####################################################################

for(k in 1:5){
  result3=com.ls[[k]]
  r=r.ls[[k]]
  interaction_matrix=intM.ls[[k]]
  nsp=nrow(interaction_matrix)
  (startsp=paste('RK',nsp,sep=''))
  
  initial_nv=result3[nrow(result3),-1]
  sp_no=length(initial_nv)
  
  colnames(interaction_matrix)=rownames(interaction_matrix)=NULL;interaction_matrix=apply(interaction_matrix,2,as.numeric)
  initial_nv=matrix(unlist(c(initial_nv)),ncol=1)
  r=unlist(c(r))
  
  #Initial settings
  time <- 0
  nv <- initial_nv
  result2<-t(rbind(time, as.data.frame(nv)))
  for(time in 1:1000)
  {
    r_t=r
    nv <- Ricker_map(nv, time, interaction_matrix, r_t)
    result2_t<-t(rbind(time, as.data.frame(nv)))
    result2<-rbind(result2, result2_t)
  }
  
  colnames(result2) <- c('time',paste('V',1:sp_no,sep=''))
  
  (nve=result2[nrow(result2),-1])
  nv <- nve
  
  initialF=initialF.seq[k]
  terminalF=terminalF.seq[k]
  equilibrium=T
  if(equilibrium){noise_on=1;sigma=0.01}else{noise_on=0}
  if(equilibrium){ # get equilibrium
    (nve=result2[nrow(result2),-1])
    nv <- nve
    result2<-t(rbind(time, as.data.frame(nv)))
    r_env=r*initialF
    for(time in 1:1000)
    {
      r_t=r_env
      nv <- Ricker_map(nv, time, interaction_matrix, r_t)
      result2_t<-t(rbind(time, as.data.frame(nv)))
      result2<-rbind(result2, result2_t)
    }
    colnames(result2) <- c('time',paste('V',1:sp_no,sep=''))
    (nve2=result2[nrow(result2),-1])
    nv <- nve2
  }
  
  result3<-NULL
  tim <- 1:1000
  wd.t=50
  xt=apply(mvda(tim,wd=wd.t),2,mean,na.rm=T)
  # Mootonic increase of bifurcation parameter 
  if(equilibrium){trend <- seq(1,terminalF,length.out=max(tim)) # from equilibrium}
  }else{trend <- seq(1,trend.seq[k],length.out=max(tim))} # from chaos}
  
  for(time in tim)
  {
    r_t=r*initialF*trend[time]
    L_noise_add <- if(noise_on!=0){sigma*nv*rnorm(1)}else{0}
    nv <- Ricker_map(nv, time, interaction_matrix, r_t)+L_noise_add
    result3_t<-t(rbind(time, as.data.frame(nv)))
    result3<-rbind(result3, result3_t)
  }
  
  #Function to calculate theoretical Jacobian matrix
  Jt=list()
  for(time in 1:nrow(result3)){
    inv=result3[time,-1]
    r_t=r*initialF*trend[time]
    Jt[[time]]=Ricker_map_diff(invec=inv, t=time, int_matrix=interaction_matrix, rr=r_t)
  }
  
  dev=lapply(Jt,function(x){
    if(any(is.na(x))){
      y=NA
    }else{
      a.e=eigen(x)$value
      ind.a1=which.max(abs(a.e))
      y=c(abs(a.e[ind.a1]),Re(a.e[ind.a1]),Im(a.e[ind.a1]))
    }
    return(y)
  })
  
  dev=matrix(unlist(dev),length(dev),ncol=3,byrow=T)
  colnames(dev)=c('abs(DEV)','Re(DEV)','Im(DEV)')
  dev[1,]
  
  win.graph(50,60);par(mfcol=c(3,1),mar=c(4,4,1,1))
  if(equilibrium){ind.t=1:1000}else{ind.t=1:792} #792
  
  xt2=apply(mvda(tim[ind.t],wd=wd.t),2,mean,na.rm=T)
  y=dev[ind.t,'abs(DEV)']
  ymv=apply(mvda(y,wd=wd.t),2,mean,na.rm=T)
  plot(y~tim[ind.t],type='l',xlab='Time',ylab='abs(DEV)')
  ind.tip=xt2[which(ymv>1)[1]]
  abline(h=1,col='grey');abline(v=ind.tip,col='red')
  lines(ymv~xt2,col='blue',lwd=1)
  abline(h=ncol(Jt[[1]]),col='grey')
  
  y=dev[ind.t,'Re(DEV)']
  ymv=apply(mvda(y,wd=wd.t),2,mean,na.rm=T)
  plot(y~tim[ind.t],type='l',xlab='Time',ylab='Re(DEV)')
  abline(h=1,col='grey');abline(v=ind.tip,col='red')
  lines(ymv~xt2,col='blue',lwd=1)
  
  y=dev[ind.t,'Im(DEV)']
  ymv=apply(mvda(y,wd=wd.t),2,mean,na.rm=T)
  plot(y~tim[ind.t],type='l',xlab='Time',ylab='Im(DEV)')
  abline(h=1,col='grey');abline(v=ind.tip,col='red')
  lines(ymv~xt2,col='blue',lwd=1)
  
  (ind.tip) # The change point (|J|=1)
  
  # graphical parameters
  if(k==1){
    win.graph(120,90);par(mfcol=c(5,5),mar=c(4,4,1,1))
  }else if(k==2){
    win.graph(120,90);par(mfcol=c(4,5),mar=c(4,4,1,1))
  }else if(k==3){
    win.graph(120,90);par(mfcol=c(3,4),mar=c(4,4,1,1))
  }else if(k==4){
    win.graph(90,90);par(mfcol=c(3,3),mar=c(4,4,1,1))
  }else if(k==5){
    win.graph(90,90);par(mfcol=c(2,2),mar=c(4,4,1,1))
  }
  
  for(i in 2:ncol(result3)){
    y=result3[,i]
    ymv=apply(mvda(y,wd=wd.t),2,mean,na.rm=T)
    plot(y~tim,type='l',xlab='Time',ylab=paste('Abundance of SP',i-1))
    abline(v=ind.tip,col='red')
    lines(ymv~xt,col='blue',lwd=1)
  }
  
  y=apply(result3[,-1],1,sum)
  ymv=apply(mvda(y,wd=wd.t),2,mean,na.rm=T)
  plot(y~tim,type='l',main='total abundance',xlab='Time',ylab='Total abundance')
  abline(v=ind.tip,col='red')
  lines(ymv~xt,col='blue',lwd=1)
  
  eqi <- if(equilibrium){'equi'}else{'chaos'}
  colnames(result3)=c('time',paste('SP',1:sp_no,sep='_'))
  if(SaveFile){write.csv(result3[(ceiling(ind.tip)-tslength+1):(ceiling(ind.tip)),],
            paste(tslength,'RemoveSeq',remove.seq,'_',startsp,'.csv',sep=''),row.names=F)}
}
