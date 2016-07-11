library(RColorBrewer)

cols <- brewer.pal(9,"Paired")

setwd("~/Documents/Onderzoek/Stanford/2016_figurePower/HCP_group_effect_sizes/SingleSubject/tfMRI_EMOTION_LR_hp200_s4.feat/")

design <- read.table("design.mat",skip=5)
TR <- 0.72

#manually make regressiors
ev1 <- c(32.053,74.196,116.338)/TR*1000
ev1c <- rep(0,length=dim(design)[1]*1000+1)
ev1c[c(ev1[1]:(ev1[1]+18000),
       ev1[2]:(ev1[2]+18000),
       ev1[3]:length(ev1c))] <- 1

ev2 <- c(10.982,53.125,95.267)/TR*1000
ev2c <- rep(0,length=dim(design)[1]*1000+1)
ev2c[c(ev2[1]:(ev2[1]+18000),
       ev2[2]:(ev2[2]+18000),
       ev2[3]:(ev2[3]+18000))] <- 1

evt <- seq(from=0,to=dim(design)[1],by=0.001)


plot(design$V1,col=cols[2],type="l",lwd=2,ylim=c(-1,1),xlab="timepoint",ylab="design matrix value")
lines(evt,ev1c,col=cols[1],lwd=2)
lines(design$V3,col=cols[4],lwd=2)
lines(evt,ev2c,col=cols[3],lwd=2)
legend(0,-0.5,c("blocked design neutral","convolved design neutral",'blocked design fear','convolved design fear'),col=cols[1:4],lwd=2)
