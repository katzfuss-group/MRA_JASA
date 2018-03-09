rm(list = ls())
setwd("S:/work/distributedComputing/code")
library(RandomFields)

reps=5
r=49
J=4
M.max=8
n=r*J^M.max
n1=sqrt(n)

# covariance function
cov.model <- RMexp(var=1, scale=.3) + RMnugget(var=.05) # exponential 


# simulate datasets and save them to disk
for(i in 1:reps) {
  set.seed(100+i)
  x.seq = y.seq = seq(0,1,length=n1)
  simu = as.matrix(RFsimulate(cov.model,x=x.seq,y=y.seq,grid=TRUE,spConform=FALSE))
  simresults=write.table(simu,paste("../results/2Dsim",i,".csv",sep=''),
                         row.names=FALSE, col.names=FALSE, sep=",")
}