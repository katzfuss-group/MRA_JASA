rm(list = ls())
setwd("S:/work/distributedComputing/code")
source('MRAfunctions.R')
library(ltsa) # for fast simulation of truth
library(fields) # needed for rdist function


### covariance function
theta=NA
cov.fun=function(locs1,locs2,theta=NA) {
  if( (length(locs1)==0) || (length(locs2)==0) ){ 0 } else { 
    x=rdist(locs1,locs2)/.05
    exp(-x) #+ (abs(x)<1e-8)*.05
    # (1 + sqrt(3)*x)*exp(-sqrt(3)*x)  # Matern with smoothness 1.5
  }
}


### set parameters
J0=3; M0=3; r0=2
nper=r0
n.all=nper*J0^M0
domain=c(-1e-5,1)

# simulate data
set.seed(10)
reg.length=diff(domain)/J0^M0
locs.all=c()
for(i in 1:J0^M0) locs.all=c(locs.all,domain[1]+runif(r0,(i-1)*reg.length,i*reg.length))
cov.mat.true=cov.fun(c(domain[1],locs.all),c(domain[1],locs.all))
cov.chol=chol(cov.mat.true)
set.seed(10)
y.all=(cov.chol%*%rnorm(n.all+1))[2:(n.all+1)]
pchs.all=c(2:(M0+1),16)
cols.all=c('blue','orange','green3','red')
ltys.all=1:(M0+1)

### 1-RAs for varying data size
M=1; r=r0*M0;
cols=cols.all[c(1,M0+1)]; pchs=pchs.all[c(1,M0+1)]; ltys=ltys.all[c(1,M0+1)]
n=n.all
J=n/r
locs=locs.all
z=y.all
  
# calculate MRA quantities
part.illus=partition.1D(J,M,r,domain,locs,z)
cov.mats=vector("list",M+1)
for(m in 0:M){
  cov.mats[[m+1]]=MRA.illus(NA,cov.fun,part.illus$data,part.illus$knots,part.illus$indices,M.cov.plot=m)[[2]]
}
u=10; xgrid=(1:(u*n.all))/(n.all*u); part.illus=partition.1D(J,M,r,domain,xgrid,1:(n.all*u))
bfs=MRA.illus(NA,cov.fun,part.illus$data,part.illus$knots,part.illus$indices)[[1]]

# plot covariance function and approximations
distmat=rdist(locs,locs)
samp=1:n^2 #sample(1:n^2,5000)
  
#pdf(file='../plots/illustration2_cov.pdf',width=6,height=3)
par(mar=c(2,3.8,.3,.3)+0.1)  # bottom, left, top and right 
curve(cov.fun(0,x),0,.15,lwd=2,xlim=c(0,.15),ylim=0:1,xlab='distance',ylab='covariance')  
for(m in 0:M){
  points(distmat[samp],cov.mats[[m+1]][samp],col=cols[m+1],pch=pchs[m+1],cex=.7+.3*m)
}
legend('topright',legend=c('truth',paste('m=',0:M,sep='')),col=c('black',cols),lty=c(1,rep(NA,M+1)),pch=c(NA,pchs),
       pt.cex=c(1,rep(1.2,M+1)),lwd=c(2,rep(NA,M+1)),cex=1.2)    
dev.off()
  
#pdf(file='../plots/illustration2_bfs_data.pdf',width=6,height=4)
par(mar=c(2,3.8,.3,.5)+0.1)   # bottom, left, top and right  
m=1; yb=c(c(min(y.all)-1.5,min(y.all))-.1,max(y.all)+.1); mp=(yb[1:(m+1)]+yb[2:(m+2)])/2
plot(locs,z,ylab='knots & y',ylim=c(yb[1],yb[m+2]),xaxs="i",yaxs='i',
     xlim=c(-.005,1.005),cex=.7,pch=pchs[M+1],yaxt='n',col=cols[M+1])
axis(side=2,at=ceiling(min(y.all)):floor(max(y.all)))
axis(side=2,at=mp[1:m],labels=expression(Q^{(0)}))
segments(c(0,part.illus$bounds[,2]),rep(yb[2]+.05,J+1),y1=rep(yb[3]-.05,J+1),col=cols[M+1],lty=ltys[M+1])
#points(part.illus$knots[[1]],rep(mp[1],times=r),cex=3,pch=pchs[1],col=cols[1])
for(j in 1:r) lines(xgrid,yb[1]+bfs[[1]][,j]/max(as.numeric(bfs[[1]]))*1.3,col=cols[1])
segments(x0=0:1,y0=yb[1]+.05,x1=0:1,y1=yb[2]-.05,col=cols[1],lwd=2,lty=ltys[1])
abline(h=yb[m+1])
dev.off()




### M-RA with M=M0 for increasing m
part.illus=partition.1D(J0,M0,r0,domain,locs.all,y.all)
cov.mats=vector("list",M0+1)
for(m in 0:M0){
  cov.mats[[m+1]]=MRA.illus(NA,cov.fun,part.illus$data,part.illus$knots,part.illus$indices,M.cov.plot=m)[[2]]
}
distmat=rdist(locs.all,locs.all)
samp=1:n.all^2 #sample(1:n^2,5000)
widths=c(2,2,2,1)
indres=vector("list",M0+1) 
for(i in 0:M0) indres[[i+1]]=if(i==0) 1 else (indres[[i]][length(indres[[i]])]+1):(indres[[i]][length(indres[[i]])]+J0^i)  
yb=c(c(min(y.all)-seq(M0,1,by=-1)*1.5,min(y.all))-.1,max(y.all)+.1); mp=(yb[1:(M0+1)]+yb[2:(M0+2)])/2
u=10; xgrid=(1:(u*n.all))/(n.all*u); part.illus=partition.1D(J0,M0,r0,domain,xgrid,1:(n.all*u))
bfs=MRA.illus(NA,cov.fun,part.illus$data,part.illus$knots,part.illus$indices)[[1]]

M=M0 #for(M in -1:M0) {
  # pdf(file=paste('../plots/illustration2_MRA_cov',M+1,'.pdf',sep=''),width=6,height=3)
  par(mar=c(2,3.8,.3,.3)+0.1)  # bottom, left, top and right 
  curve(cov.fun(0,x),0,.15,lwd=2,xlim=c(0,.15),ylim=0:1,xlab='distance',ylab='covariance')
  if(M>=0){ for(m in 0:M){
    points(distmat[samp],cov.mats[[m+1]][samp],col=cols.all[m+1],pch=pchs.all[m+1],cex=.7)
  } }
  legend('topright',legend=c('truth',paste('m=',0:M0,sep='')),col=c('black',cols.all),
         lty=c(1,rep(NA,M0+1)),pch=c(NA,pchs.all),pt.cex=c(1,rep(1.2,M0+1)),
         lwd=c(2,rep(NA,M0+1)),cex=1.2)  
  dev.off()
#}

for(M in -1:M0) {
#  pdf(file=paste('../plots/illustration2_MRA_bfs_data',M+1,'.pdf',sep=''),width=6,height=4)
  par(mar=c(2,3.8,.3,.5)+0.1)   # bottom, left, top and right  
  datacol= if(M<M0) 'black' else cols.all[M0+1]
  plot(locs.all,y.all,ylab='knots & y',ylim=c(yb[1],yb[m+2]),xaxs="i",yaxs='i',
       xlim=c(-.005,1.005),cex=.7,pch=pchs.all[M0+1],yaxt='n',col=datacol)
  axis(side=2,at=ceiling(min(y.all)):floor(max(y.all)))
  abline(h=yb[M0+1])
  if(M>=0){
    axis(side=2,at=mp[1:M0],labels=c(expression(Q^{(0)}),expression(Q^{(1)}),expression(Q^{(2)})))
    abline(h=yb,col='black')
    for(i in 0:M){
      segments(c(0,part.illus$bounds[indres[[i+1]],2]),rep(yb[i+1]+.05,J0^(i+1)+1),
               y1=rep(yb[i+2]-.05,J0^(i+1)+1),col=cols.all[i+1],lty=ltys.all[i+1],lwd=widths[i+1])
      if(i<M0){ 
#        points(unlist(part.illus$knots[indres[[i+1]]]),rep(mp[i+1],times=r0*length(indres[[i+1]])),
#             cex=M0-.5-i/2,pch=pchs.all[i+1],col=cols.all[i+1])   
        for(j in 1:ncol(bfs[[i+1]])) { lines(xgrid,
          yb[i+1]+bfs[[i+1]][,j]/max(as.numeric(bfs[[i+1]]))*1.3,
          lwd=1,col=cols.all[i+1])
        }
      }
    }  
  }
  dev.off()
}



### prediction
test.locs.ind=c(31:36) #,143:n.all)
lastfirstobs=c(locs.all[test.locs.ind[1]-1],locs.all[test.locs.ind[length(test.locs.ind)]+1])
predlocs.vec=seq(0,1,length=100)
#predlocs.vec=seq(lastfirstobs[1],lastfirstobs[2],length=30)
#test.locs.ind=c(); 
locs.ind=setdiff(1:n.all,test.locs.ind)
z=y.all[locs.ind]  # +rnorm(n)*.1
locs=locs.all[locs.ind]
JMr=rbind(c(J0,M0,r0),c(J0^M0,1,r0*M0),c(1,0,1))
preds.all=vector("list",nrow(JMr))
for(model in 1:nrow(JMr)) {  
  J=JMr[model,1]; M=JMr[model,2]; rj=JMr[model,3]
  part.pred=partition.1D(J,M,rj,domain,locs,z,predlocs.vec)
  preds=MRA.old(NA,cov.fun,part.pred$data,part.pred$knots,part.pred$indices,part.pred$pred.locs)
  post.mean= post.sd = numeric(length=0)
  for(ind in 1:length(part.pred$indices)) { if(length(part.pred$indices[[ind]])==M) {
    post.mean=c(post.mean,preds[[ind]][,1])
    post.sd=c(post.sd,ifelse(preds[[ind]][,2]>0,sqrt(preds[[ind]][,2]),0))
  }}  
  preds.all[[model]]=cbind(post.mean,post.sd)  
}
alpha=.05
q=qnorm(1-alpha/2)
range.p=range(c(z,preds.all[[1]][,1]-q*preds.all[[1]][,2],preds.all[[1]][,1]+q*preds.all[[1]][,2],
                preds.all[[2]][,1]-q*preds.all[[2]][,2],preds.all[[2]][,1]+q*preds.all[[2]][,2]))

#pdf(file='../plots/predillus.pdf',width=6,height=3)
par(mar=c(2,1.8,1,1.5)+0.1) # bottom, left, top and right  
plot(locs,z,cex=.5,xlab='',xlim=0:1,ylab='y',ylim=range.p,xaxs="i")
cols=c(2,4,3)
lwds=c(1.5,1.5,3)
for(model in seq(length(preds.all),1,-1)) {
  lines(predlocs.vec,preds.all[[model]][,1],col=cols[model],lwd=lwds[model])
  lines(predlocs.vec,preds.all[[model]][,1]-q*preds.all[[model]][,2],col=cols[model],lty=2,lwd=lwds[model])
  lines(predlocs.vec,preds.all[[model]][,1]+q*preds.all[[model]][,2],col=cols[model],lty=2,lwd=lwds[model])
}
abline(v=lastfirstobs,col='grey')
legend('bottomright',legend=c('data','0-RA','1-RA','3-RA'),col=c(1,3,4,2),lty=c(NA,1,1,1),pch=c(1,NA,NA,NA),
       pt.cex=c(.7,1,1,1),lwd=c(NA,2,1.5,1.5))
#dev.off()



# ### predictions for individual resolutions
# test.locs.ind=c(31:36) #,143:n.all)
# lastfirstobs=c(locs.all[test.locs.ind[1]-1],locs.all[test.locs.ind[length(test.locs.ind)]+1])
# predlocs.vec=seq(0,1,length=100)
# locs.ind=setdiff(1:n.all,test.locs.ind)
# z=y.all[locs.ind]  # +rnorm(n)*.1
# locs=locs.all[locs.ind]
# J=J0; M=M0; rj=r0
# part.pred=partition.1D(J,M,rj,domain,locs,z,predlocs.vec)
# predres.mean=matrix(nrow=length(predlocs.vec),ncol=M+1)
# predres.sd=matrix(nrow=length(predlocs.vec),ncol=M+1)
# for(res in 0:M) {  
#   preds=MRA(NA,cov.fun,part.pred$data,part.pred$knots,part.pred$indices,part.pred$pred.locs,predres=res)
#   post.mean= post.sd = numeric(length=0)
#   for(ind in 1:length(part.pred$indices)) { if(length(part.pred$indices[[ind]])==M) {
#     post.mean=c(post.mean,preds[[ind]][,1])
#     post.sd=c(post.sd,ifelse(preds[[ind]][,2]>0,sqrt(preds[[ind]][,2]),0))
#   }}  
#   predres.mean[,res+1]=post.mean; predres.sd[,res+1]=post.sd
# }
# 
# cols=c(3:5)
# plot(locs,z,cex=.5,xlab='',xlim=0:1,ylab='y',ylim=range.p,xaxs="i")
# lines(predlocs.vec,rowSums(predres.mean),col=2)
# for(m in M:M){ #0:M){
#   lines(predlocs.vec,predres.mean[,m+1],type='l',col=m+3)  
# }
# legend('bottomright',legend=c('data','0-RA','1-RA','3-RA'),col=c(1,3,4,2),lty=c(NA,1,1,1),pch=c(1,NA,NA,NA),
#        pt.cex=c(.7,1,1,1),lwd=c(NA,2,1.5,1.5))



### data plot for illustration of distributed computation
J=2; M=2; r=5
nper=r
n.all=nper*J^M
domain=c(-1e-5,1)
set.seed(10)
reg.length=diff(domain)/J^M
locs.all=c()
for(i in 1:J^M) locs.all=c(locs.all,domain[1]+runif(r,(i-1)*reg.length,i*reg.length))
cov.mat.true=cov.fun(c(domain[1],locs.all),c(domain[1],locs.all))
cov.chol=chol(cov.mat.true)
set.seed(100)
y.all=(cov.chol%*%rnorm(n.all+1))[2:(n.all+1)]
indres=vector("list",M+1) 
for(i in 0:M) indres[[i+1]]=if(i==0) 1 else 
  (indres[[i]][length(indres[[i]])]+1):(indres[[i]][length(indres[[i]])]+J^i)  
part.illus=partition.1D(J,M,r,domain,locs.all,y.all)


# #pdf(file='../plots/distributed_data.pdf',width=9,height=2.5)
# par(mar=c(2,3.8,.2,3.8)+0.1) # bottom, left, top and right  
# yb=c(c(min(y.all)-seq(M,1,by=-1),min(y.all))-.1,max(y.all)+.1); mp=(yb[1:(M+1)]+yb[2:(M+2)])/2
# plot(locs.all,y.all,ylab='knots & y',ylim=c(yb[1],yb[M+2]),xaxs="i",yaxs='i',
#      xlim=c(-.005,1.005),cex=.7,pch=16,yaxt='n',col='blue')
# axis(side=2,at=seq(ceiling(min(y.all*2))/2,floor(max(y.all)),by=.5))
# axis(side=2,at=mp[1:M],labels=c(expression(Q^{(0)}),expression(Q^{(1)})))
# segments(x0=0:1,y0=yb[1]+.05,x1=0:1,y1=yb[2]-.05,col=2,lwd=2)
# for(i in 0:(M-1)){
#   segments(c(0,part.illus$bounds[indres[[i+2]],2]),rep(yb[i+2]+.05,J^(i+2)+1),
#            y1=rep(yb[i+3]-.05,J^(i+2)+1),col=i+3,lty=i+1,lwd=2)
#   points(unlist(part.illus$knots[indres[[i+1]]]),rep(mp[i+1],times=r*length(indres[[i+1]])),
#          cex=M+1-i,pch=i+1,col=i+2)
# }  
# abline(h=yb,col='black'); abline(h=yb[M+1])
# dev.off()


cov.fun=function(locs1,locs2,theta=NA) {
  if( (length(locs1)==0) || (length(locs2)==0) ){ 0 } else { exp(-(rdist(locs1,locs2)/.05)) }
}
J0=J; M0=m=M; r0=r
part.illus=partition.1D(J0,M0,r0,domain,locs.all,y.all)
widths=c(2,2,2,1,1)
indres=vector("list",M0+1) 
for(i in 0:M0) indres[[i+1]]=if(i==0) 1 else (indres[[i]][length(indres[[i]])]+1):(indres[[i]][length(indres[[i]])]+J0^i)  
yb=c(c(min(y.all)-seq(M0,1,by=-1)*1.5,min(y.all))-.1,max(y.all)+.1); mp=(yb[1:(M0+1)]+yb[2:(M0+2)])/2
u=11; xgrid=(1:(u*n.all))/(n.all*u); part.illus=partition.1D(J0,M0,r0,domain,xgrid,1:(n.all*u))
bfs=MRA.illus(NA,cov.fun,part.illus$data,part.illus$knots,part.illus$indices)[[1]]
pchs.all=c(2:(M0+1),16)
cols.all=c('red','green3','blue')
ltys.all=1:(M0+1)

#pdf(file='../plots/distributed_data.pdf',width=9,height=2.5)
par(mar=c(2,3.8,.2,3.8)+0.1) # bottom, left, top and right  
plot(locs.all,y.all,ylab='knots & y',ylim=c(yb[1],yb[m+2]),xaxs="i",yaxs='i',
     xlim=c(-.005,1.005),cex=.7,pch=pchs.all[M0+1],yaxt='n',col='blue')
axis(side=2,at=ceiling(min(y.all)):floor(max(y.all)))
abline(h=yb[M0+1])
segments(c(0,part.illus$bounds[indres[[M0+1]],2]),rep(yb[M0+1]+.05,J0^(M0-1)+1),
         y1=rep(yb[M0+2]-.05,J0^(M0-1)+1),col='blue',lty=ltys.all[M0+1],lwd=widths[M0+1])
if(M>=0){
  axis(side=2,at=mp[1:M0],labels=c(expression(Q^{(0)}),expression(Q^{(1)})))
  abline(h=yb,col='black')
  for(i in 0:(M-1)){
    segments(c(0,part.illus$bounds[indres[[i+1]],2]),rep(yb[i+1]+.05,J0^(i+1)+1),
             y1=rep(yb[i+2]-.05,J0^(i+1)+1),col=cols.all[i+1],lty=ltys.all[i+1],lwd=widths[i+1])
    if(i<M0){ 
      #        points(unlist(part.illus$knots[indres[[i+1]]]),rep(mp[i+1],times=r0*length(indres[[i+1]])),
      #             cex=M0-.5-i/2,pch=pchs.all[i+1],col=cols.all[i+1])   
      for(j in 1:ncol(bfs[[i+1]])) lines(
        xgrid,yb[i+1]+bfs[[i+1]][,j]/max(as.numeric(bfs[[i+1]]))*1.3,lwd=1,col=cols.all[i+1])
    }
  }  
  abline(h=yb)
}
dev.off()



# #####################   only for talk   ##############################################
# 
# ### partitioning details for the M-RA
# indres=vector("list",M0+1) 
# for(i in 0:M0) indres[[i+1]]=if(i==0) 1 else (indres[[i]][length(indres[[i]])]+1):(indres[[i]][length(indres[[i]])]+J0^i)  
# part.illus=partition.1D(J0,M0,r,domain,locs.all,y.all)
# mi=min(y.all)-.2
# for(j in 0:(M0-1)){
#   pdf(file=paste('../plots/MRA_partitioning',j+1,'.pdf',sep=''),width=9,height=4)
#   par(mar=c(3,3.8,1,1.5)+0.1)   # bottom, left, top and right    
#   plot(locs.all,y.all,ylab='y',xlab='',ylim=c(mi+.2,max(y.all)),xaxs="i",xlim=c(-.005,1.005))
#   for(i in 0:j){
#     segments(c(0,part.illus$bounds[indres[[i+2]],2]),rep(mi+i+1,J0^(i+2)+1),y1=rep(mi+i+2,J0^(i+2)+1),col=i+2,lty=i+1,lwd=M0)
#     points(unlist(part.illus$knots[indres[[i+1]]]),rep(mi+.5+i,times=r*length(indres[[i+1]])),cex=M0-.5-i,pch=i+1,col=i+2)
#   }
#   dev.off()  
# }


##########################################################################################




###################   basis functions for nonstationary covariance   ######################

matern=function (u, phi, kappa) {
  uphi <- u/phi
  mat.cor=ifelse(u==0, 1,
                 2^(-(kappa - 1)) / gamma(kappa) * (uphi^kappa) * besselK(x = uphi, nu = kappa) )
  return(ifelse(u>600*phi,0,mat.cor))
}
cov.fun.ns=function(locs1,locs2,theta) {
  if( (length(locs1)==0) || (length(locs2)==0) ){ matrix(nrow=length(locs1),ncol=length(locs2)) } else {
    smooth1=.7*locs1+.2 # 3*(locs1-.5)^2+.2
    smooth2=.7*locs2+.2 # 3*(locs2-.5)^2+.2
    smooth.mean=(matrix(smooth1,nrow=length(smooth1),ncol=length(smooth2))+
                   matrix(smooth2,nrow=length(smooth1),ncol=length(smooth2),byrow=TRUE))/2
    matern(rdist(locs1,locs2),phi=.15,kappa=smooth.mean)
  }
}
### set parameters
J0=2; M0=4; r0=3
nper=40; n.all=nper*J0^M0
locs.all=((1:n.all)-.5)/n.all
domain=c(-1e-10,1)
cov.mat.true=cov.fun.ns(c(locs.all,1+1/n.all),c(locs.all,1+1/n.all))
cov.chol=chol(cov.mat.true)
set.seed(100)
y=(cov.chol%*%rnorm(n.all+1))[2:(n.all+1)]
test.locs.ind=c() # 396:n.all # 1:1000 #1:1000 #451:650
locs.ind=setdiff(1:n.all,test.locs.ind)
z=y[locs.ind]  # +rnorm(n)*.1
locs=locs.all[locs.ind]
part.illus=partition.1D(J0,M0,r0,domain,locs.all,y)
bfs=MRA.illus(NA,cov.fun.ns,part.illus$data,part.illus$knots,part.illus$indices)[[1]]
bounds=part.illus$bounds
indres=vector("list",M0+1)
for(m in 0:M0) indres[[m+1]]=if(m==0) 1 else (indres[[m]][length(indres[[m]])]+1):(indres[[m]][length(indres[[m]])]+J0^m)
M.plot=M0-1
range=0:1 #range(do.call(rbind,cov.mats))
cols=rep(1:6,5)

#pdf(file='../plots/nonstat_data.pdf',width=9,height=3)
par(mar=c(2,1.8,1,1)+0.1) # bottom, left, top and right  
plot(locs,z,ylab='y',xlab='s',xlim=0:1,ylim=range(z))
dev.off()

#png(file='../plots/nonstat_cov.png',width=480,height=400)
par(mar=c(1,1,1,0.6)+0.1) # bottom, left, top and right  
image.plot(apply(apply(apply(cov.mat.true,1, rev),1,rev),1,rev),
           legend.shrink=1,xlab='',ylab='',axes=FALSE,zlim=range,horizontal=FALSE)
dev.off()

#pdf(file='../plots/nonstat_bf_r3.pdf',width=6,height=3)
par(mar=c(2,4,1,1.5)+0.1) # bottom, left, top and right  
plot(-5,xlim=domain+c(-2e-3,5e-4),ylim=c(0-.3,M.plot-.3),xlab='',ylab='resolution',xaxs="i",yaxs="i",yaxt = "n")
axis(2, at = 0:(M.plot-1),labels=0:(M.plot-1))
for(res in 1:M.plot) {
#   basisvec=as.numeric(bfs[[res]])
#   basisvec[is.na(basisvec)]=0
#   lines(rep(locs.all,r0*J0^(res-1)),basisvec/max(basisvec)*.9+res-1.3,lwd=.01,col=rep(1:(r0*J0^(res-1)),each=length(locs.all)))
  for(j in 1:(r0*J0^(res-1))) lines(locs.all,bfs[[res]][,j]/max(as.numeric(bfs[[res]]))*.9+res-1.3,lwd=.01,col=cols[j])
  boundsplot=unique(as.numeric(bounds[indres[[res]],]))
  for(i in 1:length(boundsplot)) lines(rep(boundsplot[i],2),c(res-1.3,res-.5),lty=2,lwd=2)
}
abline(h=(1:(M.plot-1))-.3)
dev.off()




################    plot results from 1-D simulation study    ######################

library(ggplot2)
simresults=read.csv("../results/simscale_2mil_means.csv",header=FALSE)
names(simresults)=c('n','model','time','likincr','likfixed')
simresults[simresults=='NaN']=NA
simresults$model=factor(simresults$model,levels=c('full','1-RA S','1-RA F','M-RA S','M-RA F'),
                        labels=c('0-RA','1-RA S','1-RA F','M-RA','M-RA 8'))
simresults$imputed=factor(ifelse(is.na(simresults$time),'estimated','observed'),
                             levels=c('observed','estimated'))
simresults=simresults[-which(is.na(simresults$time) & simresults$model=='M-RA 8'),]
for(i in 1:nrow(simresults)) {
  if(is.na(simresults$time[i])) simresults$time[i]=simresults$time[i-1]*4^3
}
simresults=simresults[simresults$n>2000,]
simresults=simresults[-which(simresults$model=='M-RA 8'),]

xbreaks=c(5*1e3,1e4,2*1e4,5*1e4,1e5,2*1e5,5*1e5,1e6,2e6)
xlabs=c('5k','10k','20k','50k','100k','200k','500k','1M','2M')
timebreaks=c(1,20,60,600,60^2,60^2*24,60^2*24*100)
timelabs=c('1s','20s','1min','10min','1h','1day','100days')

#pdf(file='../plots/simtimes_2mil.pdf',width=6,height=3)
ggplot(data = simresults, aes(x = n, y = time, color = model)) + theme_bw() +       
  geom_line(size=1,aes(group=model,linetype=model)) + geom_point(size=3,aes(shape=imputed)) + 
  scale_x_log10(breaks=xbreaks,labels=xlabs) + 
  theme(legend.title=element_blank(),legend.key.size=unit(1.04,'cm')) +
  scale_y_log10(breaks=timebreaks,labels=timelabs)
dev.off()

#pdf(file='../plots/simlikincr_2mil.pdf',width=6,height=3)
ggplot(data = simresults, aes(x = n, y = likincr, color = model)) + theme_bw() +       
  geom_line(size=1,aes(group=model,linetype=model)) + geom_point(size=3) + 
  ylab("relative loglikelihood") + 
  theme(legend.title=element_blank(),legend.key.size=unit(1.2,'cm')) + 
  scale_x_log10(breaks=xbreaks,labels=xlabs) + scale_y_continuous(breaks=seq(.9,1.05,by=.05))
dev.off()

#pdf(file='../plots/simlikfixed_2mil.pdf',width=6,height=3)
ggplot(data = simresults, aes(x = n, y = likfixed, color = model)) + theme_bw() +       
  geom_line(size=1,aes(group=model,linetype=model)) + geom_point(size=3) + 
  ylab("relative loglikelihood") + 
  theme(legend.title=element_blank(),legend.key.size=unit(1.2,'cm')) + 
  scale_x_log10(breaks=xbreaks,labels=xlabs)
dev.off()


##########    simulation for fixed n   #########
library(ggplot2)
repnums=3
filename='../results/sim1D_2mill_theta05_'
dims=dim(read.csv(paste(filename,1,".csv",sep=''),header=FALSE))
simresults.all=array(dim=c(dims[1],dims[2]+1,repnums))
for(repnum in 1:repnums) {
  simresults.all[,1:dims[2],repnum]=as.matrix(read.csv(paste(
    filename,repnum,".csv",sep=''),header=FALSE))
  simresults.all[,dims[2]+1,repnum]=simresults.all[,dims[2],repnum]/simresults.all[4,dims[2],repnum]
}
simresults.all[simresults.all=='NaN']=NA
simresults=as.data.frame(apply(simresults.all, c(1,2), mean,na.rm=TRUE))
names(simresults)=c('M','r','time','loglikelihood','likscaled')
simresults$model=factor(ifelse(simresults$M >1,'M-RA','1-RA'),
                          levels=c('M-RA','1-RA'))
simresults$model[1]='M-RA'
simresults$timemin=simresults$time/60
simresults=simresults[-c(1,5,nrow(simresults)),]
ybreaks=round(simresults$likscaled[simresults$model=='M-RA'],3)

# pdf(file='../plots/simfixedn_2mil_logtime.pdf',width=6,height=3)
timebreaks=c(60*3,60*5,60*10,60*20,60*30,60^2,2*60^2)
timelabs=c('3','5','10','20','30','60','120')
ggplot(data = simresults, aes(x = time, y = likscaled, color = model)) + theme_bw() +       
  ylab('relative loglikelihood') + xlab('time (min)') +
  geom_line(size=1,aes(group=model,linetype=model)) + geom_point(size=3) + 
  scale_x_log10(breaks=timebreaks,labels=timelabs) + 
  theme(legend.title=element_blank(),legend.key.size=unit(1.2,'cm')) +
  scale_y_continuous(breaks=c(.7,.8,.9,ybreaks))
dev.off()

# pdf(file='../plots/simfixedn_2mil_xy.pdf',width=6,height=3)
timebreaks=c(60*3,60*5,60*10,60*20,60*30,60^2,2*60^2)
timelabs=c('3','5','10','20','30','60','120')
ggplot(data = simresults, aes(x = likscaled, y = time, color = model)) + theme_bw() +       
  xlab('relative loglikelihood') + ylab('time (min)') +
  geom_line(size=1,aes(group=model,linetype=model)) + geom_point(size=3) + 
  scale_y_log10(breaks=timebreaks,labels=timelabs) + 
  theme(legend.title=element_blank(),legend.key.size=unit(1.2,'cm')) +
  scale_x_continuous(breaks=c(.7,.8,.9,ybreaks))
dev.off()

# pdf(file='../plots/simfixedn_2mil.pdf',width=6,height=3)
ggplot(data = simresults, aes(x = timemin, y = likscaled, color = model)) + theme_bw() +       
  ylab('relative loglikelihood') + xlab('time (min)') +
  geom_line(size=1,aes(group=model,linetype=model)) + geom_point(size=3) + 
  theme(legend.title=element_blank(),legend.key.size=unit(1.2,'cm')) + 
  scale_y_continuous(breaks=c(.7,.8,.9,ybreaks)) +
  scale_x_continuous(breaks=seq(0,110,by=10))
dev.off()



################    plot results from 2-D simulation study    ######################

library(ggplot2)
repnums=1
filename='../results/sim2D_withnug_rep'
dims=dim(read.csv(paste(filename,1,".csv",sep=''),header=FALSE))
simresults.all=array(dim=c(dims[1],dims[2]+1,repnums))
for(repnum in 1:repnums) {
  simresults.all[,1:dims[2],repnum]=as.matrix(read.csv(paste(
    filename,repnum,".csv",sep=''),header=FALSE))
  simresults.all[,dims[2]+1,repnum]=(simresults.all[,dims[2],repnum]/simresults.all[3,dims[2],repnum])^
    sign(simresults.all[3,dims[2],repnum])
}
simresults.all[simresults.all=='NaN']=NA
simresults=as.data.frame(apply(simresults.all, c(1,2), mean,na.rm=TRUE))
names(simresults)=c('modnum','M','r','time','loglikelihood','likscaled')
simresults$model=factor(simresults$modnum,levels=1:2,labels=c('M-RA','1-RA'))
simresults$timemin=simresults$time/60

# pdf(file='../plots/sim2D_nug.pdf',width=6,height=3)
timebreaks=c(60*2,60*5,60*10,60*30,60^2,2*60^2)
timelabs=c('2','5','10','30','60','120')
ggplot(data = simresults, aes(x = time, y = likscaled, color = model)) + theme_bw() +       
  ylab('relative loglikelihood') + xlab('time (min)') +
  geom_line(size=1,aes(group=model,linetype=model)) + geom_point(size=3) + 
  scale_x_log10(breaks=timebreaks,labels=timelabs) + 
  theme(legend.title=element_blank(),legend.key.size=unit(1.2,'cm'))
dev.off()


#####################################
##########    TPW results   #########

library(maps)

### read in the observations
data=as.data.frame(read.csv("../data/MIRSmra.csv",header=FALSE))
names(data)=c('longitude','latitude','observations')
data$longitude=data$longitude-180
borders <- map_data("world")

### determine color ranges for plots
range.pred=range(data$observations); range.sd=c()
for(model in 1:3){
  preds=as.data.frame(read.csv(paste('../results/TPW_preds2_model',model,'.csv',sep=''),header=FALSE))
  range.pred=range(c(range.pred,preds[,3]))
  range.sd=range(c(range.sd,preds[,4]))
}

### plot the data
png(file='../plots/TPWdata.png',width=560,height=480)
ggplot() + theme_bw() + 
  theme(legend.title=element_blank(),legend.key.width=unit(2.8,'cm'),legend.position='bottom') +
  geom_point(data=data,aes(x=longitude,y=latitude,colour=observations),size=.5,shape=16) + 
  scale_colour_gradientn(colours = terrain.colors(101),limits=range.pred) +
  geom_path(data = borders, aes(x = long, y = lat, group = group)) +
  scale_x_continuous(limits=range(data$longitude),expand=c(0,0)) + 
  scale_y_continuous(limits=range(data$latitude),expand=c(0,0))
dev.off()


### plot predictions and st.dev.s for all methods
model=3
  
  preds=as.data.frame(read.csv(paste('../results/TPW_preds2_model',model,'.csv',sep=''),header=FALSE))
  names(preds)=c('longitude','latitude','mean','sd')
  preds$longitude=preds$longitude-180
  
  png(file=paste('../plots/TPW2pred',model,'.png',sep=''),width=560,height=480)
  ggplot() + theme_bw() + 
    theme(legend.title=element_blank(),legend.key.width=unit(2.8,'cm'),legend.position='bottom') +
    geom_point(data=preds,aes(x=longitude,y=latitude,colour=mean),size=3,shape=15) + 
    scale_colour_gradientn(colours = terrain.colors(101),limits=range.pred) +
    geom_path(data = borders, aes(x = long, y = lat, group = group)) +
    scale_x_continuous(limits=range(data$longitude),expand=c(0,0)) + 
    scale_y_continuous(limits=range(data$latitude),expand=c(0,0)) 
  dev.off()
  
  png(file=paste('../plots/TPW2sd',model,'.png',sep=''),width=560,height=480)
  ggplot() + theme_bw() + 
    theme(legend.title=element_blank(),legend.key.width=unit(2.8,'cm'),legend.position='bottom') +
    geom_point(data=preds,aes(x=longitude,y=latitude,colour=sd),size=3,shape=15) + 
    scale_colour_gradientn(colours = heat.colors(101),limits=range.sd) +
    geom_path(data = borders, aes(x = long, y = lat, group = group)) +
    scale_x_continuous(limits=range(data$longitude),expand=c(0,0)) + 
    scale_y_continuous(limits=range(data$latitude),expand=c(0,0)) 
  dev.off()
  
  png(file=paste('../plots/TPW2sd',model,'_rev.png',sep=''),width=560,height=480)
  ggplot() + theme_bw() + 
    theme(legend.title=element_blank(),legend.key.width=unit(2.8,'cm'),legend.position='bottom') +
    geom_point(data=preds,aes(x=longitude,y=latitude,colour=sd),size=3,shape=15) + 
    scale_colour_gradientn(colours = rev(heat.colors(101)),limits=range.sd) +
    geom_path(data = borders, aes(x = long, y = lat, group = group)) +
    scale_x_continuous(limits=range(data$longitude),expand=c(0,0)) + 
    scale_y_continuous(limits=range(data$latitude),expand=c(0,0)) 
  dev.off()
  

  