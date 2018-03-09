source('MRAhelpfunctions.R')

#######  partitioning et al for one-dimensional domains

partition.1D=function(J,M,r=ceiling(length(z)/J^M),domain,locs,z,predlocs.vec=FALSE,placement='regular') {
  
  ## function to return parent (numbered) index from (numbered) index
  par.ind=function(ind,J){
    full.j=indices[[ind]]
    par.j=full.j[-length(full.j)]
    num.ind(par.j,J)
  }
  
  # list of indices for all nodes in the tree structure
  indices=list(list())
  ind.list=list()
  if(M>0) {
    for(m in 1:M){
      ind.list=c(ind.list,list(1:J))
      indices=c(indices,as.list(data.frame(t(expand.grid(ind.list)[,seq(m,1,by=-1)]))))
    }
  }
  n.ind=length(indices)

  # create knot tree structure
  knots=vector("list",n.ind)
  data=vector("list",n.ind)
  bounds=matrix(ncol=2,nrow=n.ind)
  for(ind in 1:n.ind) {
    if(ind==1) bounds[ind,]=domain else {
      last.ind=indices[[ind]][length(indices[[ind]])]
      pb=bounds[par.ind(ind,J),]
      bounds[ind,]=seq(pb[1],pb[2],length=J+1)[last.ind+(0:1)]
    }
    if(length(indices[[ind]])<M) {
      if(placement=='regular') {
        knot.dist=diff(bounds[ind,])/r
        #knots[[ind]]=seq(bounds[ind,1]+knot.dist/2,bounds[ind,2]-knot.dist/2,length=r)
        knots[[ind]]=bounds[ind,1]+(1:r)/(r+1)*diff(bounds[ind,])
      } else {
        t=(r-1)/2
        knotloc=seq(1/t,1,length=t)^1.5*(t-1)/t
        knots[[ind]]=mean(bounds[ind,])+c(-knotloc[seq(t,1,-1)],0,knotloc)*diff(bounds[ind,])/2
      }
      data[[ind]]=NULL
    } else {
      ind.sub=which(locs>bounds[ind,1] & locs<=bounds[ind,2])
      knots[[ind]]=locs[ind.sub]
      data[[ind]]=z[ind.sub]
    }
  }
  
  # prediction locations
  if(typeof(predlocs.vec)=='logical') pred.locs=NA  else {
    pred.locs=vector("list",n.ind)
    for(ind in 1:n.ind) {
      if(length(indices[[ind]])>=M) {
        pred.locs[[ind]]=predlocs.vec[predlocs.vec>bounds[ind,1] & predlocs.vec<=bounds[ind,2]]
      } else pred.locs[[ind]]=knots[[ind]]
    }
  }

  return(list(indices=indices,knots=knots,data=data,bounds=bounds,pred.locs=pred.locs))
  
}




###############################   pre-C algorithm  #########################

MRA=function(theta,cov.fun,data,knots,indices,pred.locs=NULL,predres=0:length(indices[[length(indices)]]),var.eps=0) {
  
  pred=(!is.null(pred.locs)) # do prediction if locs are given (o/w return likelihood)
  co=function(locs1,locs2) {cov.fun(locs1,locs2,theta)}  
  
  ## extract dimensions and other constants
  n.ind=length(indices)
  M=length(indices[[n.ind]])
  J=if(n.ind==1) 1 else indices[[n.ind]][length(indices[[n.ind]])]
  indres=vector("list",M+1)
  for(m in 0:M) indres[[m+1]]=if(m==0) 1 else (indres[[m]][length(indres[[m]])]+1):(indres[[m]][length(indres[[m]])]+J^m)
  
  ## initialize
  Kc.B=vector("list",n.ind)
  R.prior.chol=vector("list",n.ind)
  if(pred) postmean=postvar=B.tilde=preds=vector("list",n.ind)  
  A.tilde.cur=w.tilde.cur=A.tilde.prev=w.tilde.prev=vector("list",n.ind)
  loglik=numeric(length=n.ind)
  
  ## going from coarsest to finest resolution
  for(m in 0:M){  for(ind in indres[[m+1]]){
    
    # indices of parent nodes
    par=rep(1,times=m)
    if(m>1) for(k in 1:(m-1)) par[k+1]=num.ind(indices[[ind]][1:k],J)  
    
    # create prior quantities etc
    pred.locs.j=if(pred) pred.locs[[ind]] else NULL
    temp = create.prior(co,knots[c(par,ind)],R.prior.chol[par],Kc.B[par],data[[ind]],pred.locs.j,var.eps)
    # knots.b=knots[c(par,ind)]; R.prior.chol.b=R.prior.chol[par]; Kc.B.b=Kc.B[par]; data.j=data[[ind]]
    R.prior.chol[[ind]]=temp$R.prior.chol.j; Kc.B[[ind]]=temp$Kc.B.j
    if(m==M){
      A.tilde.prev[[ind]]=temp$A.tilde.j; w.tilde.prev[[ind]]=temp$w.tilde.j
      if(pred) {
        postmean[[ind]]=temp$postmean.j; postvar[[ind]]=temp$postvar.j; B.tilde[[ind]]=temp$B.tilde.j
      } else  loglik[[ind]]=temp$loglik.j
    }
    
  }}
  
  rm(Kc.B)
  
  
  ## posterior inference (from finest to coarsest resolution)
  if(pred) R.post.chol=Kc.A=Kc.w=vector("list",n.ind)
  children=numeric(length=J)
  if(M>0) { for(m in seq(M-1,0,by=-1)){  
    
    A.tilde.cur=w.tilde.cur=vector("list",n.ind)
    
    for(ind in indres[[m+1]]){
      
      # indices of the J child nodes
      for(j in 1:J) children[j]=num.ind(as.numeric(c(indices[[ind]],j)),J)
      
      # calculate posterior quantities
      temp=post.inf(R.prior.chol[[ind]],w.tilde.prev[children],A.tilde.prev[children],pred)
      # R.prior.chol.j=R.prior.chol[[ind]];w.tilde.children=w.tilde.prev[children];A.tilde.children=A.tilde.prev[children]
      w.tilde.cur[[ind]]=temp$w.tilde.cur.j
      A.tilde.cur[[ind]]=temp$A.tilde.cur.j
      if(pred){
        R.post.chol[[ind]]=temp$R.post.chol.j; Kc.w[[ind]]=temp$Kc.w.j; Kc.A[[ind]]=temp$Kc.A.j
      } else loglik[[ind]]=temp$loglik.j
      
    }
    
    A.tilde.prev=A.tilde.cur; w.tilde.prev=w.tilde.cur
    
  }}
  
  
  ## spatial prediction
  if(pred) {
    for(ind in indres[[M+1]]){ # only for finest resolution
      if(M>0){
        
        # indices of parent nodes
        par=rep(1,times=M); if(M>1) { for(k in 1:(M-1)) par[k+1]=num.ind(indices[[ind]][1:k],J)  }
        
        # obtain predictive means and variances 
        preds[[ind]]=pred.fct(postmean[[ind]],postvar[[ind]],B.tilde[[ind]],R.post.chol[par],Kc.A[par],Kc.w[par],predres=predres)
        
      } else preds[[ind]]=cbind(postmean[[ind]],postvar[[ind]])
    }  
  }  
  
  return( if(pred) preds else sum(loglik))
  
}







###########################################

library(lava) # needed for blockdiag to create cov.mat for illustration

MRA.illus=function(theta,cov.fun,data,knots,indices,M.cov.plot=length(indices[[length(indices)]])) {
  
  pred=FALSE
  
  ## extract dimensions and other constants
  n.ind=length(indices)
  M=length(indices[[n.ind]])
  J=if(n.ind==1) 1 else indices[[n.ind]][length(indices[[n.ind]])]
  indres=vector("list",M+1)
  for(m in 0:M) indres[[m+1]]=if(m==0) 1 else (indres[[m]][length(indres[[m]])]+1):(indres[[m]][length(indres[[m]])]+J^m)
  
  ## create prior quantities
  V.prior=vector("list",n.ind)
  B=vector("list",n.ind)
  if(pred) {V.p=vector("list",n.ind); B.p=vector("list",n.ind); V.op=vector("list",n.ind); L=vector("list",n.ind)}
  R.prior=vector("list",n.ind)
  R.prior.chol=vector("list",n.ind)
  for(ind in 1:n.ind) {
    inds=indices[[ind]] # full (j) index
    m=length(inds)
    V.prior[[ind]]=vector("list",m+1)
    B[[ind]]=vector("list",m+1)
    if(pred) {V.p[[ind]]=vector("list",m+1); B.p[[ind]]=vector("list",m+1);
              V.op[[ind]]=vector("list",m+1); L[[ind]]=vector("list",m+1)}
    for(l in 0:m){  
      V.prior[[ind]][[l+1]]=vector("list",m+1)
      if(pred) {V.p[[ind]][[l+1]]=vector("list",m+1); V.op[[ind]][[l+1]]=vector("list",m+1)}
      ind.lm1=if(l<2) 1 else num.ind(inds[1:(l-1)],J)
      for(k in l:m){
        ind.k=if(k==0) 1 else num.ind(inds[1:k],J)
        V.prior[[ind]][[l+1]][[k+1]]= if(l==0) cov.fun(knots[[ind]],knots[[ind.k]],theta) else
          woodbury(V.prior[[ind]][[l]][[k+1]],V.prior[[ind]][[l]][[l]],R.prior.chol[[ind.lm1]],
                   V.prior[[ind.k]][[l]][[l]]) 
        if(pred){
          V.p[[ind]][[l+1]][[k+1]]= if(l==0) cov.fun(pred.locs[[ind]],pred.locs[[ind.k]],theta) else
            woodbury(V.p[[ind]][[l]][[k+1]],V.p[[ind]][[l]][[l]],R.prior.chol[[ind.lm1]],
                     V.p[[ind.k]][[l]][[l]])          
          V.op[[ind]][[l+1]][[k+1]]= if(l==0) cov.fun(knots[[ind]],pred.locs[[ind.k]],theta) else
            woodbury(V.op[[ind]][[l]][[k+1]],V.prior[[ind]][[l]][[l]],R.prior.chol[[ind.lm1]],
                     V.p[[ind.k]][[l]][[l]])                    
        }
      }
      if(m==M) {
        B[[ind]][[l+1]]=V.prior[[ind]][[l+1]][[l+1]]
        if(pred) {
          B.p[[ind]][[l+1]]=V.p[[ind]][[l+1]][[l+1]]
          L[[ind]][[l+1]]=V.op[[ind]][[l+1]][[l+1]]
        }
      }
    }
    R.prior[[ind]]=V.prior[[ind]][[m+1]][[m+1]]
    R.prior.chol[[ind]]=cholesky(R.prior[[ind]])
  }
  
  
  basis.functions=vector("list",M)
  K=vector("list",M)
  finest.subs=indices[indres[[M+1]]]
  temp=do.call(blockdiag,R.prior[indres[[M+1]]])
  cov.mat=if(M.cov.plot==M) temp else matrix(0,ncol=ncol(temp),nrow=nrow(temp))
  if(M>0) {  for(res in 0:(M-1)) {
    K[[res+1]]=do.call(blockdiag,lapply(R.prior[indres[[res+1]]],solve))
    bf=vector("list",length(indres[[res+1]]))
    for(sub in 1:length(indres[[res+1]])) {
      tree.ind=indices[[indres[[res+1]][sub]]]
      ind.children=sapply(finest.subs,function(x) all(x[1:res]==tree.ind))
      num.ind.children=sapply(finest.subs[ind.children],num.ind,J=J)
      bf[[sub]]=do.call(rbind,sapply(B[num.ind.children],`[`,res+1))
    }
    basis.functions[[res+1]]=do.call(blockdiag,bf)
    if(res<=M.cov.plot) cov.mat = cov.mat + basis.functions[[res+1]]%*%K[[res+1]]%*%t(basis.functions[[res+1]])    
  } }

      
  return(list(basis.functions,cov.mat,R.prior))
  
}






###########   old (unused) functions   ###############





###############################   NEW, more efficient algorithm  #########################

MRA.fast=function(theta,cov.fun,data,knots,indices,pred.locs=NULL) {
  
  pred=(!is.null(pred.locs)) # do prediction if locs are given (o/w return likelihood)
  
  ## extract dimensions and other constants
  n.ind=length(indices)
  M=length(indices[[n.ind]])
  J=if(n.ind==1) 1 else indices[[n.ind]][length(indices[[n.ind]])]
  indres=vector("list",M+1)
  for(m in 0:M) indres[[m+1]]=if(m==0) 1 else (indres[[m]][length(indres[[m]])]+1):(indres[[m]][length(indres[[m]])]+J^m)
  
  ## initialize
  Kc.B=vector("list",n.ind)
  if(pred) postmean=postvar=B.tilde=preds=vector("list",n.ind)
  R.prior.chol=vector("list",n.ind)
  A.tilde.cur=w.tilde.cur=A.tilde.prev=w.tilde.prev=vector("list",n.ind)
  loglik.j=numeric(length=n.ind)
  
  ## going from coarsest to finest resolution
  for(m in 0:M){  for(ind in indres[[m+1]]){
    
    inds=indices[[ind]] # full (j) index 
    
    # create prior quantities
    V=Kc.B[[ind]]=vector("list",m+1)    
    for(l in 0:m){  
      V[[l+1]]=vector("list",m+1)
      ind.l=if(l==0) 1 else num.ind(inds[1:l],J)
      for(k in l:m){
        ind.k=if(k==0) 1 else num.ind(inds[1:k],J)
        V[[l+1]][[k+1]]= if(l==0) cov.fun(knots[[ind]],knots[[ind.k]],theta) else
          V[[l]][[k+1]]-tp(Kc.B[[ind]][[l]],Kc.B[[ind.k]][[l]])
      }
      if(l<m) {
        Kc.B[[ind]][[l+1]]=sol(R.prior.chol[[ind.l]],t(V[[l+1]][[l+1]]))
      } else {
        R.prior=V[[m+1]][[m+1]]
        R.prior.chol[[ind]]=cholesky(R.prior)
      }
    }
    
    # begin inference for regions at finest resolution M
    if(m==M) { 
      
      # pre-compute solves
      Sic.B=vector("list",M+1)
      for(l in 0:M) Sic.B[[l+1]]=sol(R.prior.chol[[ind]],V[[l+1]][[l+1]]) 
      Sic.y=sol(R.prior.chol[[ind]],data[[ind]])
      
      # inference quantities
      w.tilde.prev[[ind]]=lapply(Sic.B,function(x) tp(x,Sic.y))
      A.tilde.prev[[ind]]=vector("list",m+1)
      for(l in 0:m) {
        A.tilde.prev[[ind]][[l+1]]=vector("list",m+1)
        for(k in l:m)  A.tilde.prev[[ind]][[l+1]][[k+1]] = tp(Sic.B[[l+1]],Sic.B[[k+1]])      
      }    
      
      # quantities for prediction or likelihood evaluation
      if(pred) {
        
        # calculate B.p and L
        V.p=Kc.Bp=V.pp=vector("list",m+1)    
        for(l in 0:M){  
          V.p[[l+1]]=vector("list",m+1)
          ind.l=if(l==0) 1 else num.ind(inds[1:l],J)
          for(k in l:M){
            ind.k=if(k==0) 1 else num.ind(inds[1:k],J)
            V.p[[l+1]][[k+1]]= if(l==0) cov.fun(pred.locs[[ind]],knots[[ind.k]],theta) else
              V.p[[l]][[k+1]]-tp(Kc.Bp[[l]],Kc.B[[ind.k]][[l]])
          }
          Kc.Bp[[l+1]]=sol(R.prior.chol[[ind.l]],t(V.p[[l+1]][[l+1]])) # Sic.L for l=M
          V.pp[[l+1]]=if(l==0) cov.fun(pred.locs[[ind]],pred.locs[[ind]],theta) else
            V.pp[[l]]-tp(Kc.Bp[[l]],Kc.Bp[[l]])
        }
        
        # initialize prediction inference
        postmean[[ind]]=postvar[[ind]]=matrix(nrow=length(pred.locs[[ind]]),ncol=M+1)
        postmean[[ind]][,M+1]=tp(Kc.Bp[[M+1]],Sic.y)
        postvar[[ind]][,M+1]=diag(V.pp[[M+1]]-tp(Kc.Bp[[M+1]],Kc.Bp[[M+1]]))
        if(M>0) {
          B.tilde[[ind]]=vector("list",M+1); B.tilde[[ind]][[M+1]]=vector("list",M)
          for(k in 0:(M-1)) B.tilde[[ind]][[M+1]][[k+1]]=V.p[[k+1]][[k+1]]-tp(Kc.Bp[[M+1]],Sic.B[[k+1]])          
        }
        
      } else {
        
        loglik.j[ind] = log.det(R.prior.chol[[ind]]) + tp(Sic.y,Sic.y)
        
      }      
      
    }
    
  }}
  
  rm(V,Kc.B)
  if(pred) rm(V.p,V.pp,Kc.Bp)
  
  
  ## posterior inference (from finest to coarsest resolution)
  R.post.chol=Kc.A=Kc.w=w.mm=vector("list",n.ind)
  children=numeric(length=J)
  if(M>0) { for(m in seq(M-1,0,by=-1)){  
    
    A.tilde.cur=w.tilde.cur=vector("list",n.ind)
    
    for(ind in indres[[m+1]]){
      
      inds=indices[[ind]] # full (j) index
      
      # sum up over children tildes
      for(j in 1:J) children[j]=num.ind(as.numeric(c(inds,j)),J)            
      w=vector("list",m+1)
      A=vector("list",m+1)    
      for(l in 0:m){
        w[[l+1]]=Reduce('+',sapply(w.tilde.prev[children],`[`,l+1))
        A[[l+1]]=vector("list",m+1)
        for(k in l:m)  A[[l+1]][[k+1]]=Reduce('+',l.ex(A.tilde.prev[children],l+1,k+1))
      }
      
      # calculate cholesky of K.inv and save relevant w
      R.post = R.prior.chol[[ind]]%*%t(R.prior.chol[[ind]]) + A[[m+1]][[m+1]]
      R.post.chol[[ind]]=cholesky(R.post)
      w.mm[[ind]]=w[[m+1]]
      
      # pre-compute the solves required later
      Kc.w[[ind]]=sol(R.post.chol[[ind]],w.mm[[ind]])
      Kc.A[[ind]]=lapply(sapply(A,`[`,m+1),function(x) sol(R.post.chol[[ind]],t(x)))
      
      if(m>0) {
        # calculate w.tilde and A.tilde
        w.tilde.cur[[ind]]=mapply('-',w,lapply(Kc.A[[ind]],function(x) tp(x,Kc.w[[ind]])),SIMPLIFY=FALSE)
        A.tilde.cur[[ind]]=vector("list",m+1)
        for(l in 0:m) {
          A.tilde.cur[[ind]][[l+1]]=vector("list",m+1)
          for(k in l:m) A.tilde.cur[[ind]][[l+1]][[k+1]] = A[[l+1]][[k+1]] - tp(Kc.A[[ind]][[l+1]],Kc.A[[ind]][[k+1]])  
        }          
      }
      
      # likelihood evaluation    
      if(!pred) loglik.j[ind] = log.det(R.post.chol[[ind]]) - log.det(R.prior.chol[[ind]]) - tp(Kc.w[[ind]],Kc.w[[ind]])
      
    }
    
    A.tilde.prev=A.tilde.cur
    w.tilde.prev=w.tilde.cur
    
  }}
  
  
  # spatial prediction
  if(pred) {
    for(ind in indres[[M+1]]){ # only for finest resolution
      if(M>0){
        if(M>1) { for(k in seq(M-1,1,by=-1)){
          ind.k=num.ind(indices[[ind]][1:k],J)
          Kc.Btilde=sol(R.post.chol[[ind.k]],t(B.tilde[[ind]][[k+2]][[k+1]]))
          B.tilde[[ind]][[k+1]]=vector("list",k)        
          for(l in seq(k-1,0,by=-1)) B.tilde[[ind]][[k+1]][[l+1]]=B.tilde[[ind]][[k+2]][[l+1]]-tp(Kc.Btilde,Kc.A[[ind.k]][[l+1]])
        } }
        for(k in 0:(M-1)) {
          ind.k=if(k==0) 1 else num.ind(indices[[ind]][1:k],J)
          Kc.Btilde.cur=sol(R.post.chol[[ind.k]],t(B.tilde[[ind]][[k+2]][[k+1]]))
          postmean[[ind]][,k+1]=tp(Kc.Btilde.cur,Kc.w[[ind.k]])
          postvar[[ind]][,k+1]=diag(tp(Kc.Btilde.cur,Kc.Btilde.cur)) # colSums(Kc.Btilde.cur^2)                  
        }      
      }
      preds[[ind]]=cbind(rowSums(postmean[[ind]]),rowSums(postvar[[ind]]))
    }  
  }
  
  return( if(pred) preds else sum(loglik.j))
  
}





## calculate A-B*inv(C)*D' using cholesky of C
woodbury=function(A,B,C.chol,D) {
  A -   if(length(B)<=1 | length(C.chol)==0 | length(D)<=1)  0 else 
    t(base::forwardsolve(C.chol,t(B)))%*%base::forwardsolve(C.chol,t(D))
}

## calculate B'*inv(C)*D using cholesky of C
quadform=function(B,C.chol,D) {
  if(length(B)<=1 | length(C.chol)==0 | length(D)<=1)  0 else #matrix(0,nrow=ncol(B),ncol=ncol(D)) else
    t(base::forwardsolve(C.chol,B))%*%forwardsolve(C.chol,D)
}