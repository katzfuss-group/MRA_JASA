##### smaller help functions

## calculate log determinant from cholesky factor of a matrix
log.det=function(mat.chol) 2*sum(log(diag(mat.chol)))

## return numbered index from tree index (i.e., inverse to indices)
num.ind=function(tree.ind,J){
  if(length(tree.ind)==0) num.index=1 else {
    l=seq(length(tree.ind)-1,0,by=-1)
    num.index=sum(J^l)+sum((tree.ind-1)*J^l)+1 #tree.ind[length(tree.ind)]
  }
  return(num.index)
}

# ## return last element in a vector
# last=function(x) { return( x[length(x)] ) }

## extract (k,l)th element of each element in a nested list
l.ex=function(list,k,l) sapply(sapply(list,`[`,k),'[',l)
# temp=list(list(1,list(2.1,2.2)),list(3,list(4.1,4.2)),list(5,list(6.1,6.2))); temp; l.ex(temp[2:3],2,1)

## cholesky decomposition
cholesky=function(A) if(length(A)<2 && (length(A)<=1 || A==0)) 0 else 
  t(tryCatch(chol(A),error=function(e) chol(A+1e-2*mean(diag(A)))*diag(nrow(A))) )


## calculate inv(C.chol)*A, where C.chol is cholesky factor
sol=function(C.chol,A)  if(length(A)<=1 | length(C.chol)<=1) 0 else base::forwardsolve(C.chol,A)

## calculate t(A)*B
tp=function(A,B)  if(length(A)<=1 | length(B)<=1) 0 else t(A)%*%B

## calculate X'*X for all (relevant) combinations
outer.tp=function(X) {
  m=length(X)-1
  outer.list=vector("list",m+1)
  for(l in 0:m) {
    outer.list[[l+1]]=vector("list",m+1)
    for(k in l:m)  outer.list[[l+1]][[k+1]] = tp(X[[l+1]],X[[k+1]])      
  }   
  return(outer.list)
}


### durbin-levinson likelihood
dllik=function(r, z) {
  EPS <- .Machine$double.eps
  n <- length(z)
  if (n != length(r)) 
    stop("arguments have unequal length")
  error <- numeric(n)
  sigmasq <- numeric(n)
  error[1] <- z[1]
  sigmasq[1] <- r[1]
  phi <- r[2]/r[1]
  error[2] <- z[2] - phi * z[1]
  sigmasqkm1 <- r[1] * (1 - phi^2)
  sigmasq[2] <- sigmasqkm1
  logg <- log(r[1]) + log(sigmasqkm1)
  for (k in 2:(n - 1)) {
    if (sigmasqkm1 < 0 || abs(sigmasqkm1) < EPS) 
      stop("r is not a p.d. sequence")
    phikk <- (r[k + 1] - phi %*% rev(r[2:k]))/sigmasqkm1
    sigmasqk <- sigmasqkm1 * (1 - phikk^2)
    phinew <- phi - phikk * rev(phi)
    phi <- c(phinew, phikk)
    sigmasqkm1 <- sigmasqk
    logg <- logg + log(sigmasqk)
    error[k + 1] <- z[k + 1] - crossprod(phi, rev(z[1:k]))
    sigmasq[k + 1] <- sigmasqk
  }
  S <- sum((error * error)/sigmasq)
  LogL <- S + logg
  LogL
}



#######  create prior (going from coarse to fine resolution)
create.prior=function(co,knots.b,R.prior.chol.b,Kc.B.b,data.j=NULL,pred.locs.j=NULL,var.eps=0) {
  
  pred=(!is.null(pred.locs.j)) # do prediction if locs are given (o/w return likelihood)
  m=length(knots.b)-1
  
  # create prior quantities
  Kc.B.c=vector("list",m); V=vector("list",m+1)
  for(l in 0:m){  
    V[[l+1]]=co(knots.b[[m+1]],knots.b[[l+1]])
    if(l>0) { for(k in 0:(l-1)){
      V[[l+1]] = if(l<m) V[[l+1]]-tp(Kc.B.c[[k+1]],Kc.B.b[[l+1]][[k+1]])  else  V[[l+1]]-tp(Kc.B.c[[k+1]],Kc.B.c[[k+1]])
    }   }
    if(l<m) Kc.B.c[[l+1]]=sol(R.prior.chol.b[[l+1]],t(V[[l+1]])) else 
      R.prior.chol.j = cholesky(V[[l+1]]) + if(is.null(data.j)) 0 else var.eps*diag(nrow(V[[l+1]]))
  }  
  Kc.B.b[[m+1]]=vector("list",m); Kc.B.b[[m+1]]=Kc.B.c
  
  
  # begin inference for regions at finest resolution M
  if(!is.null(data.j)) { 
    
    # pre-compute solves
    Sic.B=vector("list",m+1)
    for(l in 0:m) Sic.B[[l+1]]=sol(R.prior.chol.j,V[[l+1]]) 
    Sic.y=sol(R.prior.chol.j,data.j)
    
    # inference quantities
    w.tilde.j=lapply(Sic.B,function(x) tp(x,Sic.y))
    A.tilde.j=outer.tp(Sic.B)
    
    # quantities for prediction or likelihood evaluation
    if(pred) {
      
      R.prior.chol.b=c(R.prior.chol.b,list(R.prior.chol.j))
      
      # calculate B.p and L
      V.p=Kc.Bp=V.pp=vector("list",m+1)   
      for(l in 0:m){  
        V.p[[l+1]]=vector("list",m+1)
        for(k in l:m){
          V.p[[l+1]][[k+1]]= if(l==0) co(pred.locs.j,knots.b[[k+1]]) else
            V.p[[l]][[k+1]]-tp(Kc.Bp[[l]],Kc.B.b[[k+1]][[l]])
        }
        Kc.Bp[[l+1]]=sol(R.prior.chol.b[[l+1]],t(V.p[[l+1]][[l+1]])) # Sic.L for l=M
        V.pp[[l+1]]=if(l==0) co(pred.locs.j,pred.locs.j) else
          V.pp[[l]]-tp(Kc.Bp[[l]],Kc.Bp[[l]])
      }
      
      # initialize prediction inference
      postmean.j=postvar.j=matrix(nrow=length(pred.locs.j),ncol=m+1)
      postmean.j[,m+1]=tp(Kc.Bp[[m+1]],Sic.y)
      postvar.j[,m+1]=diag(V.pp[[m+1]]-tp(Kc.Bp[[m+1]],Kc.Bp[[m+1]]))
      B.tilde.j=vector("list",m+1)
      if(m>0) {
        B.tilde.j[[m+1]]=vector("list",m)
        for(k in 0:(m-1)) B.tilde.j[[m+1]][[k+1]]=V.p[[k+1]][[k+1]]-tp(Kc.Bp[[m+1]],Sic.B[[k+1]])          
      }
      
    } else {
      
      loglik.j = log.det(R.prior.chol.j) + tp(Sic.y,Sic.y)
      
    }      
    
  }
  
  return.main=list(R.prior.chol.j=R.prior.chol.j,Kc.B.j=Kc.B.b[[m+1]])
  if(is.null(data.j)) return(return.main) else {
    return.main=c(return.main,list(w.tilde.j=w.tilde.j,A.tilde.j=A.tilde.j))
    if(pred) return(c(return.main,list(postmean.j=postmean.j,postvar.j=postvar.j,B.tilde.j=B.tilde.j))) else 
      return(c(return.main,loglik.j=loglik.j))
  }
  
}



#######  posterior inference (when going from fine to coarse resolution)

post.inf=function(R.prior.chol.j,w.tilde.children,A.tilde.children,pred=FALSE) {
  
  m=length(w.tilde.children[[1]])-2
  
  # sum up over children tildes                    
  w=vector("list",m+1)
  A=vector("list",m+1)    
  for(l in 0:m){
    w[[l+1]]=Reduce('+',sapply(w.tilde.children,`[`,l+1))
    A[[l+1]]=vector("list",m+1)
    for(k in l:m)  A[[l+1]][[k+1]]=Reduce('+',l.ex(A.tilde.children,l+1,k+1))
  }
  
  # calculate cholesky of K.inv and save relevant w
  R.post = R.prior.chol.j%*%t(R.prior.chol.j) + A[[m+1]][[m+1]]
  R.post.chol.j = cholesky(R.post)
  
  # pre-compute the solves required later
  Kc.w.j=sol(R.post.chol.j,w[[m+1]])
  Kc.A.j=lapply(sapply(A,`[`,m+1),function(x) sol(R.post.chol.j,t(x)))
  
  # calculate w.tilde and A.tilde
  if(m>0) {    
    w.tilde.cur.j=mapply('-',w,lapply(Kc.A.j,function(x) tp(x,Kc.w.j)),SIMPLIFY=FALSE)
    A.tilde.cur.j=vector("list",m+1)
    for(l in 0:m) {
      A.tilde.cur.j[[l+1]]=vector("list",m+1)
      for(k in l:m) A.tilde.cur.j[[l+1]][[k+1]] = A[[l+1]][[k+1]] - tp(Kc.A.j[[l+1]],Kc.A.j[[k+1]])  
    }          
  }
  
  return.main=if(m>0) list(w.tilde.cur.j=w.tilde.cur.j,A.tilde.cur.j=A.tilde.cur.j) else list()
  if(pred) {
    return(c(return.main,list(R.post.chol.j=R.post.chol.j,Kc.w.j=Kc.w.j,Kc.A.j=Kc.A.j)))
  } else return(c(return.main,list(loglik.j=log.det(R.post.chol.j)-log.det(R.prior.chol.j)-tp(Kc.w.j,Kc.w.j))))
  
}


#######  prediction at finest resolution at the end
pred.fct=function(postmean.j,postvar.j,B.tilde.j,R.post.chol.b,Kc.A.b,Kc.w.b,predres=0:(ncol(postmean.j)-1)) {
  
  M=ncol(postmean.j)-1
  if(M>1){ for(k in seq(M-1,1,by=-1)){
    Kc.Btilde=sol(R.post.chol.b[[k+1]],t(B.tilde.j[[k+2]][[k+1]]))
    B.tilde.j[[k+1]]=vector("list",k)        
    for(l in seq(k-1,0,by=-1)) B.tilde.j[[k+1]][[l+1]]=B.tilde.j[[k+2]][[l+1]]-tp(Kc.Btilde,Kc.A.b[[k+1]][[l+1]])
  }}
  for(k in 0:(M-1)) {
    Kc.Btilde.cur=sol(R.post.chol.b[[k+1]],t(B.tilde.j[[k+2]][[k+1]]))
    postmean.j[,k+1]=tp(Kc.Btilde.cur,Kc.w.b[[k+1]])
    postvar.j[,k+1]=colSums(as.matrix(Kc.Btilde.cur)^2)  # diag(tp(Kc.Btilde.cur,Kc.Btilde.cur))
  }      
  preds.j=cbind(rowSums(as.matrix(postmean.j[,predres+1])),rowSums(as.matrix(postvar.j[,predres+1])))
  
  return(preds.j)
  
}