##### help functions

## ceiling function (for easy switching between julia versions)
function ceilfun(x)
  iceil(x)  # ceil(Int64,x)
end

## floor function (for easy switching between julia versions)
function floorfun(x)
  ifloor(x)  # floor(Int64,x)
end

## rounding function (for easy switching between julia versions)
function roundfun(x)
  iround(x)  # round(Int64,x)
end

## return numbered index from tree index (i.e., inverse to indices())
function numInd(treeInd,J)
  if length(treeInd)==0
    return 1
  else
    l=[length(treeInd)-1:-1:0;]
    return sum(J.^l)+sum((treeInd-1).*J.^l)+1 #tree.ind[length(tree.ind)]
  end
end

## return numbered index from tree index when J varies with m
function numIndJ(treeInd,J)
  if length(treeInd)==0
    return 1
  else
    m=length(treeInd)
    lastindPrevRes=sum([prod(J[1:x]) for x in 1:m-1])+1
    temp=[prod(J[x:m]) for x in 2:m+1]
    return lastindPrevRes + sum((treeInd-1).*temp) + 1
  end
end


## create r that decreases with m
function createR(r0,M)
  c=ones(M+1)
  for m=2:M
    c[m+1]=c[m]+1/c[m]^2
  end
  r=ifloor( r0 ./ (c.^2) )
  return r
end


## cholesky decomposition (returns L s.t. A=LL')
function cholesky(A)
  try
    return chol(A,:L) # chol(A+.01*eye(size(A,1)))'
  catch e
    #println("caught error: $e")
    succ=false; di=size(A,1); attempt=1; B=0;
    while succ==false
      A += 10.0^(-10+attempt)*eye(di)
      try
        B=chol(A,:L)
        succ=true
      end
      attempt += 1
    end
    #println("$attempt attempts")
    return B
  end
end


## calculate inv(C.chol)*A, where C.chol is cholesky factor (L)
function sol(Cchol,A)
  Cchol\A
end


## calculate t(A)*B
function tp(A,B)
  A'*B
end


## calculate log determinant from cholesky factor of a matrix
function logDet(matChol)
  2*sum(log(diag(matChol))) #2*logdet(matChol)
end


### converts list of predictions from MRA into data frame
using DataFrames
function dfpred(predlocs,preds,orderlist)
  plocs=zeros(0,size(predlocs[1],2)); postmean=zeros(0); postsd=zeros(0); ordervec=zeros(0);
  for ind=firstdefined(preds):length(preds)
    plocs=vcat(plocs,predlocs[ind])
    postmean=vcat(postmean,preds[ind][:,1])
    postsd=vcat(postsd,sqrt(preds[ind][:,2]))
    ordervec=vcat(ordervec,orderlist[ind])
  end
  ord=sortperm(ordervec)
  return DataFrame(predlon=plocs[ord,1],predlat=plocs[ord,2],postmean=postmean[ord],postsd=postsd[ord])
end


### returns 2-column matrix containing all combinations of two vectors
function expandgrid(vec1,vec2)
  combs=Array(typeof(vec1[1]),length(vec1)*length(vec2),2);
  counter=1;
  for x=1:length(vec1)
    for y=1:length(vec2)
      combs[counter,:]=[vec1[x];vec2[y];];
      counter += 1;
    end
  end
  return combs
end


### finds first defined element in a list (array of type Any)
function firstdefined(list)
  for i = 1:length(list)
    if isdefined(list,i)
      return i
    end
  end
  return 0
end

### round to nearest odd number (or zero)
function oddround(x)
  abs(x)<.5 ? 0 : roundfun(roundfun(x/2)*2+sign(x/2-roundfun(x/2))+1*(x/2==roundfun(x/2)))
end

## simulate stationary data on 1-D grid (function DLSimulate in R package ltsa)
function gridsim(r)
  n=length(r)
  a = randn(n) #[1:n;]  # for testing
  z = zeros(n)
  sigmasqk = r[1]
  error = a[1] * sqrt(sigmasqk)
  phi = Array(Float64,1)
  phi[1]=r[2]/r[1]
  sigmasqk = r[1] * (1 - phi[1]^2)
  error = a[2] * sqrt(sigmasqk)
  z[2] = error + phi[1] * z[1]
  sigmasqkm1 = sigmasqk
  for k = 2:(n - 1)
    phikk = (r[k + 1] .- phi' * r[k:-1:2])./sigmasqkm1
    sigmasqk = sigmasqkm1 .* (1 - phikk.^2)
    phi -= phikk .* reverse(phi)
    append!(phi,phikk)
    sigmasqkm1 = sigmasqk
    z[k + 1] = (a[k + 1]*sqrt(sigmasqk) + phi'*z[k:-1:1])[1]
  end
  z
end


## simulate stationary data on 1-D grid (based on function DLLoglikelihood in R package ltsa)
function dllik(r,z)
  n = length(z)
  error = zeros(n)
  sigmasq = zeros(n)
  error[1] = z[1]
  sigmasq[1] = r[1]
  phi=zeros(1)
  phi[1] = r[2]/r[1]
  error[2] = z[2] - phi[1] * z[1]
  sigmasqkm1 = r[1] * (1 - phi[1]^2)
  sigmasq[2] = sigmasqkm1
  logg = log(r[1]) + log(sigmasqkm1)
  for k=2:(n-1)
    phikk = (r[k + 1] .- phi' * r[k:-1:2])./sigmasqkm1
    sigmasqk = sigmasqkm1 .* (1 - phikk.^2)
    phi -= phikk .* reverse(phi)
    append!(phi,phikk)
    sigmasqkm1 = sigmasqk
    logg += log(sigmasqk)
    error[k + 1] = (z[k + 1] - phi'*z[k:-1:1])[1]
    sigmasq[k + 1] = sigmasqk[1]
  end
  S = sum((error .* error)./sigmasq)
  LogL = S + logg
  return LogL[1]
end


## simulate stationary data on 1-D grid (function DHSimulate in R package ltsa)
function DHsim(r)
  m = length(r)
  N = 2^ceilfun(log2(m-1)) #2^ceil(Int64,log2(m-1))
  acvf = [r; zeros(N-m+1)]
  append!(acvf,acvf[end-1:-1:2])
  g = real(fft(acvf))
  if (any(g .< 0.0))
    println("Davies-Harte nonnegativity condition not valid")
    return NaN
  else
    Z = complex(randn(N-1),randn(N-1))
    Z2 = 2 + sqrt(2) * randn(2)
    Z = [Z2[1]; Z; Z2[2]; conj(reverse(Z))]
    X = real(bfft(sqrt(g).*Z))/sqrt(2*N)
    z = X[1:m]/sqrt(2)
    return z
  end
end



##  function to create prior quantities
function createPrior(co,knotsb,RpriorCholb,KcBb,dataj,predlocsj=[NaN;],varEps=0)

  m=length(knotsb)-1
  mlM = ((length(dataj)==1) && (isnan(dataj[1])))  # is m less than M?

  # create prior quantities
  RpriorCholj=0
  KcBc=Array(Matrix,m)
  V=Array(Matrix,m+1)
  for l in 0:m
    V[l+1]=co(knotsb[m+1],knotsb[l+1])
    for k in 0:(l-1)
      if l<m
        V[l+1] -= tp(KcBc[k+1],KcBb[l+1][k+1]) #V[l+1] = V[l+1]-tp(KcBc[k+1],KcBb[l+1][k+1])
      else
        V[l+1] -= tp(KcBc[k+1],KcBc[k+1])
      end
    end
    if l<m
      KcBc[l+1]=sol(RpriorCholb[l+1],V[l+1]')
    else
      RpriorCholj = cholesky(V[l+1]  + (mlM ? 0 : varEps*eye(size(V[l+1],1)))  )
    end
  end

  if mlM

    return ((RpriorCholj,KcBc),(),())

  else  # begin inference for regions at finest resolution M

    # pre-compute solves
    Sicy=sol(RpriorCholj,dataj)
    SicB=Array(Matrix,m)
    for l in 0:m-1
      SicB[l+1]=sol(RpriorCholj,V[l+1])
    end

    # inference quantities
    wtj=Array(Vector,m)
    Atj=Array(Matrix,m,m)
    for l in 0:m-1
      wtj[l+1] = tp(SicB[l+1],Sicy)
      for k in l:m-1
        Atj[l+1,k+1] = tp(SicB[l+1],SicB[k+1])
      end
    end

    if (length(predlocsj)==1) && (isnan(predlocsj[1]))

      loglikj = logDet(RpriorCholj) + tp(Sicy,Sicy)
      retlikpred=loglikj[1]

    else

      Rpriorchol=[RpriorCholb;{RpriorCholj}];
      KcB=[KcBb;{KcBc}];

      # calculate Bp and L
      KcBp=Array(Matrix,m+1); Vp=Array(Matrix,m+1)
      for l=0:m
        Vp[l+1]=co(predlocsj,knotsb[l+1])
        for k=0:(l-1)
          Vp[l+1] -= tp(KcBp[k+1],KcB[l+1][k+1])
        end
        KcBp[l+1]=sol(Rpriorchol[l+1],Vp[l+1]')
      end
      Vpp=co(predlocsj,predlocsj)
      for l=1:m; Vpp -= tp(KcBp[l],KcBp[l]); end

      # initialize prediction inference
      postmeanj=Array(Float64,size(predlocsj,1),m+1); postvarj=Array(Float64,size(predlocsj,1),m+1)
      postmeanj[:,m+1]=tp(KcBp[m+1],Sicy)
      postvarj[:,m+1]=diag(Vpp-tp(KcBp[m+1],KcBp[m+1]))
      Btildej=Array(Any,m+1); Btildej[m+1]=Array(Matrix,m)
      for k=0:(m-1)
        Btildej[m+1][k+1]=Vp[k+1]-tp(KcBp[m+1],SicB[k+1])
      end

      retlikpred=(postmeanj,postvarj,Btildej)

    end

    return ((RpriorCholj,KcBc),(Atj,wtj),retlikpred)

  end

end



#######  posterior inference (when going from fine to coarse resolution)
function postInf(Rpriorcholj,wtildechildren,Atildechildren,pred=false)

  m=length(wtildechildren[1])-1

  # sum up over children tildes
  w=Array(Vector,m+1)
  A=Array(Matrix,m+1,m+1)
  for l in 0:m
    w[l+1] =sum( [ x[l+1] for x in wtildechildren ] )
    for k in l:m
      A[l+1,k+1] = sum( [ x[l+1,k+1] for x in Atildechildren] )
    end
  end

  # calculate cholesky of K.inv
  Rpost = Rpriorcholj*Rpriorcholj' + A[m+1,m+1]
  Rpostcholj = cholesky(Rpost)

  # pre-compute the solves required later
  Kcwj=sol(Rpostcholj,w[m+1])
  KcAj=Array(Matrix,m)
  for l in 0:m-1
    KcAj[l+1]=sol(Rpostcholj,A[l+1,m+1]')
  end

  # calculate w.tilde and A.tilde
  if m==0
    returnmain=(nothing,nothing)
  else
    wtildecurj=Array(Vector,m)
    Atildecurj=Array(Matrix,m,m)
    for l in 0:m-1
      wtildecurj[l+1] = w[l+1] - tp(KcAj[l+1],Kcwj)
      for k in l:m-1
        Atildecurj[l+1,k+1] = A[l+1,k+1] - tp(KcAj[l+1],KcAj[k+1])
      end
    end
    returnmain=(wtildecurj,Atildecurj)
  end

  if pred
    return (returnmain,(Rpostcholj,Kcwj,KcAj))
  else
    loglikj=logDet(Rpostcholj)-logDet(Rpriorcholj)-tp(Kcwj,Kcwj)
    return (returnmain,loglikj)
  end

end




#######  prediction at finest resolution at the end
function predfct(postmeanj,postvarj,Btildej,Rpostcholb,KcAb,Kcwb)

  M=size(postmeanj,2)-1
  for k=(M-1):-1:1
    KcBtilde=sol(Rpostcholb[k+1],Btildej[k+2][k+1]')
    Btildej[k+1]=Array(Matrix,k)
    for l=(k-1):-1:0
      Btildej[k+1][l+1]=Btildej[k+2][l+1]-tp(KcBtilde,KcAb[k+1][l+1])
    end
  end
  for k=0:(M-1)
    KcBtildecur=sol(Rpostcholb[k+1],Btildej[k+2][k+1]')
    postmeanj[:,k+1]=tp(KcBtildecur,Kcwb[k+1])
    postvarj[:,k+1]=sum(KcBtildecur.^2,1)  # diag(tp(Kc.Btilde.cur,Kc.Btilde.cur))
  end
  predsj=hcat(sum(postmeanj,2),sum(postvarj,2))
  return(predsj)

end

