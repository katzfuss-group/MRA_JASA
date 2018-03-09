include("helpfunctions.jl")


### main MRA function

function MRA(co,data,knots,indices;predlocs=zeros(0),varEpsilon=0)

  ## extract dimensions and other constants
  nInd=length(indices)
  M=length(indices[end])
  J = (nInd==1 ? 1 : indices[end])
  pred=(length(predlocs)>0)

  # indices for each resolution
  indres=Array(Vector,M+1)
  indres[1]=Int[1]
  for m in 1:M
    indres[m+1]=[(indres[m][end]+1):(indres[m][end]+prod(J[1:m]));]
  end

  ## initialize
  RpriorChol=Array(Any,nInd); KcB=Array(Any,nInd)
  Atildeprev=Array(Any,nInd); wtildeprev=Array(Any,nInd)
  loglik=Array(Float64,nInd)
  if pred
    postmean=Array(Matrix,nInd); postvar=Array(Matrix,nInd); Btilde=Array(Any,nInd)
    preds=Array(Matrix,nInd); Rpostchol=Array(Matrix,nInd); KcA=Array(Vector,nInd); Kcw=Array(Vector,nInd)
  end

  ## going from coarsest to finest resolution
  for m=0:M
    for ind in indres[m+1]

      # indices of parent nodes
      par=ones(Int64,m)
      for k in 1:(m-1)
          par[k+1]=numIndJ(indices[ind][1:k],J)
      end

      # create prior quantities etc
      predlocsj = (pred ? predlocs[ind] : [NaN;])
      # knotsb=knots[[par;ind]]; RpriorCholb=RpriorChol[par]; KcBb=KcB[par]; dataj=data[ind]; varEps=varEpsilon
      temp=createPrior(co,knots[[par;ind;]],RpriorChol[par],KcB[par],data[ind],predlocsj,varEpsilon)
      (RpriorChol[ind],KcB[ind])=temp[1]
      if m==M
        (Atildeprev[ind],wtildeprev[ind])=temp[2]
        if pred
          (postmean[ind],postvar[ind],Btilde[ind])=temp[3]
        else
          loglik[ind]=temp[3]
        end
      end

    end
  end

  KcB=0  # clear KcB from memory


  ## posterior inference (from finest to coarsest resolution)
  for m=M-1:-1:0

    children=zeros(Int64,J[m+1]); Atildecur=Array(Any,nInd); wtildecur=Array(Any,nInd)

    for ind=indres[m+1]

      # indices of the J child nodes
      for j=1:J[m+1]
        children[j]=numIndJ([indices[ind];j],J)
      end

      # calculate posterior quantities
      # Rpriorcholj=RpriorChol[ind];wtildechildren=wtildeprev[children];Atildechildren=Atildeprev[children]
      temp=postInf(RpriorChol[ind],wtildeprev[children],Atildeprev[children],pred)
      (wtildecur[ind],Atildecur[ind])=temp[1]
      if pred
        (Rpostchol[ind],Kcw[ind],KcA[ind])=temp[2]
      else
        loglik[ind]=temp[2][1]
      end

    end

    Atildeprev=Atildecur; wtildeprev=wtildecur

  end


  ## spatial prediction
  if pred

    for ind in indres[M+1]  # only for finest resolution
      if M>0
        # indices of parent nodes
        par=ones(Int64,M)
        for k in 1:(M-1)
          par[k+1]=numIndJ(indices[ind][1:k],J)
        end
        # obtain predictive means and variances
        # (postmeanj,postvarj,Btildej,Rpostcholb,KcAb,Kcwb)=(postmean[ind],postvar[ind],Btilde[ind],Rpostchol[par],KcA[par],Kcw[par])
        preds[ind]=predfct(postmean[ind],postvar[ind],Btilde[ind],Rpostchol[par],KcA[par],Kcw[par])
      else
        preds[ind]=hcat(postmean[ind],postvar[ind])
      end
    end
    return preds

  else

    return sum(loglik)

  end

end







#######################    partitioning with varying J     #######################
function partition1D(J,domain,locs,z,r=-1,predvec=zeros(0);regknots=true)

  deleteat!(J,find(J.<2)) # J should always be at least 2
  M=length(J)
  if r[1]<0
    if M==0
      r=length(z)
    else
      r=ceilfun(length(z)/prod(J))*ones(Int64,M+1)
    end
  end

  ## function to return parent (numbered) index from (numbered) index
  function parInd(ind,J)
    fullj=indices[ind] # indices is a list (Array of type Any)
    parj=fullj[1:end-1]
    return numIndJ(parj,J)
  end

  # function creates array (list) of indices for all nodes in the tree structure
  function createIndices(J)
    M=length(J)
    prevInd={Int64[]}
    indices=prevInd
    for m in 1:M
      curInd=Any[]
      for jprev in 1:length(prevInd)
        for jnew in 1:J[m]
          append!(curInd,{[prevInd[jprev]; jnew]})
        end
      end
      append!(indices,curInd)
      prevInd=curInd
    end
    return indices
  end

  indices=createIndices(J)
  nInd=length(indices)

  # create knot tree structure
  locsc=copy(locs)
  zc=copy(z)
  knots=Array(Vector,nInd)
  data=Array(Vector,nInd)
  bounds=Array(Float64,nInd,2)
  for ind=1:nInd
    m=length(indices[ind])
    if ind==1
      bounds[ind,:]=domain
    else
      lastInd=indices[ind][end]
      pb=bounds[parInd(ind,J),:]
      bounds[ind,:]=linspace(pb[1],pb[2],J[m]+1)[lastInd+(0:1)]
    end
    if m<M
      if regknots==true
        knotdist=(bounds[ind,2]-bounds[ind,1])/r[m+1]
        knots[ind]=linspace(bounds[ind,1]+knotdist/2,bounds[ind,2]-knotdist/2,r[m+1])
      else
        t=floorfun((r[m+1]-1)/2)
        knotloc=linspace(1/t,1,t).^1.5*(t-1)/t
        knots[ind]=mean(bounds[ind,:])+[-knotloc[t:-1:1];0;knotloc]*(bounds[ind,2]-bounds[ind,1])/2
      end
      data[ind]=[NaN;]
    else
      indsub=find( bounds[ind,1] .< locsc .<= bounds[ind,2] )
      knots[ind]=locsc[indsub]
      data[ind]=zc[indsub]
      deleteat!(locsc,indsub)
      deleteat!(zc,indsub)
    end
  end

  if length(predvec)==0
    predlocs=zeros(0)
  else
    predlocsvec=copy(predvec)
    predlocs=Array(Vector,nInd)
    for ind=1:nInd
      if length(indices[ind])==M
        indsub=find( bounds[ind,1] .< predlocsvec .<= bounds[ind,2] )
        predlocs[ind]=predlocsvec[indsub]
        deleteat!(predlocsvec,indsub)
      else
        predlocs[ind]=knots[ind]
      end
    end
  end

  return (indices,knots,data,bounds,predlocs)

end




####################    partitioning of 2-dimensional rectangular domain     ####################
function partition2D(J,domain,locs,z,rknots=-1,predvec=zeros(0),regknots=true)

  ## split rectangle with Jlon bins along horizontal axis, Jlat bins along vertical axis
  function splitRectangle(rectbounds,Jlon=2,Jlat=2)
    blon=linspace(rectbounds[1],rectbounds[3],Jlon+1)
    blat=linspace(rectbounds[2],rectbounds[4],Jlat+1)
    return hcat(expandgrid(blon[1:end-1],blat[1:end-1]),expandgrid(blon[2:end],blat[2:end]))
  end

  ## split rectangle in 4 or 2 subregions
  function split2d(rectbounds,Jm=4)
    horlver= (rectbounds[3]-rectbounds[1] >= rectbounds[4]-rectbounds[2]) # horizontal dim longer than vertical?
    if Jm==2 # split in 2 regions along longer axis
      if horlver
        return splitRectangle(rectbounds,2,1) # split along horizontal axis
      else
        return splitRectangle(rectbounds,1,2) # split along vertical axis
      end
    elseif Jm==8
      if horlver
        return splitRectangle(rectbounds,4,2) # split along horizontal axis
      else
        return splitRectangle(rectbounds,2,4) # split along vertical axis
      end
    elseif sqrt(Jm)==round(sqrt(Jm)) # if Jm is quadratic number
      return splitRectangle(rectbounds,roundfun(sqrt(Jm)),roundfun(sqrt(Jm)))
    elseif Jm==2^11
      return splitRectangle(rectbounds,2^6,2^5)
    else
      println("Error: Only J=2 and quadratic numbers are currently supported")
    end
  end


  ## function to return parent (numbered) index from (numbered) index
  function parInd(ind,J)
    fullj=indices[ind] # indices is a list (Array of type Any)
    parj=fullj[1:end-1]
    return numIndJ(parj,J)
  end

  # function creates array (list) of indices for all nodes in the tree structure
  function createIndices(J)
    M=length(J)
    prevInd={Int64[]}
    indices=prevInd
    for m in 1:M
      curInd=Any[]
      for jprev in 1:length(prevInd)
        for jnew in 1:J[m]
          append!(curInd,{[prevInd[jprev]; jnew]})
        end
      end
      append!(indices,curInd)
      prevInd=curInd
    end
    return indices
  end

  # function for irregular knots spacing
  function irknots(bou,nknot)
    t=floor(Int64,(nknot-1)/2) # no. of knots on each side of partition
    knotloc=linspace(1/t,1,t).^1.5*(t-1)/t
    return mean(bou)+[-knotloc[t:-1:1];0;knotloc]*(bou[2]-bou[1])/2
  end

  # initial set-up work
  deleteat!(J,find(J.<2)) # J should always be at least 2
  M=length(J)
  if rknots[1]<0
    if M==0
      rknots=length(z)
    else
      rknots=ceilfun(length(z)/prod(J))*ones(Int64,M+1)
    end
  end
  dims=2

  # create indices
  indices=createIndices(J)
  nInd=length(indices)
  indres=Array(Vector,M+1)
  indres[1]=Int[1]
  for m in 1:M
    indres[m+1]=[(indres[m][end]+1):(indres[m][end]+prod(J[1:m]));]
  end

  # determine boundaries of the regions
  bounds=Array(Float64,nInd,2*dims)
  bounds[1,:]=domain
  bounds[1,1:2] -= 1e-8  # so that all locs are *inside* the domain
  counter=1
  for m=1:M
    for ind in indres[m]
      pb=bounds[ind,:] #bounds[parInd(ind,J),:]
      bounds[counter+(1:J[m]),:]=split2d(pb,J[m])
      counter += J[m]
    end
  end

  # create knot tree structure
  knots=Array(Matrix,nInd)
  data=Array(Vector,nInd)
  for ind=1:nInd
    m=length(indices[ind])
    (l,b,r,t)=bounds[ind,:]  # left, bottom, right, top of subregion
    if m<M
      if rknots[m+1]==0
        rlon=0; rlat=0;
      else
        rlon=oddround((r-l)*sqrt(rknots[m+1]/((r-l)*(t-b))))
        rlat=oddround(rknots[m+1]/rlon)
      end
      if J[m+1]>2  # regular knots
        londist=(r-l)/rlon
        latdist=(t-b)/rlat
        lonknots=linspace(l+londist/2,r-londist/2,rlon)
        latknots=linspace(b+latdist/2,t-latdist/2,rlat)
      else  # if J=2, more knots closer to boundary
        if rlon>=rlat
          lonknots=irknots([l;r;],rlon)
          latdist=(t-b)/rlat
          latknots=linspace(b+latdist/2,t-latdist/2,rlat)
        else
          latknots=irknots([b;t;],rlat)
          londist=(r-l)/rlon
          lonknots=linspace(l+londist/2,r-londist/2,rlon)
        end
      end
      temp=[[x, y] for x in lonknots, y in latknots]
      knots[ind]=hcat([x[1] for x in temp],[x[2] for x in temp])
      data[ind]=[NaN;]
    else
      indsub=find( (l .< locs[:,1] .<= r) .* (b .< locs[:,2] .<= t) )
      knots[ind]=locs[indsub,:]
      data[ind]=z[indsub]
    end
  end

  if length(predvec)==0
    predlocs=zeros(0)
    orderlist=zeros(0)
  else
    predlocs=Array(Matrix,nInd)
    orderlist=Array(Vector,nInd)
    predorder=1:size(predvec,1)
    for ind=1:nInd
      if length(indices[ind])==M
        (l,b,r,t)=bounds[ind,:]
        indsub=find( (l .< predvec[:,1] .<= r) .* (b .< predvec[:,2] .<= t) )
        predlocs[ind]=predvec[indsub,:]
        orderlist[ind]=predorder[indsub]
      else
        predlocs[ind]=knots[ind]
      end
    end
  end

  return (indices,knots,data,bounds,predlocs,orderlist)

end
