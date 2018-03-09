# cd("S:/work/distributedComputing/code")
include("MRAfunctions.jl")
using Distances  # for distance computation (pairwise())

blas_set_num_threads(1)

### general settings
domain=[-1e-6;-1e-6;1;1]; # l,b,r,t

function co(locs1,locs2)
  x=pairwise(Euclidean(),locs1',locs2')/.3
  return exp(-x) + .05*(abs(x).<1e-10)
end;


### sizes/dimensions
r0=49;
J0=4;
M0=8;
n=r0*J0^M0;


### competitors
Mexp=[1:3;];
MRAM=2.^Mexp;
r1=[r0*J0.^[0:2;];];
modind=[rep(1,length(MRAM)); rep(2,length(r1));];
Ms=[MRAM; ones(Int64,length(r1));];
rall=[r0*ones(Int64,length(MRAM)); r1;];
nmodels=length(Ms);
rs=Array(Any,nmodels); Js=Array(Any,nmodels);
for i=1:nmodels
  if modind[i]==1
    rs[i]=r0*ones(Ms[i]+1); Js[i]=roundfun(J0^(M0/MRAM[i]))*ones(Int64,Ms[i]);
  else
    rs[i]=rall[i]*ones(2); Js[i]=[roundfun(n/rall[i]);];
  end
end


### data locations on grid
nAll=r0*J0^M0;
n1=roundfun(sqrt(nAll)); # length of grid in each dimension
gridside=linspace(0,1,n1);
xmat=Array(Float64,n1,n1);
ymat=Array(Float64,n1,n1);
for i=1:n1
  for j=1:n1
    xmat[i,j]=gridside[i];
    ymat[i,j]=gridside[j];
  end
end


### start running the algorithm
for repnum=1:5

  results=Array(Float64,nmodels,5)*NaN;

  ### sample data
  datamat=readcsv("../results/2Dsim$repnum.csv");
  n=r0*J0^M0;
  n1=roundfun(sqrt(n));
  inds=1:n1;
  z=datamat[inds,inds][:];
  locs=hcat(xmat[inds,inds][:],ymat[inds,inds][:]);

  for mod=6:nmodels

    @show (repnum,mod)
    results[mod,1]=modind[mod]
    results[mod,2]=Ms[mod]
    results[mod,3]=rs[mod][end]
    (indices,knots,data,bounds)=partition2D(Js[mod],domain,locs,z,rs[mod]);
    try
      tic();
      temp=MRA(co,data,knots,indices);
      results[mod,4]=toc();
      results[mod,5]=-0.5*temp-n/2*log(2*pi)
    catch e
      println("for model $mod, caught error: $e")
    end

    writecsv("../results/sim2D_nug_rep$repnum.csv",results)

  end

end
