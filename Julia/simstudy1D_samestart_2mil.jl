# cd("S:/work/distributedComputing/code")
include("MRAfunctions.jl")
using Distances  # for distance computation (pairwise())


### general settings

domain=[-1e-6;1];

covFunBasic=function(x)  # matern with nu=1.5
  return  .95*(1 + sqrt(3)*x).*exp(-sqrt(3)*x) + .05*(abs(x).<1e-10)
end;

function co(locs1,locs2,theta=.05)
    return covFunBasic(pairwise(Euclidean(),locs1',locs2')/theta) # pairwise calculates dists between columns
end;


### data size
r0=30;
J0=2;
M0=16;
n=r0*J0^M0;

### competitors
MRAM=2.^[0:3;];
r1=[r0*J0.^[0:5;];r0*J0^5/3;];
Ms=[MRAM; ones(Int64,length(r1));];
nmodels=length(Ms);
rs=[r0*ones(length(MRAM)); r1];
Js=Array(Any,nmodels);
for i=1:nmodels
  Js[i]= Ms[i]==1 ? [roundfun(n/rs[i]);] : roundfun(J0^(M0/Ms[i]))*ones(Int64,Ms[i]);
end

for repnum=1:5

  results=Array(Float64,nmodels,4)*NaN;

  ### create "full" dataset, to be subsampled later
  nAll=r0*J0^M0;
  locsAll=[1:nAll;]/nAll*domain[2];
  srand(99+repnum);
  if nAll<10000
    y=vec(chol(co(locsAll,locsAll))*randn(nAll,1));
  else
    y=DHsim(vec(co([locsAll[1];],locsAll)));
  end

  asym=2 # increasing- vs. fixed-domain asymptotics

  ### sample data
  if asym==1
    inds=1:n;
  else
    spacing=roundfun(nAll/n);
    inds=spacing:spacing:nAll;
  end
  locs=locsAll[inds];
  z=y[inds]; # +rnorm(n)*.1
  dom=[domain[1]; locs[end]];

  for mod=1:nmodels

    @show (repnum,mod)
    results[mod,1]=Ms[mod]
    results[mod,2]=rs[mod]
    (indices,knots,data,bounds)=partition1D(Js[mod],dom,locs,z);
    try
      tic();
      temp=MRA(co,data,knots,indices);
      results[mod,3]=toc();
      results[mod,4]=-0.5*temp-n/2*log(2*pi)
    catch e
      println("for model $mod, caught error: $e")
    end

    writecsv("../results/sim1D_2mill_theta05_$repnum.csv",results)

  end

end
