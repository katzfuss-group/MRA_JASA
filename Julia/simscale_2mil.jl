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


### output quantities
r0=30;
J0=4;
Mrange=[3:8;]; #[3:8;]  # n=r0*J0.^Mrange
nmodels=5;


### create "full" dataset, to be subsampled later
nAll=r0*J0^Mrange[end];
locsAll=[1:nAll;]/nAll*domain[2];

for repnum=1:5

  times=Array(Float64,length(Mrange),nmodels*2+1)*NaN;
  logliks=Array(Float64,length(Mrange),nmodels*2+1)*NaN;

  srand(99+repnum);
  if nAll<10000
    y=vec(chol(co(locsAll,locsAll))*randn(nAll,1));
  else
    y=DHsim(vec(co([locsAll[1];],locsAll)));
  end


  ### start loop over data sizes

  for i=1:length(Mrange)

    for asym=1:2 # increasing- vs. fixed-domain asymptotics

      ### sample data
      n=r0*J0^Mrange[i];
      if asym==1
        inds=1:n;
      else
        spacing=roundfun(nAll/n);
        inds=spacing:spacing:nAll;
      end
      locs=locsAll[inds];
      z=y[inds]; # +rnorm(n)*.1
      dom=[domain[1]; locs[end]];
      times[i,1]=n; logliks[i,1]=n;

      for mod=1:nmodels

        @show (i,asym,n,mod)

        run=1;
        if mod==1  # full model
          J=Int64[];
          if n<40000
            run=1 # only do full model if n small enough
          elseif n<500000
            run=2
          else
            run=3
          end
        elseif mod==2 # fast M-RA (J=8)
          Mfast=Mrange[i]/2;
          J=J0^2*ones(Int64,roundfun(Mfast));
          Mfast==roundfun(Mfast) ? run=true : run=false
        elseif mod==3  # slow M-RA (J=4)
          J=J0.*ones(Int64,Mrange[i]);
        elseif mod==4  # fast 1-RA
          J=Int64[J0^(Mrange[i])/8];  # r=8*r0
        else  # slow 1-RA
          J=Int64[J0^3]; # r=r0*J0.^(Mrange-3)
          (n/J[1])<4000 ? run=1 : run=3 # only run if r small enough
        end
        if n==nAll && asym==2  # don't need to run fixed-domain for largest n
          run=3;
        end

        if run==1
          (indices,knots,data,bounds)=partition1D(J,dom,locs,z);
          @show (n,length(knots[end]),J)
          try
            tic();
            temp=MRA(co,data,knots,indices);
            times[i,mod+1+(asym-1)*nmodels]=toc();
            logliks[i,mod+1+(asym-1)*nmodels]=-0.5*temp-n/2*log(2*pi)
          catch e
            println("for model $mod, caught error: $e")
          end
        elseif run==2
          temp=dllik(vec(co([locs[1];],locs)),z)
          logliks[i,mod+1+(asym-1)*nmodels]=-0.5*temp-n/2*log(2*pi)
        end

        writecsv("../results/simscale_2mil_full_times_rep$repnum.csv",times)
        writecsv("../results/simscale_2mil_full_logliks_rep$repnum.csv",logliks)

      end

    end

  end

end


#=

nm=6; nmodels=5; reps=1;
lognormIncrAll=Array(Float64,nm,nmodels,reps); lognormFixedAll=Array(Float64,nm,nmodels,reps);
timesAll=Array(Float64,nm,nmodels*2+1,reps);
for repnum=1:reps
  times=readcsv("../results/simscale_2mil_times_rep$repnum.csv");
  logliks=readcsv("../results/simscale_2mil_logliks_rep$repnum.csv");
  logliks[end,(1:nmodels)+1+nmodels]=logliks[end,(1:nmodels)+1];
  lognormIncrAll[:,:,repnum]=(logliks[:,(1:nmodels)+1]./logliks[:,3+1]).^sign(logliks[:,3+1]);
  lognormFixedAll[:,:,repnum]=(logliks[:,(1:nmodels)+1+nmodels]./logliks[:,3+1+nmodels]).^sign(logliks[:,3+1+nmodels]);
  times[end,(1:nmodels)+1+nmodels]=times[end,(1:nmodels)+1];
  timesAll[:,:,repnum]=times;
end
times=mean(timesAll,3); lognormIncr=mean(lognormIncrAll,3); lognormFixed=mean(lognormFixedAll,3);
using DataFrames
simdata=DataFrame(
  n = vcat(times[:,1],times[:,1],times[:,1],times[:,1],times[:,1]),
  model = vcat(rep("full",nm),rep("M-RA F",nm),rep("M-RA S",nm),rep("1-RA F",nm),rep("1-RA S",nm)),
  time = vec(mean(hcat(vec(times[:,(1:nmodels)+1]),vec(times[:,(1:nmodels)+1+nmodels])),2)),
  loglikIncr = vec(lognormIncr),
  loglikFixed = vec(lognormFixed)
);
# writecsv("../results/simscale_2mil_means.csv",array(simdata))


=#
