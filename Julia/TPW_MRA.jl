# cd("S:/work/distributedComputing/code")
include("MRAfunctions.jl")
using Distances  # for distance computation (pairwise())


### read in the TPW data (t=9 MIRS data)
data = readcsv("../data/MIRSmra.csv");
locs=hcat(data[:,1],data[:,2]);
meandat=mean(data[:,3]);
z=data[:,3]-meandat;
n=length(z); # 271k
domain=[[minimum(locs[:,j]) for j=1:2]; [maximum(locs[:,j]) for j=1:2]];


### covariance function
function covFun(locs1,locs2,theta)
  x=pairwise(Euclidean(),locs1',locs2')/theta[2]
  return theta[1]*(1 + sqrt(3)*x).*exp(-sqrt(3)*x)
end;


using Optim
optimsummary=Array(Float64,3,7);
modelsummary=Array(Float64,3,3);

# competitors
Js=Any[[2^8;],[2; 2; 4; 8; 8; 16],[2^10;]];
rs=Any[[0; 1000],[34; 34; 34; 15; 16; 15],[187; 187]];

for model=1:3

  ### values for partitioning ###
  J=Js[model];
  M=length(J);
  rknots=rs[model];

  ### partition the data
  predvec=expandgrid([55:0.25:106;]*1.0,[20:.25:50;]*1.0) #prediction grid
  (indices,knots,data,bounds,predlocs,orderlist)=partition2D(J,domain,locs,z,rknots,predvec);
  # using Gadfly; rs=[size(x,1)*1.0::Float64 for x in knots]; plot(x=1:length(rs), y=rs)

  ### record summaries for paper
  rall=zeros(length(knots)); for i=1:length(knots); rall[i]=size(knots[i],1); end
  function co(locs1,locs2)
    return covFun(locs1,locs2,[40.0, 2.3])
  end;
  tic(); MRA(co,data,knots,indices,varEpsilon=(4.5)^2); liktime=toc();
  modelsummary[model,:]=[prod(J),mean(rall),liktime]
  writecsv("../results/TPW_summary2.csv",modelsummary)

  ### maximize the likelihood ###
  negloglik=function(theta)
    function co(locs1,locs2)
      return covFun(locs1,locs2,theta)
    end;
    return MRA(co,data,knots,indices,varEpsilon=theta[3])  #(4.5)^2)
  end
  tic();
  optimres=optimize(negloglik,[40.0, 2.3, 4.5^2],iterations=200,store_trace=true,show_trace=true)
  totaltime=toc();
  loglikelihood= -0.5*optimres.f_minimum-n/2*log(2*pi)

  optimsummary[model,:]=[M;mean(rknots);totaltime/optimres.iterations;loglikelihood;optimres.minimum;];
  writecsv("../results/TPW_optimsummary2_model$(model).csv",optimsummary)


  ### make predictions ###
  theta=readcsv("../results/TPW_optimsummary2_model$(model).csv")[model,5:end];
  function co(locs1,locs2); covFun(locs1,locs2,theta); end;
  @time preds=MRA(co,data,knots,indices,predlocs=predlocs,varEpsilon=theta[3]); #(4.5)^2)

  ### rearrange predictions into vector/matrix
  predframe=dfpred(predlocs,preds,orderlist)
  predframe[:postmean]=predframe[:postmean]+meandat
  writecsv("../results/TPW_preds2_model$model.csv",array(predframe))


  ### prediction comparison ###

  for j=1:3
    for selection=1:2

      if selection==1
        srand(99+j);
        testind=sample(1:n,5000,replace=false) # sample randomly
      else
        srand(99+j);
        lb=domain[1:2]+rand(2).*(domain[3:4]-5-domain[1:2]) # bottom left corner of sample region
        testind=find( (lb[1].<locs[:,1].<lb[1]+5) .* (lb[2].<locs[:,2].<lb[2]+5) )
      end
      obsind=setdiff(1:n,testind)
      (indices,knots,data,bounds,predlocs,orderlist)=partition2D(J,domain,locs[obsind,:],z[obsind],rknots,locs[testind,:]);
      theta=readcsv("../results/TPW_optimsummary2_model$(model).csv")[model,5:end]
      function co(locs1,locs2); covFun(locs1,locs2,theta); end;
      preds=MRA(co,data,knots,indices,predlocs=predlocs,varEpsilon=theta[3]); #(4.5)^2)
      predframe=dfpred(predlocs,preds,orderlist)
      predframe[:postsd] = sqrt(predframe[:postsd].^2 + theta[3]) # 4.5^2)
      writecsv("../results/TPW_predcomp2_model$(model)_selection$(selection)_rep$(j).csv",
               hcat(array(predframe),z[testind]))

    end
  end

end



#######   local kriging   ###########

# use parameter estimates from block likelihood
theta=readcsv("../results/TPW_optimsummary2_model1.csv")[1,5:end]
function co(locs1,locs2); covFun(locs1,locs2,theta); end;

k=20; # number of nearest neighbors used

for j=1:3
  for selection=1:2

    # split data into test and training data
    if selection==1
      srand(99+j);
      testind=sample(1:n,5000,replace=false) # sample randomly
    else
      srand(99+j);
      lb=domain[1:2]+rand(2).*(domain[3:4]-5-domain[1:2]) # bottom left corner of sample region
      testind=find( (lb[1].<locs[:,1].<lb[1]+5) .* (lb[2].<locs[:,2].<lb[2]+5) )
    end
    obsind=setdiff(1:n,testind);
    z_obs=z[obsind];
    locs_obs=locs[obsind,:];
    locs_pred=locs[testind,:];

    # find nearest neighbors and carry out local kriging for each pred. location
    predmean=zeros(size(locs_pred,1)); predsd=zeros(size(locs_pred,1));
    tic();
    for i=1:size(locs_pred,1)

      # find k nearest neighbors
      dists=pairwise(Euclidean(),locs_pred[i,:]',locs_obs');
      nn=sortperm(squeeze(dists,1))[1:k];
      locs_nn=locs_obs[nn,:];

      # make local predictions
      nn_cov=co(locs_nn,locs_nn)+theta[3]*eye(k);
      crosscovar=co(locs_pred[i,:],locs_nn);
      predmean[i]=(crosscovar*(nn_cov\z_obs[nn]))[1];
      predsd[i]=sqrt(co(locs_pred[i,:],locs_pred[i,:])-crosscovar*(nn_cov\crosscovar') + theta[3])[1];

    end
    predtime=toc();
    print(predtime)
    writecsv("../results/TPW_predcomp2_model4_selection$(selection)_rep$(j).csv",
             hcat(locs_pred,predmean,predsd,z[testind]))

  end
end



#=
### computing the comparison results

# crps function
using Distributions
function crps(x,mu,sigma)
  xz = (x-mu)./sigma;
  sigma .* (xz.*(2*cdf(Normal(),xz)-1) + 2*pdf(Normal(),xz) - 1/sqrt(pi))
end

# interval score function
function intscore(alpha,x,int)
  (int[:,2]-int[:,1])+2/alpha*((int[:,1]-x).*(x.<int[:,1])+(x-int[:,2]).*(x.>int[:,2]))
end

# compute summary tables
reps=3
mspeAll=Array(Float64,4,2,reps); crpsAll=Array(Float64,4,2,reps); isAll=Array(Float64,4,2,reps)
alpha=.1; qu=-quantile(Normal(),alpha/2)
for selection=1:2
  for model=1:4
    for j=1:reps
      testdata=readcsv("../results/TPW_predcomp2_model$(model)_selection$(selection)_rep$(j).csv")[:,3:5]
      mspeAll[model,selection,j]=mean((testdata[:,1]-testdata[:,3]).^2)
      interval=hcat(testdata[:,1]-qu*testdata[:,2],testdata[:,1]+qu*testdata[:,2])
      isAll[model,selection,j]=mean(intscore(alpha,testdata[:,3],interval))
      crpsAll[model,selection,j]=mean(crps(testdata[:,3],testdata[:,1],testdata[:,2]))
    end
  end
end

optimsummary=readcsv("../results/TPW_optimsummary2_model3.csv")[:,[1,5,6,7,4]]
rellik=optimsummary[:,end]/optimsummary[1,end];
modelsummary=readcsv("../results/TPW_summary2.csv")[:,[2,3]]
temp=hcat(optimsummary[:,1:4],rellik,modelsummary,sqrt(squeeze(mean(mspeAll[1:3,:,:],3),3)),
  squeeze(mean(crpsAll[1:3,:,:],3),3))
summtable=temp[[2,3,1],[1,6,7,2:4,5,8,10,9,11]]
writecsv("../results/TPWtable2.csv",summtable)

# results for local kriging
hcat(sqrt(squeeze(mean(mspeAll[4,:,:],3),3)),squeeze(mean(crpsAll[4,:,:],3),3))[:,[1,3,2,4]]

=#
