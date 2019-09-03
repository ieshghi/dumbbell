module runs
using JLD
include("sim.jl");

#Parameters are to be fed into the pars array in the following order (already unitless)
# - T_1 (T_1 / U_0)
# - T_2 (T_2 / U_0)
# - k (k*l^2/U_0)
# - lambda (lambda/l)
# - gravity (g*(gamma*u0)^2/(lambda)^5)
# - potential choice (1 = sawtooth, 2 = parabolas, ...)

struct Run
  vval
  parvals
end

function linfit(xdat,ydat)
  X = zeros(size(xdat,1),2);
  X[:,1] = xdat;  X[:,2] = 1.0*ones(size(xdat,1),1);
  coeff_pred = X\ydat;

  slope = coeff_pred[1];
  intercept = coeff_pred[2];

  return slope, intercept
end

#function any_var_run(pars,runpar,parmin,parmax,n,npart = 10,nsteps = 10^7,dt = 0.001)
#  parvals = LinRange.(parmin,parmax,n);
#
#  vvals = zeros(n*ones(length(runpar)));
#
#  for i = 1:n
#    t,x = parallelcall(npart,[1,0])
#
#  end
#end

function singlevarrun(saveadd, pars, runpar, parmin, parmax, n, npart = 10, nsteps = 10^7, dt = .0001,undersamp = 100)
    parvals = LinRange(parmin,parmax,n);
    vvals = zeros(n);
    for i = 1:n
      pars[runpar] = parvals[i];
      t,x = sim.parallelcall(npart, [1,0], dt, nsteps, pars,undersamp);
      s,b = linfit(t[100:end],x[100:end]);
      vvals[i] = s;
    end

    dir = string("/home/data/ie355/Documents/code/dumbbell/data/",saveadd,".jld");
    save(dir,"parnum",runpar,"parvals",parvals,"vvals",vvals);
end

function savetrajs(name, pars, runpar, parmin, parmax, n, npart = 10, nsteps = 10^7, dt = .001, undersamp = 100)
    parvals = LinRange(parmin,parmax,n);
    npoints = Int(nsteps / undersamp);
    xvals = zeros(npoints,n);
    t = zeros(npoints);
    for i = 1:n
      pars[runpar] = parvals[i];
      t,x = sim.parallelcall(npart,[1,0],dt,nsteps,pars,undersamp);
      xvals[:,i] = x;
    end

    dir = string("/home/data/ie355/Documents/code/dumbbell/data/individual_paths/",name,".jld");
    save(dir,"t",t,"xvals",xvals,"parvals",parvals);
end


function currentset()
  #print("First, get some trajectories")
  #savetrajs("basic",[0.1,0.1,1.0,0.1,0,2],2,0.1,3,5);

  #print("Next, get some t-runs. First, set cold particle at very cold and vary hot")
  #singlevarrun("temp_runs/colder",[0.1,0.1,1.0,0.1,0,2],2,0.1,1.0,40);
  #print("Then, a t-run with a warmer cold particle")
  #singlevarrun("temp_runs/warmer",[0.5,0.5,1.0,0.1,0,2],2,0.5,1.5,40);
  #print("And then a t-run with a stiffer spring, where our theory should be better")
  #singlevarrun("temp_runs/stiffer",[0.1,0.1,10.0,0.1,0,2],2,0.1,1.0,40);
  #print("And then a t-run with a stifferer spring, where our theory should be better")
  #singlevarrun("temp_runs/stifferer",[0.1,0.1,20.0,0.1,0,2],2,0.1,1.0,40);

  #print("Then, a soft-k run with colder cold particle")
  #singlevarrun("k_runs/soft",[0.1,1.0,1.0,0.1,0,2],3,.05,1,100);
  #print("Then, a stiff-k run with colder cold particle")
  #singlevarrun("k_runs/stiff",[0.1,1.0,1.0,0.1,0,2],3,1.0,20,100);

  print("Then, a gravity run with small range and cold particle")
  singlevarrun("grav_runs/small_cold_dt0",[0.1,0.1,1.0,0.1,0,2],5,0,1.5,100);
  singlevarrun("grav_runs/small_cold_dt1",[0.1,1.1,1.0,0.1,0,2],5,0,1.5,100);
  singlevarrun("grav_runs/small_cold_dt2",[0.1,2.1,1.0,0.1,0,2],5,0,1.5,100);
  singlevarrun("grav_runs/small_cold_dt3",[0.1,3.1,1.0,0.1,0,2],5,0,1.5,100);
  #print("Then, a gravity run with small range and hot particle")
  #singlevarrun("grav_runs/small_hot_dt0",[0.5,0.5,1.0,0.1,0,2],5,0,.01,100);
  #singlevarrun("grav_runs/small_hot_dt1",[0.5,1.5,1.0,0.1,0,2],5,0,.01,100);
  #singlevarrun("grav_runs/small_hot_dt2",[0.5,2.5,1.0,0.1,0,2],5,0,.01,100);
  #singlevarrun("grav_runs/small_hot_dt3",[0.5,3.5,1.0,0.1,0,2],5,0,.01,100);
  #print("Then, a gravity run with large range and cold particle")
  #singlevarrun("grav_runs/large_cold_dt0",[0.1,.1,1.0,0.1,0,2],5,0,8,100);
  #singlevarrun("grav_runs/large_cold_dt1",[0.1,1.1,1.0,0.1,0,2],5,0,8,100);
  #singlevarrun("grav_runs/large_cold_dt2",[0.1,2.1,1.0,0.1,0,2],5,0,8,100);
  #singlevarrun("grav_runs/large_cold_dt3",[0.1,3.1,1.0,0.1,0,2],5,0,8,100);

  #print("A few t-runs with large k-values")
  #singlevarrun("temp_runs/k30",[0.1,0.1,30,0.1,0,2],2,0.1,1.0,40);
  #singlevarrun("temp_runs/k50",[0.1,0.1,50,0.1,0,2],2,0.1,1.0,40); #take smaller timesteps and undersample... 
  #singlevarrun("temp_runs/k70",[0.1,0.1,70,0.1,0,2],2,0.1,1.0,40);
  #singlevarrun("temp_runs/k100",[0.1,0.1,100,0.1,0,2],2,0.1,1.0,40);
  #
  #print("\n one t-run which goes much further in k")
  #singlevarrun("temp_runs/k30_bigt",[0.1,0.1,30,0.1,0,2],2,0.1,5.0,40);
end

end
