n = ARGS[1] #number of timesteps

include("../base_sim.jl")
using HDF5

pars = [0.1,1.0,1.0,0.7,0];
l,x,y = simbase.run_and_bin([0,0].*1.0,.01,pars,n,100,100,5.0);

