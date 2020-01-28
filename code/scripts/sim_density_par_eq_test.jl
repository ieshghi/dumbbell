n = parse(Int,ARGS[1]) #number of timesteps
id = ARGS[2]

include("../base_sim.jl")
using HDF5

pars = [0.1,0.1,1.0,0.7,0];
l,x,y = simbase.run_and_bin([0,0].*1.0,.001,pars,n,1000,1000,5.0,0); #tells the simulation not to consider the ratchet force

h5open(string("/home/ie355/code/dumbbell/data/hists/run",id,".h5"),"w") do file
	write(file,"l",l)
	write(file,"x",x)
	write(file,"y",y)
end
