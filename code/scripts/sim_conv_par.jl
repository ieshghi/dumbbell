n = parse(Int,ARGS[1]) #number of timesteps
nx = parse(Int,ARGS[2]) #grid size
id = ARGS[3]

include("../base_sim.jl")
using HDF5

pars = [0.1,0.1,1.0,0.7,0];
l,x,y = simbase.run_and_bin([0,0].*1.0,.01,pars,n,1000,1000,5.0,0);

h5open(string("/home/ie355/code/dumbbell/data/hists/",nx,"run",id,"ifext",0,".h5"),"w") do file
	write(file,"l",l)
	write(file,"x",x)
	write(file,"y",y)
end
