module current_runs
include("current.jl")
using HDF5

function calcnet(x,y,p,jx,jy,pars)
	dx = unique(x)[2]-unique(x)[1];
	dy = unique(y)[2]-unique(y)[1];
	jleft = jx[x.== dx];
	jright = jx[x.== 1-dx];

	return dy*sum(jleft[jleft.<0]+jright[jright.>0])
end

function prep_vals()

	tc = 0
	thvals = [0.1,0.3,0.5,0.8,1.1];
	for i = thvals
		pars = [tc,th,1.0,0.3,0,2]
		x,y,p,jx,jy = current.getrun(10,pars,100,100,3,.01,10^7)
		current.savestuff(x,y,p,jx,jy,pars,string("extra/tc",tc,"th",i))
	end
	tc = 0.01
	for i = thvals
		pars = [tc,i,1.0,0.3,0,2]
		x,y,p,jx,jy = current.getrun(10,pars,100,100,3,.01,10^7)
		current.savestuff(x,y,p,jx,jy,pars,string("extra/tc",tc,"th",i))
	end
	h5write("/home/ieshghi/Documents/code/dumbbell/data/hists/extra/params.h5","pars",[10,pars,100,100,3,.01,10^7])
end
end
