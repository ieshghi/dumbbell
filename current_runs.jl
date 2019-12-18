module current_runs
include("current.jl")
using Plots
pyplot()
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
		pars = [tc,i,1.0,0.3,0,2]
		x,y,p,jx,jy = current.getrun(10,pars,100,100,3,.01,10^7)
		current.savestuff(x,y,p,jx,jy,pars,string("extra/tc",tc,"th",i))
	end
	tc = 0.01
	for i = thvals
		pars = [tc,i,1.0,0.3,0,2]
		x,y,p,jx,jy = current.getrun(10,pars,100,100,3,.01,10^7)
		current.savestuff(x,y,p,jx,jy,pars,string("extra/tc",tc,"th",i))
	end
end

function compare_vals()
    
    thvals = LinRange(0.1,1,40)
    v0 = h5read["/home/data/ie355/Documents/code/dumbbell/data/temp_runs/decoup_k1_t0.h5", "vvals"]
    v1 = h5read["/home/data/ie355/Documents/code/dumbbell/data/temp_runs/decoup_k1_t01.h5", "vvals"]
    

    thvals_curr = [0.1,0.3,0.5,0.8,1.1]
    ncurr = length(thvals_curr)
    v0_curr = zeros(ncurr);
    v1_curr = zeros(ncurr);
    
    for i = 1:ncurr
        x,y,p,jx,jy = current.loadstuff(string("extra/tc",0,"th",th[i]))
        v0_curr[i] = calcnet(x,y,p,jx,jy,pars)
        x,y,p,jx,jy = current.loadstuff(string("extra/tc",0.01,"th",th[i]))
        v1_curr[i] = calcnet(x,y,p,jx,jy,pars)
    end

    plot(thvals,v0,"-")
    plot(thvals_curr,v0_curr,".")

end


end
