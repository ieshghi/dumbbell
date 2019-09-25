module current
include("sim.jl")
include("histutils.jl")
using JLD
using DistributedArrays

# pars = [t1,t2,k,lambda,g,potchoice]

#calculates probability density heat maps from runs
function getprob(npart,pars,nx,ny,yl,dt,nsteps,coupling)
	xbins = linrang(0,1,nx);
	ybins = linrang(-yl,yl,ny);
	pairs = histutils.genmesh(xbins,ybins)[1];
	hists = @DArray[sim.bintrajs([0,0],dt,nsteps,pars,xbins,ybins,coupling) for j = 1:npart];
	hist = sum(hists);
	prob = hist./sum(hist); #convert histogram to probability distribution
	return hists
end

function getrun(npart,pars,nx,ny,yl,dt = 0.001,nsteps = 10^7, coupling = [1,0])
end

function linrang(a,b,n) #like linrange, but without the end. Get it?
	return range(a, stop=b, length=n+1)[1:end-1]
end


end
