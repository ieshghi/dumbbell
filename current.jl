module current
include("sim.jl")
include("histutils.jl")
using JLD
using DistributedArrays
using Plots
pyplot()

# pars = [t1,t2,k,lambda,g,potchoice]

#calculates probability density heat maps from runs
function getprob(npart,pars,nx,ny,yl,dt,nsteps,coupling)
	xbins = linrang(0,1,nx);
	ybins = linrang(-yl,yl,ny);
	pairs = histutils.genmesh(xbins,ybins)[1];
	hists = @DArray[sim.bintrajs([rand(),rand()],dt,nsteps,pars,xbins,ybins,coupling) for j = 1:npart];
	hist = sum(hists);
	prob = hist./sum(hist); #convert histogram to probability distribution
	return prob,pairs[:,1],pairs[:,2]
end

function calc_curr(x,y,p,pars,dx,dy)
	xm,ym,pm = histutils.shape_prob(x,y,p);
	f1 = -sim.gpot.(xm,pars[4],2)+pars[3].*ym;
	f2 = sim.gpot.(xm,pars[4],2)-2*pars[3].*ym;
	dxp = findif(pm,dx,1);
	dyp = findif(pm,dy,2);
	tx = pars[1];
	ty = pars[2];
	D = [[tx,tx],[tx,tx+ty]];
	jx = f1.*pm-(D[1][1]).*dxp-(D[1][2]).*dyp;
	jy = f2.*pm-(D[2][1]).*dxp-(D[2][2]).*dyp;
	
	return xm,ym,pm,jx,jy
end

function getrun(npart,pars,nx,ny,yl,dt = 0.01,nsteps = 10^7, coupling = [1,0])
	p,x,y = getprob(npart,pars,nx,ny,yl,dt,nsteps,coupling);
	dx = unique(x)[2]-unique(x)[1];
	dy = unique(y)[2]-unique(y)[1];
	xm,ym,pm,jx,jy = calc_curr(x,y,p,pars,dx,dy);

	return xm,ym,pm,jx,jy
end

function linrang(a,b,n) #like linrange, but without the end. Get it?
	return range(a, stop=b, length=n+1)[1:end-1]
end

function findif(f,d,dir)
	if ndims(f)>1
	if dir == 1
		a = [0,1];
	elseif dir == 2
		a = [1,0];
	end
	else
		a = 1;
	end
	return 1/(2*d)*(circshift(f,-a)-circshift(f,a));
end

function plotstuff(x,y,p,jx,jy,pars)
	rat = length(unique(x))/length(unique(y));
	prplot = heatmap(p,aspect_ratio=rat);
	xl = linrang(pars[4],1+pars[4],1000);
	potplot = plot(xl,sim.potential.(xl,pars[4],2));
	boltzplot = heatmap(boltz(x,y,pars),aspect_ratio=rat);
	plot(prplot,boltzplot,layout=(2,1),legend=false,size=(600,600))
end

function boltz(x,y,pars)
	t1 = pars[1];
	t2 = pars[2];
	k = pars[3];
	l = pars[4];
	xsc = mod.(x.+pars[4],1);
	pe = sim.potential.(xsc,pars[4],2);
	spr = (y.^2).*1/2*k;
	E = pe+spr;
	p = exp.(-2*E/(t1+t2));
	return p./(sum(p))
end

end
