module current
include("sim.jl")
include("histutils.jl")
using DistributedArrays
using Plots
pyplot()
using HDF5
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

function calc_curr(x,y,p,pars,dx,dy,ifshape=true)
	if ifshape
		xm,ym,pm = histutils.shape_prob(x,y,p);
	else
		xm = x;
		ym = y;
		pm = p;
	end
	f1 = -sim.gpot.(xm.+pars[4],pars[4],2)+pars[3].*ym;
	f2 = sim.gpot.(xm.+pars[4],pars[4],2)-2*pars[3].*ym;
	dxp = findif(pm,dx,1);
	dyp = findif(pm,dy,2);
	tx = pars[1];
	ty = pars[2];
	D = [[tx,tx],[tx,tx+ty]];
	jx = f1.*pm-(D[1][1]).*dxp-(D[1][2]).*dyp;
	jy = f2.*pm+(D[2][1]).*dxp-(D[2][2]).*dyp;

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

function curl(fx,fy,dx,dy);
	return findif(fy,dx,1)-findif(fx,dy,2);
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
	theme(:default);
	col = (:plasma);
	usamp = 2;
	mult = 5;
	jxl = LinRange(minimum(jx),maximum(jx),50);
	jyl = LinRange(minimum(jy),maximum(jy),50);
	rat = (maximum(x)-minimum(x))/(maximum(y)-minimum(y));
	xl = linrang(pars[4],1+pars[4],1000);
	xs = linrang(0,1,length(p[:,1]));
	ys = LinRange(minimum(y),maximum(y),length(p[1,:]));
	dx = xs[2]-xs[1];
	dy = ys[2]-ys[1];
	cf = curl(jx,jy,dx,dy);
	cfl = LinRange(minimum(cf),maximum(cf),50);

	prplot = contourf(xs,ys,p,aspect_ratio=rat,c=col,colorbar=true,title="Potential");
	potplot = plot(xl,sim.potential.(xl,pars[4],2));
	boltzplot = contourf(xs,ys,boltz(x,y,pars),c=col,aspect_ratio=rat);
	jxplot = contourf(xs,ys,jx,aspect_ratio = rat,colorbar=true,levels=collect(jxl),title="X-current");
	jyplot = contourf(xs,ys,jy,aspect_ratio = rat,colorbar=true,levels=collect(jyl),title="S-current");
	qplot = quiv_under(x,y,jx,jy,4,100,rat);	
	curlplot = contourf(xs,ys,cf,aspect_ratio=rat,colorbar=true,levels=collect(cfl),title="curl");
	sliceplot = plot(unique(y),[jx[x.==0.05][2:end],jx[x.==0.95][2:end]],title="Current at s=0",label=["jx","jy"])
	
	print(sum(jx[x.==0.05]))

#	l = @layout [ grid(2,2)
#				  b{0.2h} ];
	plot(prplot,qplot,jxplot,jyplot,curlplot,sliceplot,layout = (3,2),legend=false,size=(800,800),reuse=false)
end

function plotgauss(x,y)
	p = exp.(-x.^2-y.^2-x.*y);
	pnorm = sum(p);
	p = p./pnorm;
	pars = [.01,1,1,0.3,0,2];
	dx = unique(x)[2]-unique(x)[1];
	dy = unique(y)[2]-unique(y)[1];
	rat = (maximum(x)-minimum(x))/(maximum(y)-minimum(y));
	theme(:default);
	col = (:plasma);
	xs = linrang(0.01,1,length(p[:,1]));
	ys = LinRange(minimum(y),maximum(y),length(p[1,:]));

	x,y,p,jx,jy = calc_curr(x,y,p,pars,dx,dy,false);	
	jxl = LinRange(minimum(jx),maximum(jx),50);
	jyl = LinRange(minimum(jy),maximum(jy),50);

	prplot = contourf(xs,ys,p,aspect_ratio=rat,c=col,colorbar=true,title="potential");
	jxplot = contourf(xs,ys,jx,aspect_ratio = rat,colorbar=true,levels=collect(jxl),title="X-current");
	jyplot = contourf(xs,ys,jy,aspect_ratio = rat,colorbar=true,levels=collect(jyl),title="S-current");
	qplot = quiv_under(x,y,jx,jy,5,200,rat);	
	sliceplot = plot(unique(x),[jx[y.==0.],jy[y.==0.]],title="current at s=0")
	l = @layout [ grid(2,2)
				  b{0.2h} ];
	plot(prplot,qplot,jxplot,jyplot,sliceplot,layout=l,legend=false,size=(700,700))
end

function quiv_under(x,y,jx,jy,u,m,rat = 1)
	x = x[1:u:end,1:u:end];
	y = y[1:u:end,1:u:end];
	jx = m*jx[1:u:end,1:u:end];
	jy = m*jy[1:u:end,1:u:end];
	return quiver(x,y,quiver=maketuples(jx,jy),aspect_ratio = rat,arrow=arrow(.3,.1),linecolor=:steelblue);
end

function savestuff(x,y,p,jx,jy,pars,name)
	h5write(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"x",x);
	h5write(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"y",y);
	h5write(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"p",p);
	h5write(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"jx",jx);
	h5write(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"jy",jy);
	h5write(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"pars",pars);
end
function loadstuff(name)
	x = h5read(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"x");
	y = h5read(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"y");
	p = h5read(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"p");
	jx = h5read(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"jx");
	jy = h5read(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"jy");
	pars = h5read(string("/home/ieshghi/Documents/code/dumbbell/data/hists/",name,".h5"),"pars");
	return x,y,p,jx,jy,pars
end

function maketuples(a1,a2)
	s = size(a1);
	a = Array{Tuple{Float64,Float64},2}(undef,s[1],s[2]);
	for i = 1:s[1]
		for j = 1:s[2]
			a[i,j] = (a1[i,j],a2[i,j]);
		end
	end
	return a
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
