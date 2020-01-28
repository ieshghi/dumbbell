#
#	file: base_sim
#	Description: runs the brownian dynamics simulation for a two temperature dumbbell on a parabolic ratchet
#	
#	Author: Iraj Eshghi
#
#

module simbase

function ext_pot(x::Float64,l::Float64)
	xm = x%1-(sign(x)-1)/2

	if xm<=l
		return (x-l)^2/(l^2)
	else
		return (x-l)^2/((1-l)^2)
	end
end

function ext_force(x::Float64,l::Float64)
	if x==0
		xm = 0;
	else
		xm = x%1-(sign(x)-1)/2
	end
	if xm<l
		return 2*(xm-l)./(l^2)
	else
		return 2*(xm-l)./((1-l)^2)
	end
end

function timestep!(x::Array{Float64,1},dt::Float64,pars::Array{Float64,1},noise::Array{Float64,1},ifext::Int)
	b = sqrt.(2*dt*[pars[1],pars[2]])
	a = rhs(x,pars,ifext)
	hop = a[1]*dt+b[1]*noise[1]
	x[1] += hop
	x[2] += a[2]*dt+b[2]*noise[2]-b[1]*noise[1]
	return hop
end

function rhs(x::Array{Float64,1},pars::Array{Float64,1},ifext::Int)
	l = pars[4]
	k = pars[3]
	g = pars[5]	
	fex = ifext*ext_force(x[1],l)*[-1,1]
	fspr = k*x[2].*[1,-2]
	fg = g*[1,0]
	return fex+fspr+fg
end

function run_and_bin(x0::Array{Float64,1},dt::Float64,pars::Array{Float64,1},nsteps::Int,nx::Int,ny::Int,ymax::Float64,ifext::Int=1)
	pars[4] = 0.1*round(pars[4]*10) #Rounds to nearest 1/10th. we don't need better than this precision on this point. Helps lattice construction
	lat,xvals,yvals = gen_c_lattice(nx,ny,ymax,pars[4])
	noise = randn(nsteps,2)
	x = x0
	maxhop = 0
	@inbounds for i = 1:nsteps
			hop = timestep!(x,dt,pars,noise[i,:],ifext)
			place_in_lattice!(x,lat,xvals,yvals)
			maxhop += abs(hop)
	end
	nlat = norm_lat(lat,xvals,yvals)
	return nlat,xvals,yvals,maxhop/nsteps
end

function place_in_lattice!(x::Array{Float64,1},lat::Array{Float64,2},xvals::Array{Float64,1},yvals::Array{Float64,1})
	xm = x[1]%1 - (sign(x[1])-1)/2
	
	dx = xvals[2]-xvals[1]
	dy = yvals[2]-yvals[1]
	xind = (abs.(xm.-xvals).<(dx/2))
	yind = (abs.(x[2].-yvals).<(dy/2))
	
	lat[xind,yind] .+= 1
end

function gen_c_lattice(nx::Int,ny::Int,ymax::Float64,l::Float64)
	xvals_1 = range(0.0,stop = 1.0,length = nx+1)[1:end-1]
	dx = xvals_1[2]-xvals_1[1]
	xvals_2 = xvals_1.+dx/2
	nshift = sum(xvals_2.>1)
	xvals = Array(circshift(xvals_2.%1,nshift))
	yvals = Array(range(-ymax,stop = ymax,length = ny))
	lat = zeros(nx,ny)
	
	if sum(xvals.==l)>0
		print("c lattice coincides with bad point!\n")
	end
	return lat,xvals,yvals
end

function norm_lat(lat::Array{Float64,2},x::Array{Float64,1},y::Array{Float64,1})
	dx = x[2]-x[1]
	dy = y[2]-y[1]
	ntot = sum(lat)
	return lat./(ntot*dy*dx)
end

end
