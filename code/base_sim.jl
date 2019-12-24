#
#	file: base_sim
#	Description: runs the brownian dynamics simulation for a two temperature dumbbell on a parabolic ratchet
#	
#	Author: Iraj Eshghi
#
#

module simbase

function ext_force(x::Float64,l::Float64)
	if x<=l
		return 2*(l.-x[x.<=l])./(l^2)
	else
		return 2*(l.-x[x.>l])./((1-l)^2)
	end
end

function timestep(x0::Array{Float64,2},dt::Float64,pars::Array{Float64,5},coupling::Array{Bool,2},noise::Float64)
	b = sqrt.(2*dt*[pars[1],pars[2]])
	x = x0;
	a = rhs(x,pars,coupling)
	x += a*dt+b.*noise
	return x
end

function rhs(x::Array{Float64,2},pars::Array{Float64,5},coupling::Array{Bool,2})
	l = pars[4]
	k = pars[3]
	g = pars[5]	
	fex = coupling.*ext_force(x,l)
	fspr = k*(x[1]-x[2]).*[-1,1]
	fg = g*ones(size(x))
	return fex+fspr+fg
end

function run_and_bin(x0::Array{Float64,2},dt::Float64,pars::Array{Float64,5},nsteps::Int,nx::Int,ny::Int,ymax::Float64)
	xvals,yvals,lat = gen_c_lattice(nx,ny,ymax)
	@inbounds for i = 1:nsteps
			
	end
end

function place_in_lattice!(x::Array{Float64,2},lat::Array{Float64,2},)
	
end

function gen_c_lattice(nx::Int,ny::Int,ymax::Float64,l::Float64)
	nx1 = Int(l*nx)
	nx2 = nx-nx1
	xvals1 = range(0,stop=l,length=nx1+1)[1:end-1]
	xvals2 = range(l,stop=1,length=nx2+2)[2:end-1]
	xvals = vcat(xvals1,xvals2)
	yvals = range(-ymax,stop=ymax,length=ny)
	lat = zeros(nx,ny)
	
	return lat,xvals,yvals
end
end
