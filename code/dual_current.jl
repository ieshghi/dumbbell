#######################################################################
#
#	File: dual_current
#	Description: from the probability distributions calculated by sim_base.jl, 
#				 calculates currents on a dual lattice which sits on pathological points
#				 where J is well defined. This lattice may have a slightly different number
#				 of pixels from the original lattice
#
#	Author: Iraj Eshghi
#
#######################################################################

module dual_current
include("base_sim.jl")
using Interpolations

function interp(f::Array{Float64,1})
	return (f+circshift(f,1))/2
end

function newforce(xo::Array{Float64,1},l::Float64)
	fold = simbase.ext_force.(xo,l)
	return interp(fold)
end

function gen_dual_lat(lat::Array{Float64,2},x::Array{Float64,1},y::Array{Float64,1})
	n = length(x)
	xn = Array(range(0.,stop = 1.0,length = n+1)[1:end-1])
 	lat_n = zeros(size(xn*y'))
	@inbounds for j = 1:length(y)
		lat_n[:,j] = interp(lat[:,j])
	end
	
	dy = y[2]-y[1]
	dx = xn[2]-xn[1]
	lat_dy = 1/(2*dy)*(circshift(lat_n,[0,-1])-circshift(lat_n,[0,1]))
	lat_dx = 1/dx*(lat-circshift(lat,[1,0]))
	
	return lat_n,lat_dx,lat_dy,xn,y
end

function force(x::Array{Float64,1},xo::Array{Float64,1},y::Array{Float64,1},l::Float64,k::Float64)
	newf = newforce(xo,l)
	xf = -newf*ones(size(y')) + k*ones(size(x))*y'	
	yf = newf*ones(size(y')) - 2*k*ones(size(x))*y'	

	return xf,yf
end

function calc_curr(c::Array{Float64,2},pars::Array{Float64,1},xo::Array{Float64,1},y::Array{Float64,1})
	cn,cdx,cdy,xn,y = gen_dual_lat(c,xo,y)
	l = pars[4]
	k = pars[3]
	t1 = pars[1]
	t2 = pars[2]
	fx,fy = force(xn,xo,y,l,k)
	D = [[t1,t1],[t1,t1+t2]]
	
	jx = fx.*cn.-D[1][1].*cdx.+D[1][2].*cdy
	jy = fy.*cn.+D[2][1].*cdx.-D[2][2].*cdy

	return jx,jy,cn,xn,y,cdx,cdy	
end


end
