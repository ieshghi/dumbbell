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

function gen_dual_lat(lat::Array{Float64,2},x::Array{Float64,1},y::Array{Float64,1})
	n = length(x)
	xn = Array(range(0.,stop = 1.0,length = n+1)[1:end-1]) #This works when n is a multiple of 10 and the patching point is at a multiple of 0.1. No need to generalize further. These constraints are already enforced in the base_sim module.

	lat_n = 1/2*(lat + circshift(lat,[1,0]))#Interpolate the probability distribution onto the new lattice. At every point, the new lattice values are the average of the old lattice neighbour values.
	
	dy = y[2]-y[1]
	dx = xn[2]-xn[1]
	
	lat_dy = 1/(2*dy)*(circshift(lat_n,[0,-1])-circshift(lat_n,[0,1])) #Calculate derivatives with second-order centered finite difference. Assume periodicity (technically untrue in y direction but the probability distribution is basically 0 at large y values so this shouldn't matter.
	lat_dx = 1/(dx)*(lat-circshift(lat,[1,0])) #This is the slope at the midpoint between two old lattice points.
	
	return lat_n,lat_dx,lat_dy,xn,y
end

function force(x::Array{Float64,1},xo::Array{Float64,1},y::Array{Float64,1},l::Float64,k::Float64)
	oldf = simbase.ext_force.(xo,l) #External force on old lattice
	newf = simbase.ext_force.(x,l) #External force on new lattice (has one discontinuity at 0!)
	newf[1] = 1/2*(oldf[end]+oldf[1]) #the only problematic point for the force, the discontinuity at 0, is smoothed by taking the value of the force there to be the average of the neighbours on the old lattice 

	xf = -newf*ones(size(y')) + k*ones(size(x))*y'
	yf = newf*ones(size(y')) - 2*k*ones(size(x))*y' #Force fields on new lattice as given by Langevin equations

	return xf,yf
end

function calc_curr(c::Array{Float64,2},pars::Array{Float64,1},xo::Array{Float64,1},y::Array{Float64,1})
	l = pars[4]
	k = pars[3]
	t1 = pars[1]
	t2 = pars[2]
	D = [[t1,-t1],[-t1,t1+t2]]

	cn,cdx,cdy,xn,y = gen_dual_lat(c,xo,y)
	fx,fy = force(xn,xo,y,l,k)
	
	jx = fx.*cn.-D[1][1].*cdx.-D[1][2].*cdy #See Fokker-Planck equation for the langevin equations.
	jy = fy.*cn.-D[2][1].*cdx.-D[2][2].*cdy

	return jx,jy,cn,xn,y,fx,fy,cdx,cdy	
end

end
