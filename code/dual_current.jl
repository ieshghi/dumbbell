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

function edges(xe::Float64,xi::Float64,fe::Float64,fi::Float64,xm::Array{Float64,1})
	xm[xm.<xe] .+= 1	
	xi += 1
	itp = LinearInterpolation([xe,xi],[fe,fi],extrapolation_bc=Flat())
	return itp.(xm)
end

function interp(xo::Array{Float64,1},xn::Array{Float64,1},f::Array{Float64,1})
	fn = zeros(size(xn))
	itp = LinearInterpolation(xo,f,extrapolation_bc=Flat())
	q = (xn.>xo[1]).&(xn.<xo[end])
	fn[q] = itp.(xn[q])
	fn[.!q] = edges(xo[end],xo[1],f[end],f[1],xn[.!q])		
	return fn
end

function dual_lat_grad(lat::Array{Float64,2},lat_n::Array{Float64,2},xo::Array{Float64,1},xn::Array{Float64,1},y::Array{Float64,1})
	dy = y[2]-y[1]
	dxo = xo[2]-xo[1]
	dxn = xn[2]-xn[1]
	lat_dy = 1/(2*dy)*(circshift(lat_n,[0,-1])-circshift(lat_n,[0,1]))
	lat_dx = zeros(size(lat_n))
	@inbounds for i = 1:length(y)
		lat_dx[:,i] = 5.0/dxo*(interp(xo,xn.+dxo/10,lat[:,i])-interp(xo,xn.-dxo/10,lat[:,i]))
	end
	return lat_dx,lat_dy
end

function gen_dual_lat(lat::Array{Float64,2},x::Array{Float64,1},y::Array{Float64,1})
	nnew = Int(10*round(length(x)*1.0/10)) #round to the nearest multiple of 10. Combined with l being rounded to the nearest 1/10th, this guarantees
# that the j-lattice will hit the pathological points

	xj = Array(range(0.,stop = 1.0,length = nnew+1)[1:end-1])
 	newlat = zeros(size(xj*y'))
	@inbounds for j = 1:length(y)
		newlat[:,j] = interp(x,xj,lat[:,j])
	end
	
	return newlat,xj,y
end

function up(xj::Array{Float64,1},xo::Array{Float64,1},l::Float64)
	uold = simbase.ext_force.(xo,l)
	return interp(xo,xj,uold)
end

function boltzmann(x::Array{Float64,1},y::Array{Float64,1},pars::Array{Float64,1}) #boltzmann distribution for a set of parameters. Can be useful for tests
	l = pars[4]
	k = pars[3]
	t1 = pars[1]
	t2 = pars[2]
	tb = mean([t1,t2])
	if t1!=t2
		print("Warning: boltzmann temperature is not well defined! Taking average.\n\n")
	end
	xarr = x*ones(size(y'))
	yarr = ones(size(x))*y'

	pot = k/2*yarr.^2 + simbase.ext_pot.(xarr,l)
	return exp.(-pot./tb) #We can always normalize later, shouldn't matter much
end

function force(x::Array{Float64,1},xo::Array{Float64,1},y::Array{Float64,1},l::Float64,k::Float64)
	xf = zeros(size(x*y'))
	yf = zeros(size(x*y'))
	
	uprime = up(x,xo,l)

	#@inbounds for i = 1:length(y)
	#	xf[:,i] = -uprime.+k*y[i]
	#	yf[:,i] = uprime.-2*k*y[i]
	#end
	@inbounds for i = 1:length(x)
		xf[i,:] = -k*x[i]
		yf[i,:] = -2*k*x[i]
	end
	
	return xf,yf
end

function calc_curr(c::Array{Float64,2},pars::Array{Float64,1},xo::Array{Float64,1},y::Array{Float64,1})
	cn,xn,y = gen_dual_lat(c,xo,y)
	cdx,cdy = dual_lat_grad(c,cn,xo,xn,y)
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
